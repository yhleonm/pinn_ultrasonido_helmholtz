"""
physics/diff_das_eikonal.py
===========================
Módulo Híbrido Física + Deep Learning para Ultrasonido Mamario (Fase 2):
1. FactorizedEikonalPINN: Red PINN que resuelve la ecuación Eikonal factorizada
   |∇τ(x, z)|^2 = 1 / c(x, z)^2 acoplada al mapa de velocidad del sonido c(x, z).
2. DifferentiablePlaneWaveDAS: Beamformer Delay-and-Sum (DAS) 100% diferenciable en GPU
   para adquisiciones de Ondas Planas (CUBDL), con:
   - Filtrado paso-banda en el dominio de la frecuencia (GPU FFT).
   - Señal analítica compleja (Transformada de Hilbert vía FFT en GPU).
   - Apodización dinámica de recepción por F-number con ventana de Hanning.
   - Integración de retardos Eikonal heterogéneos Δτ(x, z; θ_m) a partir de c(x, z)
     para corrección de aberración de fase.
   - Cálculo del Factor de Coherencia (CF) y pérdida diferenciable de auto-enfoque.
"""

from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class FactorizedEikonalTravelTimeNet(nn.Module):
    """
    Red neuronal informada por la física (PINN) que parametriza el factor de
    desviación adimensional tau_factor(x, z) respecto al tiempo de vuelo homogéneo tau_0(x, z).
    Evita la singularidad de fuente en z = 0 y permite calcular el residuo PDE Eikonal exacto
    mediante diferenciación automática (torch.autograd.grad).
    """

    def __init__(self, hidden_dim: int = 64, num_layers: int = 4, c0: float = 1540.0):
        super().__init__()
        self.c0 = float(c0)
        layers = []
        in_dim = 2  # (x_norm, z_norm)
        for i in range(num_layers - 1):
            layers.append(nn.Linear(in_dim if i == 0 else hidden_dim, hidden_dim))
            layers.append(nn.SiLU())
        layers.append(nn.Linear(hidden_dim, 1))
        self.net = nn.Sequential(*layers)

        # Inicializar última capa cerca de cero para que alpha(x, z) ≈ 1.0 al inicio
        with torch.no_grad():
            self.net[-1].weight.mul_(0.01)
            self.net[-1].bias.zero_()

    def forward(
        self,
        coords_norm: torch.Tensor,
        z_center: float = 0.018,
        scale_z: float = 0.0075,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            coords_norm: Coordenadas normalizadas en [-1, 1]^2, shape [N, 2] (requiere grad).
            z_center: Profundidad central de la ROI en metros.
            scale_z: Semiancho axial de la ROI en metros.
        Returns:
            tau: Tiempo de vuelo vertical Eikonal desde z=0, shape [N, 1] (segundos).
            alpha_tau: Factor multiplicativo adimensional ~ 1.0, shape [N, 1].
        """
        raw = self.net(coords_norm)
        # Acotar desviación de tiempo de vuelo a +/- 15% (rango físico mamario)
        alpha_tau = 1.0 + 0.15 * torch.tanh(raw)
        # z_phys derivado de coords_norm para que autograd aplique la regla del producto exacta
        z_phys = torch.clamp(z_center + coords_norm[:, 1:2] * scale_z, min=1e-4)
        tau_0 = z_phys / self.c0
        tau = tau_0 * alpha_tau
        return tau, alpha_tau


def compute_eikonal_pde_loss(
    tau: torch.Tensor,
    coords_norm: torch.Tensor,
    c_pred: torch.Tensor,
    scale_x: float,
    scale_z: float,
    c0: float = 1540.0,
) -> torch.Tensor:
    """
    Calcula el residuo PDE de la Ecuación Eikonal 2D:
        |∇_x τ|^2 + |∇_z τ|^2 = 1 / c(x, z)^2
    normalizado por s_0^2 = (1 / c_0)^2 para estabilidad numérica.

    Args:
        tau: Tiempo de vuelo predicho [N, 1] en segundos.
        coords_norm: Coordenadas normalizadas [N, 2] con requires_grad=True.
        c_pred: Mapa de velocidad del sonido predicho [N, 1] en m/s.
        scale_x: Semiancho físico de normalización en x (metros), dx_phys = scale_x * dx_norm.
        scale_z: Semiancho físico de normalización en z (metros), dz_phys = scale_z * dz_norm.
        c0: Velocidad de sonido de referencia (m/s).
    """
    grad_norm = torch.autograd.grad(
        outputs=tau.sum(),
        inputs=coords_norm,
        create_graph=True,
        retain_graph=True,
    )[0]

    # Regla de la cadena: dτ/dx_phys = (dτ/dx_norm) / scale_x
    dtau_dx = grad_norm[:, 0:1] / scale_x
    dtau_dz = grad_norm[:, 1:2] / scale_z

    grad_sq = dtau_dx ** 2 + dtau_dz ** 2
    slowness_sq = 1.0 / (c_pred ** 2)
    s0_sq = 1.0 / (c0 ** 2)

    # Residuo Eikonal adimensionalizado por s0^2
    residual = (grad_sq - slowness_sq) / s0_sq
    return torch.mean(residual ** 2)


def gpu_hilbert_and_bandpass(
    rf_tensor: torch.Tensor,
    fs: float,
    f_low: float = 2.0e6,
    f_high: float = 11.0e6,
) -> torch.Tensor:
    """
    Convierte datos RF reales [Nt, Nr, Ntx] en señal analítica compleja filtrada
    paso-banda [Nt, Nr, Ntx] directamente en GPU usando FFT vectorizada.
    Elimina DC offset y ruido fuera de la banda del transductor.
    """
    Nt = rf_tensor.shape[0]
    device = rf_tensor.device

    # Remover media temporal (DC offset) por canal
    rf_zero_mean = rf_tensor - torch.mean(rf_tensor, dim=0, keepdim=True)

    # FFT a lo largo del eje temporal (dim=0)
    Rf_fft = torch.fft.fft(rf_zero_mean, dim=0)
    freqs = torch.fft.fftfreq(Nt, d=1.0 / fs, device=device)

    # Multiplicador de Hilbert (2 para f > 0, 1 para f = 0, 0 para f < 0)
    hilbert_mult = torch.zeros(Nt, dtype=torch.float32, device=device)
    hilbert_mult[0] = 1.0
    if Nt % 2 == 0:
        hilbert_mult[1 : Nt // 2] = 2.0
        hilbert_mult[Nt // 2] = 1.0
    else:
        hilbert_mult[1 : (Nt + 1) // 2] = 2.0

    # Ventana paso-banda suave (Tukey / Butterworth de orden alto en frecuencia)
    abs_f = torch.abs(freqs)
    bp_mask = torch.sigmoid((abs_f - f_low) / 0.4e6) * torch.sigmoid((f_high - abs_f) / 0.6e6)

    weight = (hilbert_mult * bp_mask).view(Nt, 1, 1)
    analytic = torch.fft.ifft(Rf_fft * weight, dim=0)
    return analytic.to(torch.complex64)


class DifferentiablePlaneWaveDAS(nn.Module):
    """
    Beamformer Delay-and-Sum (DAS) Diferenciable para Ondas Planas (PWI).
    Corrige los dos bugs del beamformer previo:
      1. Usa la geometría exacta de frente de onda plano en transmisión (sin sumar
         tx_delay[r] espurio al canal receptor r).
      2. Aplica apodización dinámica por recepción (F-number con ventana Hanning) y
         corrección de aberración de fase Eikonal integrada desde el mapa c(x, z).
    """

    def __init__(
        self,
        fs: float,
        c0: float = 1540.0,
        f_number: float = 1.4,
        fc: float = 5.5e6,
    ):
        super().__init__()
        self.fs = float(fs)
        self.c0 = float(c0)
        self.f_number = float(f_number)
        self.fc = float(fc)

    @staticmethod
    def compute_eikonal_delay_perturbation(
        c_grid: torch.Tensor,
        x_grid: torch.Tensor,
        z_grid: torch.Tensor,
        angles: torch.Tensor,
        c0: float,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Integra de forma diferenciable la perturbación de lentitud Δs(x, z) = 1/c(x, z) - 1/c0
        a lo largo de los rayos de transmisión (para cada ángulo θ_m) y de recepción:
            Δτ(x, z; θ_m) = ∫_0^z [1/c(x - (z - z') tan θ_m, z') - 1/c0] dz' / cos θ_m

        Args:
            c_grid: Mapa 2D de velocidad del sonido [Nz, Nx] en m/s.
            x_grid: Vector lateral [Nx] en metros.
            z_grid: Vector axial [Nz] en metros.
            angles: Ángulos de transmisión [Ntx] en radianes.
            c0: Velocidad homogénea de fondo (m/s).
        Returns:
            delta_tau_tx: Perturbación de retardo en transmisión [Ntx, Nz, Nx] (segundos).
            delta_tau_rx_vert: Perturbación de retardo vertical en recepción [Nz, Nx] (segundos).
        """
        Nz, Nx = c_grid.shape
        dz = float(z_grid[1] - z_grid[0]) if Nz > 1 else 1e-4

        # Perturbación de lentitud Δs(z, x) = 1/c(z, x) - 1/c0
        delta_s = (1.0 / c_grid) - (1.0 / c0)  # [Nz, Nx]

        # Integral acumulativa vertical (para recepción y θ = 0)
        delta_tau_rx_vert = torch.cumsum(delta_s, dim=0) * dz  # [Nz, Nx]

        # Para cada ángulo θ_m, integrar a lo largo de la dirección inclinada θ_m
        # Usando integración por capas (marcha de rayos diferenciable)
        cos_theta = torch.cos(angles).view(-1, 1, 1)  # [Ntx, 1, 1]
        # Aproximación de primer orden + corrección angular para rapidez y diferenciabilidad exacta
        delta_tau_tx = delta_tau_rx_vert.unsqueeze(0) / torch.clamp(cos_theta, min=0.85)
        return delta_tau_tx, delta_tau_rx_vert

    def forward(
        self,
        analytic_rf: torch.Tensor,
        x_grid: torch.Tensor,
        z_grid: torch.Tensor,
        rx_x: torch.Tensor,
        angles: torch.Tensor,
        c_grid: torch.Tensor | None = None,
        t0: float = 0.0,
        tx_origin_mode: str = "edge",
        return_coherence: bool = True,
    ) -> dict[str, torch.Tensor]:
        """
        Ejecuta el Beamforming DAS diferenciable sobre la malla (z_grid, x_grid).

        Args:
            analytic_rf: Señal analítica compleja [Nt, Nr, Ntx] en GPU.
            x_grid: Coordenadas laterales [Nx] en metros.
            z_grid: Coordenadas axiales [Nz] en metros.
            rx_x: Posiciones laterales de los elementos del transductor [Nr] en metros.
            angles: Ángulos de onda plana [Ntx] en radianes.
            c_grid: Mapa opcional 2D de velocidad del sonido [Nz, Nx] predicho por la PINN.
            t0: Tiempo inicial de adquisición (segundos).
            tx_origin_mode: 'edge' si el disparo inicia en t=0 en el borde del array (convención
                            min(tx_delay)=0 de CUBDL/Verasonics), o 'center' si inicia en x=0.
            return_coherence: Si True, calcula el mapa de Factor de Coherencia CF(z, x).
        Returns:
            Diccionario con:
              - 'complex_image': Imagen compleja compuesta [Nz, Nx]
              - 'envelope': Envolvente lineal [Nz, Nx]
              - 'bmode_db': Imagen Modo-B en dB normalizada [Nz, Nx]
              - 'coherence_factor': Mapa CF [Nz, Nx] en [0, 1]
              - 'cf_weighted_bmode_db': Imagen Modo-B ponderada por coherencia (suprime lóbulos/ruido)
        """
        device = analytic_rf.device
        Nt, Nr, Ntx = analytic_rf.shape
        Nx = x_grid.shape[0]
        Nz = z_grid.shape[0]

        # Malla 2D en orden [Nz, Nx] (estándar de imagen: filas = z, columnas = x)
        Z, X = torch.meshgrid(z_grid, x_grid, indexing="ij")  # [Nz, Nx]

        # 1. Distancias y tiempos de recepción base: [Nr, Nz, Nx]
        xr = rx_x.view(Nr, 1, 1)
        dx_rx = X.unsqueeze(0) - xr  # [Nr, Nz, Nx]
        dist_rx = torch.sqrt(dx_rx ** 2 + Z.unsqueeze(0) ** 2 + 1e-12)
        t_rx_base = dist_rx / self.c0  # [Nr, Nz, Nx]

        # 2. Apodización dinámica por recepción (F-number con ventana de Hanning suave)
        # semi-apertura activa a profundidad Z: a(Z) = Z / (2 * F#) + apertura mínima (1.5 mm)
        half_aperture = torch.clamp(Z.unsqueeze(0) / (2.0 * self.f_number), min=1.5e-3)
        u_apod = torch.abs(dx_rx) / half_aperture
        apod_mask = (u_apod <= 1.0).to(torch.float32)
        apod_weights = 0.5 * (1.0 + torch.cos(np.pi * torch.clamp(u_apod, 0.0, 1.0))) * apod_mask
        # Normalizar pesos por píxel para mantener ganancia uniforme con la profundidad
        apod_norm = torch.sum(apod_weights, dim=0, keepdim=True) + 1e-6
        apod_weights = apod_weights / apod_norm  # [Nr, Nz, Nx]

        # 3. Perturbaciones de tiempo de vuelo Eikonal si se provee c_grid(z, x)
        if c_grid is not None:
            delta_tau_tx, delta_tau_rx_vert = self.compute_eikonal_delay_perturbation(
                c_grid, x_grid, z_grid, angles, self.c0
            )
            # Oblicuidad de recepción para escalar la integral vertical a cada receptor r
            cos_rx = torch.clamp(Z.unsqueeze(0) / (dist_rx + 1e-9), min=0.5)
            delta_tau_rx = delta_tau_rx_vert.unsqueeze(0) / cos_rx  # [Nr, Nz, Nx]
        else:
            delta_tau_tx = None
            delta_tau_rx = None

        # 4. Bucle sobre ángulos de transmisión (acumulación coherente e incoherente)
        coherent_sum = torch.zeros((Nz, Nx), dtype=torch.complex64, device=device)
        incoherent_energy = torch.zeros((Nz, Nx), dtype=torch.float32, device=device)

        x_min_arr = torch.min(rx_x)
        x_max_arr = torch.max(rx_x)
        r_indices = torch.arange(Nr, device=device).view(Nr, 1, 1).expand(Nr, Nz, Nx)
        omega_c = 2.0 * np.pi * self.fc

        for m in range(Ntx):
            theta = angles[m]
            sin_t = torch.sin(theta)
            cos_t = torch.cos(theta)

            if tx_origin_mode == "edge":
                # El primer elemento en disparar (en t=0) es x_min si theta >= 0, o x_max si theta < 0
                x_orig = torch.where(theta >= 0, x_min_arr, x_max_arr)
                d_tx = (X - x_orig) * sin_t + Z * cos_t
            else:
                d_tx = X * sin_t + Z * cos_t

            t_tx = d_tx / self.c0  # [Nz, Nx]

            if delta_tau_tx is not None and delta_tau_rx is not None:
                dt_aberr = delta_tau_tx[m].unsqueeze(0) + delta_tau_rx  # [Nr, Nz, Nx]
            else:
                dt_aberr = 0.0

            t_total = t_tx.unsqueeze(0) + t_rx_base + dt_aberr  # [Nr, Nz, Nx]

            # Índice fraccionario de muestra temporal
            sample_idx = (t_total - t0) * self.fs
            valid_mask = ((sample_idx >= 1.0) & (sample_idx < float(Nt - 2))).to(torch.float32)
            sample_idx_clipped = torch.clamp(sample_idx, 0.0, float(Nt - 2))

            idx_low = sample_idx_clipped.long()
            idx_high = idx_low + 1
            w_high = (sample_idx_clipped - idx_low.to(torch.float32)).to(torch.complex64)
            w_low = (1.0 - (sample_idx_clipped - idx_low.to(torch.float32))).to(torch.complex64)

            rf_m = analytic_rf[:, :, m]  # [Nt, Nr]
            val_low = rf_m[idx_low, r_indices]
            val_high = rf_m[idx_high, r_indices]
            interp_val = (val_low * w_low + val_high * w_high) * valid_mask.to(torch.complex64)

            # Si c_grid tiene gradiente activo, aplicar rotación de fase explícita suave
            # para garantizar flujo de gradiente limpio hacia c_grid sin bloqueo por cuantización de índices
            if c_grid is not None and c_grid.requires_grad:
                phase_corr = torch.exp(-1j * omega_c * dt_aberr)
                interp_val = interp_val * phase_corr

            # Suma ponderada por apodización dinámica sobre los receptores
            weighted_rx = interp_val * apod_weights.to(torch.complex64)
            img_angle = torch.sum(weighted_rx, dim=0)  # [Nz, Nx]

            coherent_sum = coherent_sum + img_angle
            if return_coherence:
                incoherent_energy = incoherent_energy + torch.abs(img_angle) ** 2

        # 5. Envolvente, factor de coherencia inter-angular y compresión logarítmica
        envelope = torch.abs(coherent_sum)
        env_max = torch.max(envelope) + 1e-12
        env_norm = envelope / env_max
        bmode_db = 20.0 * torch.log10(env_norm + 1e-4)  # Rango dinámico [-80, 0] dB

        if return_coherence:
            # Factor de Coherencia Multi-Ángulo: mide qué tan bien alineadas en fase
            # están las reconstrucciones de los Ntx ángulos de onda plana
            cf = (envelope ** 2) / (float(Ntx) * incoherent_energy + 1e-12)
            cf = torch.clamp(cf, 0.0, 1.0)
            # Ponderación parcial por coherencia (γ = 0.5) para realzar bordes y contraste de lesión
            env_cf = envelope * torch.sqrt(cf + 1e-6)
            env_cf_norm = env_cf / (torch.max(env_cf) + 1e-12)
            cf_bmode_db = 20.0 * torch.log10(env_cf_norm + 1e-4)
        else:
            cf = torch.ones_like(envelope)
            cf_bmode_db = bmode_db

        return {
            "complex_image": coherent_sum,
            "envelope": envelope,
            "bmode_db": bmode_db,
            "coherence_factor": cf,
            "cf_weighted_bmode_db": cf_bmode_db,
        }
