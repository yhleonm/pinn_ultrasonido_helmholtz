import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def solve_helmholtz_1d_fd(z, frequencies, c_profile, alpha0_profile):
    """
    Solves the 1D Helmholtz equation exactly on a grid z using finite differences.
    d^2 P/dz^2 + k^2 P = 0
    BCs:
      - Dirichlet at z=0: P(0) = 1.0
      - Sommerfeld radiation at z=L: dP/dz - i * k(L) * P(L) = 0
    """
    N = len(z)
    dz = z[1] - z[0]
    P_dict = {}
    
    for f in frequencies:
        omega = 2.0 * np.pi * f
        k = omega / c_profile + 1j * (alpha0_profile * f)
        
        A = np.zeros((N, N), dtype=complex)
        B = np.zeros(N, dtype=complex)
        
        A[0, 0] = 1.0
        B[0] = 1.0 + 0j
        
        for j in range(1, N - 1):
            A[j, j - 1] = 1.0
            A[j, j] = (k[j]**2) * (dz**2) - 2.0
            A[j, j + 1] = 1.0
            B[j] = 0.0
            
        A[N - 1, N - 2] = -1.0
        A[N - 1, N - 1] = 1.0 - 1j * k[N - 1] * dz
        B[N - 1] = 0.0
        
        P = np.linalg.solve(A, B)
        P_dict[f] = P
        
    return P_dict

class SyntheticHelmholtz1D:
    """
    Generator for synthetic 1D ultrasound propagation data.
    """
    def __init__(self, z_max=0.01, dz=1e-5, frequencies=None):
        self.z_max = z_max
        self.dz = dz
        self.z = np.arange(0, z_max + dz, dz)
        self.frequencies = frequencies if frequencies is not None else [1e6, 1.5e6, 2e6]
        
        self.c_gt = np.ones_like(self.z) * 1500.0
        self.alpha0_gt = np.ones_like(self.z) * (0.5 / 8.686e4)
        
        idx_inclusion = (self.z >= 0.004) & (self.z <= 0.007)
        self.c_gt[idx_inclusion] = 1580.0
        self.alpha0_gt[idx_inclusion] = 1.5 / 8.686e4
        
        self.P_gt = solve_helmholtz_1d_fd(self.z, self.frequencies, self.c_gt, self.alpha0_gt)
        
    def get_data(self):
        return {
            'z': self.z,
            'frequencies': self.frequencies,
            'c_gt': self.c_gt,
            'alpha0_gt': self.alpha0_gt,
            'P_gt': self.P_gt
        }

def solve_helmholtz_2d_fd(x, z, frequencies, c_map, alpha0_map):
    """
    Solves the 2D Helmholtz equation exactly on a grid (x, z) using finite differences.
    d^2 P/dx^2 + d^2 P/dz^2 + k^2 P = 0
    BCs:
      - Dirichlet at z=0: P(x, 0) = HannWindow(x) (Transducer aperture)
      - Sommerfeld radiation on left, right, and bottom borders:
        dP/dn - i * k * P = 0
    """
    Nx = len(x)
    Nz = len(z)
    h = x[1] - x[0]  # Assumes dx == dz
    P_dict = {}
    
    for f in frequencies:
        omega = 2.0 * np.pi * f
        # Wavenumber k = omega/c + i * alpha0 * f
        k = omega / c_map + 1j * (alpha0_map * f)
        
        row_indices = []
        col_indices = []
        data_values = []
        B = np.zeros(Nz * Nx, dtype=complex)
        
        for j in range(Nz):
            for i in range(Nx):
                p = j * Nx + i
                
                # 1. Boundary j=0 (Dirichlet transducer source)
                if j == 0:
                    row_indices.append(p)
                    col_indices.append(p)
                    data_values.append(1.0)
                    
                    # Transducer aperture: active width from -1.5 cm to +1.5 cm
                    if -0.015 <= x[i] <= 0.015:
                        B[p] = np.sin(np.pi * (x[i] + 0.015) / 0.030)**2 + 0.0j
                    else:
                        B[p] = 0.0j
                
                # 2. Boundary i=0 (Left border Sommerfeld BC)
                elif i == 0:
                    row_indices.append(p)
                    col_indices.append(p)
                    data_values.append(1.0 + 1j * k[j, i] * h)
                    
                    row_indices.append(p)
                    col_indices.append(p + 1)
                    data_values.append(-1.0)
                    B[p] = 0.0j
                
                # 3. Boundary i=Nx-1 (Right border Sommerfeld BC)
                elif i == Nx - 1:
                    row_indices.append(p)
                    col_indices.append(p - 1)
                    data_values.append(-1.0)
                    
                    row_indices.append(p)
                    col_indices.append(p)
                    data_values.append(1.0 - 1j * k[j, i] * h)
                    B[p] = 0.0j
                
                # 4. Boundary j=Nz-1 (Bottom border Sommerfeld BC)
                elif j == Nz - 1:
                    row_indices.append(p)
                    col_indices.append(p - Nx)
                    data_values.append(-1.0)
                    
                    row_indices.append(p)
                    col_indices.append(p)
                    data_values.append(1.0 - 1j * k[j, i] * h)
                    B[p] = 0.0j
                
                # 5. Interior nodes (5-point Laplcian stencil)
                else:
                    row_indices.append(p)
                    col_indices.append(p - 1)
                    data_values.append(1.0)
                    
                    row_indices.append(p)
                    col_indices.append(p + 1)
                    data_values.append(1.0)
                    
                    row_indices.append(p)
                    col_indices.append(p - Nx)
                    data_values.append(1.0)
                    
                    row_indices.append(p)
                    col_indices.append(p + Nx)
                    data_values.append(1.0)
                    
                    row_indices.append(p)
                    col_indices.append(p)
                    data_values.append(k[j, i]**2 * h**2 - 4.0)
                    B[p] = 0.0j
                    
        # Solve sparse system
        A = sp.coo_matrix((data_values, (row_indices, col_indices)), shape=(Nz*Nx, Nz*Nx)).tocsr()
        P_vec = spla.spsolve(A, B)
        P_dict[f] = P_vec.reshape((Nz, Nx))
        
    return P_dict

class SyntheticHelmholtz2D:
    """
    Generator for synthetic 2D ultrasound propagation data.
    Allows generic phantom definition.
    """
    def __init__(self, x_max=0.02, z_max=0.04, dx=1e-4, dz=1e-4, frequencies=None,
                 c_bg=1500.0, c_in=1580.0, alpha0_bg_db=0.5, alpha0_in_db=1.5,
                 inclusion_radius=0.006, inclusion_z=0.02):
        self.x_max = x_max
        self.z_max = z_max
        self.dx = dx
        self.dz = dz
        
        self.x = np.arange(-x_max, x_max + dx, dx)
        self.z = np.arange(0, z_max + dz, dz)
        self.frequencies = frequencies if frequencies is not None else [1e6, 1.5e6, 2e6]
        
        Nx, Nz = len(self.x), len(self.z)
        
        self.c_gt = np.ones((Nz, Nx)) * c_bg
        self.alpha0_gt = np.ones((Nz, Nx)) * (alpha0_bg_db / 8.686e4)
        
        # Inclusion: a circle of radius inclusion_radius centered at (x=0, z=inclusion_z)
        X, Z = np.meshgrid(self.x, self.z)
        dist_sq = (X - 0.0)**2 + (Z - inclusion_z)**2
        idx_inclusion = dist_sq <= inclusion_radius**2
        
        self.c_gt[idx_inclusion] = c_in
        self.alpha0_gt[idx_inclusion] = alpha0_in_db / 8.686e4
        
        self.P_gt = solve_helmholtz_2d_fd(self.x, self.z, self.frequencies, self.c_gt, self.alpha0_gt)
        
    def get_data(self):
        return {
            'x': self.x,
            'z': self.z,
            'frequencies': self.frequencies,
            'c_gt': self.c_gt,
            'alpha0_gt': self.alpha0_gt,
            'P_gt': self.P_gt
        }

class CUBDLLoader:
    """
    Loader for raw channel RF data in CUBDL format (Johns Hopkins University challenge).
    Supports Scipy .mat loading and H5PY loading (MATLAB v7.3).
    """
    def __init__(self, file_path):
        self.file_path = file_path
        self.fs = None
        self.c0 = 1540.0
        self.rf = None
        self.rx_positions = None
        self.tx_delays = None
        self.angles = None
        self.t0 = 0.0
        self._load_data()
        
    def _load_data(self):
        try:
            # Intentar carga estándar con scipy (MATLAB v7)
            import scipy.io as sio
            mat_data = sio.loadmat(self.file_path)
            
            self.fs = float(mat_data.get('fs', [40e6])[0])
            self.c0 = float(mat_data.get('c', [1540.0])[0])
            self.t0 = float(mat_data.get('t0', [0.0])[0])
            self.rf = mat_data.get('rf')  # shape: [samples, channels, transmits]
            
            if 'angles' in mat_data:
                self.angles = np.array(mat_data['angles']).flatten()
            elif 'transmit_direction' in mat_data:
                self.angles = np.array(mat_data['transmit_direction'])[0, :].flatten()
                
            if 'rx_pos' in mat_data:
                self.rx_positions = mat_data['rx_pos']
            elif 'rx_position' in mat_data:
                self.rx_positions = mat_data['rx_position']
            else:
                n_chan = self.rf.shape[1] if self.rf is not None else 128
                pitch = 0.3e-3
                x_pos = (np.arange(n_chan) - (n_chan - 1) / 2.0) * pitch
                self.rx_positions = np.stack([x_pos, np.zeros_like(x_pos)], axis=1)
                
            self.tx_delays = mat_data.get('tx_delay', np.zeros((1, self.rx_positions.shape[0])))
            
        except Exception as e:
            # Fallback a h5py (MATLAB v7.3 HDF5 o HDF5 genérico)
            print(f"Scipy loadmat falló: {e}. Intentando carga v7.3/HDF5 con h5py...")
            import h5py
            with h5py.File(self.file_path, 'r') as f:
                # 1. sampling frequency
                self.fs = None
                for k in ['fs', 'sampling_frequency']:
                    if k in f:
                        val = f[k][()]
                        self.fs = float(val.flatten()[0]) if hasattr(val, 'flatten') else float(val)
                        break
                if self.fs is None:
                    self.fs = 40e6
                    
                # 2. sound speed
                self.c0 = None
                for k in ['c', 'sound_speed']:
                    if k in f:
                        val = f[k][()]
                        self.c0 = float(val.flatten()[0]) if hasattr(val, 'flatten') else float(val)
                        break
                if self.c0 is None:
                    self.c0 = 1540.0
                    
                # 3. start_time (t0)
                self.t0 = 0.0
                for k in ['t0', 'start_time']:
                    if k in f:
                        val = f[k][()]
                        self.t0 = float(val.flatten()[0]) if hasattr(val, 'flatten') else float(val)
                        break
                        
                # 4. angles/transmit direction
                angles = None
                if 'angles' in f:
                    angles = np.array(f['angles']).flatten()
                elif 'transmit_direction' in f:
                    angles = np.array(f['transmit_direction'])[0, :].flatten()
                    
                # 5. rf channel data
                if 'rf' in f:
                    rf_ref = f['rf']
                    self.rf = np.array(rf_ref).T
                elif 'channel_data' in f:
                    raw_data = np.array(f['channel_data'])
                    if len(raw_data.shape) == 3:
                        # shape: [transmits, channels, samples] -> transpose to [samples, channels, transmits]
                        self.rf = np.transpose(raw_data, (2, 1, 0))
                    elif len(raw_data.shape) == 2:
                        # shape: [channels, transmits * samples] -> reshape & transpose
                        num_channels = raw_data.shape[0]
                        if angles is not None:
                            num_transmits = len(angles)
                            num_samples = raw_data.shape[1] // num_transmits
                            reshaped = raw_data.reshape((num_channels, num_transmits, num_samples))
                            self.rf = np.transpose(reshaped, (2, 0, 1))
                        else:
                            self.rf = raw_data
                
                # 6. receiver positions
                if 'rx_pos' in f:
                    self.rx_positions = np.array(f['rx_pos']).T
                elif 'rx_position' in f:
                    self.rx_positions = np.array(f['rx_position']).T
                elif 'element_positions' in f:
                    self.rx_positions = np.array(f['element_positions'])[:2, :].T
                else:
                    n_chan = self.rf.shape[1] if self.rf is not None else 128
                    pitch = 0.3e-3
                    x_pos = (np.arange(n_chan) - (n_chan - 1) / 2.0) * pitch
                    self.rx_positions = np.stack([x_pos, np.zeros_like(x_pos)], axis=1)
                    
                # 7. transmit delays
                if 'tx_delay' in f:
                    self.tx_delays = np.array(f['tx_delay']).T
                elif 'tx_delays' in f:
                    self.tx_delays = np.array(f['tx_delays']).T
                elif angles is not None:
                    n_trans = len(angles)
                    n_chan = self.rx_positions.shape[0]
                    self.tx_delays = np.zeros((n_trans, n_chan))
                    for t_idx, theta in enumerate(angles):
                        delays = self.rx_positions[:, 0] * np.sin(theta) / self.c0
                        delays = delays - np.min(delays)
                        self.tx_delays[t_idx, :] = delays
                else:
                    n_chan = self.rx_positions.shape[0]
                    self.tx_delays = np.zeros((1, n_chan))
                    
                self.angles = angles
                
    def extract_boundary_pressure(self, selected_frequencies, transmit_index=0):
        """
        Extrae la presión compleja de contorno en los sensores aplicando
        Transformada de Hilbert para obtener la fase analítica y compensando
        los retardos de disparo (tx delays).
        
        Returns:
            sensor_x: Array con posiciones x de los sensores [N_channels]
            P_meas: Dict mapeando frecuencia (Hz) -> presión compleja [N_channels]
        """
        from scipy.signal import hilbert
        
        if len(self.rf.shape) == 3:
            rf_tx = self.rf[:, :, transmit_index]
        else:
            rf_tx = self.rf
            
        N_samples, N_channels = rf_tx.shape
        
        # 1. Señal analítica por transformada de Hilbert
        rf_complex = hilbert(rf_tx, axis=0)
        
        # 2. FFT temporal
        rf_fft = np.fft.fft(rf_complex, axis=0)
        freqs = np.fft.fftfreq(N_samples, d=1.0/self.fs)
        
        P_meas = {}
        sensor_x = self.rx_positions[:, 0]
        
        # Obtener los delays de disparo para esta transmisión
        if len(self.tx_delays.shape) == 2:
            tx_delay = self.tx_delays[transmit_index, :]
        else:
            tx_delay = self.tx_delays
            
        for f in selected_frequencies:
            idx_freq = np.argmin(np.abs(freqs - f))
            P_raw = rf_fft[idx_freq, :]
            
            # 3. Compensación de retardo: e^(i * omega * (tau_tx - t0))
            omega = 2.0 * np.pi * f
            P_calib = P_raw * np.exp(1j * omega * (tx_delay - self.t0))
            
            P_meas[f] = P_calib
            
        return sensor_x, P_meas


class OASBUDLoader:
    """
    Cargador clínico para el dataset OASBUD (Open Access Series of Breast Ultrasound Data,
    Zenodo 545928 - Piotrzkowska-Wróblewska et al., Medical Physics 2017).

    Gestiona señales RF 2D post-beamformed (líneas A-scan focalizadas de un transductor
    lineal Ultrasonix L14-5/38 a fs = 40 MHz) y las máscaras binarias (roi1/roi2)
    delineadas por radiólogos expertos para 100 lesiones mamarias (48 benignas, 52 malignas).
    """

    _MAT_CACHE: dict[str, np.ndarray] = {}

    def __init__(
        self,
        mat_path: str = "data_raw/OASBUD.mat",
        case_id: str | int = "1mm",
        view: int = 1,
    ):
        self.mat_path = str(mat_path)
        self.fs = 40.0e6          # Frecuencia de muestreo Ultrasonix SonixTouch (40 MHz)
        self.c0 = 1540.0          # Velocidad de referencia en tejido mamario (m/s)
        self.aperture_m = 0.038   # Ancho lateral activo del transductor L14-5/38 (38 mm)
        self.t0 = 0.0

        self.case_id: str = ""
        self.view: int = int(view)
        self.birads: str = ""
        self.label: int = 0       # 0 = Benigno, 1 = Maligno
        self.clinical_label: str = ""
        self.rf_2d: np.ndarray | None = None
        self.roi_mask: np.ndarray | None = None
        self.x_grid: np.ndarray | None = None
        self.z_grid: np.ndarray | None = None
        self.dx: float = 0.0
        self.dz: float = 0.0

        self._load_case(case_id=case_id, view=view)

    def _load_case(self, case_id: str | int, view: int) -> None:
        import scipy.io as sio

        if self.mat_path not in OASBUDLoader._MAT_CACHE:
            mat = sio.loadmat(self.mat_path, squeeze_me=True)
            OASBUDLoader._MAT_CACHE[self.mat_path] = mat["data"]

        data = OASBUDLoader._MAT_CACHE[self.mat_path]
        record = None
        if isinstance(case_id, int):
            record = data[case_id]
        else:
            target = str(case_id).replace("OASBUD_", "").strip()
            for item in data:
                if str(item["id"]).strip() == target:
                    record = item
                    break
        if record is None:
            raise ValueError(f"Caso '{case_id}' no encontrado en {self.mat_path}.")

        rf_key = "rf1" if view == 1 else "rf2"
        roi_key = "roi1" if view == 1 else "roi2"

        self.case_id = str(record["id"]).strip()
        self.birads = str(record["birads"]).strip()
        self.label = int(record["class"])
        pathology_str = "MALIGNO" if self.label == 1 else "BENIGNO"
        self.clinical_label = f"{pathology_str} - BI-RADS {self.birads} (OASBUD)"

        self.rf_2d = np.asarray(record[rf_key], dtype=np.float32)
        self.roi_mask = np.asarray(record[roi_key] > 0, dtype=bool)

        Nz, Nx = self.rf_2d.shape
        self.dz = float(self.c0 / (2.0 * self.fs))  # 19.25 um
        self.dx = float(self.aperture_m / Nx)       # ~74.51 um para 510 líneas
        self.x_grid = ((np.arange(Nx, dtype=np.float32) - (Nx - 1) / 2.0) * self.dx)
        self.z_grid = (np.arange(Nz, dtype=np.float32) * self.dz)

    def get_roi_centroid_and_bounds(self, margin_factor: float = 0.82) -> dict[str, float]:
        """
        Obtiene el centro geométrico y la ventana de búsqueda alrededor de la lesión
        para centrar el dominio computacional de la PINN (sin filtrar la forma de roi_mask
        al entrenamiento).
        """
        z_idx, x_idx = np.where(self.roi_mask)
        x_center = float(self.x_grid[x_idx].mean())
        z_center = float(self.z_grid[z_idx].mean())
        span_x = float(self.x_grid[x_idx.max()] - self.x_grid[x_idx.min()])
        span_z = float(self.z_grid[z_idx.max()] - self.z_grid[z_idx.min()])
        half_x = max(0.0045, margin_factor * span_x)
        half_z = max(0.0045, margin_factor * span_z)
        half_iso = max(half_x, half_z)
        # Radio equivalente área-efectiva solo para referencia escalar
        area_mm2 = float(self.roi_mask.sum()) * (self.dx * 1e3) * (self.dz * 1e3)
        equiv_radius_mm = float(np.sqrt(area_mm2 / np.pi))
        return {
            "x_center_m": x_center,
            "z_center_m": z_center,
            "half_x_m": half_x,
            "half_z_m": half_z,
            "roi_half_mm": round(float(half_iso * 1e3), 2),
            "equiv_radius_mm": round(equiv_radius_mm, 2),
        }

    def extract_roi_contour_mm(self, max_points: int = 140) -> list[list[float]]:
        """
        Extrae el polígono cerrado real [x_mm, z_mm] del contorno de la máscara
        delineada por el radiólogo (roi_mask) ordenando los píxeles de borde angularmente
        o mediante matplotlib contour, para renderizado en gráficas y Medical Insight 3D.
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        X_mm, Z_mm = np.meshgrid(self.x_grid * 1e3, self.z_grid * 1e3)
        cs = ax.contour(X_mm, Z_mm, self.roi_mask.astype(np.float32), levels=[0.5])
        segs = cs.allsegs[0] if len(cs.allsegs) > 0 else []
        plt.close(fig)

        if not segs:
            return []
        longest = max(segs, key=lambda s: len(s))
        if len(longest) > max_points:
            idx = np.linspace(0, len(longest) - 1, max_points, dtype=int)
            longest = longest[idx]
        return [[round(float(pt[0]), 3), round(float(pt[1]), 3)] for pt in longest]

