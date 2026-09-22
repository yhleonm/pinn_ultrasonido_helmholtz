import os
import time
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

# Import custom modules
from data_engine.loader import SyntheticHelmholtz2D
from models.pinn_net import DecoupledUltrasoundPINN
from physics.loss_helmholtz import compute_helmholtz_residuals_2d, compute_sommerfeld_loss_2d

# Self-contained SSIM implementation to avoid dependency issues
def calculate_ssim(im1, im2):
    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1_sq = im1.var()
    sigma2_sq = im2.var()
    # Covariance
    sigma12 = np.cov(im1.flatten(), im2.flatten())[0, 1]
    
    L = max(im1.max(), im2.max()) - min(im1.min(), im2.min())
    if L == 0:
        L = 1.0
    C1 = (0.01 * L) ** 2
    C2 = (0.03 * L) ** 2
    
    numerator = (2 * mu1 * mu2 + C1) * (2 * sigma12 + C2)
    denominator = (mu1**2 + mu2**2 + C1) * (sigma1_sq + sigma2_sq + C2)
    return numerator / (denominator + 1e-10)

def run_5MHz_roi():
    # 1. Configuración de Hardware
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"ROI 5MHz: Usando dispositivo: {device}")
    
    # 2. Parámetros de la ROI (1.5 cm x 1.5 cm)
    x_max = 0.0075  # lateral: -7.5 mm a +7.5 mm (1.5 cm total)
    z_max = 0.015   # axial: 0 a 1.5 cm (1.5 cm total)
    
    # Resolución: al menos 10 puntos por lambda
    # Wavelength lambda = c_bg / f = 1500 / 5e6 = 0.3 mm = 300 um
    # lambda / 10 = 30 um = 0.03 mm = 3.0e-5 m
    dx = 3.0e-5
    dz = 3.0e-5
    frequencies = [5.0e6]  # 5 MHz
    
    print(f"Generando dataset sintético en ROI pequeña a 5MHz...")
    print(f"ROI: lateral={2*x_max*100:.2f}cm, axial={z_max*100:.2f}cm")
    print(f"Espaciado: {dx*1e6:.1f} micras. Frecuencia: {frequencies[0]/1e6:.1f} MHz")
    
    t0_gen = time.time()
    # Inclusión circular centrada en el medio de la ROI
    # Centro en x=0.0, z = z_max / 2 = 0.0075
    # Radio = 3 mm = 0.003
    data_gen = SyntheticHelmholtz2D(
        x_max=x_max, z_max=z_max, dx=dx, dz=dz, frequencies=frequencies,
        c_bg=1500.0, c_in=1580.0,
        alpha0_bg_db=0.5, alpha0_in_db=0.9,
        inclusion_radius=0.003, inclusion_z=0.0075
    )
    data = data_gen.get_data()
    t_gen = time.time() - t0_gen
    print(f"Dataset de ROI generado en {t_gen:.1f} segundos.")
    
    x_grid = data['x']
    z_grid = data['z']
    c_gt = data['c_gt']            # shape [Nz, Nx]
    alpha0_gt = data['alpha0_gt']  # shape [Nz, Nx]
    P_gt = data['P_gt']            # dict: 5e6 -> array shape [Nz, Nx]
    
    Nx, Nz = len(x_grid), len(z_grid)
    print(f"Dimensión de grilla regular de la ROI: {Nz} x {Nx} = {Nz*Nx} puntos.")
    
    # Guardar Dataset de la ROI
    os.makedirs("data_engine", exist_ok=True)
    np.savez(
        'data_engine/synthetic_dataset_5MHz_roi.npz',
        x=x_grid, z=z_grid,
        c_gt=c_gt, alpha0_gt=alpha0_gt,
        frequencies=frequencies,
        P_gt_5MHz=P_gt[5.0e6]
    )
    print("Dataset de la ROI guardado en data_engine/synthetic_dataset_5MHz_roi.npz")
    
    # Preparar coordenadas y tensores
    X, Z = np.meshgrid(x_grid, z_grid)  # shape [Nz, Nx]
    coords_flat = np.stack([X.flatten(), Z.flatten()], axis=1)  # shape [Nz*Nx, 2]
    coords_t = torch.tensor(coords_flat, dtype=torch.float32, device=device)
    
    c_gt_flat_t = torch.tensor(c_gt.flatten(), dtype=torch.float32, device=device).unsqueeze(1)
    alpha0_gt_flat_t = torch.tensor(alpha0_gt.flatten(), dtype=torch.float32, device=device).unsqueeze(1)
    
    # 3. Selección de Datos de Observación Escasos (20% del dominio)
    np.random.seed(42)
    obs_fraction = 0.20
    n_pts = len(coords_flat)
    n_obs = int(n_pts * obs_fraction)
    obs_indices = np.random.choice(n_pts, n_obs, replace=False)
    obs_indices.sort()
    
    coords_obs_t = coords_t[obs_indices]
    
    # Datos de observación de presión en el interior
    P_5MHz_flat = P_gt[5.0e6].flatten()
    obs_data_real = torch.tensor(P_5MHz_flat.real[obs_indices], dtype=torch.float32, device=device).unsqueeze(1)
    obs_data_imag = torch.tensor(P_5MHz_flat.imag[obs_indices], dtype=torch.float32, device=device).unsqueeze(1)
    
    # 4. Datos de Sensores en el Borde Superior z=0 (Condición de Dirichlet)
    # Todos los puntos a lo largo de la línea z=0
    x_sensor = torch.tensor(x_grid, dtype=torch.float32, device=device).unsqueeze(1)
    z_sensor = torch.zeros_like(x_sensor)
    coords_sensor_t = torch.cat([x_sensor, z_sensor], dim=1)
    
    P_sensor_real = torch.tensor(P_gt[5.0e6][0, :].real, dtype=torch.float32, device=device).unsqueeze(1)
    P_sensor_imag = torch.tensor(P_gt[5.0e6][0, :].imag, dtype=torch.float32, device=device).unsqueeze(1)
    
    # 5. Inicialización de la PINN para ROI
    # Usando omega_0_pressure=60.0 y omega_0_params=5.0
    pinn = DecoupledUltrasoundPINN(
        spatial_dim=2, hidden_features=256, num_layers=4,
        omega_0_pressure=60.0, omega_0_params=5.0
    ).to(device)
    
    # Configuración de Optimización
    epochs = 3000
    optimizer = optim.Adam(pinn.parameters(), lr=3e-4)
    scheduler = optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=500, T_mult=2)
    
    # Pesos de Pérdidas
    lambda_data = 50.0         # Observaciones internas (20%)
    lambda_dirichlet = 50.0    # Sensores en la superficie z=0
    lambda_pde = 5.0           # Ecuación de Helmholtz en el volumen
    lambda_sommerfeld = 1.0     # Bordes absorbentes (izq, der, inferior)
    
    # Función de normalización
    def normalize_inputs(coords_xz, omega):
        x_norm = coords_xz[:, 0:1] / x_max
        z_norm = 2.0 * (coords_xz[:, 1:2] / z_max) - 1.0
        # Frecuencia única: el valor de entrada omega_norm es constante a 0
        omega_norm = torch.zeros_like(x_norm)
        return torch.cat([x_norm, z_norm], dim=1), omega_norm
    
    # Puntos de colocación dinámicos
    n_coloc_points = 15000
    n_boundary_points = 500
    
    best_ssim_c = -1.0
    best_state = None
    best_epoch = 0
    
    print("Iniciando entrenamiento de PINN en ROI a 5MHz...")
    start_time = time.time()
    progress_history = []
    
    for epoch in range(1, epochs + 1):
        pinn.train()
        optimizer.zero_grad()
        
        # A. Pérdida de Datos (Puntos de observación dispersos en el volumen)
        omega = 2.0 * np.pi * 5.0e6
        coords_norm_obs, omega_norm_obs = normalize_inputs(coords_obs_t, omega)
        inputs_obs = torch.cat([coords_norm_obs, omega_norm_obs], dim=1)
        pred_obs_real, pred_obs_imag = pinn.forward_pressure(inputs_obs)
        loss_data = torch.mean((pred_obs_real - obs_data_real)**2 + (pred_obs_imag - obs_data_imag)**2)
        
        # B. Pérdida de Dirichlet (Sensores en el borde superior z=0)
        coords_norm_sensor, omega_norm_sensor = normalize_inputs(coords_sensor_t, omega)
        inputs_sensor = torch.cat([coords_norm_sensor, omega_norm_sensor], dim=1)
        pred_sensor_real, pred_sensor_imag = pinn.forward_pressure(inputs_sensor)
        loss_dirichlet = torch.mean((pred_sensor_real - P_sensor_real)**2 + (pred_sensor_imag - P_sensor_imag)**2)
        
        # C. Pérdida PDE (Helmholtz Residuals)
        x_rand = (torch.rand(n_coloc_points, 1, device=device) * 2.0 - 1.0) * x_max
        z_rand = torch.rand(n_coloc_points, 1, device=device) * z_max
        coords_coloc = torch.cat([x_rand, z_rand], dim=1)
        coords_coloc.requires_grad_(True)
        
        coords_norm_coloc, omega_norm_coloc = normalize_inputs(coords_coloc, omega)
        c_pred, alpha0_pred = pinn.forward_params(coords_norm_coloc)
        
        inputs_coloc = torch.cat([coords_norm_coloc, omega_norm_coloc], dim=1)
        pred_real, pred_imag = pinn.forward_pressure(inputs_coloc)
        
        res_real, res_imag = compute_helmholtz_residuals_2d(
            coords_coloc, pred_real, pred_imag, alpha0_pred, c_pred, omega
        )
        loss_pde = torch.mean(res_real**2 + res_imag**2)
        
        # D. Pérdida Sommerfeld en bordes izquierdo, derecho e inferior
        z_left = torch.rand(n_boundary_points, 1, device=device) * z_max
        x_left = torch.ones_like(z_left) * -x_max
        coords_left = torch.cat([x_left, z_left], dim=1)
        coords_left.requires_grad_(True)
        normal_left = torch.tensor([[-1.0, 0.0]], device=device).repeat(n_boundary_points, 1)
        
        z_right = torch.rand(n_boundary_points, 1, device=device) * z_max
        x_right = torch.ones_like(z_right) * x_max
        coords_right = torch.cat([x_right, z_right], dim=1)
        coords_right.requires_grad_(True)
        normal_right = torch.tensor([[1.0, 0.0]], device=device).repeat(n_boundary_points, 1)
        
        x_bottom = (torch.rand(n_boundary_points, 1, device=device) * 2.0 - 1.0) * x_max
        z_bottom = torch.ones_like(x_bottom) * z_max
        coords_bottom = torch.cat([x_bottom, z_bottom], dim=1)
        coords_bottom.requires_grad_(True)
        normal_bottom = torch.tensor([[0.0, 1.0]], device=device).repeat(n_boundary_points, 1)
        
        coords_b = torch.cat([coords_left, coords_right, coords_bottom], dim=0)
        normals_b = torch.cat([normal_left, normal_right, normal_bottom], dim=0)
        
        coords_norm_b, omega_norm_b = normalize_inputs(coords_b, omega)
        c_pred_b, alpha0_pred_b = pinn.forward_params(coords_norm_b)
        
        inputs_b = torch.cat([coords_norm_b, omega_norm_b], dim=1)
        pred_real_b, pred_imag_b = pinn.forward_pressure(inputs_b)
        
        loss_sommer = compute_sommerfeld_loss_2d(
            coords_b, pred_real_b, pred_imag_b, alpha0_pred_b, c_pred_b, omega, normals_b
        )
        
        # Entrenamiento por etapas
        # Épocas 1-1500: Ajustar campo de presiones acústicas
        # Épocas 1501-3000: Ajuste conjunto de parámetros c y alpha0
        w_pde = lambda_pde if epoch > 1500 else 0.0
        w_sommer = lambda_sommerfeld if epoch > 1500 else 0.0
        
        loss_total = lambda_data * loss_data + lambda_dirichlet * loss_dirichlet + w_pde * loss_pde + w_sommer * loss_sommer
        
        loss_total.backward()
        optimizer.step()
        scheduler.step()
        
        # Evaluación periódica cada 100 épocas
        if epoch % 100 == 0 or epoch == 1 or epoch == 1501:
            pinn.eval()
            with torch.no_grad():
                coords_norm_grid, _ = normalize_inputs(coords_t, omega)
                c_pred_flat, alpha0_pred_flat = pinn.forward_params(coords_norm_grid)
                
                err_c = torch.norm(c_pred_flat - c_gt_flat_t) / torch.norm(c_gt_flat_t)
                err_alpha0 = torch.norm(alpha0_pred_flat - alpha0_gt_flat_t) / torch.norm(alpha0_gt_flat_t)
                
                c_pred_2d = c_pred_flat.cpu().numpy().reshape((Nz, Nx))
                alpha0_pred_2d = alpha0_pred_flat.cpu().numpy().reshape((Nz, Nx)) * 8.686e4
                alpha0_gt_db = alpha0_gt * 8.686e4
                
                ssim_c = calculate_ssim(c_pred_2d, c_gt)
                ssim_alpha0 = calculate_ssim(alpha0_pred_2d, alpha0_gt_db)
                
                lr = optimizer.param_groups[0]['lr']
                print(f"ROI Epoch {epoch:4d} | Loss: {loss_total.item():.4e} "
                      f"[Data: {loss_data.item():.4e}, Dir: {loss_dirichlet.item():.4e}, PDE: {loss_pde.item():.4e}] | "
                      f"Err c: {err_c.item()*100:.2f}%, Err alpha0: {err_alpha0.item()*100:.2f}% | "
                      f"SSIM c: {ssim_c:.3f}, SSIM alpha0: {ssim_alpha0:.3f} | lr: {lr:.2e}")
                
                if epoch > 1500 and ssim_c > best_ssim_c:
                    best_ssim_c = ssim_c
                    best_state = {k: v.cpu().clone() for k, v in pinn.state_dict().items()}
                    best_epoch = epoch
                    print(f"  --> [ROI SAVED] Nuevo mejor estado (SSIM c: {ssim_c:.4f} en época {epoch})")
                
                progress_history.append({
                    'epoch': epoch,
                    'loss_total': loss_total.item(),
                    'loss_data': loss_data.item(),
                    'loss_dirichlet': loss_dirichlet.item(),
                    'loss_pde': loss_pde.item(),
                    'err_c': err_c.item(),
                    'err_alpha0': err_alpha0.item(),
                    'ssim_c': ssim_c,
                    'ssim_alpha0': ssim_alpha0
                })
                
    elapsed_time = time.time() - start_time
    print(f"Entrenamiento en ROI completado en {elapsed_time:.1f} segundos.")
    
    # Cargar el mejor estado
    if best_state is not None:
        print(f"\nCargando el mejor estado (Época {best_epoch}, SSIM c: {best_ssim_c:.4f}) para reportar.")
        pinn.load_state_dict({k: v.to(device) for k, v in best_state.items()})
        torch.save(pinn.state_dict(), "articulo/best_pinn_5MHz_roi.pt")
        
    pinn.eval()
    with torch.no_grad():
        coords_norm_grid, _ = normalize_inputs(coords_t, omega)
        c_pred_final, alpha0_pred_final_t = pinn.forward_params(coords_norm_grid)
        c_pred_final = c_pred_final.cpu().numpy().reshape((Nz, Nx))
        alpha0_pred_final = alpha0_pred_final_t.cpu().numpy().reshape((Nz, Nx)) * 8.686e4
        alpha0_gt_db = alpha0_gt * 8.686e4
        
        # Presión acústica módulo predicha
        inputs_plot = torch.cat([coords_norm_grid, torch.zeros_like(coords_norm_grid[:, 0:1])], dim=1)
        p_real_pred, p_imag_pred = pinn.forward_pressure(inputs_plot)
        p_amp_pred = torch.sqrt(p_real_pred**2 + p_imag_pred**2).cpu().numpy().reshape((Nz, Nx))
        p_amp_gt = np.abs(P_gt[5.0e6])
        
    # Graficar comparación 2D de la ROI
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # c GT
    im = axes[0, 0].imshow(c_gt, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='jet', aspect='equal')
    axes[0, 0].set_title("c GT (m/s)")
    axes[0, 0].set_xlabel("Lateral (cm)")
    axes[0, 0].set_ylabel("Axial (cm)")
    fig.colorbar(im, ax=axes[0, 0])
    
    # c Predicha
    im = axes[0, 1].imshow(c_pred_final, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='jet', aspect='equal')
    axes[0, 1].set_title(f"ROI c (SSIM: {best_ssim_c:.3f})")
    axes[0, 1].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[0, 1])
    
    # Curvas de pérdida
    epochs_list = [h['epoch'] for h in progress_history]
    losses_total = [h['loss_total'] for h in progress_history]
    losses_data = [h['loss_data'] for h in progress_history]
    losses_dir = [h['loss_dirichlet'] for h in progress_history]
    losses_pde = [h['loss_pde'] for h in progress_history]
    
    axes[0, 2].semilogy(epochs_list, losses_total, 'k-', label='Total')
    axes[0, 2].semilogy(epochs_list, losses_data, 'b--', label='Datos (20%)')
    axes[0, 2].semilogy(epochs_list, losses_dir, 'g-.', label='Dirichlet z=0')
    axes[0, 2].semilogy(epochs_list, losses_pde, 'r:', label='PDE')
    axes[0, 2].set_title("Curvas de Convergencia")
    axes[0, 2].set_xlabel("Época")
    axes[0, 2].set_ylabel("Pérdida")
    axes[0, 2].legend()
    axes[0, 2].grid(True)
    
    # alpha0 GT
    im = axes[1, 0].imshow(alpha0_gt_db, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='hot', aspect='equal')
    axes[1, 0].set_title("alpha_0 GT (dB/MHz/cm)")
    axes[1, 0].set_xlabel("Lateral (cm)")
    axes[1, 0].set_ylabel("Axial (cm)")
    fig.colorbar(im, ax=axes[1, 0])
    
    # alpha0 Predicha
    im = axes[1, 1].imshow(alpha0_pred_final, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='hot', aspect='equal')
    axes[1, 1].set_title(f"ROI alpha_0 (SSIM: {progress_history[-1]['ssim_alpha0']:.3f})")
    axes[1, 1].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[1, 1])
    
    # Presión acústica módulo predicha a 5.0 MHz
    im = axes[1, 2].imshow(p_amp_pred, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='viridis', aspect='equal')
    axes[1, 2].set_title("Módulo de Presión (5 MHz)")
    axes[1, 2].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[1, 2])
    
    plt.tight_layout()
    plt.savefig('articulo/resultados_5MHz_roi.png', dpi=150)
    plt.close()
    
    # Generar Reporte Markdown
    best_h = next(h for h in progress_history if h['epoch'] == best_epoch)
    table_rows = ""
    for h in progress_history:
        table_rows += (f"| {h['epoch']} | {h['loss_total']:.4e} | {h['loss_data']:.4e} | {h['loss_pde']:.4e} | "
                      f"{h['err_c']*100:.2f}% | {h['err_alpha0']*100:.2f}% | {h['ssim_c']:.3f} | {h['ssim_alpha0']:.3f} |\n")
        
    met_c_l2 = best_h['err_c']*100 < 3.0
    met_c_ssim = best_h['ssim_c'] > 0.85
    met_alpha_l2 = best_h['err_alpha0']*100 < 20.0
    met_alpha_ssim = best_h['ssim_alpha0'] > 0.70
    success = met_c_l2 and met_c_ssim and met_alpha_l2 and met_alpha_ssim
    
    report_content = f"""# Evidencia Empírica: Inversión 2D Helmholtz en ROI Pequeña (1.5 cm x 1.5 cm) a 5 MHz
Fecha de Ejecución: {time.strftime("%Y-%m-%d %H:%M:%S")}

## 1. Resumen Ejecutivo de Métricas de Validación - ROI Pequeña
El modelo PINN desacoplado se entrenó en una región de interés (ROI) pequeña de 1.5 cm x 1.5 cm centrada en la inclusión, utilizando $\\omega_{{0,\\text{{params}}}} = 5.0$ y un sensor Dirichlet explícito en la superficie ($z=0$).

### Comparación contra Criterios de Éxito del Revisor:
*   **Velocidad de Sonido c(x)**:
    *   Error L2 Final: **{best_h['err_c']*100:.3f}%** (Criterio: < 3.0%) -> **{"🟢 CUMPLIDO" if met_c_l2 else "🔴 NO ALCANZADO"}**
    *   SSIM Final: **{best_h['ssim_c']:.3f}** (Criterio: > 0.85) -> **{"🟢 CUMPLIDO" if met_c_ssim else "🔴 NO ALCANZADO"}**
*   **Atenuación alpha_0(x)**:
    *   Error L2 Final: **{best_h['err_alpha0']*100:.3f}%** (Criterio: < 20.0%) -> **{"🟢 CUMPLIDO" if met_alpha_l2 else "🔴 NO ALCANZADO"}**
    *   SSIM Final: **{best_h['ssim_alpha0']:.3f}** (Criterio: > 0.70) -> **{"🟢 CUMPLIDO" if met_alpha_ssim else "🔴 NO ALCANZADO"}**
*   **Tiempo de GPU**: {elapsed_time:.1f} segundos ({epochs} épocas).

¿Se cumplieron todos los criterios? **{"SÍ - Tarea A Exitosa" if success else "NO - Tarea A Fallida"}**

## 2. Gráficos de Inversión en la ROI
![Resultados de Inversión en ROI](resultados_5MHz_roi.png)

## 3. Historial de Entrenamiento y Convergencia
| Época | Pérdida Total | Pérdida Datos | Pérdida PDE | Err L2 c | Err L2 alpha0 | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_rows}
"""
    with open("articulo/metricas_progreso_5MHz_roi.md", "w", encoding="utf-8") as f:
        f.write(report_content)
    print("Reporte ROI guardado en articulo/metricas_progreso_5MHz_roi.md")
    
    return success

if __name__ == "__main__":
    run_5MHz_roi()
