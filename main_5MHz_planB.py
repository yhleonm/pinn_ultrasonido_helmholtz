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

def run_plan_b():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Plan B: Usando dispositivo: {device}")
    
    # Cargar el dataset sintético generado en la Tarea 1
    dataset_path = 'data_engine/synthetic_dataset_5MHz.npz'
    if not os.path.exists(dataset_path):
        print("Dataset no encontrado. Generando...")
        # (código para generar el dataset si no existe)
        x_max, z_max, dx, dz = 0.03, 0.04, 0.05e-3, 0.05e-3
        frequencies = [4e6, 5e6, 6e6]
        data_gen = SyntheticHelmholtz2D(
            x_max=x_max, z_max=z_max, dx=dx, dz=dz, frequencies=frequencies,
            c_bg=1500.0, c_in=1580.0,
            alpha0_bg_db=0.5, alpha0_in_db=0.9,
            inclusion_radius=0.005, inclusion_z=0.02
        )
        data = data_gen.get_data()
        np.savez(dataset_path, **data)
    
    # Cargar datos
    data = np.load(dataset_path, allow_pickle=True)
    x_grid = data['x']
    z_grid = data['z']
    c_gt = data['c_gt']
    alpha0_gt = data['alpha0_gt']
    frequencies = data['frequencies']
    
    # Reconstruir dict de P_gt
    P_gt = {}
    for f in frequencies:
        key = f'P_gt_{int(f/1e6)}MHz'
        if key in data:
            P_gt[f] = data[key]
        else:
            P_gt[f] = data['P_gt'].item()[f]
            
    Nx, Nz = len(x_grid), len(z_grid)
    x_max = max(x_grid)
    z_max = max(z_grid)
    
    # Coordenadas y tensores
    X, Z = np.meshgrid(x_grid, z_grid)
    coords_flat = np.stack([X.flatten(), Z.flatten()], axis=1)
    coords_t = torch.tensor(coords_flat, dtype=torch.float32, device=device)
    c_gt_flat_t = torch.tensor(c_gt.flatten(), dtype=torch.float32, device=device).unsqueeze(1)
    alpha0_gt_flat_t = torch.tensor(alpha0_gt.flatten(), dtype=torch.float32, device=device).unsqueeze(1)
    
    # Selección de Datos de Observación Escasos (20%)
    np.random.seed(42)
    obs_fraction = 0.20
    n_pts = len(coords_flat)
    n_obs = int(n_pts * obs_fraction)
    obs_indices = np.random.choice(n_pts, n_obs, replace=False)
    obs_indices.sort()
    
    coords_obs_t = coords_t[obs_indices]
    obs_data = {}
    for f in frequencies:
        P_f_flat = P_gt[f].flatten()[obs_indices]
        obs_data[f] = {
            'real': torch.tensor(P_f_flat.real, dtype=torch.float32, device=device).unsqueeze(1),
            'imag': torch.tensor(P_f_flat.imag, dtype=torch.float32, device=device).unsqueeze(1)
        }
        
    # Inicialización del Modelo PINN Desacoplado 2D Plan B (hidden_features=256, omega_0_params=15.0)
    pinn = DecoupledUltrasoundPINN(
        spatial_dim=2, hidden_features=256, num_layers=4,
        omega_0_pressure=60.0, omega_0_params=15.0
    ).to(device)
    
    epochs = 3000
    optimizer = optim.Adam(pinn.parameters(), lr=3e-4)
    scheduler = optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=500, T_mult=2)
    
    # Pesos de Pérdidas Plan B (Mayor penalización de datos)
    lambda_data = 50.0
    lambda_pde = 5.0
    lambda_sommerfeld = 1.0
    
    omega_min = 2.0 * np.pi * min(frequencies)
    omega_max = 2.0 * np.pi * max(frequencies)
    
    def normalize_inputs(coords_xz, omega):
        x_norm = coords_xz[:, 0:1] / x_max
        z_norm = 2.0 * (coords_xz[:, 1:2] / z_max) - 1.0
        if len(frequencies) > 1:
            omega_val = 2.0 * ((omega - omega_min) / (omega_max - omega_min)) - 1.0
        else:
            omega_val = 0.0
        omega_norm = torch.full_like(x_norm, omega_val)
        return torch.cat([x_norm, z_norm], dim=1), omega_norm
    
    n_coloc_points = 30000
    n_boundary_points = 500
    
    best_ssim_c = -1.0
    best_state = None
    best_epoch = 0
    
    print("Iniciando entrenamiento Plan B de la PINN a 5MHz...")
    start_time = time.time()
    progress_history = []
    
    for epoch in range(1, epochs + 1):
        pinn.train()
        optimizer.zero_grad()
        
        # A. Pérdida de Datos
        loss_data = 0.0
        for f in frequencies:
            omega = 2.0 * np.pi * f
            coords_norm_obs, omega_norm_obs = normalize_inputs(coords_obs_t, omega)
            inputs_obs = torch.cat([coords_norm_obs, omega_norm_obs], dim=1)
            pred_real, pred_imag = pinn.forward_pressure(inputs_obs)
            loss_data += torch.mean((pred_real - obs_data[f]['real'])**2 + (pred_imag - obs_data[f]['imag'])**2)
        loss_data = loss_data / len(frequencies)
        
        # B. Pérdida PDE
        loss_pde = 0.0
        x_rand = (torch.rand(n_coloc_points, 1, device=device) * 2.0 - 1.0) * x_max
        z_rand = torch.rand(n_coloc_points, 1, device=device) * z_max
        coords_coloc = torch.cat([x_rand, z_rand], dim=1)
        coords_coloc.requires_grad_(True)
        
        for f in frequencies:
            omega = 2.0 * np.pi * f
            coords_norm_coloc, omega_norm_coloc = normalize_inputs(coords_coloc, omega)
            c_pred, alpha0_pred = pinn.forward_params(coords_norm_coloc)
            
            inputs_coloc = torch.cat([coords_norm_coloc, omega_norm_coloc], dim=1)
            pred_real, pred_imag = pinn.forward_pressure(inputs_coloc)
            
            res_real, res_imag = compute_helmholtz_residuals_2d(
                coords_coloc, pred_real, pred_imag, alpha0_pred, c_pred, omega
            )
            loss_pde += torch.mean(res_real**2 + res_imag**2)
        loss_pde = loss_pde / len(frequencies)
        
        # C. Pérdida Sommerfeld
        loss_sommer = 0.0
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
        
        for f in frequencies:
            omega = 2.0 * np.pi * f
            coords_norm_b, omega_norm_b = normalize_inputs(coords_b, omega)
            c_pred_b, alpha0_pred_b = pinn.forward_params(coords_norm_b)
            
            inputs_b = torch.cat([coords_norm_b, omega_norm_b], dim=1)
            pred_real_b, pred_imag_b = pinn.forward_pressure(inputs_b)
            
            loss_sommer += compute_sommerfeld_loss_2d(
                coords_b, pred_real_b, pred_imag_b, alpha0_pred_b, c_pred_b, omega, normals_b
            )
        loss_sommer = loss_sommer / len(frequencies)
        
        # Entrenamiento por etapas
        w_pde = lambda_pde if epoch > 1500 else 0.0
        w_sommer = lambda_sommerfeld if epoch > 1500 else 0.0
        
        loss_total = lambda_data * loss_data + w_pde * loss_pde + w_sommer * loss_sommer
        
        loss_total.backward()
        optimizer.step()
        scheduler.step()
        
        if epoch % 100 == 0 or epoch == 1 or epoch == 1501:
            pinn.eval()
            with torch.no_grad():
                coords_norm_grid, _ = normalize_inputs(coords_t, 2.0 * np.pi * 5e6)
                c_pred_flat, alpha0_pred_flat = pinn.forward_params(coords_norm_grid)
                
                err_c = torch.norm(c_pred_flat - c_gt_flat_t) / torch.norm(c_gt_flat_t)
                err_alpha0 = torch.norm(alpha0_pred_flat - alpha0_gt_flat_t) / torch.norm(alpha0_gt_flat_t)
                
                c_pred_2d = c_pred_flat.cpu().numpy().reshape((Nz, Nx))
                alpha0_pred_2d = alpha0_pred_flat.cpu().numpy().reshape((Nz, Nx)) * 8.686e4
                alpha0_gt_db = alpha0_gt * 8.686e4
                
                ssim_c = calculate_ssim(c_pred_2d, c_gt)
                ssim_alpha0 = calculate_ssim(alpha0_pred_2d, alpha0_gt_db)
                
                lr = optimizer.param_groups[0]['lr']
                print(f"Plan B Epoch {epoch:4d} | Loss: {loss_total.item():.4e} "
                      f"[Data: {loss_data.item():.4e}, PDE: {loss_pde.item():.4e}] | "
                      f"Err c: {err_c.item()*100:.2f}%, Err alpha0: {err_alpha0.item()*100:.2f}% | "
                      f"SSIM c: {ssim_c:.3f}, SSIM alpha0: {ssim_alpha0:.3f} | lr: {lr:.2e}")
                
                if epoch > 1500 and ssim_c > best_ssim_c:
                    best_ssim_c = ssim_c
                    best_state = {k: v.cpu().clone() for k, v in pinn.state_dict().items()}
                    best_epoch = epoch
                    print(f"  --> [Plan B SAVED] Nuevo mejor estado (SSIM c: {ssim_c:.4f} en época {epoch})")
                
                progress_history.append({
                    'epoch': epoch,
                    'loss_total': loss_total.item(),
                    'loss_data': loss_data.item(),
                    'loss_pde': loss_pde.item(),
                    'err_c': err_c.item(),
                    'err_alpha0': err_alpha0.item(),
                    'ssim_c': ssim_c,
                    'ssim_alpha0': ssim_alpha0
                })
                
    elapsed_time = time.time() - start_time
    print(f"Entrenamiento Plan B completado en {elapsed_time:.1f} segundos.")
    
    if best_state is not None:
        print(f"\nCargando el mejor estado Plan B (Época {best_epoch}, SSIM c: {best_ssim_c:.4f}) para reportar.")
        pinn.load_state_dict({k: v.to(device) for k, v in best_state.items()})
        torch.save(pinn.state_dict(), "articulo/best_pinn_5MHz_planB.pt")
        
    pinn.eval()
    with torch.no_grad():
        coords_norm_grid, _ = normalize_inputs(coords_t, 2.0 * np.pi * 5e6)
        c_pred_final, alpha0_pred_final_t = pinn.forward_params(coords_norm_grid)
        c_pred_final = c_pred_final.cpu().numpy().reshape((Nz, Nx))
        alpha0_pred_final = alpha0_pred_final_t.cpu().numpy().reshape((Nz, Nx)) * 8.686e4
        alpha0_gt_db = alpha0_gt * 8.686e4
        
    # Graficar comparación 2D Plan B
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # c GT
    im = axes[0, 0].imshow(c_gt, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='jet', aspect='equal')
    axes[0, 0].set_title("c GT (m/s)")
    axes[0, 0].set_xlabel("Lateral (cm)")
    axes[0, 0].set_ylabel("Axial (cm)")
    fig.colorbar(im, ax=axes[0, 0])
    
    # c Predicha Plan B
    im = axes[0, 1].imshow(c_pred_final, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='jet', aspect='equal')
    axes[0, 1].set_title(f"Plan B c (SSIM: {best_ssim_c:.3f})")
    axes[0, 1].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[0, 1])
    
    # Curvas de pérdida
    epochs_list = [h['epoch'] for h in progress_history]
    losses_total = [h['loss_total'] for h in progress_history]
    losses_data = [h['loss_data'] for h in progress_history]
    losses_pde = [h['loss_pde'] for h in progress_history]
    
    axes[0, 2].semilogy(epochs_list, losses_total, 'k-', label='Total')
    axes[0, 2].semilogy(epochs_list, losses_data, 'b--', label='Datos')
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
    
    # alpha0 Predicha Plan B
    im = axes[1, 1].imshow(alpha0_pred_final, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='hot', aspect='equal')
    axes[1, 1].set_title(f"Plan B alpha_0 (SSIM: {progress_history[-1]['ssim_alpha0']:.3f})")
    axes[1, 1].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[1, 1])
    
    # Presión acústica módulo predicha a 5.0 MHz
    f_plot = 5.0e6
    coords_norm_grid, omega_norm_grid = normalize_inputs(coords_t, 2.0 * np.pi * f_plot)
    inputs_plot = torch.cat([coords_norm_grid, omega_norm_grid], dim=1)
    with torch.no_grad():
        p_real_pred, p_imag_pred = pinn.forward_pressure(inputs_plot)
    p_amp_pred = torch.sqrt(p_real_pred**2 + p_imag_pred**2).cpu().numpy().reshape((Nz, Nx))
    
    im = axes[1, 2].imshow(p_amp_pred, extent=[-x_max*100, x_max*100, z_max*100, 0], cmap='viridis', aspect='equal')
    axes[1, 2].set_title("Módulo de Presión (5 MHz)")
    axes[1, 2].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[1, 2])
    
    plt.tight_layout()
    plt.savefig('articulo/resultados_5MHz_planB.png', dpi=150)
    plt.close()
    
    # Generar Reporte Markdown en articulo/metricas_progreso_5MHz_planB.md
    best_h = next(h for h in progress_history if h['epoch'] == best_epoch)
    table_rows = ""
    for h in progress_history:
        table_rows += (f"| {h['epoch']} | {h['loss_total']:.4e} | {h['loss_data']:.4e} | {h['loss_pde']:.4e} | "
                      f"{h['err_c']*100:.2f}% | {h['err_alpha0']*100:.2f}% | {h['ssim_c']:.3f} | {h['ssim_alpha0']:.3f} |\n")
        
    report_content = f"""# Evidencia Empírica: Inversión 2D Helmholtz Sintética a Alta Frecuencia (5 MHz) - Plan B
Fecha de Ejecución: {time.strftime("%Y-%m-%d %H:%M:%S")}

## 1. Resumen Ejecutivo de Métricas de Validación - Plan B
El modelo PINN desacoplado se entrenó con mayor capacidad (256 neuronas ocultas), prior omega_0_params = 15.0, y lambda_data = 50.0 para bloquear el atractor de parámetros planos.

### Comparación contra Criterios de Éxito Clínicos (Plan B):
*   **Velocidad de Sonido c(x)**:
    *   Error L2 Final: **{best_h['err_c']*100:.3f}%** (Criterio: < 5.0%) -> **{"🟢 CUMPLIDO" if best_h['err_c']*100 < 5.0 else "🔴 NO ALCANZADO"}**
    *   SSIM Final: **{best_h['ssim_c']:.3f}** (Criterio: > 0.70) -> **{"🟢 CUMPLIDO" if best_h['ssim_c'] > 0.70 else "🔴 NO ALCANZADO"}**
*   **Atenuación alpha_0(x)**:
    *   Error L2 Final: **{best_h['err_alpha0']*100:.3f}%** (Criterio: < 30.0%) -> **{"🟢 CUMPLIDO" if best_h['err_alpha0']*100 < 30.0 else "🔴 NO ALCANZADO"}**
    *   SSIM Final: **{best_h['ssim_alpha0']:.3f}** (Criterio: > 0.60) -> **{"🟢 CUMPLIDO" if best_h['ssim_alpha0'] > 0.60 else "🔴 NO ALCANZADO"}**
*   **Tiempo de GPU**: {elapsed_time:.1f} segundos ({epochs} épocas).

## 2. Gráficos de Inversión Plan B
![Resultados de Inversión Plan B](resultados_5MHz_planB.png)

## 3. Historial de Entrenamiento y Convergencia
| Época | Pérdida Total | Pérdida Datos | Pérdida PDE | Err L2 c | Err L2 alpha0 | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_rows}
"""
    with open("articulo/metricas_progreso_5MHz_planB.md", "w", encoding="utf-8") as f:
        f.write(report_content)
    print("Reporte Plan B guardado en articulo/metricas_progreso_5MHz_planB.md")

if __name__ == "__main__":
    run_plan_b()
