import numpy as np
import json
from pathlib import Path

# Rutas
articulo_dir = Path(r"C:\Users\yelem\workspace\proyecto_pinn_ultrasonido_mama_fase_2\articulo")

# Buscar todos los JSON de métricas
metric_files = sorted(articulo_dir.glob("metrics_*.json"))

print(f"Encontrados {len(metric_files)} archivos de métricas.\n")
print(f"{'Dataset':<10} | {'Tumor X (mm)':<12} | {'Tumor Z (mm)':<12} | {'Peak Dist c (mm)':<16} | {'Peak Dist alpha (mm)':<20}")
print("-" * 80)

roi_half = 0.0075  # 7.5 mm en metros
MARGIN_FRAC = 0.15  # 15% de margen

for mf in metric_files:
    dataset_id = mf.stem.replace("metrics_", "")
    
    c_path = articulo_dir / f"c_map_{dataset_id}.npy"
    alpha_path = articulo_dir / f"alpha_map_{dataset_id}.npy"
    
    if not c_path.exists() or not alpha_path.exists():
        continue
        
    # Cargar datos
    c_map = np.load(c_path)
    alpha_map = np.load(alpha_path)
    
    with open(mf, 'r') as f:
        metrics = json.load(f)
        
    tumor_x_mm = metrics.get("tumor_x_mm")
    tumor_z_mm = metrics.get("tumor_z_mm")
    
    if tumor_x_mm is None or tumor_z_mm is None:
        print(f"{dataset_id:<10} | Sin coordenadas de tumor en el JSON.")
        continue
        
    tumor_x = tumor_x_mm / 1e3
    tumor_z = tumor_z_mm / 1e3
        
    Nz, Nx = c_map.shape
    x_grid = np.linspace(-0.01905, 0.01905, Nx)
    z_grid = np.linspace(0.0, 0.032, Nz)
    X2D, Z2D = np.meshgrid(x_grid, z_grid)
    
    # Extraer sub-array ROI 2D
    mask_x = (x_grid >= tumor_x - roi_half) & (x_grid <= tumor_x + roi_half)
    mask_z = (z_grid >= tumor_z - roi_half) & (z_grid <= tumor_z + roi_half)
    z_idx = np.where(mask_z)[0]
    x_idx = np.where(mask_x)[0]
    
    z0, z1 = z_idx.min(), z_idx.max() + 1
    x0, x1 = x_idx.min(), x_idx.max() + 1
    
    X2D_roi = X2D[z0:z1, x0:x1]
    Z2D_roi = Z2D[z0:z1, x0:x1]
    
    roi_c = c_map[z0:z1, x0:x1]
    roi_alpha = alpha_map[z0:z1, x0:x1]
    
    Nz_roi, Nx_roi = roi_c.shape
    mz = int(MARGIN_FRAC * Nz_roi)
    mx = int(MARGIN_FRAC * Nx_roi)
    
    interior_mask = np.zeros_like(roi_c, dtype=bool)
    interior_mask[mz:-mz, mx:-mx] = True
    
    # Calcular para c
    roi_c_int = np.where(interior_mask, roi_c, -1e9)
    max_c_idx = np.unravel_index(np.argmax(roi_c_int), roi_c_int.shape)
    dist_c = np.sqrt((X2D_roi[max_c_idx] - tumor_x)**2 + (Z2D_roi[max_c_idx] - tumor_z)**2) * 1e3
    
    # Calcular para alpha
    roi_alpha_int = np.where(interior_mask, roi_alpha, -1e9)
    max_alpha_idx = np.unravel_index(np.argmax(roi_alpha_int), roi_alpha_int.shape)
    dist_alpha = np.sqrt((X2D_roi[max_alpha_idx] - tumor_x)**2 + (Z2D_roi[max_alpha_idx] - tumor_z)**2) * 1e3
    
    print(f"{dataset_id:<10} | {tumor_x_mm:<12.2f} | {tumor_z_mm:<12.2f} | {dist_c:<16.2f} | {dist_alpha:<20.2f}")
    
    # Actualizar el JSON agregando peak_distance_alpha_mm y restaurando peak_distance_mm (c)
    metrics["peak_distance_alpha_mm"] = float(dist_alpha)
    metrics["peak_distance_mm"] = float(dist_c)
    
    with open(mf, 'w') as f:
        json.dump(metrics, f, indent=2)
        
print("\nLos archivos metrics_*.json han sido actualizados con los peak_distance correctos (excluyendo el 15% de margen).")
