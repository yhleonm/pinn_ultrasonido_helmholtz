import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scratch.evaluate_factorial import inspect_config
from data_engine.loader import CUBDLLoader

def inspect_config_details(c_map_path, dataset_path="data_cubdl_prueba_1/MYO001.hdf5"):
    if not os.path.exists(c_map_path):
        return None
        
    loader = CUBDLLoader(dataset_path)
    rx_x = loader.rx_positions[:, 0]
    arr = np.load(c_map_path)
    Nz, Nx = arr.shape
    
    x_grid = np.linspace(min(rx_x), max(rx_x), Nx)
    z_grid = np.linspace(0.0, 0.032, Nz)
        
    x_tumor, z_tumor = 2.00e-3, 18.50e-3
    roi_half = 0.0075
    
    X2D, Z2D = np.meshgrid(x_grid - x_tumor, z_grid - z_tumor)
    mask = (x_grid >= x_tumor - roi_half) & (x_grid <= x_tumor + roi_half)
    mask_z = (z_grid >= z_tumor - roi_half) & (z_grid <= z_tumor + roi_half)
    
    z_idx = np.where(mask_z)[0]
    x_idx = np.where(mask)[0]
    
    z0, z1 = z_idx.min(), z_idx.max() + 1
    x0, x1 = x_idx.min(), x_idx.max() + 1
    
    roi_2d = arr[z0:z1, x0:x1]
    X2D_roi = X2D[z0:z1, x0:x1]
    Z2D_roi = Z2D[z0:z1, x0:x1]
    
    Nz_roi, Nx_roi = roi_2d.shape
    mz = int(0.15 * Nz_roi)
    mx = int(0.15 * Nx_roi)
    
    interior_mask = np.zeros_like(roi_2d, dtype=bool)
    interior_mask[mz:-mz, mx:-mx] = True
    
    roi_int = np.where(interior_mask, roi_2d, -1e9)
    max_idx = np.unravel_index(np.argmax(roi_int), roi_int.shape)
    
    dx_mm = X2D_roi[max_idx] * 1e3
    dz_mm = Z2D_roi[max_idx] * 1e3
    dist_mm = np.sqrt(dx_mm**2 + dz_mm**2)
    
    dist_from_center = np.sqrt(X2D_roi**2 + Z2D_roi**2) * 1e3 # in mm
    inclusion_mask = dist_from_center <= 3.0
    background_mask = (dist_from_center > 4.5) & (dist_from_center <= 7.0)
    
    inclusion_mask_int = inclusion_mask & interior_mask
    background_mask_int = background_mask & interior_mask
    
    val_inc = roi_2d[inclusion_mask_int]
    val_bg = roi_2d[background_mask_int]
    
    mean_inc = val_inc.mean()
    mean_bg = val_bg.mean()
    var_bg = val_bg.var()
    std_bg = np.sqrt(var_bg)
    
    cnr = np.abs(mean_inc - mean_bg) / (std_bg + 1e-12)
    cnr_db = 20 * np.log10(cnr) if cnr > 0 else -99.0
    
    return {
        "max_val": roi_2d[max_idx],
        "dist_mm": dist_mm,
        "mean_inc": mean_inc,
        "mean_bg": mean_bg,
        "std_bg": std_bg,
        "cnr_db": cnr_db,
    }

configs = [
    {"name": "Config A (PML 15, Sharp)", "path": "articulo/c_map_MYO001_lu_verify.npy"},
    {"name": "Config B (PML 35, Sharp)", "path": "articulo/c_map_MYO001_pml35_sharp.npy"},
    {"name": "Config C (PML 15, Taper)", "path": "articulo/c_map_MYO001_pml15_taper.npy"},
    {"name": "Config D (PML 35, Taper)", "path": "articulo/c_map_MYO001_pml35_taper.npy"},
]

print("=" * 95)
print(f"{'Config':<25} | {'Mean Inc':<10} | {'Mean BG':<10} | {'Std BG':<10} | {'CNR (dBr)':<10} | {'Peak Dist':<10}")
print("=" * 95)
for cfg in configs:
    res = inspect_config_details(cfg["path"])
    if res:
        print(f"{cfg['name']:<25} | {res['mean_inc']:>10.1f} | {res['mean_bg']:>10.1f} | {res['std_bg']:>10.2f} | {res['cnr_db']:>10.2f} | {res['dist_mm']:>10.2f}")
print("=" * 95)
