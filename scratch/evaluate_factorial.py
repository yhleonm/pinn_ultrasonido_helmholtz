"""
Script to evaluate and compare the results of the 4 factorial configurations on MYO001.
Computes:
- Peak velocity location in the ROI interior.
- Distance of peak velocity to tumor center (mm).
- Diagnosis (LIMPIO if distance <= 3.0mm, else FUGA).
- CNR (Contrast-to-Noise Ratio) in dBr for the reconstructed c map.
"""
import numpy as np
import os
import sys
# Add project root to sys.path to allow absolute imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from data_engine.loader import CUBDLLoader

def inspect_config(c_map_path, dataset_path="data_cubdl_prueba_1/MYO001.hdf5"):
    if not os.path.exists(c_map_path):
        return None
        
    loader = CUBDLLoader(dataset_path)
    rx_x = loader.rx_positions[:, 0]
    
    arr = np.load(c_map_path)
    Nz, Nx = arr.shape
    
    # Reconstruct grids
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
    
    # Calculate CNR for c map
    # Target inclusion: within 3mm from tumor center
    # Background: outside 4.5mm from tumor center, within ROI
    dist_from_center = np.sqrt(X2D_roi**2 + Z2D_roi**2) * 1e3 # in mm
    inclusion_mask = dist_from_center <= 3.0
    background_mask = (dist_from_center > 4.5) & (dist_from_center <= 7.0)
    
    inclusion_mask_int = inclusion_mask & interior_mask
    background_mask_int = background_mask & interior_mask
    
    mean_inc = roi_2d[inclusion_mask_int].mean()
    mean_bg = roi_2d[background_mask_int].mean()
    var_bg = roi_2d[background_mask_int].var() + 1e-12
    
    cnr = np.abs(mean_inc - mean_bg) / np.sqrt(var_bg)
    cnr_db = 20 * np.log10(cnr) if cnr > 0 else -99.0
    
    return {
        "max_val": roi_2d[max_idx],
        "dist_mm": dist_mm,
        "is_clean": dist_mm <= 3.0,
        "cnr_db": cnr_db,
    }

def main():
    configs = [
        {"name": "Config A (PML 15, Sharp)", "path": "articulo/c_map_MYO001_lu_verify.npy"},
        {"name": "Config B (PML 35, Sharp)", "path": "articulo/c_map_MYO001_pml35_sharp.npy"},
        {"name": "Config C (PML 15, Taper)", "path": "articulo/c_map_MYO001_pml15_taper.npy"},
        {"name": "Config D (PML 35, Taper)", "path": "articulo/c_map_MYO001_pml35_taper.npy"},
    ]
    
    print("=" * 75)
    print(f"{'Configuration':<25} | {'Peak c (m/s)':<12} | {'Dist (mm)':<10} | {'Diagnosis':<10} | {'CNR (dBr)':<10}")
    print("=" * 75)
    
    for cfg in configs:
        res = inspect_config(cfg["path"])
        if res is None:
            print(f"{cfg['name']:<25} | {'PENDING':<12} | {'-':<10} | {'-':<10} | {'-':<10}")
        else:
            diag = "LIMPIO" if res["is_clean"] else "FUGA"
            print(f"{cfg['name']:<25} | {res['max_val']:>12.1f} | {res['dist_mm']:>10.2f} | {diag:<10} | {res['cnr_db']:>10.2f}")
    print("=" * 75)

if __name__ == "__main__":
    main()
