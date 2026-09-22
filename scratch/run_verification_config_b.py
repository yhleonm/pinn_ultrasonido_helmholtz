"""
Verification runner and evaluator for Config B on MYO002, MYO004, and MYO006.
Executes the FWI pipeline with npml=35 and sharp ROI for the target datasets,
then performs a comparative analysis against the baseline Config A (npml=15).
"""
import os
import sys
import time
import numpy as np

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from main_5MHz_fwi_real_data import run_fwi_real_data, delay_and_sum
from data_engine.loader import CUBDLLoader

def get_tumor_location_das(loader):
    x_bmode = np.arange(-0.02, 0.02 + 0.5e-3, 0.5e-3)
    z_bmode = np.arange(0.002, 0.035 + 0.5e-3, 0.5e-3)
    bmode_img = delay_and_sum(loader, x_bmode, z_bmode, compound=True)
    
    depth_mask = z_bmode > 0.010
    img_valid = bmode_img[depth_mask, :]
    max_idx = np.unravel_index(np.argmax(img_valid), img_valid.shape)
    z_tumor = z_bmode[depth_mask][max_idx[0]]
    x_tumor = x_bmode[max_idx[1]]
    return x_tumor, z_tumor

def evaluate_run(c_map_path, dataset_path):
    if not os.path.exists(c_map_path):
        return None
        
    loader = CUBDLLoader(dataset_path)
    rx_x = loader.rx_positions[:, 0]
    
    arr = np.load(c_map_path)
    Nz, Nx = arr.shape
    
    x_grid = np.linspace(min(rx_x), max(rx_x), Nx)
    z_grid = np.linspace(0.0, 0.032, Nz)
    
    x_tumor, z_tumor = get_tumor_location_das(loader)
    
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
    
    dist_from_center = np.sqrt(X2D_roi**2 + Z2D_roi**2) * 1e3
    inclusion_mask = dist_from_center <= 3.0
    background_mask = (dist_from_center > 4.5) & (dist_from_center <= 7.0)
    
    inclusion_mask_int = inclusion_mask & interior_mask
    background_mask_int = background_mask & interior_mask
    
    val_inc = roi_2d[inclusion_mask_int]
    val_bg = roi_2d[background_mask_int]
    
    mean_inc = val_inc.mean()
    mean_bg = val_bg.mean()
    std_bg = val_bg.std()
    
    cnr = np.abs(mean_inc - mean_bg) / (std_bg + 1e-12)
    cnr_db = 20 * np.log10(cnr) if cnr > 0 else -99.0
    
    return {
        "max_val": roi_2d[max_idx],
        "dist_mm": dist_mm,
        "mean_inc": mean_inc,
        "mean_bg": mean_bg,
        "std_bg": std_bg,
        "cnr_db": cnr_db,
        "x_tumor_das": x_tumor,
        "z_tumor_das": z_tumor
    }

def main():
    datasets = {
        "MYO002": "data_cubdl_prueba_1/MYO002.hdf5",
        "MYO004": "data_cubdl_prueba_1/MYO004.hdf5",
        "MYO006": "data_cubdl_prueba_1/MYO006.hdf5",
    }
    
    print("=" * 80)
    print("RUNNING VERIFICATION SCRIPT FOR CONFIG B ON MYO002, MYO004, MYO006")
    print("=" * 80)
    
    results = {}
    
    for name, path in datasets.items():
        if not os.path.exists(path):
            print(f"Error: {path} not found.")
            continue
            
        print(f"\n>>> Checking Config B for {name}...")
        c_map_b = f"articulo/c_map_{name}_pml35_sharp.npy"
        
        if not os.path.exists(c_map_b):
            print(f"Config B results not found. Running training on {name} with Config B parameters...")
            start_t = time.time()
            try:
                run_fwi_real_data(
                    file_path=path,
                    seed=42,
                    output_suffix="_pml35_sharp",
                    npml=35,
                    roi_taper=False
                )
                elapsed = time.time() - start_t
                print(f"Finished training for {name} with Config B in {elapsed/60:.1f} minutes.")
            except Exception as e:
                print(f"Error training {name} with Config B: {e}")
                import traceback
                traceback.print_exc()
                continue
        else:
            print(f"Config B results found for {name}. Skipping training.")
            
        # Also ensure Config A exists (from baseline, which has no suffix)
        c_map_a = f"articulo/c_map_{name}.npy"
        if not os.path.exists(c_map_a):
            print(f"Warning: Baseline Config A (no suffix) results not found for {name}. Running Config A first...")
            try:
                run_fwi_real_data(
                    file_path=path,
                    seed=42,
                    output_suffix="",
                    npml=15,
                    roi_taper=False
                )
            except Exception as e:
                print(f"Error training baseline for {name}: {e}")
                continue
                
        # Evaluate both configs
        res_a = evaluate_run(c_map_a, path)
        res_b = evaluate_run(c_map_b, path)
        
        results[name] = {"A": res_a, "B": res_b}
        
    # Print comparison table
    print("\n" + "=" * 90)
    print("COMPARATIVE EVALUATION: CONFIG A (PML 15, Sharp) vs CONFIG B (PML 35, Sharp)")
    print("=" * 90)
    print(f"{'Dataset':<10} | {'Config':<10} | {'Peak c (m/s)':<12} | {'Dist (mm)':<10} | {'Verdict':<10} | {'Mean Inc':<10} | {'Mean BG':<10} | {'CNR (dBr)':<10}")
    print("-" * 90)
    
    for name in datasets.keys():
        if name not in results or results[name]["A"] is None or results[name]["B"] is None:
            print(f"{name:<10} | PENDING OR FAILED")
            continue
            
        res_a = results[name]["A"]
        res_b = results[name]["B"]
        
        verdict_a = "LIMPIO" if res_a["dist_mm"] <= 3.0 else "FUGA"
        verdict_b = "LIMPIO" if res_b["dist_mm"] <= 3.0 else "FUGA"
        
        print(f"{name:<10} | {'Config A':<10} | {res_a['max_val']:>12.1f} | {res_a['dist_mm']:>10.2f} | {verdict_a:<10} | {res_a['mean_inc']:>10.1f} | {res_a['mean_bg']:>10.1f} | {res_a['cnr_db']:>10.2f}")
        print(f"{'':<10} | {'Config B':<10} | {res_b['max_val']:>12.1f} | {res_b['dist_mm']:>10.2f} | {verdict_b:<10} | {res_b['mean_inc']:>10.1f} | {res_b['mean_bg']:>10.1f} | {res_b['cnr_db']:>10.2f}")
        print("-" * 90)
    print("=" * 90)

if __name__ == "__main__":
    main()
