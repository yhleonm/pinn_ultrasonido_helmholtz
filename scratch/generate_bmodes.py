import sys
from pathlib import Path
import json
import numpy as np
import torch

# Add workspace to path
sys.path.append(str(Path(__file__).parent.parent))

from data_engine.loader import CUBDLLoader

def delay_and_sum_gpu(loader, x_grid, z_grid, num_tx_to_use=11):
    """
    GPU-accelerated and fully-vectorized Delay-and-Sum (DAS) beamformer.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device} for beamforming")
    
    fs = float(loader.fs)
    c0 = float(loader.c0)
    t0 = float(loader.t0)
    rx_x = torch.from_numpy(loader.rx_positions[:, 0]).to(device, dtype=torch.float32) # (Nr,)
    
    # 2D Grid
    x_grid_t = torch.from_numpy(x_grid).to(device, dtype=torch.float32) # (Nx,)
    z_grid_t = torch.from_numpy(z_grid).to(device, dtype=torch.float32) # (Nz,)
    
    # Create meshgrid (Nx, Nz)
    X, Z = torch.meshgrid(x_grid_t, z_grid_t, indexing='ij')
    Nx, Nz = X.shape
    
    rf_data = loader.rf # (Nt, Nr, Ntx) or (Nt, Nr)
    if len(rf_data.shape) == 2:
        rf_data = rf_data[:, :, np.newaxis]
    
    Nt, Nr, Ntx = rf_data.shape
    
    # Select subset of transmitters to compound to speed up even more
    if Ntx > num_tx_to_use:
        tx_indices = np.linspace(0, Ntx - 1, num_tx_to_use, dtype=int)
    else:
        tx_indices = np.arange(Ntx)
        
    compounded_image = torch.zeros((Nx, Nz), dtype=torch.complex64, device=device)
    
    from scipy.signal import hilbert
    
    # Pre-calculate rx distances on GPU to save memory/computation
    xr = rx_x.view(Nr, 1, 1)
    dist_rx = torch.sqrt((X.unsqueeze(0) - xr)**2 + Z.unsqueeze(0)**2) # (Nr, Nx, Nz)
    t_rx = dist_rx / c0 # (Nr, Nx, Nz)
    
    for tx_idx in tx_indices:
        rf = rf_data[:, :, tx_idx] # (Nt, Nr)
        tx_delay = loader.tx_delays[tx_idx, :] # (Nr,)
        tx_center_idx = np.argmin(tx_delay)
        x_tx = rx_x[tx_center_idx]
        
        # Distance from transmitter to pixels
        dist_tx = torch.sqrt((X - x_tx)**2 + Z**2) # (Nx, Nz)
        t_tx = dist_tx / c0
        
        # Hilbert transform (done on CPU, fast enough)
        rf_complex_np = hilbert(rf, axis=0) # (Nt, Nr)
        rf_complex = torch.from_numpy(rf_complex_np).to(device, dtype=torch.complex64) # (Nt, Nr)
        
        # Total time of flight: (Nr, Nx, Nz)
        tx_del_t = torch.from_numpy(tx_delay).to(device, dtype=torch.float32).view(Nr, 1, 1)
        t_total = t_tx.unsqueeze(0) + t_rx + tx_del_t
        
        # Sample index: (Nr, Nx, Nz)
        sample_idx = (t_total - t0) * fs
        
        # Interpolate
        sample_idx_clipped = torch.clamp(sample_idx, 0.0, float(Nt - 2))
        idx_low = sample_idx_clipped.long()
        idx_high = idx_low + 1
        weight_high = sample_idx_clipped - idx_low.float()
        weight_low = 1.0 - weight_high
        
        r_indices = torch.arange(Nr, device=device).view(Nr, 1, 1).expand(Nr, Nx, Nz)
        
        val_low = rf_complex[idx_low, r_indices] # (Nr, Nx, Nz)
        val_high = rf_complex[idx_high, r_indices] # (Nr, Nx, Nz)
        
        interpolated = val_low * weight_low + val_high * weight_high
        image_tx = torch.sum(interpolated, dim=0) # (Nx, Nz)
        compounded_image += image_tx
        
    # Return as [Nz, Nx] NumPy array
    compounded_np = compounded_image.cpu().numpy()
    return np.abs(compounded_np).T

def main():
    workspace_dir = Path("C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2")
    articulo_dir = workspace_dir / "articulo"
    data_dir = workspace_dir / "data_cubdl_prueba_1"
    
    datasets = ["MYO001", "MYO002", "MYO003", "MYO004", "MYO005", "MYO006", "TSH002"]
    
    for ds in datasets:
        print(f"\nProcessing {ds}...")
        
        # Load maps to get shape
        c_map_path = articulo_dir / f"c_map_{ds}.npy"
        if not c_map_path.exists():
            print(f"c_map not found for {ds}, skipping.")
            continue
            
        c_arr = np.load(c_map_path)
        Nz, Nx = c_arr.shape
        
        # Load metrics sidecar
        sidecar_path = articulo_dir / f"metrics_{ds}.json"
        if sidecar_path.exists():
            with open(sidecar_path, "r") as f:
                metrics = json.load(f)
            x_min_mm = metrics.get("x_min_mm", -19.05)
            x_max_mm = metrics.get("x_max_mm", 19.05)
            z_min_mm = metrics.get("z_min_mm", 0.0)
            z_max_mm = metrics.get("z_max_mm", 32.0)
        else:
            print(f"Metrics not found for {ds}, using defaults.")
            x_min_mm = -19.05
            x_max_mm = 19.05
            z_min_mm = 0.0
            z_max_mm = 32.0
            
        # Find HDF5 file
        h5_path = data_dir / f"{ds}.hdf5"
        if not h5_path.exists():
            # Check TSH002 subdirectory
            h5_path = data_dir / ds / f"{ds}.hdf5"
            
        if not h5_path.exists():
            print(f"HDF5 file not found for {ds}, skipping.")
            continue
            
        print(f"Loading {h5_path.name}...")
        loader = CUBDLLoader(str(h5_path))
        
        # Grid for B-mode (matching reconstruction grid exactly)
        x_grid = np.linspace(x_min_mm * 1e-3, x_max_mm * 1e-3, Nx)
        z_grid = np.linspace(z_min_mm * 1e-3, z_max_mm * 1e-3, Nz)
        
        print(f"Reconstructing B-mode ({Nz}x{Nx})...")
        bmode_img = delay_and_sum_gpu(loader, x_grid, z_grid, num_tx_to_use=11)
        
        # Normalize and log compress
        bmode_normalized = bmode_img / np.max(bmode_img)
        bmode_log = 20 * np.log10(bmode_normalized + 1e-5)
        
        out_path = articulo_dir / f"bmode_{ds}.npy"
        np.save(out_path, bmode_log)
        print(f"Saved B-mode to {out_path}")

if __name__ == "__main__":
    main()
