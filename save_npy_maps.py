import os
import torch
import numpy as np
from data_engine.loader import CUBDLLoader
from models.pinn_net import DecoupledUltrasoundPINN
from main_5MHz_fwi_real_data import delay_and_sum

def generate_maps(file_path):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    file_prefix = os.path.splitext(os.path.basename(file_path))[0]
    model_path = f"articulo/best_pinn_5MHz_fwi_real_{file_prefix}.pt"
    
    if not os.path.exists(model_path):
        print(f"Model path {model_path} does not exist.")
        return
        
    print(f"\nGenerating maps for {file_prefix} using {model_path}...")
    loader = CUBDLLoader(file_path)
    
    # 2. Reconstruct B-mode image using DAS
    x_bmode = np.arange(-0.02, 0.02 + 0.5e-3, 0.5e-3)
    z_bmode = np.arange(0.002, 0.035 + 0.5e-3, 0.5e-3)
    bmode_img = delay_and_sum(loader, x_bmode, z_bmode, compound=True)
    
    # Auto-localize tumor center
    depth_mask = z_bmode > 0.010
    img_valid = bmode_img[depth_mask, :]
    max_idx = np.unravel_index(np.argmax(img_valid), img_valid.shape)
    z_tumor = z_bmode[depth_mask][max_idx[0]]
    x_tumor = x_bmode[max_idx[1]]
    
    print(f"--> Localización del tumor por Beamforming B-mode: x = {x_tumor*1e3:.2f} mm, z = {z_tumor*1e3:.2f} mm")
    
    roi_half_w = 0.0075
    roi_half_h = 0.0075
    x_min_roi, x_max_roi = x_tumor - roi_half_w, x_tumor + roi_half_w
    z_min_roi, z_max_roi = z_tumor - roi_half_h, z_tumor + roi_half_h
    
    def normalize_roi_coords(coords):
        x_norm = (coords[:, 0:1] - x_tumor) / roi_half_w
        z_norm = (coords[:, 1:2] - z_tumor) / roi_half_h
        return torch.cat([x_norm, z_norm], dim=1)
        
    pinn = DecoupledUltrasoundPINN(
        spatial_dim=2, hidden_features=128, num_layers=4,
        omega_0_params=5.0
    ).to(device)
    
    pinn.load_state_dict(torch.load(model_path, map_location=device))
    pinn.eval()
    
    # Grid parameters for Stage 3
    dx = 4.5e-5
    dz = dx
    rx_x = loader.rx_positions[:, 0]
    x_grid = np.arange(min(rx_x), max(rx_x) + dx, dx)
    z_grid = np.arange(0.0, 0.032 + dz, dz)
    Nx, Nz = len(x_grid), len(z_grid)
    
    X, Z = np.meshgrid(x_grid, z_grid)
    coords_flat = np.stack([X.flatten(), Z.flatten()], axis=1)
    
    mask_roi_flat = (coords_flat[:, 0] >= x_min_roi) & (coords_flat[:, 0] <= x_max_roi) & \
                    (coords_flat[:, 1] >= z_min_roi) & (coords_flat[:, 1] <= z_max_roi)
    mask_roi_2d = mask_roi_flat.reshape((Nz, Nx))
    
    coords_roi_t = torch.tensor(coords_flat[mask_roi_flat], dtype=torch.float32, device=device)
    coords_roi_norm = normalize_roi_coords(coords_roi_t)
    
    with torch.no_grad():
        c_roi, alpha0_roi = pinn.forward_params(coords_roi_norm)
        c_grid = torch.ones((Nz, Nx), dtype=torch.float32, device=device) * 1540.0
        alpha0_grid = torch.ones((Nz, Nx), dtype=torch.float32, device=device) * (0.5 / 8.686e4)
        c_grid[mask_roi_2d] = c_roi.squeeze()
        alpha0_grid[mask_roi_2d] = alpha0_roi.squeeze()
        
        c_final_np = c_grid.cpu().numpy()
        alpha0_final_np = alpha0_grid.cpu().numpy() * 8.686e4
        
    c_out_path = f"articulo/c_map_{file_prefix}.npy"
    alpha_out_path = f"articulo/alpha_map_{file_prefix}.npy"
    np.save(c_out_path, c_final_np)
    np.save(alpha_out_path, alpha0_final_np)
    print(f"Saved {c_out_path}")
    print(f"Saved {alpha_out_path}")

if __name__ == "__main__":
    generate_maps("data_cubdl_prueba_1/MYO001.hdf5")
    generate_maps("data_cubdl_prueba_1/TSH002/TSH002.hdf5")
