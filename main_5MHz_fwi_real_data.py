import os
import time
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

# Set seeds for reproducibility
torch.manual_seed(42)
np.random.seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

# Import custom modules
from data_engine.loader import CUBDLLoader
from models.pinn_net import DecoupledUltrasoundPINN

# 1. Vectorized Helmholtz matrix assembly with Sommerfeld BC at all boundaries
# (Top, Left, Right, Bottom).
# Top boundary is modeled as Sommerfeld with outward normal pointing upwards (-z).
# dP/dn - i*k*P = 0 => -dP/dz - i*k*P = 0 => dP/dz + i*k*P = 0
# Discretized: (P_{1, i} - P_{0, i})/h + i*k_i*P_{0, i} = S_i => (1 - i*k_i*h)*P_{0, i} - P_{1, i} = -h*S_i
def assemble_helmholtz_matrix_sommerfeld_top(x, z, f, c_map, alpha0_map, Npml=15):
    import scipy.sparse as sp
    Nx = len(x)
    Nz = len(z)
    h = x[1] - x[0]
    omega = 2.0 * np.pi * f
    k = omega / c_map + 1j * (alpha0_map * f)
    N = Nz * Nx
    
    # Define PML parameters (Left, Right, Bottom) — dynamic Npml
    Lpml = Npml * h
    # Target reflection R = 1e-4
    sigma_max = (3.0 * 1540.0 / (2.0 * Lpml)) * np.log(1.0 / 1e-4)
    
    # 1D PML profiles
    sigma_x_1d = np.zeros(Nx)
    for i in range(Nx):
        if i < Npml:
            d = (Npml - i) * h
            sigma_x_1d[i] = sigma_max * (d / Lpml)**2
        elif i >= Nx - Npml:
            d = (i - (Nx - Npml - 1)) * h
            sigma_x_1d[i] = sigma_max * (d / Lpml)**2
            
    sigma_z_1d = np.zeros(Nz)
    for j in range(Nz):
        if j >= Nz - Npml:
            d = (j - (Nz - Npml - 1)) * h
            sigma_z_1d[j] = sigma_max * (d / Lpml)**2
            
    # Complex coordinate stretching profiles (positive sign for e^-iwt convention)
    gamma_x_1d = 1.0 + 1j * sigma_x_1d / omega
    gamma_z_1d = 1.0 + 1j * sigma_z_1d / omega
    
    # Mid-point values (averages)
    gamma_x_plus_1d = np.zeros(Nx, dtype=complex)
    gamma_x_minus_1d = np.zeros(Nx, dtype=complex)
    for i in range(Nx):
        ip = min(i + 1, Nx - 1)
        im = max(i - 1, 0)
        gamma_x_plus_1d[i] = 1.0 + 1j * 0.5 * (sigma_x_1d[i] + sigma_x_1d[ip]) / omega
        gamma_x_minus_1d[i] = 1.0 + 1j * 0.5 * (sigma_x_1d[i] + sigma_x_1d[im]) / omega
        
    gamma_z_plus_1d = np.zeros(Nz, dtype=complex)
    gamma_z_minus_1d = np.zeros(Nz, dtype=complex)
    for j in range(Nz):
        jp = min(j + 1, Nz - 1)
        jm = max(j - 1, 0)
        gamma_z_plus_1d[j] = 1.0 + 1j * 0.5 * (sigma_z_1d[j] + sigma_z_1d[jp]) / omega
        gamma_z_minus_1d[j] = 1.0 + 1j * 0.5 * (sigma_z_1d[j] + sigma_z_1d[jm]) / omega
        
    # Reconstruct 2D meshes
    gamma_x = np.tile(gamma_x_1d, (Nz, 1))
    gamma_x_plus = np.tile(gamma_x_plus_1d, (Nz, 1))
    gamma_x_minus = np.tile(gamma_x_minus_1d, (Nz, 1))
    
    gamma_z = np.tile(gamma_z_1d[:, np.newaxis], (1, Nx))
    gamma_z_plus = np.tile(gamma_z_plus_1d[:, np.newaxis], (1, Nx))
    gamma_z_minus = np.tile(gamma_z_minus_1d[:, np.newaxis], (1, Nx))
    
    J, I = np.meshgrid(np.arange(Nz), np.arange(Nx), indexing='ij')
    p = J * Nx + I
    
    mask_top = (J == 0)
    mask_left = (I == 0) & (~mask_top)
    mask_right = (I == Nx - 1) & (~mask_top)
    mask_bottom = (J == Nz - 1) & (~mask_top) & (~mask_left) & (~mask_right)
    mask_interior = (~mask_top) & (~mask_left) & (~mask_right) & (~mask_bottom)
    
    rows = []
    cols = []
    data = []
    
    # Top boundary (Sommerfeld: (1 - i * k_i * h) * P_0 - P_1 = -h * S)
    p_top = p[mask_top]
    rows.append(p_top)
    cols.append(p_top)
    data.append(1.0 - 1j * k[mask_top] * h)
    rows.append(p_top)
    cols.append(p_top + Nx)
    data.append(-np.ones_like(p_top, dtype=complex))
    
    # Left boundary
    p_left = p[mask_left]
    rows.append(p_left)
    cols.append(p_left)
    data.append(1.0 - 1j * k[mask_left] * h)
    rows.append(p_left)
    cols.append(p_left + 1)
    data.append(-np.ones_like(p_left, dtype=complex))
    
    # Right boundary
    p_right = p[mask_right]
    rows.append(p_right)
    cols.append(p_right - 1)
    data.append(-np.ones_like(p_right, dtype=complex))
    rows.append(p_right)
    cols.append(p_right)
    data.append(1.0 - 1j * k[mask_right] * h)
    
    # Bottom boundary
    p_bottom = p[mask_bottom]
    rows.append(p_bottom)
    cols.append(p_bottom - Nx)
    data.append(-np.ones_like(p_bottom, dtype=complex))
    rows.append(p_bottom)
    cols.append(p_bottom)
    data.append(1.0 - 1j * k[mask_bottom] * h)
    
    # Interior nodes (PML stencil)
    p_int = p[mask_interior]
    
    c_right = 1.0 / (gamma_x[mask_interior] * gamma_x_plus[mask_interior])
    c_left = 1.0 / (gamma_x[mask_interior] * gamma_x_minus[mask_interior])
    c_bottom = 1.0 / (gamma_z[mask_interior] * gamma_z_plus[mask_interior])
    c_top = 1.0 / (gamma_z[mask_interior] * gamma_z_minus[mask_interior])
    c_center = k[mask_interior]**2 * h**2 - (c_right + c_left + c_bottom + c_top)
    
    rows.append(p_int)
    cols.append(p_int - Nx)
    data.append(c_top)
    
    rows.append(p_int)
    cols.append(p_int - 1)
    data.append(c_left)
    
    rows.append(p_int)
    cols.append(p_int + 1)
    data.append(c_right)
    
    rows.append(p_int)
    cols.append(p_int + Nx)
    data.append(c_bottom)
    
    rows.append(p_int)
    cols.append(p_int)
    data.append(c_center)
    
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.concatenate(data)
    
    A = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    return A

# 2. Differentiable PyTorch Autograd Helmholtz Solver with Sommerfeld Top
class DifferentiableHelmholtzSolverSommerfeldTop(torch.autograd.Function):
    # Class-level timing storage for instrumentation
    _last_timing = {}
    
    @staticmethod
    def forward(ctx, c, alpha0, x, z, f, B_vec, Npml=15):
        from scipy.sparse.linalg import splu
        
        ctx.x = x
        ctx.z = z
        ctx.f = f
        ctx.B_vec = B_vec
        
        c_np = c.detach().cpu().numpy()
        alpha0_np = alpha0.detach().cpu().numpy()
        
        t0 = time.perf_counter()
        A = assemble_helmholtz_matrix_sommerfeld_top(x, z, f, c_np, alpha0_np, Npml)
        t1 = time.perf_counter()
        LU = splu(A.tocsc())      # Factorize once, reuse in backward
        t2 = time.perf_counter()
        P_vec = LU.solve(B_vec)   # Forward solve using cached factorization
        t3 = time.perf_counter()
        
        # Store timing breakdown and LU stats
        DifferentiableHelmholtzSolverSommerfeldTop._last_timing.update({
            'assembly_s': t1 - t0,
            'splu_s': t2 - t1,
            'fwd_solve_s': t3 - t2,
            'grid_shape': f"{len(z)}x{len(x)}",
            'nnz_A': A.nnz,
            'nnz_LU': LU.nnz,
            'fillin_ratio': LU.nnz / A.nnz,
        })
        
        ctx.save_for_backward(c, alpha0)
        ctx.P_vec = P_vec
        ctx.LU = LU               # Cache factorization for adjoint solve
        ctx.A = A                  # Keep A for dA/dc gradient computation
        
        P = P_vec.reshape((len(z), len(x)))
        dtype_out = torch.complex128 if c.dtype == torch.float64 else torch.complex64
        return torch.tensor(P, dtype=dtype_out, device=c.device)
        
    @staticmethod
    def backward(ctx, grad_output):
        c, alpha0 = ctx.saved_tensors
        P_vec = ctx.P_vec
        LU = ctx.LU
        A = ctx.A
        x = ctx.x
        z = ctx.z
        f = ctx.f
        
        Nx = len(x)
        Nz = len(z)
        h = x[1] - x[0]
        omega = 2.0 * np.pi * f
        
        c_np = c.detach().cpu().numpy()
        alpha0_np = alpha0.detach().cpu().numpy()
        
        grad_output_vec = grad_output.detach().cpu().numpy().flatten()
        rhs_adjoint = grad_output_vec.conj()
        
        # Solve adjoint equation A^T * lambda = rhs using cached LU factorization
        # trans='T' = simple transpose (NOT 'H' = conjugate transpose)
        t_bwd_0 = time.perf_counter()
        lambda_vec = LU.solve(rhs_adjoint, trans='T')
        t_bwd_1 = time.perf_counter()
        
        J, I = np.meshgrid(np.arange(Nz), np.arange(Nx), indexing='ij')
        mask_top = (J == 0)
        mask_left = (I == 0) & (~mask_top)
        mask_right = (I == Nx - 1) & (~mask_top)
        mask_bottom = (J == Nz - 1) & (~mask_top) & (~mask_left) & (~mask_right)
        mask_interior = (~mask_top) & (~mask_left) & (~mask_right) & (~mask_bottom)
        
        c_flat = c_np.flatten()
        alpha0_flat = alpha0_np.flatten()
        
        # Derivatives of k^2 (interior)
        dk2_dc = -2.0 * omega**2 / (c_flat**3) - 2j * omega * f * alpha0_flat / (c_flat**2)
        dk2_dalpha0 = -2.0 * alpha0_flat * (f**2) + 2j * omega * f / c_flat
        
        # Derivatives of k (Sommerfeld boundaries)
        dk_dc = -omega / (c_flat**2)
        dk_dalpha0 = np.full_like(c_flat, 1j * f, dtype=complex)
        
        dApp_dc = np.zeros_like(c_flat, dtype=complex)
        dApp_dalpha0 = np.zeros_like(alpha0_flat, dtype=complex)
        
        # Interior nodes
        flat_mask_interior = mask_interior.flatten()
        dApp_dc[flat_mask_interior] = (h**2) * dk2_dc[flat_mask_interior]
        dApp_dalpha0[flat_mask_interior] = (h**2) * dk2_dalpha0[flat_mask_interior]
        
        # Top boundary
        flat_mask_top = mask_top.flatten()
        dApp_dc[flat_mask_top] = -1j * h * dk_dc[flat_mask_top]
        dApp_dalpha0[flat_mask_top] = -1j * h * dk_dalpha0[flat_mask_top]
        
        # Left boundary
        flat_mask_left = mask_left.flatten()
        dApp_dc[flat_mask_left] = -1j * h * dk_dc[flat_mask_left]
        dApp_dalpha0[flat_mask_left] = -1j * h * dk_dalpha0[flat_mask_left]
        
        # Right boundary
        flat_mask_right = mask_right.flatten()
        dApp_dc[flat_mask_right] = -1j * h * dk_dc[flat_mask_right]
        dApp_dalpha0[flat_mask_right] = -1j * h * dk_dalpha0[flat_mask_right]
        
        # Bottom boundary
        flat_mask_bottom = mask_bottom.flatten()
        dApp_dc[flat_mask_bottom] = -1j * h * dk_dc[flat_mask_bottom]
        dApp_dalpha0[flat_mask_bottom] = -1j * h * dk_dalpha0[flat_mask_bottom]
        
        # Gradient formula: grad = -Re( lambda * dA/dtheta * P )
        grad_c_vec = -np.real(lambda_vec * dApp_dc * P_vec)
        grad_alpha0_vec = -np.real(lambda_vec * dApp_dalpha0 * P_vec)
        
        grad_c_t = torch.tensor(grad_c_vec.reshape((Nz, Nx)), dtype=c.dtype, device=c.device)
        grad_alpha0_t = torch.tensor(grad_alpha0_vec.reshape((Nz, Nx)), dtype=alpha0.dtype, device=alpha0.device)
        t_bwd_2 = time.perf_counter()
        
        # Store backward timing
        DifferentiableHelmholtzSolverSommerfeldTop._last_timing.update({
            'bwd_adjoint_solve_s': t_bwd_1 - t_bwd_0,
            'bwd_grad_compute_s': t_bwd_2 - t_bwd_1,
        })
        
        return grad_c_t, grad_alpha0_t, None, None, None, None, None

# 3. Delay-and-Sum (DAS) Beamformer
def delay_and_sum(loader, x_grid, z_grid, transmit_index=0, compound=False):
    """
    Standard or Compounded Delay-and-Sum (DAS) Beamformer to reconstruct B-mode envelope image.
    Used to locate the inclusion/tumor center automatically.
    """
    if compound:
        print("Running Compounded Delay-and-Sum (DAS) beamforming to localize tumor...")
    else:
        print("Running Delay-and-Sum (DAS) beamforming to localize tumor...")
        
    fs = loader.fs
    c0 = loader.c0
    t0 = loader.t0
    rx_x = loader.rx_positions[:, 0]
    
    # 2D Meshgrid for reconstruction
    X, Z = np.meshgrid(x_grid, z_grid, indexing='ij')  # shape [Nx, Nz]
    Nx, Nz = X.shape
    
    from scipy.signal import hilbert
    
    if compound:
        compounded_image = np.zeros((Nx, Nz), dtype=complex)
        num_transmits = loader.rf.shape[2] if len(loader.rf.shape) == 3 else 1
        for tx_idx in range(num_transmits):
            rf = loader.rf[:, :, tx_idx] if len(loader.rf.shape) == 3 else loader.rf
            tx_delay = loader.tx_delays[tx_idx, :] if len(loader.tx_delays.shape) == 2 else loader.tx_delays
            
            # Estimate transmit source center (element with minimum delay)
            tx_center_idx = np.argmin(tx_delay)
            x_tx = rx_x[tx_center_idx]
            
            # Distance from transmitter center to pixel
            dist_tx = np.sqrt((X - x_tx)**2 + Z**2)
            t_tx = dist_tx / c0
            
            # Hilbert transform for analytical complex signal (envelope detection)
            rf_complex = hilbert(rf, axis=0)
            
            # Sum backscattered signals from all receivers
            image_tx = np.zeros((Nx, Nz), dtype=complex)
            for r_idx, xr in enumerate(rx_x):
                # Distance from pixel to receiver
                dist_rx = np.sqrt((X - xr)**2 + Z**2)
                t_rx = dist_rx / c0
                
                # Total time of flight
                t_total = t_tx + t_rx + tx_delay[r_idx]
                
                # Interpolate RF signal at travel time
                sample_idx = (t_total - t0) * fs
                flat_idx = sample_idx.flatten()
                flat_val = np.interp(flat_idx, np.arange(len(rf_complex)), rf_complex[:, r_idx], left=0, right=0)
                image_tx += flat_val.reshape((Nx, Nz))
            compounded_image += image_tx
            
        # Transpose back to [Nz, Nx] to align with standard coordinates
        return np.abs(compounded_image).T
    else:
        rf = loader.rf[:, :, transmit_index] if len(loader.rf.shape) == 3 else loader.rf
        tx_delay = loader.tx_delays[transmit_index, :] if len(loader.tx_delays.shape) == 2 else loader.tx_delays
        tx_center_idx = np.argmin(tx_delay)
        x_tx = rx_x[tx_center_idx]
        
        # Distance from transmitter center to pixel
        dist_tx = np.sqrt((X - x_tx)**2 + Z**2)
        t_tx = dist_tx / c0
        
        # Hilbert transform for analytical complex signal (envelope detection)
        rf_complex = hilbert(rf, axis=0)
        
        image = np.zeros((Nx, Nz), dtype=complex)
        # Sum backscattered signals from all receivers
        for r_idx, xr in enumerate(rx_x):
            # Distance from pixel to receiver
            dist_rx = np.sqrt((X - xr)**2 + Z**2)
            t_rx = dist_rx / c0
            
            # Total time of flight
            t_total = t_tx + t_rx + tx_delay[r_idx]
            
            # Interpolate RF signal at travel time
            sample_idx = (t_total - t0) * fs
            flat_idx = sample_idx.flatten()
            flat_val = np.interp(flat_idx, np.arange(len(rf_complex)), rf_complex[:, r_idx], left=0, right=0)
            image += flat_val.reshape((Nx, Nz))
            
        # Transpose back to [Nz, Nx] to align with standard coordinates
        return np.abs(image).T

# Self-contained SSIM implementation
def calculate_ssim(im1, im2):
    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1_sq = im1.var()
    sigma2_sq = im2.var()
    sigma12 = np.cov(im1.flatten(), im2.flatten())[0, 1]
    L = max(im1.max(), im2.max()) - min(im1.min(), im2.min())
    if L == 0: L = 1.0
    C1 = (0.01 * L) ** 2
    C2 = (0.03 * L) ** 2
    return ((2 * mu1 * mu2 + C1) * (2 * sigma12 + C2)) / ((mu1**2 + mu2**2 + C1) * (sigma1_sq + sigma2_sq + C2) + 1e-10)

def tv_loss_alpha(alpha):
    diff_x = torch.abs(alpha[:, 1:] - alpha[:, :-1])
    diff_z = torch.abs(alpha[1:, :] - alpha[:-1, :])
    return diff_x.mean() + diff_z.mean()

def run_fwi_real_data(file_path="data_engine/cubdl_sample_dataset.mat", seed=42, output_suffix="", npml=15, roi_taper=False, taper_cells=5, lr_decay=False):
    # Set seeds for reproducibility
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        
    file_prefix = os.path.splitext(os.path.basename(file_path))[0] + output_suffix
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n==================================================")
    print(f"CUBDL FWI Pipeline: Running on {file_path}")
    print(f"Device: {device}")
    print(f"==================================================")
    
    if not os.path.exists(file_path):
        print(f"Archivo {file_path} no encontrado. Ejecutando la prueba de sanity para generar el mock dataset...")
        from data_engine.sanity_check_cubdl import run_sanity_check
        run_sanity_check()
        
    # 1. Load CUBDL data
    loader = CUBDLLoader(file_path)
    print(f"Dataset cargado. Frecuencia de muestreo (fs): {loader.fs/1e6:.1f} MHz")
    print(f"Velocidad de fondo (c0): {loader.c0:.1f} m/s")
    print(f"RF data shape: {loader.rf.shape}")
    
    # Determine central transmit angle index
    tx_idx = 0
    if loader.angles is not None:
        tx_idx = int(np.argmin(np.abs(loader.angles)))
        print(f"Using central transmit index {tx_idx} (angle: {loader.angles[tx_idx]:.4f} rad) for the pipeline.")
    else:
        print(f"Using default transmit index {tx_idx} for the pipeline.")

    # 2. Reconstruct B-mode image using DAS
    # Define a grid for beamforming (lateral -2 cm to 2 cm, axial 0.2 cm to 3.5 cm)
    x_bmode = np.arange(-0.02, 0.02 + 0.5e-3, 0.5e-3)
    z_bmode = np.arange(0.002, 0.035 + 0.5e-3, 0.5e-3)
    
    bmode_img = delay_and_sum(loader, x_bmode, z_bmode, compound=True)
    
    # Normalize and log compress for display
    bmode_normalized = bmode_img / np.max(bmode_img)
    bmode_log = 20 * np.log10(bmode_normalized + 1e-5)
    
    # Auto-localize tumor center: find max envelope below 1.0 cm to avoid near-field clutter
    depth_mask = z_bmode > 0.010
    img_valid = bmode_img[depth_mask, :]
    max_idx = np.unravel_index(np.argmax(img_valid), img_valid.shape)
    z_tumor = z_bmode[depth_mask][max_idx[0]]
    x_tumor = x_bmode[max_idx[1]]
    
    print(f"--> Localización del tumor por Beamforming B-mode: x = {x_tumor*1e3:.2f} mm, z = {z_tumor*1e3:.2f} mm")
    
    # Define ROI (1.5 cm x 1.5 cm) centered on the localized tumor
    roi_half_w = 0.0075
    roi_half_h = 0.0075
    
    x_min_roi, x_max_roi = x_tumor - roi_half_w, x_tumor + roi_half_w
    z_min_roi, z_max_roi = z_tumor - roi_half_h, z_tumor + roi_half_h
    
    def normalize_roi_coords(coords):
        # Normalize relative to tumor center and ROI width/height
        x_norm = (coords[:, 0:1] - x_tumor) / roi_half_w
        z_norm = (coords[:, 1:2] - z_tumor) / roi_half_h
        return torch.cat([x_norm, z_norm], dim=1)
        
    pinn = DecoupledUltrasoundPINN(
        spatial_dim=2, hidden_features=128, num_layers=4,
        omega_0_params=5.0
    ).to(device)
    
    # Initialize biases for homogeneous background: c = 1540 m/s, alpha0 = 0.5 dB/MHz/cm
    with torch.no_grad():
        pinn.param_out.bias[0] = -0.4200  # maps to c = 1540.0
        pinn.param_out.bias[1] = -1.3863  # maps to alpha0 = 0.5
        pinn.param_out.weight.fill_(0.0)
        
    # 5. Multiscale FWI loop
    f_stages = [1.0e6, 2.0e6, 3.5e6]
    stage_epochs = [30, 30, 80]
    lr_stages = [1.5e-3, 1.2e-3, 1.2e-3]
    dx_stages = [1.0e-4, 7.5e-5, 4.5e-5]  # Adaptive dx to satisfy stability criteria kh < pi/4
    
    optimizer = optim.Adam([
        {'params': pinn.parameters(), 'lr': lr_stages[0]},
    ], lr=lr_stages[0])
    
    # Añadir prior físico de α por institución
    institution = os.path.basename(file_path)[:3].upper()  # 'MYO', 'OSL', etc.
    alpha_ref_db = {'OSL': 0.5, 'MYO': 0.7, 'EUT': 0.6, 'INS': 0.5}.get(institution, 0.5)
    alpha_ref = torch.tensor(alpha_ref_db / 8.686e4, dtype=torch.float32, device=device)
    
    print("\nIniciando optimización FWI diferenciable por etapas de frecuencia (Marching)...")
    start_time = time.time()
    
    progress_history = []
    global_epoch = 0
    
    for stage_idx, f in enumerate(f_stages):
        lr_stage = lr_stages[stage_idx]
        epochs_stage = stage_epochs[stage_idx]
        dx = dx_stages[stage_idx]
        dz = dx
        omega = 2.0 * np.pi * f
        
        # 3. Setup FWI Grid for this stage
        # Cover transducer range laterally, and up to 3.2 cm axially
        rx_x = loader.rx_positions[:, 0]
        x_grid = np.arange(min(rx_x), max(rx_x) + dx, dx)
        z_grid = np.arange(0.0, 0.032 + dz, dz)
        Nx, Nz = len(x_grid), len(z_grid)
        h = x_grid[1] - x_grid[0]
        
        # Assert stability condition: kh < pi/4 (using minimum sound speed c0 = 1400.0)
        k_max = omega / 1400.0
        kh = k_max * h
        print(f"\n--- Etapa {stage_idx+1}/{len(f_stages)}: Frecuencia {f/1e6:.1f} MHz (Grid: {Nz}x{Nx}, dx: {dx*1e6:.1f} um, kh: {kh:.4f}, LR: {lr_stage:.2e}) ---")
        assert kh < np.pi / 4, f"Inestabilidad Helmholtz detectada: kh = {kh:.4f} >= pi/4 ({np.pi/4:.4f}) a f = {f/1e6:.1f} MHz y dx = {dx*1e6:.1f} um"
        
        # Map transducer receiver coordinates to grid indices
        rx_indices = []
        for rx in rx_x:
            idx = np.argmin(np.abs(x_grid - rx))
            rx_indices.append(idx)
        rx_indices = np.array(rx_indices)
        
        # Coordinates for meshgrid
        X, Z = np.meshgrid(x_grid, z_grid)
        coords_flat = np.stack([X.flatten(), Z.flatten()], axis=1)
        
        # Mask of points inside the ROI
        mask_roi_flat = (coords_flat[:, 0] >= x_min_roi) & (coords_flat[:, 0] <= x_max_roi) & \
                        (coords_flat[:, 1] >= z_min_roi) & (coords_flat[:, 1] <= z_max_roi)
        mask_roi_2d = mask_roi_flat.reshape((Nz, Nx))
        
        # Neural parameterizer coordinates for the new grid
        coords_roi_t = torch.tensor(coords_flat[mask_roi_flat], dtype=torch.float32, device=device)
        coords_roi_norm = normalize_roi_coords(coords_roi_t)
        
        # Define 2D weight matrix W for ROI transition
        if roi_taper:
            i_indices = np.where((x_grid >= x_min_roi) & (x_grid <= x_max_roi))[0]
            j_indices = np.where((z_grid >= z_min_roi) & (z_grid <= z_max_roi))[0]
            if len(i_indices) > 0 and len(j_indices) > 0:
                i_min, i_max = i_indices[0], i_indices[-1]
                j_min, j_max = j_indices[0], j_indices[-1]
                
                weight_x = np.zeros(Nx)
                weight_z = np.zeros(Nz)
                
                for i in range(Nx):
                    if i_min <= i <= i_max:
                        dist_x = min(i - i_min, i_max - i)
                        if dist_x >= taper_cells:
                            weight_x[i] = 1.0
                        else:
                            weight_x[i] = 0.5 * (1.0 - np.cos(np.pi * dist_x / taper_cells))
                for j in range(Nz):
                    if j_min <= j <= j_max:
                        dist_z = min(j - j_min, j_max - j)
                        if dist_z >= taper_cells:
                            weight_z[j] = 1.0
                        else:
                            weight_z[j] = 0.5 * (1.0 - np.cos(np.pi * dist_z / taper_cells))
                W_np = np.outer(weight_z, weight_x)
            else:
                W_np = mask_roi_2d.astype(np.float32)
            W = torch.tensor(W_np, dtype=torch.float32, device=device)
        else:
            W = torch.tensor(mask_roi_2d, dtype=torch.float32, device=device)
            
        # Reconstruct GT on the new grid
        dist_sq = (X - x_tumor)**2 + (Z - z_tumor)**2
        c_gt = np.ones((Nz, Nx)) * 1540.0
        c_gt[dist_sq <= 0.003**2] = 1620.0  # +80 m/s tumor speed
        
        alpha0_gt = np.ones((Nz, Nx)) * (0.5 / 8.686e4)
        alpha0_gt[dist_sq <= 0.003**2] = 1.2 / 8.686e4 # high attenuation in tumor
        alpha0_gt_db = alpha0_gt * 8.686e4
        
        # Resetear estado del optimizador sin perder los pesos del PINN (borra el momentum anterior)
        optimizer = optim.Adam([
            {'params': pinn.parameters(), 'lr': lr_stage}
        ], lr=lr_stage)
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.95) if lr_decay else None
            
        # 5a. Extract boundary measured data from loader
        _, P_meas_dict = loader.extract_boundary_pressure([f], transmit_index=tx_idx)
        P_meas = P_meas_dict[f]  # Complex array [128]
        P_meas_t = torch.tensor(P_meas, dtype=torch.complex64, device=device)
        
        # 5b. Construct the source term B_vec at the top boundary
        tx_delay = loader.tx_delays[tx_idx, :] if len(loader.tx_delays.shape) == 2 else loader.tx_delays
        
        # Source term: S = exp(-1j * omega * tx_delay)
        S = np.exp(-1j * omega * tx_delay)
        
        # Map source to grid B_vec (Sommerfeld top is (1 - i*k*h)*P0 - P1 = -h*S)
        B_vec = np.zeros(Nz * Nx, dtype=complex)
        B_vec[rx_indices] = -h * S
        
        # 5c. Precompute Incident Field in background medium (homogenous)
        # c_bg = 1540, alpha0_bg = 0.5
        c_bg_np = np.ones((Nz, Nx)) * 1540.0
        alpha0_bg_np = np.ones((Nz, Nx)) * (0.5 / 8.686e4)
        
        A_bg = assemble_helmholtz_matrix_sommerfeld_top(x_grid, z_grid, f, c_bg_np, alpha0_bg_np, npml)
        import scipy.sparse.linalg as spla
        P_inc_vec = spla.spsolve(A_bg, B_vec)
        P_inc_rec = P_inc_vec.reshape((Nz, Nx))[0, rx_indices]
        P_inc_rec_t = torch.tensor(P_inc_rec, dtype=torch.complex64, device=device)
        
        for epoch in range(1, epochs_stage + 1):
            global_epoch += 1
            t_epoch_start = time.perf_counter()
            optimizer.zero_grad()
            
            # Predict parameters in the ROI
            c_roi, alpha0_roi = pinn.forward_params(coords_roi_norm)
            
            # Reconstruct the full 2D map: ROI values + homogeneous background elsewhere
            c_pinn_grid = torch.ones((Nz, Nx), dtype=torch.float32, device=device) * 1540.0
            alpha0_pinn_grid = torch.ones((Nz, Nx), dtype=torch.float32, device=device) * (0.5 / 8.686e4)
            
            # Set the ROI indices
            c_pinn_grid[mask_roi_2d] = c_roi.squeeze()
            alpha0_pinn_grid[mask_roi_2d] = alpha0_roi.squeeze()
            
            # Apply transition taper weighting
            c_grid = 1540.0 + W * (c_pinn_grid - 1540.0)
            alpha0_grid = (0.5 / 8.686e4) + W * (alpha0_pinn_grid - (0.5 / 8.686e4))
            
            # Solve Helmholtz wave equation
            P_pred = DifferentiableHelmholtzSolverSommerfeldTop.apply(
                c_grid, alpha0_grid, x_grid, z_grid, f, B_vec, npml
            )
            
            # Extract predicted pressure at receivers (top row j=0, indices rx_indices)
            P_pred_rec = P_pred[0, rx_indices]
            
            # Calculate scattered field: P_scat = P_pred - P_inc
            P_scat_rec = P_pred_rec - P_inc_rec_t
            
            # Calculate complex scaling factor beta to align amplitudes and phases
            num = torch.sum(P_meas_t * torch.conj(P_scat_rec))
            den = torch.sum(torch.abs(P_scat_rec)**2) + 1e-12
            beta_raw = num / den
            
            # Restricción de fase en beta para evitar la inversión de contraste (ambigüedad fase-velocidad)
            beta_phase = torch.angle(beta_raw)
            beta_phase_clipped = torch.clamp(beta_phase, -np.pi/4, np.pi/4)
            beta = torch.abs(beta_raw) * torch.exp(1j * beta_phase_clipped)
            
            # Calibrated prediction
            P_scat_calib = beta * P_scat_rec
            
            # L2 Loss on the boundary receivers
            loss_data = torch.mean(torch.abs(P_meas_t - P_scat_calib)**2)
            
            # Loss compuesta con prior, regularización TV y penalización de fase
            alpha_roi_mean = alpha0_roi.mean()
            loss_alpha_prior = (alpha_roi_mean - alpha_ref) ** 2
            loss_tv = tv_loss_alpha(alpha0_grid)
            loss_phase_reg = 0.1 * beta_phase ** 2
            
            loss = loss_data + loss_phase_reg + 0.1 * loss_alpha_prior + 1e-3 * loss_tv
            
            loss.backward()
            optimizer.step()
            if scheduler is not None:
                scheduler.step()
            t_epoch_end = time.perf_counter()
            
            # Print detailed timing breakdown for first 3 epochs of each stage
            if epoch <= 3:
                tm = DifferentiableHelmholtzSolverSommerfeldTop._last_timing
                print(f"  [TIMING] Stage {stage_idx+1} Epoch {epoch}/{epochs_stage} | "
                      f"Grid: {tm.get('grid_shape','?')} | "
                      f"Assembly: {tm.get('assembly_s',0)*1e3:.1f}ms | "
                      f"spLU: {tm.get('splu_s',0)*1e3:.1f}ms | "
                      f"Fwd solve: {tm.get('fwd_solve_s',0)*1e3:.1f}ms | "
                      f"Bwd adjoint: {tm.get('bwd_adjoint_solve_s',0)*1e3:.1f}ms | "
                      f"Bwd grad: {tm.get('bwd_grad_compute_s',0)*1e3:.1f}ms | "
                      f"Epoch total: {(t_epoch_end - t_epoch_start)*1e3:.1f}ms | "
                      f"nnz(A): {tm.get('nnz_A',0):,} | nnz(LU): {tm.get('nnz_LU',0):,} | "
                      f"Fill-in: {tm.get('fillin_ratio',0):.2f}x")
            
            if epoch % 10 == 0 or epoch == 1:
                with torch.no_grad():
                    c_np = c_grid.cpu().numpy()
                    alpha0_np = alpha0_grid.cpu().numpy() * 8.686e4
                    alpha0_gt_db = alpha0_gt * 8.686e4
                    
                    err_c = np.linalg.norm(c_np[mask_roi_2d] - c_gt[mask_roi_2d]) / np.linalg.norm(c_gt[mask_roi_2d])
                    err_alpha0 = np.linalg.norm(alpha0_np[mask_roi_2d] - alpha0_gt_db[mask_roi_2d]) / np.linalg.norm(alpha0_gt_db[mask_roi_2d])
                    
                    ssim_c = calculate_ssim(c_np, c_gt)
                    ssim_alpha0 = calculate_ssim(alpha0_np, alpha0_gt_db)
                    
                    print(f"Epoca {epoch:2d}/{epochs_stage:2d} | Loss: {loss.item():.4e} | "
                          f"Err c (ROI): {err_c*100:.2f}%, Err alpha0 (ROI): {err_alpha0*100:.2f}% | "
                          f"SSIM c: {ssim_c:.3f}, SSIM alpha0: {ssim_alpha0:.3f}")
                    
                    progress_history.append({
                        'stage': stage_idx + 1,
                        'frequency': f,
                        'epoch': global_epoch,
                        'loss': loss.item(),
                        'err_c': err_c,
                        'err_alpha0': err_alpha0,
                        'ssim_c': ssim_c,
                        'ssim_alpha0': ssim_alpha0
                    })
                    
    elapsed_time = time.time() - start_time
    print(f"\nInversión FWI en datos reales completada en {elapsed_time:.1f} segundos.")
    
    # 6. Save and plot results
    with torch.no_grad():
        c_roi, alpha0_roi = pinn.forward_params(coords_roi_norm)
        c_pinn_grid = torch.ones((Nz, Nx), dtype=torch.float32, device=device) * 1540.0
        alpha0_pinn_grid = torch.ones((Nz, Nx), dtype=torch.float32, device=device) * (0.5 / 8.686e4)
        c_pinn_grid[mask_roi_2d] = c_roi.squeeze()
        alpha0_pinn_grid[mask_roi_2d] = alpha0_roi.squeeze()
        
        c_grid = 1540.0 + W * (c_pinn_grid - 1540.0)
        alpha0_grid = (0.5 / 8.686e4) + W * (alpha0_pinn_grid - (0.5 / 8.686e4))
        
        c_final_np = c_grid.cpu().numpy()
        alpha0_final_np = alpha0_grid.cpu().numpy() * 8.686e4
        alpha0_gt_db = alpha0_gt * 8.686e4
        
        # Save optimized model
        torch.save(pinn.state_dict(), f"articulo/best_pinn_5MHz_fwi_real_{file_prefix}.pt")
        
        # Guardar arrays numpy para CNR
        np.save(f"articulo/c_map_{file_prefix}.npy", c_final_np)
        np.save(f"articulo/alpha_map_{file_prefix}.npy", alpha0_final_np)
        
        # Solve wave field at the last stage frequency for plotting
        P_final_np = DifferentiableHelmholtzSolverSommerfeldTop.apply(
            c_grid, alpha0_grid, x_grid, z_grid, f_stages[-1], B_vec, npml
        ).cpu().numpy()
        
    # Plot results
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # 1. B-mode Image showing auto-localization
    im = axes[0, 0].imshow(bmode_log, extent=[x_bmode[0]*100, x_bmode[-1]*100, z_bmode[-1]*100, z_bmode[0]*100], cmap='gray', aspect='auto')
    axes[0, 0].axvline(x_tumor*100, color='r', linestyle='--', label='Centro Tumor')
    axes[0, 0].axhline(z_tumor*100, color='r', linestyle='--')
    # Draw ROI box
    roi_box = plt.Rectangle((x_min_roi*100, z_min_roi*100), roi_half_w*200, roi_half_h*200, edgecolor='cyan', facecolor='none', linewidth=2, label='ROI FWI')
    axes[0, 0].add_patch(roi_box)
    axes[0, 0].set_title("Beamforming DAS (B-mode)")
    axes[0, 0].set_xlabel("Lateral (cm)")
    axes[0, 0].set_ylabel("Axial (cm)")
    axes[0, 0].legend()
    fig.colorbar(im, ax=axes[0, 0], label='Decibelios (dB)')
    
    # 2. Velocity GT
    im = axes[0, 1].imshow(c_gt, extent=[x_grid[0]*100, x_grid[-1]*100, z_grid[-1]*100, z_grid[0]*100], cmap='jet', aspect='auto')
    axes[0, 1].set_title("c GT (m/s)")
    axes[0, 1].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[0, 1])
    
    # 3. Velocity FWI predicted
    im = axes[0, 2].imshow(c_final_np, extent=[x_grid[0]*100, x_grid[-1]*100, z_grid[-1]*100, z_grid[0]*100], cmap='jet', aspect='auto')
    # Draw ROI box
    roi_box = plt.Rectangle((x_min_roi*100, z_min_roi*100), roi_half_w*200, roi_half_h*200, edgecolor='cyan', facecolor='none', linewidth=1.5)
    axes[0, 2].add_patch(roi_box)
    axes[0, 2].set_title(f"FWI c (ROI L2: {progress_history[-1]['err_c']*100:.2f}%)")
    axes[0, 2].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[0, 2])
    
    # 4. Convergence
    epochs_list = [h['epoch'] for h in progress_history]
    losses = [h['loss'] for h in progress_history]
    axes[1, 0].semilogy(epochs_list, losses, 'k-o')
    axes[1, 0].set_title("Convergencia FWI en Datos de Canal")
    axes[1, 0].set_xlabel("Época")
    axes[1, 0].set_ylabel("Misfit Receivers (MSE)")
    axes[1, 0].grid(True)
    
    # 5. Attenuation GT
    im = axes[1, 1].imshow(alpha0_gt_db, extent=[x_grid[0]*100, x_grid[-1]*100, z_grid[-1]*100, z_grid[0]*100], cmap='hot', aspect='auto')
    axes[1, 1].set_title("alpha_0 GT (dB/MHz/cm)")
    axes[1, 1].set_xlabel("Lateral (cm)")
    axes[1, 1].set_ylabel("Axial (cm)")
    fig.colorbar(im, ax=axes[1, 1])
    
    # 6. Attenuation FWI predicted
    im = axes[1, 2].imshow(alpha0_final_np, extent=[x_grid[0]*100, x_grid[-1]*100, z_grid[-1]*100, z_grid[0]*100], cmap='hot', aspect='auto')
    # Draw ROI box
    roi_box = plt.Rectangle((x_min_roi*100, z_min_roi*100), roi_half_w*200, roi_half_h*200, edgecolor='cyan', facecolor='none', linewidth=1.5)
    axes[1, 2].add_patch(roi_box)
    axes[1, 2].set_title(f"FWI alpha_0 (ROI L2: {progress_history[-1]['err_alpha0']*100:.2f}%)")
    axes[1, 2].set_xlabel("Lateral (cm)")
    fig.colorbar(im, ax=axes[1, 2])
    
    plt.tight_layout()
    plt.savefig(f'articulo/resultados_fwi_real_cubdl_{file_prefix}.png', dpi=150)
    plt.close()
    
    print(f"Resultados guardados en articulo/resultados_fwi_real_cubdl_{file_prefix}.png")
    
    # Build a Markdown report for Tarea C
    best_h = progress_history[-1]
    table_rows = ""
    for h in progress_history:
        table_rows += (f"| {h['epoch']} | {h['stage']} | {h['frequency']/1e6:.1f} MHz | {h['loss']:.4e} | "
                      f"{h['err_c']*100:.2f}% | {h['err_alpha0']*100:.2f}% | {h['ssim_c']:.3f} | {h['ssim_alpha0']:.3f} |\n")
                      
    report_content = f"""# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: {time.strftime("%Y-%m-%d %H:%M:%S")}

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{{\\text{{inc}}}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{{\\text{{scat}}}} = P_{{\\text{{pred}}}} - P_{{\\text{{inc}}}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = {x_tumor*1e3:.2f} mm, z = {z_tumor*1e3:.2f} mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: {best_h['err_c']*100:.2f}%
*   **Error L2 alpha0 final en ROI**: {best_h['err_alpha0']*100:.2f}%
*   **SSIM Final c**: {best_h['ssim_c']:.3f}
*   **SSIM Final alpha0**: {best_h['ssim_alpha0']:.3f}

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_{file_prefix}.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_{file_prefix}.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
{table_rows}
"""
    with open(f"articulo/metricas_progreso_fwi_real_cubdl_{file_prefix}.md", "w", encoding="utf-8") as f:
        f.write(report_content)
    print(f"Reporte guardado en articulo/metricas_progreso_fwi_real_cubdl_{file_prefix}.md")
    
    # Escribir sidecars JSON para la webapp automáticamente al finalizar
    import json
    
    # Parámetros para calcular CNR y peak distances en la ROI activa
    roi_half = 0.0075  # 7.5 mm en metros
    MARGIN_FRAC = 0.15  # 15% de margen
    
    try:
        # Meshgrid para la grilla final
        X_grid_2d, Z_grid_2d = np.meshgrid(x_grid, z_grid)
        
        # Extraer sub-array ROI 2D
        mask_x = (x_grid >= x_tumor - roi_half) & (x_grid <= x_tumor + roi_half)
        mask_z = (z_grid >= z_tumor - roi_half) & (z_grid <= z_tumor + roi_half)
        z_idx = np.where(mask_z)[0]
        x_idx = np.where(mask_x)[0]
        
        z0, z1 = z_idx.min(), z_idx.max() + 1
        x0, x1 = x_idx.min(), x_idx.max() + 1
        
        X2D_roi = X_grid_2d[z0:z1, x0:x1]
        Z2D_roi = Z_grid_2d[z0:z1, x0:x1]
        
        roi_c = c_final_np[z0:z1, x0:x1]
        roi_alpha = alpha0_final_np[z0:z1, x0:x1]
        
        Nz_roi, Nx_roi = roi_c.shape
        mz = int(MARGIN_FRAC * Nz_roi)
        mx = int(MARGIN_FRAC * Nx_roi)
        
        # Resguardo para grillas muy pequeñas
        if mz == 0: mz = 1
        if mx == 0: mx = 1
        
        interior_mask = np.zeros_like(roi_c, dtype=bool)
        interior_mask[mz:-mz, mx:-mx] = True
        
        # CNR c robusto con margen
        vals_roi_c = roi_c[interior_mask]
        p90_c = np.percentile(vals_roi_c, 90)
        p50_c = np.percentile(vals_roi_c, 50)
        bg_vals_c = vals_roi_c[vals_roi_c < p90_c]
        sigma_c = bg_vals_c.std() if bg_vals_c.size > 0 else vals_roi_c.std()
        cnr_c_linear = abs(p90_c - p50_c) / (sigma_c + 1e-12)
        cnr_c_db = float(20 * np.log10(cnr_c_linear + 1e-12))
        
        # CNR alpha robusto con margen
        vals_roi_a = roi_alpha[interior_mask]
        p90_a = np.percentile(vals_roi_a, 90)
        p50_a = np.percentile(vals_roi_a, 50)
        bg_vals_a = vals_roi_a[vals_roi_a < p90_a]
        sigma_a = bg_vals_a.std() if bg_vals_a.size > 0 else vals_roi_a.std()
        cnr_a_linear = abs(p90_a - p50_a) / (sigma_a + 1e-12)
        cnr_a_db = float(20 * np.log10(cnr_a_linear + 1e-12))
        
        # Distancia del pico para c
        roi_c_int = np.where(interior_mask, roi_c, -1e9)
        max_c_idx = np.unravel_index(np.argmax(roi_c_int), roi_c_int.shape)
        dist_c = float(np.sqrt((X2D_roi[max_c_idx] - x_tumor)**2 + (Z2D_roi[max_c_idx] - z_tumor)**2) * 1e3)
        
        # Distancia del pico para alpha
        roi_alpha_int = np.where(interior_mask, roi_alpha, -1e9)
        max_alpha_idx = np.unravel_index(np.argmax(roi_alpha_int), roi_alpha_int.shape)
        dist_alpha = float(np.sqrt((X2D_roi[max_alpha_idx] - x_tumor)**2 + (Z2D_roi[max_alpha_idx] - z_tumor)**2) * 1e3)
        
        # Guardar sidecar de métricas
        metrics_sidecar = {
            "tumor_x_mm": float(x_tumor * 1e3),
            "tumor_z_mm": float(z_tumor * 1e3),
            "cnr_c_dbr": cnr_c_db,
            "cnr_alpha_dbr": cnr_a_db,
            "peak_distance_mm": dist_c,
            "x_min_mm": float(x_grid.min() * 1e3),
            "x_max_mm": float(x_grid.max() * 1e3),
            "z_min_mm": float(z_grid.min() * 1e3),
            "z_max_mm": float(z_grid.max() * 1e3),
            "alpha_unit": "dB/MHz/cm",
            "peak_distance_alpha_mm": dist_alpha
        }
        
        metrics_path = f"articulo/metrics_{file_prefix}.json"
        with open(metrics_path, "w", encoding="utf-8") as f_metrics:
            json.dump(metrics_sidecar, f_metrics, indent=2)
        print(f"Sidecar de métricas guardado en {metrics_path}")
        
        # Guardar sidecar de convergencia
        epochs_list = [int(h['epoch']) for h in progress_history]
        loss_list = [float(h['loss']) for h in progress_history]
        err_c_list = [float(h['err_c'] * 100) for h in progress_history]
        err_alpha_list = [float(h['err_alpha0'] * 100) for h in progress_history]
        
        stage_boundaries = [
            int(progress_history[i]['epoch']) for i in range(1, len(progress_history))
            if progress_history[i]['stage'] != progress_history[i-1]['stage']
        ]
        
        convergence_sidecar = {
            "epochs": epochs_list,
            "loss": loss_list,
            "err_c": err_c_list,
            "err_alpha": err_alpha_list,
            "stage_boundaries": stage_boundaries
        }
        
        convergence_path = f"articulo/convergence_{file_prefix}.json"
        with open(convergence_path, "w", encoding="utf-8") as f_conv:
            json.dump(convergence_sidecar, f_conv, indent=2)
        print(f"Sidecar de convergencia guardado en {convergence_path}")
        
    except Exception as e:
        print(f"Error al generar los sidecars de la webapp: {e}")
        
    return True

if __name__ == "__main__":
    import argparse
    import os
    
    default_path = "data_engine/cubdl_sample_dataset.mat"
    myo_path = "data_cubdl_prueba_1/MYO001.hdf5"
    default_target_path = myo_path if os.path.exists(myo_path) else default_path
    
    parser = argparse.ArgumentParser(description="FWI Real Data Reconstruction")
    parser.add_argument("path", nargs="?", default=default_target_path, help="Path to CUBDL file")
    parser.add_argument("seed", nargs="?", type=int, default=42, help="Seed for reproducibility")
    parser.add_argument("suffix", nargs="?", default="", help="Suffix for output files")
    parser.add_argument("--npml", type=int, default=15, help="Number of PML cells")
    parser.add_argument("--roi-taper", action="store_true", help="Use smooth raised cosine taper on ROI boundary")
    parser.add_argument("--taper-cells", type=int, default=5, help="Taper cells thickness")
    parser.add_argument("--lr-decay", action="store_true", help="Enable Exponential learning rate decay")
    
    args = parser.parse_args()
    
    run_fwi_real_data(
        file_path=args.path,
        seed=args.seed,
        output_suffix=args.suffix,
        npml=args.npml,
        roi_taper=args.roi_taper,
        taper_cells=args.taper_cells,
        lr_decay=args.lr_decay
    )
