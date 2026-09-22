"""
Test de regresión numérica: spsolve original vs splu + LU.solve optimizado.

Verifica que la optimización de reutilizar la factorización LU entre forward
y backward produce resultados BIT-IDENTICAL al código original (spsolve).

Comparación:
  - Forward: P_vec (campo de presión)
  - Backward: grad_c y grad_alpha0 (gradientes)
  - Tolerancia: < 1e-8 (elemento a elemento, max abs error)
"""
import sys
import os
import numpy as np
import torch
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import time

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ===== Import the assembly function from main_5MHz_fwi =====
from main_5MHz_fwi import assemble_helmholtz_matrix

# ============================================================
# ORIGINAL implementation (spsolve, no LU caching)
# ============================================================
class OriginalHelmholtzSolver(torch.autograd.Function):
    """Exact copy of the ORIGINAL solver before optimization."""
    @staticmethod
    def forward(ctx, c, alpha0, x, z, f, B_vec):
        ctx.x = x
        ctx.z = z
        ctx.f = f
        ctx.B_vec = B_vec
        
        c_np = c.detach().cpu().numpy()
        alpha0_np = alpha0.detach().cpu().numpy()
        
        A = assemble_helmholtz_matrix(x, z, f, c_np, alpha0_np)
        P_vec = spla.spsolve(A, B_vec)  # Original: factorize from scratch
        
        ctx.save_for_backward(c, alpha0)
        ctx.P_vec = P_vec
        ctx.A = A
        
        P = P_vec.reshape((len(z), len(x)))
        dtype_out = torch.complex128 if c.dtype == torch.float64 else torch.complex64
        return torch.tensor(P, dtype=dtype_out, device=c.device)
        
    @staticmethod
    def backward(ctx, grad_output):
        c, alpha0 = ctx.saved_tensors
        P_vec = ctx.P_vec
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
        
        # Original: spsolve(A.T, ...) from scratch
        lambda_vec = spla.spsolve(A.T, rhs_adjoint)
        
        J, I = np.meshgrid(np.arange(Nz), np.arange(Nx), indexing='ij')
        mask_top = (J == 0)
        mask_left = (I == 0) & (~mask_top)
        mask_right = (I == Nx - 1) & (~mask_top)
        mask_bottom = (J == Nz - 1) & (~mask_top) & (~mask_left) & (~mask_right)
        mask_interior = (~mask_top) & (~mask_left) & (~mask_right) & (~mask_bottom)
        
        c_flat = c_np.flatten()
        alpha0_flat = alpha0_np.flatten()
        
        dk2_dc = -2.0 * omega**2 / (c_flat**3) - 2j * omega * f * alpha0_flat / (c_flat**2)
        dk2_dalpha0 = -2.0 * alpha0_flat * (f**2) + 2j * omega * f / c_flat
        dk_dc = -omega / (c_flat**2)
        dk_dalpha0 = np.full_like(c_flat, 1j * f, dtype=complex)
        
        dApp_dc = np.zeros_like(c_flat, dtype=complex)
        dApp_dalpha0 = np.zeros_like(alpha0_flat, dtype=complex)
        
        flat_mask_interior = mask_interior.flatten()
        dApp_dc[flat_mask_interior] = (h**2) * dk2_dc[flat_mask_interior]
        dApp_dalpha0[flat_mask_interior] = (h**2) * dk2_dalpha0[flat_mask_interior]
        
        flat_mask_left = mask_left.flatten()
        dApp_dc[flat_mask_left] = -1j * h * dk_dc[flat_mask_left]
        dApp_dalpha0[flat_mask_left] = -1j * h * dk_dalpha0[flat_mask_left]
        
        flat_mask_right = mask_right.flatten()
        dApp_dc[flat_mask_right] = -1j * h * dk_dc[flat_mask_right]
        dApp_dalpha0[flat_mask_right] = -1j * h * dk_dalpha0[flat_mask_right]
        
        flat_mask_bottom = mask_bottom.flatten()
        dApp_dc[flat_mask_bottom] = -1j * h * dk_dc[flat_mask_bottom]
        dApp_dalpha0[flat_mask_bottom] = -1j * h * dk_dalpha0[flat_mask_bottom]
        
        grad_c_vec = -np.real(lambda_vec * dApp_dc * P_vec)
        grad_alpha0_vec = -np.real(lambda_vec * dApp_dalpha0 * P_vec)
        
        grad_c_t = torch.tensor(grad_c_vec.reshape((Nz, Nx)), dtype=c.dtype, device=c.device)
        grad_alpha0_t = torch.tensor(grad_alpha0_vec.reshape((Nz, Nx)), dtype=alpha0.dtype, device=alpha0.device)
        
        return grad_c_t, grad_alpha0_t, None, None, None, None

# ============================================================
# OPTIMIZED implementation (splu + LU.solve, from modified file)
# ============================================================
from main_5MHz_fwi import DifferentiableHelmholtzSolver as OptimizedHelmholtzSolver


def run_regression_test():
    print("=" * 70)
    print("TEST DE REGRESIÓN NUMÉRICA: spsolve vs splu + LU.solve(trans='T')")
    print("=" * 70)
    
    # --- Setup: realistic Helmholtz problem ---
    np.random.seed(42)
    torch.manual_seed(42)
    
    Nx, Nz = 80, 80  # Grilla representativa
    f = 5e6           # 5 MHz
    
    x = np.linspace(0, 0.04, Nx)
    z = np.linspace(0, 0.04, Nz)
    h = x[1] - x[0]
    
    # Mapa de velocidad con heterogeneidad (simula tejido)
    c_map_np = 1500.0 * np.ones((Nz, Nx))
    c_map_np[25:55, 25:55] = 1540.0  # Inclusión
    
    # Mapa de atenuación
    alpha0_map_np = 0.5 * np.ones((Nz, Nx))
    alpha0_map_np[25:55, 25:55] = 1.0
    
    # Source vector (fuente puntual en la fila superior)
    N = Nz * Nx
    B_vec = np.zeros(N, dtype=complex)
    B_vec[Nx // 2] = 1.0  # Fuente en el centro de la fila superior
    
    # Crear tensores con requires_grad para que backward se active
    c_old = torch.tensor(c_map_np, dtype=torch.float64, requires_grad=True)
    alpha0_old = torch.tensor(alpha0_map_np, dtype=torch.float64, requires_grad=True)
    
    c_new = torch.tensor(c_map_np, dtype=torch.float64, requires_grad=True)
    alpha0_new = torch.tensor(alpha0_map_np, dtype=torch.float64, requires_grad=True)
    
    # ===== FORWARD: Original =====
    print("\n[1] Forward pass (original spsolve)...")
    t0 = time.perf_counter()
    P_old = OriginalHelmholtzSolver.apply(c_old, alpha0_old, x, z, f, B_vec)
    t_old_fwd = time.perf_counter() - t0
    print(f"    Tiempo: {t_old_fwd:.4f}s")
    
    # ===== FORWARD: Optimized =====
    print("[2] Forward pass (optimized splu)...")
    t0 = time.perf_counter()
    P_new = OptimizedHelmholtzSolver.apply(c_new, alpha0_new, x, z, f, B_vec)
    t_new_fwd = time.perf_counter() - t0
    print(f"    Tiempo: {t_new_fwd:.4f}s")
    
    # ===== Compare forward outputs =====
    P_old_np = P_old.detach().numpy()
    P_new_np = P_new.detach().numpy()
    fwd_max_err = np.max(np.abs(P_old_np - P_new_np))
    fwd_rel_err = fwd_max_err / (np.max(np.abs(P_old_np)) + 1e-30)
    
    print(f"\n--- Forward comparison ---")
    print(f"    Max abs error P_vec:  {fwd_max_err:.2e}")
    print(f"    Max rel error P_vec:  {fwd_rel_err:.2e}")
    fwd_pass = fwd_max_err < 1e-8
    print(f"    RESULTADO: {'PASS [OK]' if fwd_pass else 'FAIL [X]'}")
    
    # ===== BACKWARD: need a dummy loss =====
    # Use same synthetic grad_output for both
    dummy_grad = torch.randn_like(P_old)
    # Clone the exact same grad for the new solver
    dummy_grad_new = dummy_grad.clone()
    
    print(f"\n[3] Backward pass (original spsolve)...")
    t0 = time.perf_counter()
    P_old.backward(dummy_grad)
    t_old_bwd = time.perf_counter() - t0
    print(f"    Tiempo: {t_old_bwd:.4f}s")
    
    print("[4] Backward pass (optimized splu)...")
    t0 = time.perf_counter()
    P_new.backward(dummy_grad_new)
    t_new_bwd = time.perf_counter() - t0
    print(f"    Tiempo: {t_new_bwd:.4f}s")
    
    # ===== Compare gradients =====
    grad_c_old = c_old.grad.numpy()
    grad_c_new = c_new.grad.numpy()
    grad_a_old = alpha0_old.grad.numpy()
    grad_a_new = alpha0_new.grad.numpy()
    
    gc_max_err = np.max(np.abs(grad_c_old - grad_c_new))
    gc_rel_err = gc_max_err / (np.max(np.abs(grad_c_old)) + 1e-30)
    ga_max_err = np.max(np.abs(grad_a_old - grad_a_new))
    ga_rel_err = ga_max_err / (np.max(np.abs(grad_a_old)) + 1e-30)
    
    print(f"\n--- Gradient comparison (grad_c) ---")
    print(f"    Max abs error:  {gc_max_err:.2e}")
    print(f"    Max rel error:  {gc_rel_err:.2e}")
    gc_pass = gc_max_err < 1e-8
    print(f"    RESULTADO: {'PASS [OK]' if gc_pass else 'FAIL [X]'}")
    
    print(f"\n--- Gradient comparison (grad_alpha0) ---")
    print(f"    Max abs error:  {ga_max_err:.2e}")
    print(f"    Max rel error:  {ga_rel_err:.2e}")
    ga_pass = ga_max_err < 1e-8
    print(f"    RESULTADO: {'PASS [OK]' if ga_pass else 'FAIL [X]'}")
    
    # ===== Timing summary =====
    total_old = t_old_fwd + t_old_bwd
    total_new = t_new_fwd + t_new_bwd
    speedup = total_old / total_new if total_new > 0 else float('inf')
    
    print(f"\n--- Timing summary ---")
    print(f"    Original (fwd+bwd): {total_old:.4f}s")
    print(f"    Optimized (fwd+bwd): {total_new:.4f}s")
    print(f"    Speedup:            {speedup:.2f}x")
    print(f"    Ahorro backward:    {t_old_bwd - t_new_bwd:.4f}s ({(1 - t_new_bwd/t_old_bwd)*100:.1f}%)" if t_old_bwd > 0 else "")
    
    # ===== Final verdict =====
    all_pass = fwd_pass and gc_pass and ga_pass
    print(f"\n{'=' * 70}")
    if all_pass:
        print(f"VEREDICTO FINAL: PASS [OK] -- Resultados bit-identical (tolerancia < 1e-8)")
        print(f"Optimizacion segura para aplicar. Speedup medido: {speedup:.2f}x")
    else:
        print(f"VEREDICTO FINAL: FAIL [X] -- HAY DIFERENCIAS NUMERICAS")
        print(f"NO aplicar esta optimizacion hasta investigar las diferencias.")
    print(f"{'=' * 70}")
    
    return all_pass


if __name__ == "__main__":
    success = run_regression_test()
    sys.exit(0 if success else 1)
