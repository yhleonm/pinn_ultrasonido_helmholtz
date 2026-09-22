import torch

def compute_helmholtz_residuals_1d(coords, p_real, p_imag, alpha0, c, omega, c0=1500.0):
    """
    Computes normalized 1D Helmholtz PDE residuals (scaled by k0^2).
    """
    coords.requires_grad_(True)
    
    # First derivative
    dp_real_dz = torch.autograd.grad(p_real.sum(), coords, create_graph=True)[0]
    dp_imag_dz = torch.autograd.grad(p_imag.sum(), coords, create_graph=True)[0]
    
    # Second derivative
    d2p_real_dz2 = torch.autograd.grad(dp_real_dz.sum(), coords, create_graph=True)[0]
    d2p_imag_dz2 = torch.autograd.grad(dp_imag_dz.sum(), coords, create_graph=True)[0]
    
    k0 = omega / c0
    f = omega / (2.0 * 3.14159265)
    alpha = alpha0 * f
    
    w_over_c = omega / c
    k2_real = w_over_c**2 - alpha**2
    k2_imag = 2.0 * w_over_c * alpha
    
    # Normalize by dividing the entire Helmholtz equation by k0^2
    res_real = (d2p_real_dz2 + p_real * k2_real - p_imag * k2_imag) / (k0**2)
    res_imag = (d2p_imag_dz2 + p_imag * k2_real + p_real * k2_imag) / (k0**2)
    
    return res_real, res_imag

def compute_sommerfeld_loss_1d(boundary_coords, p_real, p_imag, alpha0, c, omega, c0=1500.0):
    """
    Computes normalized 1D Sommerfeld loss (scaled by k0).
    """
    boundary_coords.requires_grad_(True)
    
    dp_real_dz = torch.autograd.grad(p_real.sum(), boundary_coords, create_graph=True)[0]
    dp_imag_dz = torch.autograd.grad(p_imag.sum(), boundary_coords, create_graph=True)[0]
    
    k0 = omega / c0
    f = omega / (2.0 * 3.14159265)
    alpha = alpha0 * f
    
    k_real = omega / c
    k_imag = alpha
    
    # Normalize by dividing by k0
    res_sommer_real = (dp_real_dz + k_imag * p_real + k_real * p_imag) / k0
    res_sommer_imag = (dp_imag_dz + k_imag * p_imag - k_real * p_real) / k0
    
    return torch.mean(res_sommer_real**2 + res_sommer_imag**2)

def compute_helmholtz_residuals_2d(coords, p_real, p_imag, alpha0, c, omega, c0=1500.0):
    """
    Computes normalized 2D Helmholtz PDE residuals (scaled by k0^2).
    """
    coords.requires_grad_(True)
    
    # First derivatives
    dp_real_dcoords = torch.autograd.grad(p_real.sum(), coords, create_graph=True)[0]
    dp_imag_dcoords = torch.autograd.grad(p_imag.sum(), coords, create_graph=True)[0]
    dp_real_dx, dp_real_dz = dp_real_dcoords[:, 0:1], dp_real_dcoords[:, 1:2]
    dp_imag_dx, dp_imag_dz = dp_imag_dcoords[:, 0:1], dp_imag_dcoords[:, 1:2]
    
    # Second derivatives
    d2p_real_dx2 = torch.autograd.grad(dp_real_dx.sum(), coords, create_graph=True)[0][:, 0:1]
    d2p_real_dz2 = torch.autograd.grad(dp_real_dz.sum(), coords, create_graph=True)[0][:, 1:2]
    d2p_imag_dx2 = torch.autograd.grad(dp_imag_dx.sum(), coords, create_graph=True)[0][:, 0:1]
    d2p_imag_dz2 = torch.autograd.grad(dp_imag_dz.sum(), coords, create_graph=True)[0][:, 1:2]
    
    laplacian_real = d2p_real_dx2 + d2p_real_dz2
    laplacian_imag = d2p_imag_dx2 + d2p_imag_dz2
    
    k0 = omega / c0
    f = omega / (2.0 * 3.14159265)
    alpha = alpha0 * f
    
    w_over_c = omega / c
    k2_real = w_over_c**2 - alpha**2
    k2_imag = 2.0 * w_over_c * alpha
    
    # Normalize by dividing the entire Helmholtz equation by k0^2
    res_real = (laplacian_real + p_real * k2_real - p_imag * k2_imag) / (k0**2)
    res_imag = (laplacian_imag + p_imag * k2_real + p_real * k2_imag) / (k0**2)
    
    return res_real, res_imag

def compute_sommerfeld_loss_2d(boundary_coords, p_real, p_imag, alpha0, c, omega, normals, c0=1500.0):
    """
    Computes normalized 2D Sommerfeld loss (scaled by k0).
    """
    boundary_coords.requires_grad_(True)
    
    dp_real_dcoords = torch.autograd.grad(p_real.sum(), boundary_coords, create_graph=True)[0]
    dp_imag_dcoords = torch.autograd.grad(p_imag.sum(), boundary_coords, create_graph=True)[0]
    
    dp_real_dn = (dp_real_dcoords * normals).sum(dim=1, keepdim=True)
    dp_imag_dn = (dp_imag_dcoords * normals).sum(dim=1, keepdim=True)
    
    k0 = omega / c0
    f = omega / (2.0 * 3.14159265)
    alpha = alpha0 * f
    
    k_real = omega / c
    k_imag = alpha
    
    # Normalize by dividing by k0
    res_sommer_real = (dp_real_dn + k_imag * p_real + k_real * p_imag) / k0
    res_sommer_imag = (dp_imag_dn + k_imag * p_imag - k_real * p_real) / k0
    
    return torch.mean(res_sommer_real**2 + res_sommer_imag**2)
