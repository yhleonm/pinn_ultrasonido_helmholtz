import torch
import torch.nn as nn
import numpy as np

class SineLayer(nn.Module):
    """
    Sine activation layer (SIREN) with specialized weight initialization
    to handle high-frequency spatial representations.
    """
    def __init__(self, in_features, out_features, bias=True, is_first=False, omega_0=60.0):
        super().__init__()
        self.omega_0 = omega_0
        self.is_first = is_first
        self.linear = nn.Linear(in_features, out_features, bias=bias)
        self.init_weights()
        
    def init_weights(self):
        with torch.no_grad():
            if self.is_first:
                self.linear.weight.uniform_(-1 / self.linear.in_features, 1 / self.linear.in_features)
            else:
                self.linear.weight.uniform_(-np.sqrt(6 / self.linear.in_features) / self.omega_0, 
                                             np.sqrt(6 / self.linear.in_features) / self.omega_0)
                
    def forward(self, x):
        return torch.sin(self.omega_0 * self.linear(x))

class DecoupledUltrasoundPINN(nn.Module):
    """
    Decoupled architecture:
    1. Independent Parameter Subnetworks:
       - VelocityNet (c_net): maps (x, z) -> c(x, z) using low omega_0_c (e.g., 3.8)
         matching the Fresnel-zone resolution of wave-equation travel-time tomography.
       - AttenuationNet (alpha_net): maps (x, z) -> alpha_0(x, z) using omega_0_alpha (e.g., 5.2)
         matching log-spectral amplitude decay dynamics without sharing hidden latents with c(x, z).
    2. PressureNet: maps (x, z, omega) -> (P_real, P_imag)
       Enforces wave field propagation (using a higher omega_0, e.g., 60.0).
    """
    def __init__(
        self,
        spatial_dim=1,
        hidden_features=64,
        num_layers=4,
        omega_0_pressure=60.0,
        omega_0_params=15.0,
        omega_0_alpha=None,
        decouple_c_alpha=True,
    ):
        super().__init__()
        self.spatial_dim = spatial_dim
        self.decouple_c_alpha = decouple_c_alpha
        omega_c = omega_0_params
        omega_a = omega_0_alpha if omega_0_alpha is not None else (omega_0_params * 1.35)

        if self.decouple_c_alpha:
            # Independent Subnetwork 1A: Phase Velocity Network c(x, z)
            c_layers = [SineLayer(spatial_dim, hidden_features, is_first=True, omega_0=omega_c)]
            for _ in range(num_layers - 2):
                c_layers.append(SineLayer(hidden_features, hidden_features, omega_0=omega_c))
            self.c_net = nn.Sequential(*c_layers)
            self.c_out = nn.Linear(hidden_features, 1)

            # Independent Subnetwork 1B: Acoustic Attenuation Network alpha_0(x, z)
            a_layers = [SineLayer(spatial_dim, hidden_features, is_first=True, omega_0=omega_a)]
            for _ in range(num_layers - 2):
                a_layers.append(SineLayer(hidden_features, hidden_features, omega_0=omega_a))
            self.alpha_net = nn.Sequential(*a_layers)
            self.alpha_out = nn.Linear(hidden_features, 1)
        else:
            param_layers = [SineLayer(spatial_dim, hidden_features, is_first=True, omega_0=omega_0_params)]
            for _ in range(num_layers - 2):
                param_layers.append(SineLayer(hidden_features, hidden_features, omega_0=omega_0_params))
            self.param_net = nn.Sequential(*param_layers)
            self.param_out = nn.Linear(hidden_features, 2)

        # Subnetwork 2: Pressure Network (for P_real and P_imag)
        pressure_layers = [SineLayer(spatial_dim + 1, hidden_features, is_first=True, omega_0=omega_0_pressure)]
        for _ in range(num_layers - 2):
            pressure_layers.append(SineLayer(hidden_features, hidden_features, omega_0=omega_0_pressure))
        self.pressure_net = nn.Sequential(*pressure_layers)
        self.pressure_out = nn.Linear(hidden_features, 2)

    def forward_params(self, coords):
        """
        Predict physical parameters from spatial coordinates using independent latent representations.
        coords: [N, spatial_dim]
        Returns:
            c: Speed of sound map [N, 1] in m/s
            alpha0: Attenuation map [N, 1] in nepers/(Hz*m)
        """
        if self.decouple_c_alpha:
            out_c = self.c_out(self.c_net(coords))
            out_a = self.alpha_out(self.alpha_net(coords))
        else:
            h = self.param_net(coords)
            out = self.param_out(h)
            out_c = out[:, 0:1]
            out_a = out[:, 1:2]

        # Physical bounds
        # c in [1400, 1700] m/s
        c = 1400.0 + torch.sigmoid(out_c) * 300.0

        # alpha0 in [0, 2.5] dB/(MHz*cm) -> converted to nepers/(Hz*m)
        alpha0_db = torch.sigmoid(out_a) * 2.5
        alpha0 = alpha0_db / (8.686e4)

        return c, alpha0
        
    def forward_pressure(self, coords_and_freq):
        """
        Predict complex pressure field.
        coords_and_freq: [N, spatial_dim + 1]
        Returns:
            p_real: [N, 1]
            p_imag: [N, 1]
        """
        h = self.pressure_net(coords_and_freq)
        out = self.pressure_out(h)
        p_real = out[:, 0:1]
        p_imag = out[:, 1:2]
        return p_real, p_imag
        
    def forward(self, coords_and_freq):
        """
        Full forward pass.
        coords_and_freq: [N, spatial_dim + 1]
        """
        coords = coords_and_freq[:, :self.spatial_dim]
        c, alpha0 = self.forward_params(coords)
        p_real, p_imag = self.forward_pressure(coords_and_freq)
        return p_real, p_imag, alpha0, c
