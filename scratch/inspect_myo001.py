import numpy as np
import json
from pathlib import Path

# Rutas
articulo_dir = Path(r"C:\Users\yelem\workspace\proyecto_pinn_ultrasonido_mama_fase_2\articulo")
c_path = articulo_dir / "c_map_MYO001.npy"
alpha_path = articulo_dir / "alpha_map_MYO001.npy"
metrics_path = articulo_dir / "metrics_MYO001.json"

if not c_path.exists() or not alpha_path.exists() or not metrics_path.exists():
    print("Faltan archivos de MYO001 en articulo/")
    exit(1)

# Cargar datos
c_map = np.load(c_path)
alpha_map = np.load(alpha_path)
with open(metrics_path, 'r') as f:
    metrics = json.load(f)

print("--- METADATOS METRICS_MYO001.JSON ---")
for k, v in metrics.items():
    print(f"  {k}: {v}")

Nz, Nx = c_map.shape
print(f"\nDimensiones de los mapas: {Nz} x {Nx}")

# Reconstruir coordenadas físicas
x_min, x_max = metrics["x_min_mm"], metrics["x_max_mm"]
z_min, z_max = metrics["z_min_mm"], metrics["z_max_mm"]

x_vals = np.linspace(x_min, x_max, Nx)
z_vals = np.linspace(z_min, z_max, Nz)

# Ubicar picos de reconstrucción
max_c_idx = np.unravel_index(np.argmax(c_map), c_map.shape)
max_alpha_idx = np.unravel_index(np.argmax(alpha_map), alpha_map.shape)

peak_c_x = x_vals[max_c_idx[1]]
peak_c_z = z_vals[max_c_idx[0]]

peak_alpha_x = x_vals[max_alpha_idx[1]]
peak_alpha_z = z_vals[max_alpha_idx[0]]

print("\n--- PICOS ENCONTRADOS ---")
print(f"Pico de velocidad c(x,z) max: {c_map[max_c_idx]:.2f} m/s")
print(f"  en index: {max_c_idx} -> físico: x = {peak_c_x:.2f} mm, z = {peak_c_z:.2f} mm")
print(f"Pico de atenuación alpha(x,z) max: {alpha_map[max_alpha_idx]:.4e}")
print(f"  en index: {max_alpha_idx} -> físico: x = {peak_alpha_x:.2f} mm, z = {peak_alpha_z:.2f} mm")

# Tumor/inclusión esperada
tumor_x = metrics["tumor_x_mm"]
tumor_z = metrics["tumor_z_mm"]
print(f"\nCentro esperado del tumor (DAS): x = {tumor_x:.2f} mm, z = {tumor_z:.2f} mm")

dist_c = np.sqrt((peak_c_x - tumor_x)**2 + (peak_c_z - tumor_z)**2)
dist_alpha = np.sqrt((peak_alpha_x - tumor_x)**2 + (peak_alpha_z - tumor_z)**2)

print(f"Distancia del pico c al tumor: {dist_c:.2f} mm")
print(f"Distancia del pico alpha al tumor: {dist_alpha:.2f} mm")
