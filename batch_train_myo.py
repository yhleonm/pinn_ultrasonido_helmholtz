"""
Batch runner: entrena FWI-PINN sobre MYO002-MYO006 secuencialmente.
MYO001 ya fue procesado previamente.
Cada run genera:
  - articulo/best_pinn_5MHz_fwi_real_MYO00X.pt
  - articulo/c_map_MYO00X.npy
  - articulo/alpha_map_MYO00X.npy
  - articulo/metricas_progreso_fwi_real_cubdl_MYO00X.md
  - articulo/resultados_fwi_real_cubdl_MYO00X.png
"""
import sys
import os
import time

# Los datasets a procesar (MYO001 ya existe)
datasets = [
    "data_cubdl_prueba_1/MYO001.hdf5",
    "data_cubdl_prueba_1/MYO002.hdf5",
    "data_cubdl_prueba_1/MYO003.hdf5",
    "data_cubdl_prueba_1/MYO004.hdf5",
    "data_cubdl_prueba_1/MYO005.hdf5",
    "data_cubdl_prueba_1/MYO006.hdf5",
]

from main_5MHz_fwi_real_data import run_fwi_real_data

total_start = time.time()

for i, ds_path in enumerate(datasets):
    ds_name = os.path.splitext(os.path.basename(ds_path))[0]
    print(f"\n{'#'*70}")
    print(f"# BATCH RUN {i+1}/{len(datasets)}: {ds_name}")
    print(f"{'#'*70}")
    
    if not os.path.exists(ds_path):
        print(f"  ERROR: {ds_path} no existe, saltando.")
        continue
    
    # Verificar si ya fue procesado
    npy_check = f"articulo/c_map_{ds_name}.npy"
    if os.path.exists(npy_check):
        print(f"  {npy_check} ya existe. Saltando re-entrenamiento.")
        print(f"  (Elimina el archivo .npy si quieres re-entrenar.)")
        continue
    
    run_start = time.time()
    try:
        run_fwi_real_data(ds_path)
    except Exception as e:
        print(f"  ERROR en {ds_name}: {e}")
        import traceback
        traceback.print_exc()
        continue
    
    run_elapsed = time.time() - run_start
    print(f"\n  {ds_name} completado en {run_elapsed:.0f}s ({run_elapsed/60:.1f} min)")

total_elapsed = time.time() - total_start
print(f"\n{'#'*70}")
print(f"# BATCH COMPLETADO: {len(datasets)} datasets en {total_elapsed:.0f}s ({total_elapsed/60:.1f} min)")
print(f"{'#'*70}")
print(f"\nAhora ejecuta: python articulo/sweep_consistencia_myo.py")
