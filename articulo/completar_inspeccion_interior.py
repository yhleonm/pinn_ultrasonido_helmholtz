"""
articulo/completar_inspeccion_interior.py
==========================================
Aplica la inspección física de distancia al centro de la ROI para clasificar
el pico de velocidad en los 6 datasets (MYO001-MYO006).

Criterio de clasificación física:
    1. Extraer la ROI de 1.5x1.5 cm centrada en el tumor.
    2. Excluir el 15% de margen del SIREN para obtener el interior de la ROI.
    3. Localizar el índice del valor máximo de velocidad (c) en el interior.
    4. Calcular la distancia física en milímetros de este pico respecto al centro del tumor (0.0, 0.0).
    5. Si distancia <= 3.0 mm (radio del tumor) -> LIMPIO (el pico está dentro de la inclusión física).
    6. Si distancia > 3.0 mm -> FUGA (el pico está en la periferia, dominado por el artefacto de borde).
"""

from __future__ import annotations
from pathlib import Path
import numpy as np

ARTICULO_DIR = Path(__file__).parent
MARGIN_FRAC = 0.15          # margen de borde ya excluido
TUMOR_RADIUS_MM = 3.0       # radio físico del tumor en el fantoma CIRS 040GSE

DATASETS = ["MYO001", "MYO002", "MYO003", "MYO004", "MYO005", "MYO006"]

TUMOR_LOCATIONS = {
    "MYO001": {"x": 2.00e-3, "z": 18.50e-3},
    "MYO002": {"x": 0.00e-3, "z": 23.50e-3},
    "MYO003": {"x": 16.00e-3, "z": 18.00e-3},
    "MYO004": {"x": -2.50e-3, "z": 11.00e-3},
    "MYO005": {"x": 4.00e-3, "z": 11.50e-3},
    "MYO006": {"x": -16.50e-3, "z": 10.50e-3},
}


def load_map_safe(dataset: str, param: str) -> np.ndarray | None:
    fname = ARTICULO_DIR / f"{param}_map_{dataset}.npy"
    if not fname.exists():
        print(f"   {fname.name} no encontrado")
        return None
    return np.load(fname).astype(np.float32)


def extract_roi_subarray(arr: np.ndarray, x_tumor: float, z_tumor: float, roi_half: float = 0.0075) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extrae el sub-array 2D rectangular de la ROI activa y sus coordenadas relativas."""
    Nz, Nx = arr.shape
    x_grid = np.linspace(-0.01905, 0.01905, Nx)
    z_grid = np.linspace(0.0, 0.032, Nz)
    X2D, Z2D = np.meshgrid(x_grid - x_tumor, z_grid - z_tumor)
    
    mask = (x_grid >= x_tumor - roi_half) & (x_grid <= x_tumor + roi_half)
    mask_z = (z_grid >= z_tumor - roi_half) & (z_grid <= z_tumor + roi_half)
    
    z_idx = np.where(mask_z)[0]
    x_idx = np.where(mask)[0]
    
    z0, z1 = z_idx.min(), z_idx.max() + 1
    x0, x1 = x_idx.min(), x_idx.max() + 1
    
    return arr[z0:z1, x0:x1], X2D[z0:z1, x0:x1], Z2D[z0:z1, x0:x1]


def inspect_peak_distance(roi_2d: np.ndarray, X2D_roi: np.ndarray, Z2D_roi: np.ndarray) -> dict:
    """
    Clasifica el dataset según la distancia del pico interior al centro físico.
    """
    Nz_roi, Nx_roi = roi_2d.shape
    mz = int(MARGIN_FRAC * Nz_roi)
    mx = int(MARGIN_FRAC * Nx_roi)

    # Crear máscara de interior
    interior_mask = np.zeros_like(roi_2d, dtype=bool)
    interior_mask[mz:-mz, mx:-mx] = True

    # Buscar el máximo dentro del interior
    roi_int = np.where(interior_mask, roi_2d, -1e9)
    max_idx = np.unravel_index(np.argmax(roi_int), roi_int.shape)

    # Calcular distancia física del pico al centro
    dx_mm = X2D_roi[max_idx] * 1e3
    dz_mm = Z2D_roi[max_idx] * 1e3
    dist_mm = np.sqrt(dx_mm**2 + dz_mm**2)

    is_clean = dist_mm <= TUMOR_RADIUS_MM

    return {
        "max_val": roi_2d[max_idx],
        "dist_mm": dist_mm,
        "is_clean": is_clean,
        "diagnosis": "LIMPIO" if is_clean else "FUGA",
        "coords": max_idx
    }


def main():
    print("="*70)
    print("Verificación completa: inspección centro-vs-periferia física sobre ROI")
    print("="*70)

    results = {}
    for ds in DATASETS:
        c_arr = load_map_safe(ds, "c")
        if c_arr is None:
            print(f"{ds}: SALTADO - archivo no encontrado")
            continue

        loc = TUMOR_LOCATIONS[ds]
        roi_2d, X2D_roi, Z2D_roi = extract_roi_subarray(c_arr, loc["x"], loc["z"])
        r = inspect_peak_distance(roi_2d, X2D_roi, Z2D_roi)
        results[ds] = r

        print(f"{ds}:")
        print(f"  Valor máximo en interior:         {r['max_val']:.1f} m/s")
        print(f"  Ubicación del pico (idx):         {r['coords']}")
        print(f"  Distancia al centro del tumor:    {r['dist_mm']:.2f} mm")
        print(f"  Diagnóstico: {r['diagnosis']}")
        print()

    # Tabla final consolidada
    print("="*70)
    print("TABLA CONSOLIDADA: los 6 datasets con criterio uniforme")
    print("="*70)

    cnr_c_values = {
        "MYO001": 1.51, "MYO002": 2.97, "MYO003": 6.26,
        "MYO004": 9.96, "MYO005": 3.53, "MYO006": 9.20,
    }
    z_values = {
        "MYO001": 18.5, "MYO002": 23.5, "MYO003": 18.0,
        "MYO004": 11.0, "MYO005": 11.5, "MYO006": 10.5,
    }

    print(f"\n{'Dataset':<8} {'z(mm)':>6} {'CNR_c':>6} {'Dist(mm)':>10} {'Diagnóstico':<10}")
    for ds in DATASETS:
        res = results.get(ds, {})
        diag = res.get("diagnosis", "?")
        dist = res.get("dist_mm", 0.0)
        print(f"{ds:<8} {z_values[ds]:>6.1f} {cnr_c_values[ds]:>6.2f} {dist:>10.2f} {diag:<10}")

    # Recalcular correlación SOLO con datasets confirmados LIMPIOS
    limpios = [ds for ds in DATASETS if results.get(ds, {}).get("is_clean")]
    print(f"\nDatasets confirmados LIMPIOS (sin fuga de borde): {limpios}")

    if len(limpios) >= 3:
        z_clean = np.array([z_values[ds] for ds in limpios])
        c_clean = np.array([cnr_c_values[ds] for ds in limpios])
        corr_clean = np.corrcoef(z_clean, c_clean)[0,1]
        print(f"\nCorrelación z vs CNR_c usando SOLO datasets verificados limpios (n={len(limpios)}):")
        print(f"  r = {corr_clean:.3f}")
    else:
        print(f"\nSolo {len(limpios)} datasets confirmados limpios - insuficiente para")
        print("  una correlación confiable (se requiere n>=4 para Pearson con N pequeño).")
        print("  Recomendación para el artículo: reportar como una LIMITACIÓN METODOLÓGICA.")
        print("  El método sin PML no permite aún medir la dependencia de la profundidad")
        print("  de forma confiable porque el artefacto de frontera contamina 5/6 de las ejecuciones,")
        print("  inflando artificialmente el CNR de forma no uniforme.")


if __name__ == "__main__":
    main()
