"""
articulo/sweep_consistencia_myo.py
====================================
Sweep de consistencia N=5 sobre fantomas MYO001-MYO006 (CIRS 040GSE).

Objetivo: verificar que el contraste CNR robusto (P90/P50 con margen
de borde excluido) es reproducible across multiples adquisiciones del
MISMO fantoma fisico, no un artefacto de una sola ejecucion afortunada.

Si los 5-6 valores de CNR caen en un rango estrecho y consistente,
es evidencia fuerte de que el metodo generaliza. Si varian erraticamente
o cambian de signo, el problema es mas profundo que el artefacto de borde
ya conocido, y PML pasa a ser prioridad inmediata.

Uso:
    1. Ejecutar main_5MHz_fwi_real_data.py sobre cada MYO00X.hdf5
       (genera c_map_MYO00X.npy y alpha_map_MYO00X.npy en articulo/)
    2. Ejecutar este script:
       python articulo/sweep_consistencia_myo.py

Salida:
    articulo/sweep_consistencia_myo.md    tabla + estadisticas de dispersion
    articulo/sweep_consistencia_myo.png   grafico de CNR por fantoma con banda de consistencia
"""

from __future__ import annotations
import json
import os
import sys
from pathlib import Path
import numpy as np

# Add project root to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

ARTICULO_DIR = Path(__file__).parent

#  Configuracion comun a todos los MYO (mismo fantoma CIRS 040GSE) 
MARGIN_FRAC = 0.15          # margen de borde excluido (igual al usado en MYO001/TSH002)
C_REF_BG    = 1540.0        # m/s
C_REF_INC   = 1580.0        # m/s  corregido: c_opt real del CIRS 040GSE (no 1620, ese fue un valor de ejemplo previo)
A_REF_BG    = 0.5           # dB/MHz/cm
A_REF_INC   = 0.7           # dB/MHz/cm  segun CSV de fabricante (no 1.2, revisar)

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
        print(f"   {fname.name} no encontrado  omitiendo {dataset}")
        return None
    return np.load(fname).astype(np.float32)


def extract_roi_subarray(arr: np.ndarray, x_tumor: float, z_tumor: float, roi_half: float = 0.0075) -> np.ndarray:
    """Extrae el sub-array 2D rectangular de la ROI activa reconstruida."""
    Nz, Nx = arr.shape
    x_grid = np.linspace(-0.01905, 0.01905, Nx)
    z_grid = np.linspace(0.0, 0.032, Nz)
    X2D, Z2D = np.meshgrid(x_grid, z_grid)
    
    mask = (X2D >= x_tumor - roi_half) & (X2D <= x_tumor + roi_half) & \
           (Z2D >= z_tumor - roi_half) & (Z2D <= z_tumor + roi_half)
    
    z_idx, x_idx = np.where(mask)
    if len(z_idx) == 0 or len(x_idx) == 0:
        return arr
        
    z0, z1 = z_idx.min(), z_idx.max() + 1
    x0, x1 = x_idx.min(), x_idx.max() + 1
    return arr[z0:z1, x0:x1]


def compute_robust_cnr(arr: np.ndarray, x_tumor: float, z_tumor: float, margin_frac: float = MARGIN_FRAC) -> dict:
    """
    CNR robusto con margen de borde excluido.
    Calcula percentiles sobre la ROI interior.
    """
    roi_2d = extract_roi_subarray(arr, x_tumor, z_tumor)
    Nz_roi, Nx_roi = roi_2d.shape
    mz = int(margin_frac * Nz_roi)
    mx = int(margin_frac * Nx_roi)

    if mz == 0 or mx == 0 or Nz_roi - 2*mz <= 0 or Nx_roi - 2*mx <= 0:
        interior = roi_2d
    else:
        interior = roi_2d[mz:-mz, mx:-mx]

    p90 = np.percentile(interior, 90)
    p50 = np.percentile(interior, 50)
    bg_vals = interior[interior < p90]
    sigma_local = bg_vals.std() if bg_vals.size > 0 else interior.std()

    cnr_linear = abs(p90 - p50) / (sigma_local + 1e-12)
    cnr_db = 20 * np.log10(cnr_linear + 1e-12)

    # Posicion del pico para verificar si sigue cayendo en el borde de la ROI completa
    iz, ix = np.unravel_index(np.argmax(roi_2d), roi_2d.shape)
    edge_frac_z = min(iz, Nz_roi - 1 - iz) / Nz_roi
    edge_frac_x = min(ix, Nx_roi - 1 - ix) / Nx_roi
    near_edge = min(edge_frac_z, edge_frac_x) < margin_frac

    return {
        "p90": p90, "p50": p50, "sigma_local": sigma_local,
        "cnr_linear": cnr_linear, "cnr_db": cnr_db,
        "peak_near_edge": near_edge,
        "edge_frac": min(edge_frac_z, edge_frac_x),
    }


def main():
    results = []

    print("="*60)
    print("Sweep de consistencia  MYO001-MYO006 (CIRS 040GSE)")
    print("="*60)

    for ds in DATASETS:
        print(f"\n{ds}:")
        c_arr = load_map_safe(ds, "c")
        a_arr = load_map_safe(ds, "alpha")

        if c_arr is None or a_arr is None:
            results.append({"dataset": ds, "available": False})
            continue

        loc = TUMOR_LOCATIONS[ds]
        cnr_c = compute_robust_cnr(c_arr, loc["x"], loc["z"])
        cnr_a = compute_robust_cnr(a_arr, loc["x"], loc["z"])

        print(f"  CNR c:  {cnr_c['cnr_db']:.2f} dBr (pico cerca de borde: {cnr_c['peak_near_edge']})")
        print(f"  CNR a:  {cnr_a['cnr_db']:.2f} dBr (pico cerca de borde: {cnr_a['peak_near_edge']})")

        # Calcular distancia del pico real dentro de la ROI
        roi_2d = extract_roi_subarray(c_arr, loc["x"], loc["z"])
        Nz_roi, Nx_roi = roi_2d.shape
        
        Nx = c_arr.shape[1]
        Nz = c_arr.shape[0]
        x_grid = np.linspace(-0.01905, 0.01905, Nx)
        z_grid = np.linspace(0.0, 0.032, Nz)
        X2D, Z2D = np.meshgrid(x_grid, z_grid)
        
        roi_half = 0.0075
        mask_x = (x_grid >= loc["x"] - roi_half) & (x_grid <= loc["x"] + roi_half)
        mask_z = (z_grid >= loc["z"] - roi_half) & (z_grid <= loc["z"] + roi_half)
        z_idx = np.where(mask_z)[0]
        x_idx = np.where(mask_x)[0]
        
        z0, z1 = z_idx.min(), z_idx.max() + 1
        x0, x1 = x_idx.min(), x_idx.max() + 1
        
        X2D_roi = X2D[z0:z1, x0:x1]
        Z2D_roi = Z2D[z0:z1, x0:x1]
        
        mz = int(MARGIN_FRAC * Nz_roi)
        mx = int(MARGIN_FRAC * Nx_roi)
        interior_mask = np.zeros_like(roi_2d, dtype=bool)
        interior_mask[mz:-mz, mx:-mx] = True
        
        roi_int = np.where(interior_mask, roi_2d, -1e9)
        max_idx = np.unravel_index(np.argmax(roi_int), roi_int.shape)
        
        dx_mm = (X2D_roi[max_idx] - loc["x"]) * 1e3
        dz_mm = (Z2D_roi[max_idx] - loc["z"]) * 1e3
        dist_mm = float(np.sqrt(dx_mm**2 + dz_mm**2))

        # Obtener los bounds reales desde el archivo HDF5
        try:
            from data_engine.loader import CUBDLLoader
            hdf5_path = Path("data_cubdl_prueba_1") / f"{ds}.hdf5"
            loader = CUBDLLoader(str(hdf5_path))
            rx_x = loader.rx_positions[:, 0]
            x_min_mm = float(min(rx_x) * 1e3)
            x_max_mm = float(max(rx_x) * 1e3)
            z_min_mm = 0.0
            z_max_mm = 32.0  # Profundidad constante del resolvedor FWI
        except Exception as e:
            print(f"  Warning: No se pudieron extraer los bounds reales de {ds}: {e}. Usando defaults.")
            x_min_mm = -19.05
            x_max_mm = 19.05
            z_min_mm = 0.0
            z_max_mm = 32.0

        # Escribir sidecar JSON para la webapp
        metrics_sidecar = {
            "tumor_x_mm": float(loc["x"] * 1e3),
            "tumor_z_mm": float(loc["z"] * 1e3),
            "cnr_c_dbr": float(cnr_c["cnr_db"]),
            "cnr_alpha_dbr": float(cnr_a["cnr_db"]),
            "peak_distance_mm": dist_mm,
            "x_min_mm": x_min_mm,
            "x_max_mm": x_max_mm,
            "z_min_mm": z_min_mm,
            "z_max_mm": z_max_mm,
            "alpha_unit": "dB/MHz/cm"
        }
        
        sidecar_path = ARTICULO_DIR / f"metrics_{ds}.json"
        with open(sidecar_path, "w", encoding="utf-8") as f_side:
            json.dump(metrics_sidecar, f_side, indent=2)
        print(f"  Sidecar guardado: {sidecar_path.name} | peak_dist: {dist_mm:.2f} mm | bounds: X[{x_min_mm:.2f}, {x_max_mm:.2f}] Z[{z_min_mm:.1f}, {z_max_mm:.1f}]")

        results.append({
            "dataset": ds, "available": True,
            "cnr_c_db": cnr_c["cnr_db"], "cnr_a_db": cnr_a["cnr_db"],
            "peak_edge_c": cnr_c["peak_near_edge"], "peak_edge_a": cnr_a["peak_near_edge"],
        })

    #  Estadisticas de consistencia 
    available = [r for r in results if r["available"]]
    n_available = len(available)

    print(f"\n{'='*60}")
    print(f"Datasets disponibles: {n_available}/{len(DATASETS)}")

    if n_available >= 2:
        cnr_c_vals = [r["cnr_c_db"] for r in available]
        cnr_a_vals = [r["cnr_a_db"] for r in available]

        print(f"\nCNR c    media: {np.mean(cnr_c_vals):.2f} dBr | std: {np.std(cnr_c_vals):.2f} | rango: [{min(cnr_c_vals):.2f}, {max(cnr_c_vals):.2f}]")
        print(f"CNR a    media: {np.mean(cnr_a_vals):.2f} dBr | std: {np.std(cnr_a_vals):.2f} | rango: [{min(cnr_a_vals):.2f}, {max(cnr_a_vals):.2f}]")

        # Criterio de consistencia: std < 30% de la media absoluta
        consist_c = np.std(cnr_c_vals) < 0.3 * abs(np.mean(cnr_c_vals)) if np.mean(cnr_c_vals) != 0 else False
        consist_a = np.std(cnr_a_vals) < 0.3 * abs(np.mean(cnr_a_vals)) if np.mean(cnr_a_vals) != 0 else False

        print(f"\nVeredicto consistencia c: {'CONSISTENTE' if consist_c else 'DISPERSO'}")
        print(f"Veredicto consistencia a: {'CONSISTENTE' if consist_a else 'DISPERSO'}")

    #  Reporte Markdown 
    md = [
        "# Sweep de consistencia - Fantoma CIRS 040GSE (MYO001-MYO006)",
        "_Medical Insight 3D - Phase 2_\n",
        f"Margen de borde excluido: {MARGIN_FRAC*100:.0f}%\n",
        "## Resultados por dataset",
        "| Dataset | Disponible | CNR c (dBr) | CNR a (dBr) | Pico cerca borde (c) | Pico cerca borde (a) |",
        "|---|---|---|---|---|---|",
    ]
    for r in results:
        if r["available"]:
            md.append(
                f"| {r['dataset']} | si | {r['cnr_c_db']:.2f} | {r['cnr_a_db']:.2f} "
                f"| {'si' if r['peak_edge_c'] else 'no'} "
                f"| {'si' if r['peak_edge_a'] else 'no'} |"
            )
        else:
            md.append(f"| {r['dataset']} | no procesado | - | - | - | - |")

    if n_available >= 2:
        md += [
            "",
            "## Estadisticas de consistencia",
            f"- CNR c:  media = {np.mean(cnr_c_vals):.2f} dBr, std = {np.std(cnr_c_vals):.2f} dBr",
            f"- CNR a:  media = {np.mean(cnr_a_vals):.2f} dBr, std = {np.std(cnr_a_vals):.2f} dBr",
            f"- Veredicto c: {'CONSISTENTE (std < 30% de la media)' if consist_c else 'DISPERSO - investigar causa'}",
            f"- Veredicto a: {'CONSISTENTE (std < 30% de la media)' if consist_a else 'DISPERSO - investigar causa'}",
            "",
            "## Interpretacion",
            "Si ambos parametros son consistentes: el contraste fisico es reproducible",
            "across multiples adquisiciones del mismo fantoma, y el artefacto de borde",
            "(documentado en runs previos) puede tratarse como limitacion conocida sin",
            "bloquear la publicacion del resultado principal.",
            "",
            "Si alguno es disperso: revisar si los archivos 'pico cerca borde' coinciden",
            "con los outliers - eso indicaria que el artefacto de borde no es uniforme",
            "y SI requiere PML antes de reportar cualquier CNR.",
        ]

    out_md = ARTICULO_DIR / "sweep_consistencia_myo.md"
    out_md.write_text("\n".join(md), encoding="utf-8")
    print(f"\nReporte guardado: {out_md}")

    #  Grafico 
    if n_available >= 2:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            names = [r["dataset"] for r in available]

            for ax, vals, title, color in zip(
                axes, [cnr_c_vals, cnr_a_vals],
                ["CNR velocidad c", "CNR atenuacion alpha"],
                ["#185FA5", "#BA7517"]
            ):
                ax.bar(names, vals, color=color, alpha=0.7)
                ax.axhline(np.mean(vals), color="black", linestyle="--", linewidth=1, label=f"media={np.mean(vals):.2f}")
                ax.axhline(0, color="gray", linewidth=0.5)
                ax.fill_between([-0.5, len(names)-0.5],
                                np.mean(vals)-np.std(vals), np.mean(vals)+np.std(vals),
                                color=color, alpha=0.15, label=f"1 sigma={np.std(vals):.2f}")
                ax.set_title(title, fontweight="bold")
                ax.set_ylabel("CNR (dBr)")
                ax.legend(fontsize=9)
                ax.set_xlim(-0.5, len(names)-0.5)

            plt.suptitle("Sweep de consistencia - CIRS 040GSE (N={})".format(n_available), fontweight="bold")
            plt.tight_layout()
            out_png = ARTICULO_DIR / "sweep_consistencia_myo.png"
            plt.savefig(out_png, dpi=150, bbox_inches="tight")
            plt.close()
            print(f"Figura guardada: {out_png}")
        except ImportError:
            print("matplotlib no disponible - figura omitida")


if __name__ == "__main__":
    main()
