"""
backend/app/services/reconstruction_store.py
==============================================
Lee los resultados ya calculados (c_map_*.npy, alpha_map_*.npy, y los
reportes .md de CNR/distancia radial) desde la carpeta articulo/ del
proyecto de investigación, y los expone en un formato listo para la API.

IMPORTANTE — decisión de arquitectura:
Esta capa NUNCA dispara un entrenamiento. Solo lee archivos que ya
existen en disco, generados por main_5MHz_fwi_real_data.py corriendo
por separado (línea de comandos, o vía el job_queue en modo background).
Esto es intencional: un entrenamiento toma 20-60 min incluso optimizado,
y ningún framework web debe bloquear un request HTTP por ese tiempo.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

import numpy as np

from app.models.schemas import (
    DatasetSummary, ReconstructionMaps, ConvergenceHistory,
    JobStatus, TissueType, ClinicalInterpretation, ReliabilityLevel,
)


# ---------------------------------------------------------------------------
# Conocimiento físico conocido por institución de fantoma
# (mismo mapeo que se usó en calcular_cnr.py y el sweep de consistencia)
# ---------------------------------------------------------------------------

INSTITUTION_MAP = {
    "MYO": "CIRS 040GSE",
    "OSL": "CIRS 054GS",
    "EUT": "CIRS 050 / ATS549",
    "INS": "Variado",
    "TSH": None,   # in vivo, sin institución de fantoma
}

TUMOR_RADIUS_MM = 3.0   # radio físico conocido de la inclusión CIRS 040GSE


def _dataset_prefix(dataset_id: str) -> str:
    match = re.match(r"([A-Za-z]+)", dataset_id)
    return match.group(1).upper() if match else ""


def _tissue_type(dataset_id: str) -> TissueType:
    prefix = _dataset_prefix(dataset_id)
    return TissueType.IN_VIVO if prefix in ("TSH", "OASBUD") else TissueType.PHANTOM


class ReconstructionStore:
    """
    Punto único de acceso a los resultados guardados en articulo/.

    Args:
        articulo_dir: carpeta donde viven los .npy y .md generados
                      por el pipeline de investigación (train.py /
                      main_5MHz_fwi_real_data.py + calcular_cnr.py).
    """

    def __init__(self, articulo_dir: Path):
        self.dir = Path(articulo_dir)
        if not self.dir.exists():
            raise FileNotFoundError(
                f"No se encontró el directorio de resultados: {self.dir}. "
                f"Verifica ARTICULO_DIR en la configuración del backend."
            )

    # ------------------------------------------------------------------
    def list_datasets(self) -> list[DatasetSummary]:
        """
        Descubre todos los datasets con al menos un c_map_*.npy disponible.
        """
        summaries = []
        for c_map_path in sorted(self.dir.glob("c_map_*.npy")):
            dataset_id = c_map_path.stem.replace("c_map_", "")
            summaries.append(self._build_summary(dataset_id))
        return summaries

    # ------------------------------------------------------------------
    def _build_summary(self, dataset_id: str) -> DatasetSummary:
        metrics = self._load_cnr_metrics(dataset_id)
        prefix = _dataset_prefix(dataset_id)
        eff_radius = float(metrics.get("equiv_radius_mm", TUMOR_RADIUS_MM))

        return DatasetSummary(
            dataset_id=dataset_id,
            tissue_type=_tissue_type(dataset_id),
            institution=metrics.get("clinical_label") or INSTITUTION_MAP.get(prefix),
            status=JobStatus.DONE,
            tumor_x_mm=metrics.get("tumor_x_mm"),
            tumor_z_mm=metrics.get("tumor_z_mm"),
            cnr_c_dbr=metrics.get("cnr_c_dbr"),
            cnr_alpha_dbr=metrics.get("cnr_alpha_dbr"),
            is_reliable_velocity=metrics.get("peak_distance_mm", 99.0) <= eff_radius,
        )

    # ------------------------------------------------------------------
    def _load_cnr_metrics(self, dataset_id: str) -> dict:
        """
        Intenta leer un sidecar JSON de métricas si existe
        (recomendado: guardar esto junto a cada .npy al final del
        entrenamiento, ver services/job_runner.py::_save_metrics_sidecar).
        Si no existe, devuelve valores vacíos sin fallar.
        """
        sidecar = self.dir / f"metrics_{dataset_id}.json"
        if sidecar.exists():
            try:
                return json.loads(sidecar.read_text())
            except json.JSONDecodeError:
                pass
        return {}

    # ------------------------------------------------------------------
    def get_reconstruction(
        self,
        dataset_id: str,
        max_resolution: int = 200,
    ) -> ReconstructionMaps:
        """
        Carga c_map y alpha_map desde disco, los sub-muestrea a una
        resolución razonable para el navegador (los mapas reales pueden
        tener 600k+ puntos, demasiado para renderizar o transferir tal cual),
        y arma el payload completo con metadatos físicos.
        """
        c_path = self.dir / f"c_map_{dataset_id}.npy"
        a_path = self.dir / f"alpha_map_{dataset_id}.npy"

        if not c_path.exists() or not a_path.exists():
            raise FileNotFoundError(f"Reconstrucción no encontrada para {dataset_id}")

        c_arr = np.load(c_path).astype(np.float64)
        a_arr = np.load(a_path).astype(np.float64)

        metrics = self._load_cnr_metrics(dataset_id)
        tumor_x = metrics.get("tumor_x_mm")
        tumor_z = metrics.get("tumor_z_mm")
        
        # Obtener los límites del grid completo
        x_min_mm = metrics.get("x_min_mm", -19.05)
        x_max_mm = metrics.get("x_max_mm", 19.05)
        z_min_mm = metrics.get("z_min_mm", 0.0)
        z_max_mm = metrics.get("z_max_mm", 32.0)
        
        # Recorte de la ROI (zoom automático solo para fantomas sintéticos/CUBDL;
        # en OASBUD in-vivo se preserva el mapa de campo completo continuo)
        is_full_field = dataset_id.startswith("OASBUD_") or bool(metrics.get("full_field_view", False))
        if not is_full_field and tumor_x is not None and tumor_z is not None:
            # ROI de 1.5 cm x 1.5 cm (± 7.5 mm)
            roi_half = float(metrics.get("roi_half_mm", 7.5))
            x_min_roi = max(x_min_mm, tumor_x - roi_half)
            x_max_roi = min(x_max_mm, tumor_x + roi_half)
            z_min_roi = max(z_min_mm, tumor_z - roi_half)
            z_max_roi = min(z_max_mm, tumor_z + roi_half)
            
            Nz, Nx = c_arr.shape
            x_vals = np.linspace(x_min_mm, x_max_mm, Nx)
            z_vals = np.linspace(z_min_mm, z_max_mm, Nz)
            
            x_idx = np.where((x_vals >= x_min_roi) & (x_vals <= x_max_roi))[0]
            z_idx = np.where((z_vals >= z_min_roi) & (z_vals <= z_max_roi))[0]
            
            if len(x_idx) > 0 and len(z_idx) > 0:
                c_arr = c_arr[z_idx.min():z_idx.max()+1, x_idx.min():x_idx.max()+1]
                a_arr = a_arr[z_idx.min():z_idx.max()+1, x_idx.min():x_idx.max()+1]
                
                # Actualizar límites de la respuesta a la ROI recortada
                x_min_mm = float(x_vals[x_idx.min()])
                x_max_mm = float(x_vals[x_idx.max()])
                z_min_mm = float(z_vals[z_idx.min()])
                z_max_mm = float(z_vals[z_idx.max()])

        alpha_unit = metrics.get("alpha_unit", "Np/m")

        # Convertir alpha de Np/m a dB/MHz/cm si el sidecar indica Np/m
        if alpha_unit == "Np/m":
            a_arr_db = a_arr * 8.686e4  # Np/m -> dB/MHz/cm
        else:
            a_arr_db = a_arr

        c_small = self._downsample(c_arr, max_resolution)
        a_small = self._downsample(a_arr_db, max_resolution)

        # Cargar B-mode if available
        b_path = self.dir / f"bmode_{dataset_id}.npy"
        bmode_small = None
        if b_path.exists():
            try:
                b_arr = np.load(b_path).astype(np.float64)
                if tumor_x is not None and tumor_z is not None and 'z_idx' in locals() and len(z_idx) > 0 and len(x_idx) > 0:
                    b_arr = b_arr[z_idx.min():z_idx.max()+1, x_idx.min():x_idx.max()+1]
                bmode_small = self._downsample(b_arr, max_resolution).tolist()
            except Exception as e:
                # Si falla por alguna razón (por ejemplo, archivo incompleto), registrar y continuar
                print(f"Error al cargar B-mode para {dataset_id}: {e}")

        peak_dist = metrics.get("peak_distance_mm")
        peak_dist_alpha = metrics.get("peak_distance_alpha_mm")

        # Generar interpretación clínica con el radio equivalente real de la lesión (OASBUD) o del fantoma (3.0 mm)
        eff_radius_mm = max(TUMOR_RADIUS_MM, float(metrics.get("equiv_radius_mm", TUMOR_RADIUS_MM)))
        from app.services.interpretation import build_interpretation
        clinical_interpretation = build_interpretation(
            tumor_x_mm=tumor_x,
            tumor_z_mm=tumor_z,
            cnr_c_dbr=metrics.get("cnr_c_dbr"),
            cnr_alpha_dbr=metrics.get("cnr_alpha_dbr"),
            peak_distance_mm=peak_dist,
            peak_distance_alpha_mm=peak_dist_alpha,
            expected_radius_mm=eff_radius_mm,
        )

        reliability_note = clinical_interpretation.reliability_reason if clinical_interpretation.reliability_level != ReliabilityLevel.HIGH else None

        return ReconstructionMaps(
            dataset_id=dataset_id,
            grid_shape=(c_small.shape[0], c_small.shape[1]),
            x_min_mm=x_min_mm,
            x_max_mm=x_max_mm,
            z_min_mm=z_min_mm,
            z_max_mm=z_max_mm,
            c_map=c_small.tolist(),
            alpha_map=a_small.tolist(),
            bmode_map=bmode_small,
            clinical_interpretation=clinical_interpretation,
            tumor_x_mm=tumor_x,
            tumor_z_mm=tumor_z,
            tumor_radius_mm=eff_radius_mm,
            cnr_c_dbr=metrics.get("cnr_c_dbr"),
            cnr_alpha_dbr=metrics.get("cnr_alpha_dbr"),
            peak_distance_mm=peak_dist,
            peak_distance_alpha_mm=peak_dist_alpha,
            roi_contour_mm=metrics.get("roi_contour_mm"),
            pred_contour_mm=metrics.get("pred_contour_mm"),
            clinical_label=metrics.get("clinical_label"),
            dice_c=metrics.get("dice_c"),
            dice_alpha=metrics.get("dice_alpha"),
            c_roi_mean_ms=metrics.get("c_roi_mean_ms"),
            alpha_roi_mean_db=metrics.get("alpha_roi_mean_db"),
            reliability_note=reliability_note,
        )

    # ------------------------------------------------------------------
    def get_convergence(self, dataset_id: str) -> ConvergenceHistory:
        """
        Lee el historial de convergencia desde metricas_progreso_*.md
        (parseo simple de la tabla Markdown) o desde un sidecar JSON
        si está disponible (preferido, más robusto).
        """
        json_path = self.dir / f"convergence_{dataset_id}.json"
        if json_path.exists():
            data = json.loads(json_path.read_text())
            return ConvergenceHistory(dataset_id=dataset_id, **data)

        md_path = self.dir / f"metricas_progreso_fwi_real_cubdl_{dataset_id}.md"
        if not md_path.exists():
            raise FileNotFoundError(f"No hay historial de convergencia para {dataset_id}")

        epochs, losses, errs_c, errs_a, stages = self._parse_convergence_md(md_path)
        return ConvergenceHistory(
            dataset_id=dataset_id,
            epochs=epochs, loss=losses, err_c=errs_c, err_alpha=errs_a,
            stage_boundaries=stages,
        )

    # ------------------------------------------------------------------
    @staticmethod
    def _parse_convergence_md(path: Path):
        """Parsea la tabla Markdown '| Época | Etapa | Frecuencia | Pérdida ...'."""
        epochs, losses, errs_c, errs_a, stage_col = [], [], [], [], []
        for line in path.read_text().splitlines():
            if not line.startswith("|") or "Época" in line or ":---" in line:
                continue
            cols = [c.strip() for c in line.strip("|").split("|")]
            if len(cols) < 6:
                continue
            try:
                epochs.append(int(cols[0]))
                stage_col.append(int(cols[1]))
                losses.append(float(cols[3].replace(",", "").replace("e", "e")))
                errs_c.append(float(cols[4].replace("%", "")))
                errs_a.append(float(cols[5].replace("%", "")))
            except (ValueError, IndexError):
                continue

        # Detectar cambios de etapa para marcar límites en el gráfico
        stage_boundaries = [
            epochs[i] for i in range(1, len(stage_col))
            if stage_col[i] != stage_col[i-1]
        ]
        return epochs, losses, errs_c, errs_a, stage_boundaries

    # ------------------------------------------------------------------
    @staticmethod
    def _downsample(arr: np.ndarray, max_res: int) -> np.ndarray:
        """Sub-muestreo simple por decimación para no transferir 600k+ puntos."""
        Nz, Nx = arr.shape
        step_z = max(1, Nz // max_res)
        step_x = max(1, Nx // max_res)
        return arr[::step_z, ::step_x]
