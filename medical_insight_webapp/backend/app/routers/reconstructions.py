"""
backend/app/routers/reconstructions.py
========================================
Endpoints de solo lectura sobre reconstrucciones ya calculadas.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from app.models.schemas import DatasetSummary, ReconstructionMaps, ConvergenceHistory
from app.services.reconstruction_store import ReconstructionStore
from app.config import get_settings

router = APIRouter(prefix="/api/reconstructions", tags=["reconstructions"])


def _get_store() -> ReconstructionStore:
    settings = get_settings()
    return ReconstructionStore(settings.articulo_dir)


@router.get("/", response_model=list[DatasetSummary])
def list_datasets():
    """Lista todas las adquisiciones con reconstrucción disponible."""
    return _get_store().list_datasets()


@router.get("/{dataset_id}/maps", response_model=ReconstructionMaps)
def get_maps(
    dataset_id: str,
    resolution: int = Query(128, ge=32, le=512, description="Resolución máxima por eje para el navegador"),
):
    """
    Devuelve los mapas c(x,z) y alpha(x,z) sub-muestreados, listos para
    dibujar como heatmap en el frontend, junto con la localización del
    tumor y las métricas de contraste (CNR).
    """
    try:
        return _get_store().get_reconstruction(dataset_id, max_resolution=resolution)
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.get("/{dataset_id}/convergence", response_model=ConvergenceHistory)
def get_convergence(dataset_id: str):
    """Historial de pérdida y error por época, para el gráfico de convergencia."""
    try:
        return _get_store().get_convergence(dataset_id)
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
