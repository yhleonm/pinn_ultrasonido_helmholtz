"""
backend/app/models/schemas.py
=============================
Esquemas Pydantic para las respuestas de la API.

Filosofía de diseño: esta API NUNCA ejecuta un entrenamiento FWI dentro
de un request HTTP (eso tomaría 20-60 minutos). Sirve resultados ya
calculados y guardados en articulo/*.npy, o encola nuevos trabajos que
corren en un proceso de fondo separado (ver services/job_queue.py).
"""

from __future__ import annotations
from datetime import datetime
from enum import Enum
from typing import Optional
from pydantic import BaseModel, Field


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"


class TissueType(str, Enum):
    PHANTOM = "phantom"       # fantoma físico (CIRS, ATS, etc.)
    IN_VIVO = "in_vivo"       # tejido humano real


class DatasetSummary(BaseModel):
    """Resumen de una adquisición disponible para visualizar."""
    dataset_id: str = Field(..., examples=["MYO001"])
    tissue_type: TissueType
    institution: Optional[str] = Field(None, examples=["CIRS 040GSE"])
    status: JobStatus
    tumor_x_mm: Optional[float] = None
    tumor_z_mm: Optional[float] = None
    cnr_c_dbr: Optional[float] = None
    cnr_alpha_dbr: Optional[float] = None
    is_reliable_velocity: bool = Field(
        True,
        description="False si el pico de velocidad está fuera del radio "
                    "físico esperado del tumor (artefacto de borde conocido)."
    )
    created_at: Optional[datetime] = None


class ReliabilityLevel(str, Enum):
    HIGH = "alta"
    MEDIUM = "media"
    LOW = "baja"


class ClinicalInterpretation(BaseModel):
    finding_summary: str
    reliability_level: ReliabilityLevel
    reliability_reason: str
    reliability_reason_velocity: Optional[str] = None
    reliability_reason_attenuation: Optional[str] = None
    velocity_reliable: bool
    attenuation_reliable: bool


class ReconstructionMaps(BaseModel):
    """
    Payload completo de una reconstrucción: los mapas c(x,z) y alpha(x,z)
    ya sub-muestreados a una resolución razonable para el navegador,
    más metadatos físicos para dibujar correctamente los ejes.
    """
    dataset_id: str
    grid_shape: tuple[int, int]                  # (Nz, Nx) tras submuestreo
    x_min_mm: float
    x_max_mm: float
    z_min_mm: float
    z_max_mm: float
    c_map: list[list[float]]                     # m/s
    alpha_map: list[list[float]]                  # dB/MHz/cm (ya convertido)
    bmode_map: Optional[list[list[float]]] = None # B-mode en escala logarítmica (dB)
    clinical_interpretation: Optional[ClinicalInterpretation] = None
    tumor_x_mm: Optional[float]
    tumor_z_mm: Optional[float]
    tumor_radius_mm: float = 3.0
    cnr_c_dbr: Optional[float]
    cnr_alpha_dbr: Optional[float]
    peak_distance_mm: Optional[float] = Field(
        None, description="Distancia del pico reconstruido al centro DAS del tumor."
    )
    peak_distance_alpha_mm: Optional[float] = Field(
        None, description="Distancia del pico de atenuacion reconstruido al centro DAS."
    )
    reliability_note: Optional[str] = Field(
        None,
        description="Nota clínica/técnica si la reconstrucción tiene "
                    "limitaciones conocidas (ej. artefacto de borde)."
    )


class ConvergenceHistory(BaseModel):
    dataset_id: str
    epochs: list[int]
    loss: list[float]
    err_c: list[float]
    err_alpha: list[float]
    stage_boundaries: list[int] = Field(
        default_factory=list,
        description="Épocas donde cambia la etapa de frecuencia (1.0/2.0/3.5 MHz)."
    )


class NewJobRequest(BaseModel):
    """Solicitud para encolar un nuevo cálculo FWI sobre un archivo subido."""
    dataset_id: str
    tissue_type: TissueType
    notes: Optional[str] = None


class NewJobResponse(BaseModel):
    job_id: str
    dataset_id: str
    status: JobStatus
    estimated_minutes: int = Field(
        30, description="Estimación conservadora basada en el tiempo histórico "
                        "por etapa de frecuencia con el solver optimizado."
    )
