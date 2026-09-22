"""
backend/app/routers/jobs.py
=============================
Endpoints para encolar un nuevo cálculo FWI a partir de un archivo
.h5/.mat subido, y consultar su progreso.

Ningún endpoint aquí es síncrono/bloqueante: enqueue_job devuelve
inmediatamente un job_id, y el cómputo real corre en un subproceso
separado (ver services/job_queue.py). El frontend hace polling sobre
GET /api/jobs/{job_id} cada pocos segundos.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from fastapi import APIRouter, HTTPException, UploadFile, File, Form

from app.models.schemas import NewJobResponse, JobStatus
from app.config import get_settings

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


@router.post("/", response_model=NewJobResponse)
async def enqueue_job(
    dataset_id: str = Form(...),
    file: UploadFile = File(...),
):
    """
    Sube un archivo de canal RF (.h5/.hdf5/.mat) y encola un nuevo
    cálculo FWI. Devuelve inmediatamente un job_id para hacer seguimiento;
    el entrenamiento real toma 20-60 minutos y corre en segundo plano.
    """
    settings = get_settings()

    if not file.filename.lower().endswith((".h5", ".hdf5", ".mat")):
        raise HTTPException(400, "Formato no soportado. Usa .h5, .hdf5 o .mat")

    uploads_dir = settings.uploads_dir
    uploads_dir.mkdir(parents=True, exist_ok=True)
    dest_path = uploads_dir / file.filename

    with open(dest_path, "wb") as f:
        shutil.copyfileobj(file.file, f)

    job = settings.job_queue.enqueue(dataset_id=dataset_id, hdf5_path=dest_path)

    return NewJobResponse(
        job_id=job.job_id,
        dataset_id=job.dataset_id,
        status=job.status,
        estimated_minutes=30,  # basado en el tiempo histórico con solver optimizado
    )


@router.get("/{job_id}")
def get_job_status(job_id: str):
    settings = get_settings()
    job = settings.job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(404, "Trabajo no encontrado")

    return {
        "job_id": job.job_id,
        "dataset_id": job.dataset_id,
        "status": job.status,
        "queue_position": settings.job_queue.queue_position(job_id),
        "error": job.error,
        "log_tail": settings.job_queue.get_log_tail(job_id, n_lines=20),
    }
