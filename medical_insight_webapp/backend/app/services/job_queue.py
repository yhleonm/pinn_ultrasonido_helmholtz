"""
backend/app/services/job_queue.py
===================================
Cola de trabajos en segundo plano para lanzar nuevas reconstrucciones FWI.

Decisión de arquitectura clave: un entrenamiento FWI real toma 20-60
minutos (verificado empíricamente con el solver optimizado LU compartido,
ver metricas_progreso_fwi_real_cubdl_*.md). Ningún endpoint HTTP puede
bloquear ese tiempo. Este módulo lanza el entrenamiento como un
subproceso independiente (main_5MHz_fwi_real_data.py) y el estado se
consulta por polling desde el frontend, no por WebSocket bloqueante
ni por response síncrona.

Para producción real (múltiples médicos, múltiples GPUs), esto debería
migrar a una cola real (Celery + Redis, o RQ) — este diseño en memoria
es apropiado para un solo investigador/laboratorio con una GPU compartida,
que es el caso de uso actual del proyecto.
"""

from __future__ import annotations

import subprocess
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

from app.models.schemas import JobStatus


@dataclass
class Job:
    job_id: str
    dataset_id: str
    status: JobStatus
    created_at: datetime
    log_path: Path
    process: Optional[subprocess.Popen] = None
    error: Optional[str] = None


class JobQueue:
    """
    Cola en memoria de un solo worker (una GPU, un entrenamiento a la vez).

    IMPORTANTE: solo permite UN trabajo corriendo simultáneamente, porque
    el solver Helmholtz en CPU (spsolve/splu) satura los núcleos disponibles
    y correr dos entrenamientos en paralelo degradaría a ambos (ver la
    discusión sobre por qué el sweep factorial corre secuencialmente,
    no en paralelo).
    """

    def __init__(self, script_path: Path, logs_dir: Path, data_dir: Path):
        self.script_path = script_path
        self.logs_dir = logs_dir
        self.data_dir = data_dir
        self.logs_dir.mkdir(parents=True, exist_ok=True)

        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()
        self._current_job_id: Optional[str] = None

    # ------------------------------------------------------------------
    def enqueue(self, dataset_id: str, hdf5_path: Path) -> Job:
        job_id = str(uuid.uuid4())[:8]
        log_path = self.logs_dir / f"job_{job_id}_{dataset_id}.log"

        job = Job(
            job_id=job_id,
            dataset_id=dataset_id,
            status=JobStatus.QUEUED,
            created_at=datetime.now(),
            log_path=log_path,
        )
        with self._lock:
            self._jobs[job_id] = job

        # Lanzar en un hilo separado que espera su turno si hay otro corriendo
        thread = threading.Thread(
            target=self._run_job, args=(job, hdf5_path), daemon=True
        )
        thread.start()
        return job

    # ------------------------------------------------------------------
    def _run_job(self, job: Job, hdf5_path: Path) -> None:
        # Esperar turno: solo un entrenamiento a la vez sobre el solver CPU
        while True:
            with self._lock:
                if self._current_job_id is None:
                    self._current_job_id = job.job_id
                    break
            threading.Event().wait(2.0)

        job.status = JobStatus.RUNNING
        try:
            with open(job.log_path, "w") as log_file:
                job.process = subprocess.Popen(
                    [
                        "python", "-u", str(self.script_path),
                        str(hdf5_path),  # Path al archivo HDF5 (posicional)
                        "42",            # Semilla (posicional)
                        "",              # Sufijo (posicional, vacío por defecto)
                    ],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=str(self.script_path.parent),
                )
                job.process.wait()

            if job.process.returncode == 0:
                job.status = JobStatus.DONE
            else:
                job.status = JobStatus.FAILED
                job.error = f"El proceso terminó con código {job.process.returncode}. Revisa {job.log_path.name}."

        except Exception as e:
            job.status = JobStatus.FAILED
            job.error = str(e)

        finally:
            with self._lock:
                if self._current_job_id == job.job_id:
                    self._current_job_id = None

    # ------------------------------------------------------------------
    def get_job(self, job_id: str) -> Optional[Job]:
        return self._jobs.get(job_id)

    def get_log_tail(self, job_id: str, n_lines: int = 50) -> list[str]:
        job = self._jobs.get(job_id)
        if job is None or not job.log_path.exists():
            return []
        lines = job.log_path.read_text(errors="ignore").splitlines()
        return lines[-n_lines:]

    def queue_position(self, job_id: str) -> int:
        """0 si está corriendo, N si hay N trabajos delante en cola."""
        with self._lock:
            if self._current_job_id == job_id:
                return 0
            queued_before = [
                j for j in self._jobs.values()
                if j.status == JobStatus.QUEUED
                and j.created_at < self._jobs[job_id].created_at
            ]
            return len(queued_before) + (1 if self._current_job_id else 0)
