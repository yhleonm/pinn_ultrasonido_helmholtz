"""
backend/app/config.py
=======================
Configuración centralizada. Ajusta ARTICULO_DIR, DATA_DIR y SCRIPT_PATH
para que apunten a tu workspace real de investigación
(ej. C:\\Users\\yelem\\workspace\\proyecto_pinn_ultrasonido_mama_fase_2).
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

from app.services.job_queue import JobQueue


class Settings:
    def __init__(self):
        # --- AJUSTA ESTAS RUTAS A TU WORKSPACE REAL ---
        base = Path(os.environ.get(
            "MEDICAL_INSIGHT_WORKSPACE",
            r"C:\Users\yelem\workspace\proyecto_pinn_ultrasonido_mama_fase_2"
        ))

        self.articulo_dir = base / "articulo"
        self.uploads_dir = base / "data_raw" / "uploads"
        self.logs_dir = base / "logs" / "webapp_jobs"
        self.script_path = base / "main_5MHz_fwi_real_data.py"

        self.job_queue = JobQueue(
            script_path=self.script_path,
            logs_dir=self.logs_dir,
            data_dir=self.articulo_dir,
        )


@lru_cache
def get_settings() -> Settings:
    return Settings()
