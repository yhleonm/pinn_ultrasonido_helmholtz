"""
backend/app/main.py
=====================
Punto de entrada de la API. Sigue el mismo patrón que ya usaste en
Fase 1 y en el simulador SPR-LMR (FastAPI + uvicorn).

Ejecutar:
    conda activate pinn_med_env   (o el env que uses para Fase 2)
    pip install fastapi uvicorn python-multipart
    export MEDICAL_INSIGHT_WORKSPACE=/ruta/a/tu/proyecto_pinn_ultrasonido_mama_fase_2
    python -m uvicorn backend.app.main:app --reload --port 8000
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.routers import reconstructions, jobs

app = FastAPI(
    title="Medical Insight 3D — API de Reconstrucciones FWI",
    description=(
        "API de solo-lectura sobre reconstrucciones FWI-PINN ya calculadas, "
        "más un endpoint de encolado para nuevos cálculos en segundo plano. "
        "Ningún endpoint bloquea esperando un entrenamiento (20-60 min)."
    ),
    version="0.1.0",
)

# CORS abierto para desarrollo local — restringir en despliegue real
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(reconstructions.router)
app.include_router(jobs.router)


@app.get("/api/health")
def health():
    return {"status": "ok", "service": "medical-insight-3d-api"}
