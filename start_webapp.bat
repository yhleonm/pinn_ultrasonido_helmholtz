@echo off
title Medical Insight 3D - WebApp Runner
cd /d "%~dp0"

echo ==========================================================
echo   Iniciando Medical Insight 3D (Frontend + Backend)
echo ==========================================================

REM Configurar ruta base del workspace
set MEDICAL_INSIGHT_WORKSPACE=%~dp0

echo [1/2] Iniciando Backend FastAPI en http://localhost:8000 ...
start "Medical Insight Backend" cmd /k "cd /d "%~dp0medical_insight_webapp\backend" && python -m uvicorn app.main:app --reload --port 8000"

echo [2/2] Iniciando Frontend Vite en http://localhost:5173 ...
start "Medical Insight Frontend" cmd /k "cd /d "%~dp0medical_insight_webapp\frontend" && npm run dev"

echo.
echo ==========================================================
echo   Servicios lanzados exitosamente:
echo   - Backend API: http://localhost:8000 (Docs: http://localhost:8000/docs)
echo   - Frontend UI: http://localhost:5173
echo ==========================================================
pause
