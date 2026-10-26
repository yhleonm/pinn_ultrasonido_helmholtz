# PowerShell Script para lanzar Frontend y Backend simultáneamente
$workspace = Split-Path -Parent $MyInvocation.MyCommand.Path
$env:MEDICAL_INSIGHT_WORKSPACE = $workspace

Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host "  Iniciando Medical Insight 3D (Frontend + Backend)" -ForegroundColor Cyan
Write-Host "==========================================================" -ForegroundColor Cyan

# Lanzar Backend FastAPI
$backendPath = Join-Path $workspace "medical_insight_webapp\backend"
Start-Process powershell -ArgumentList "-NoExit", "-Command", "`$env:MEDICAL_INSIGHT_WORKSPACE='$workspace'; cd '$backendPath'; Write-Host 'Iniciando Backend en http://localhost:8000...' -ForegroundColor Green; python -m uvicorn app.main:app --reload --port 8000"

# Lanzar Frontend React Vite
$frontendPath = Join-Path $workspace "medical_insight_webapp\frontend"
Start-Process powershell -ArgumentList "-NoExit", "-Command", "cd '$frontendPath'; Write-Host 'Iniciando Frontend en http://localhost:5173...' -ForegroundColor Green; npm run dev"

Write-Host "`nServicios en ejecucion:" -ForegroundColor Green
Write-Host "  - Backend API: http://localhost:8000 (Swagger: http://localhost:8000/docs)" -ForegroundColor Yellow
Write-Host "  - Frontend UI: http://localhost:5173" -ForegroundColor Yellow
Write-Host "==========================================================" -ForegroundColor Cyan
