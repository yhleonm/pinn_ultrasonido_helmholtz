agu

# Medical Insight 3D — Aplicación Web (Visor de Reconstrucciones FWI)

## Decisión de arquitectura central

**Ningún request HTTP ejecuta un entrenamiento FWI.** Un entrenamiento real
toma entre 20 y 60 minutos incluso con el solver optimizado (factorización
LU compartida forward/backward). Por eso la app tiene dos modos separados:

1. **Visualización (síncrona, rápida)** — lee mapas `.npy` ya calculados
   desde `articulo/` y los sirve al navegador en milisegundos.
2. **Nuevo cálculo (asíncrono, en segundo plano)** — el usuario sube un
   `.h5`/`.mat`, el backend lanza `main_5MHz_fwi_real_data.py` como
   subproceso, devuelve un `job_id` de inmediato, y el frontend hace
   polling cada 4 segundos para mostrar el progreso.

Esto es exactamente el mismo patrón que un visor PACS médico: no
recalculas la resonancia magnética cada vez que un doctor la abre — la
ves ya reconstruida, y encolas estudios nuevos por separado.

---

## Cómo conectar esto a tu proyecto real

### 1. Backend

```bash
cd backend
pip install fastapi uvicorn python-multipart numpy

# Apunta al workspace real donde vive articulo/, main_5MHz_fwi_real_data.py, etc.
export MEDICAL_INSIGHT_WORKSPACE="/ruta/a/proyecto_pinn_ultrasonido_mama_fase_2"
# En Windows (PowerShell):
#   $env:MEDICAL_INSIGHT_WORKSPACE="C:\Users\yelem\workspace\proyecto_pinn_ultrasonido_mama_fase_2"

python -m uvicorn app.main:app --reload --port 8000
```

### 2. Frontend

```bash
cd frontend
npm install
npm run dev
```

Abre `http://localhost:5173`.

---

## Piezas que debes adaptar a tu código real

### A. Sidecar de métricas (`metrics_{dataset_id}.json`)

Ahora mismo `reconstruction_store.py` intenta leer un archivo
`articulo/metrics_MYO001.json` con esta forma:

```json
{
  "tumor_x_mm": 2.0,
  "tumor_z_mm": 18.5,
  "cnr_c_dbr": 3.23,
  "cnr_alpha_dbr": 6.02,
  "peak_distance_mm": 5.56,
  "x_min_mm": -19.0, "x_max_mm": 19.0,
  "z_min_mm": 0.0, "z_max_mm": 32.0
}
```

Si no existe, la API sigue funcionando pero sin esos metadatos (mapas
sin overlay de tumor, sin badge de fiabilidad). **Recomendación:** al
final de `calcular_cnr.py` o `sweep_consistencia_myo.py`, agrega una
función que escriba este JSON junto a cada `.npy` — es la forma más
robusta de conectar la investigación con la app sin re-parsear Markdown.

### B. Interfaz de línea de comandos de `main_5MHz_fwi_real_data.py`

`job_queue.py` asume que el script acepta:

```
python main_5MHz_fwi_real_data.py --input <ruta.h5> --dataset-id <ID> --output-dir <articulo/>
```

Si tu script usa otros nombres de argumento (ya vimos en esta conversación
que en algún punto se agregó `--sufijo` y semilla por CLI), ajusta la
lista de argumentos en `job_queue.py::_run_job`.

### C. Conversión de unidades de atenuación

`reconstruction_store.py` tiene una heurística (`if a_arr.max() < 5.0`)
para detectar si alpha está en Np/m o ya en dB/MHz/cm. Esto es frágil —
mejor que el sidecar JSON incluya explícitamente `"alpha_unit": "Np/m"`
y elimines la heurística.

---

## Qué NO hace esta versión (deliberadamente, para no sobre-construir)

- No tiene autenticación — apropiado para un solo investigador/laboratorio,
  no para producción con múltiples usuarios sin control de acceso.
- La cola de trabajos es en memoria (un diccionario Python) — se pierde
  al reiniciar el servidor. Para uso real de laboratorio con múltiples
  personas, migrar a Celery + Redis o RQ.
- No valida que el archivo subido sea físicamente compatible (misma
  frecuencia, mismo probe) antes de lanzar el entrenamiento — el propio
  script de entrenamiento debe fallar con un mensaje claro si el HDF5
  no tiene el formato esperado.
- Nada de esto reemplaza la validación clínica, ética o regulatoria
  necesaria para uso diagnóstico real — es una herramienta de
  visualización de investigación.
