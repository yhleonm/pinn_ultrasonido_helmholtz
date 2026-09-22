# Prototipos y Experimentos Históricos (Legacy)

Este directorio almacena las implementaciones de etapas anteriores del proyecto, preservadas con fines de trazabilidad metodológica:

* `main_1d.py`: Prueba de concepto inicial 1D de PINN para Helmholtz con atenuación.
* `main_2d.py`: Extensión a 2D sobre dominio sintético con condición de Sommerfeld simple.
* `main_5MHz.py`: Exploración inicial a frecuencia clínica de 5 MHz.
* `main_5MHz_planB.py`: Variante alternativa de formulación de funciones de pérdida.
* `main_5MHz_roi.py`: Prototipo de acotamiento de la región de interés (ROI).

### Pipeline Actual de Producción:
El pipeline principal para datos experimentales reales de CUBDL es:
* `main_5MHz_fwi_real_data.py`: Inversión de forma de onda completa (FWI) multietapa con capas absorbentes (PML), formulación adyunta acelerada con sparse LU (`splu`) y exportación de sidecars de métricas para `medical_insight_webapp`.
