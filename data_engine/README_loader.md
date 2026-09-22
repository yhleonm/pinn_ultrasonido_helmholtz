# CUBDLLoader: Cargador y Preprocesador de Datos de Canal Reales (RF)

Este módulo provee la clase `CUBDLLoader` para la ingesta y preparación de señales de radiofrecuencia (RF) crudas provenientes del desafío **CUBDL (Challenge on Ultrasound Beamforming with Deep Learning)** de la Universidad Johns Hopkins.

---

## 1. Funcionamiento y Calibración Física

El loader transforma las señales reales registradas por los sensores en el tiempo en presiones complejas monocromáticas en el dominio de la frecuencia para alimentar las condiciones de Dirichlet del modelo PINN:

1.  **Construcción de la Señal Analítica**:
    Aplica la Transformada de Hilbert sobre los datos de RF temporales $s(t)$ de cada canal para construir la señal analítica compleja:
    $$s_a(t) = s(t) + i \mathcal{H}\{s(t)\}$$
    Esto remueve redundancias y provee una estimación instantánea de fase y envolvente.
2.  **Transformada de Fourier (FFT)**:
    Se realiza una FFT a lo largo de la dimensión temporal para extraer la respuesta en las frecuencias seleccionadas $f_k$.
3.  **Compensación de Retardo de Disparo ($\tau_{\text{tx}}$)**:
    En adquisiciones reales con apertura focalizada o dirigida (phased/linear arrays), cada elemento emisor se dispara con un retardo temporal $\tau_{\text{tx}, m}$ específico. Para alinear el frente de fase incidente y obtener la presión compleja coherente en los receptores, se aplica una compensación de fase:
    $$P_{\text{calibrated}}(x_m, 0, \omega_k) = P_{\text{meas}}(x_m, 0, \omega_k) \cdot e^{i \omega_k (\tau_{\text{tx}, m} - t_0)}$$
    Donde $t_0$ es el tiempo inicial de muestreo del digitalizador.

---

## 2. Ejemplo de Uso

```python
from data_engine.loader import CUBDLLoader

# 1. Inicializar cargador (detecta automáticamente formatos .mat y .h5 de v7.3)
loader = CUBDLLoader("ruta/al/archivo/cubdl_data.mat")

# 2. Especificar frecuencias de interés en Hz (dentro del ancho de banda)
selected_freqs = [3.5e6, 5.0e6, 6.5e6]

# 3. Extraer la presión compleja en los sensores para la transmisión indexada (ej. tx 0)
sensor_x, P_meas = loader.extract_boundary_pressure(
    selected_frequencies=selected_freqs, 
    transmit_index=0
)

# P_meas es un diccionario: { 3500000.0: [N_channels] complejos, ... }
print("Coordenadas de sensores (x):", sensor_x)
print("Presión compleja a 5.0 MHz:", P_meas[5e6])
```

---

## 3. Fuentes Consultadas y Referencias

Para la especificación del mapeo de variables y las dimensiones del tensor de canal en CUBDL, se consultaron los repositorios y estándares oficiales de beamforming:
1.  **CUBDL Johns Hopkins Repository**: [https://github.com/jhu-lumi/CUBDL](https://github.com/jhu-lumi/CUBDL) ( Johns Hopkins University, LUMI Lab).
2.  **USTB (Ultrasound Toolbox Norway)**: [https://www.ustb.no/](https://www.ustb.no/) (Software para la lectura de datos tipo L9-4 o simulados mediante Field II).
3.  **HDF5 Dataset Mapping in Python**: Mapeo de referencias de memoria MATLAB v7.3 mediante la librería `h5py` para resolver la transposición de arrays `(samples, channels, transmits)`.
