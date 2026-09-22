# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-15 21:42:33

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 6.53%
*   **Error L2 alpha0 final en ROI**: 61.13%
*   **SSIM Final c**: -0.166
*   **SSIM Final alpha0**: -0.277

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_decay_seed100.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_decay_seed100.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2837e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3052e+06 | 2.78% | 42.45% | 0.012 | 0.008 |
| 20 | 1 | 1.0 MHz | 8.9366e+06 | 2.69% | 42.05% | -0.024 | 0.036 |
| 30 | 1 | 1.0 MHz | 8.7561e+06 | 2.73% | 41.77% | -0.050 | 0.054 |
| 31 | 2 | 2.0 MHz | 8.5670e+06 | 2.74% | 41.78% | -0.048 | 0.056 |
| 40 | 2 | 2.0 MHz | 7.2630e+06 | 3.05% | 45.72% | -0.011 | 0.049 |
| 50 | 2 | 2.0 MHz | 6.9592e+06 | 3.57% | 43.62% | -0.046 | -0.029 |
| 60 | 2 | 2.0 MHz | 6.8119e+06 | 3.74% | 43.30% | -0.027 | -0.069 |
| 61 | 3 | 3.5 MHz | 1.4094e+08 | 3.74% | 43.37% | -0.022 | -0.072 |
| 70 | 3 | 3.5 MHz | 1.3648e+08 | 6.40% | 39.82% | -0.200 | 0.049 |
| 80 | 3 | 3.5 MHz | 1.3209e+08 | 6.71% | 42.65% | -0.176 | -0.051 |
| 90 | 3 | 3.5 MHz | 1.2848e+08 | 6.61% | 50.25% | -0.170 | -0.190 |
| 100 | 3 | 3.5 MHz | 1.2549e+08 | 6.62% | 55.78% | -0.167 | -0.245 |
| 110 | 3 | 3.5 MHz | 1.2266e+08 | 6.61% | 58.66% | -0.167 | -0.264 |
| 120 | 3 | 3.5 MHz | 1.2041e+08 | 6.58% | 60.04% | -0.166 | -0.270 |
| 130 | 3 | 3.5 MHz | 1.1905e+08 | 6.56% | 60.75% | -0.166 | -0.274 |
| 140 | 3 | 3.5 MHz | 1.1829e+08 | 6.53% | 61.13% | -0.166 | -0.277 |

