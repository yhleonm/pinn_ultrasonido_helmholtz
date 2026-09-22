# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-07 15:12:11

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -2.50 mm, z = 11.00 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 3.94%
*   **Error L2 alpha0 final en ROI**: 89.34%
*   **SSIM Final c**: 0.036
*   **SSIM Final alpha0**: -0.243

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO004_pml35_sharp.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO004_pml35_sharp.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.1803e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.2322e+06 | 2.13% | 44.15% | 0.262 | -0.095 |
| 20 | 1 | 1.0 MHz | 8.4271e+06 | 2.47% | 39.78% | 0.076 | 0.148 |
| 30 | 1 | 1.0 MHz | 8.1925e+06 | 2.43% | 37.61% | -0.008 | 0.337 |
| 31 | 2 | 2.0 MHz | 8.4452e+06 | 2.41% | 37.53% | 0.004 | 0.349 |
| 40 | 2 | 2.0 MHz | 6.7943e+06 | 2.84% | 43.02% | -0.209 | 0.364 |
| 50 | 2 | 2.0 MHz | 6.3250e+06 | 2.76% | 37.22% | -0.076 | 0.386 |
| 60 | 2 | 2.0 MHz | 6.0577e+06 | 2.72% | 39.23% | 0.013 | 0.218 |
| 61 | 3 | 3.5 MHz | 2.3720e+08 | 2.78% | 40.11% | -0.001 | 0.183 |
| 70 | 3 | 3.5 MHz | 2.2609e+08 | 3.75% | 54.12% | -0.130 | 0.014 |
| 80 | 3 | 3.5 MHz | 2.0496e+08 | 3.89% | 71.36% | -0.091 | 0.032 |
| 90 | 3 | 3.5 MHz | 1.8489e+08 | 3.81% | 79.32% | -0.033 | 0.011 |
| 100 | 3 | 3.5 MHz | 1.7490e+08 | 3.77% | 77.33% | -0.001 | -0.086 |
| 110 | 3 | 3.5 MHz | 1.5458e+08 | 3.90% | 80.56% | 0.017 | -0.169 |
| 120 | 3 | 3.5 MHz | 1.4626e+08 | 3.87% | 83.92% | 0.036 | -0.209 |
| 130 | 3 | 3.5 MHz | 1.3850e+08 | 3.90% | 86.44% | 0.040 | -0.246 |
| 140 | 3 | 3.5 MHz | 1.3235e+08 | 3.94% | 89.34% | 0.036 | -0.243 |

