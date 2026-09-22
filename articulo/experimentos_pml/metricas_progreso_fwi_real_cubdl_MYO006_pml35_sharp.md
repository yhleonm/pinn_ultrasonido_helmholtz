# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-07 15:45:15

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -16.50 mm, z = 10.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.90%
*   **Error L2 alpha0 final en ROI**: 79.86%
*   **SSIM Final c**: 0.004
*   **SSIM Final alpha0**: -0.226

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO006_pml35_sharp.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO006_pml35_sharp.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.6389e+07 | 3.02% | 43.64% | -0.288 | 0.040 |
| 10 | 1 | 1.0 MHz | 1.1435e+07 | 3.84% | 46.87% | -0.204 | 0.001 |
| 20 | 1 | 1.0 MHz | 1.0612e+07 | 4.34% | 49.33% | -0.191 | -0.063 |
| 30 | 1 | 1.0 MHz | 1.0372e+07 | 4.61% | 50.47% | -0.153 | -0.054 |
| 31 | 2 | 2.0 MHz | 1.1163e+07 | 4.66% | 50.63% | -0.153 | -0.054 |
| 40 | 2 | 2.0 MHz | 7.8155e+06 | 3.67% | 46.82% | 0.137 | 0.192 |
| 50 | 2 | 2.0 MHz | 7.5430e+06 | 3.83% | 45.79% | 0.185 | 0.226 |
| 60 | 2 | 2.0 MHz | 7.3617e+06 | 4.10% | 46.65% | 0.135 | 0.170 |
| 61 | 3 | 3.5 MHz | 2.4424e+08 | 4.13% | 47.01% | 0.117 | 0.152 |
| 70 | 3 | 3.5 MHz | 2.2866e+08 | 4.63% | 47.42% | -0.075 | 0.065 |
| 80 | 3 | 3.5 MHz | 2.0859e+08 | 5.32% | 54.12% | -0.102 | 0.005 |
| 90 | 3 | 3.5 MHz | 2.0133e+08 | 5.63% | 61.03% | -0.091 | -0.065 |
| 100 | 3 | 3.5 MHz | 1.9994e+08 | 5.62% | 64.04% | -0.054 | -0.114 |
| 110 | 3 | 3.5 MHz | 1.9773e+08 | 5.66% | 65.58% | -0.000 | -0.091 |
| 120 | 3 | 3.5 MHz | 1.9340e+08 | 5.83% | 69.78% | -0.002 | -0.121 |
| 130 | 3 | 3.5 MHz | 1.8904e+08 | 5.87% | 74.61% | 0.006 | -0.167 |
| 140 | 3 | 3.5 MHz | 1.8572e+08 | 5.90% | 79.86% | 0.004 | -0.226 |

