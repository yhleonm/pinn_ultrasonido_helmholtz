# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-06 18:55:28

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.79%
*   **Error L2 alpha0 final en ROI**: 91.41%
*   **SSIM Final c**: 0.197
*   **SSIM Final alpha0**: -0.254

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_pml15_taper.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_pml15_taper.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.5760e+07 | 2.57% | 39.23% | -0.268 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.1195e+07 | 3.12% | 36.38% | 0.055 | 0.282 |
| 20 | 1 | 1.0 MHz | 1.0404e+07 | 3.33% | 36.36% | 0.260 | 0.306 |
| 30 | 1 | 1.0 MHz | 9.7105e+06 | 3.10% | 38.96% | 0.354 | 0.226 |
| 31 | 2 | 2.0 MHz | 1.0794e+07 | 3.15% | 39.52% | 0.362 | 0.210 |
| 40 | 2 | 2.0 MHz | 1.0422e+07 | 4.67% | 48.12% | -0.095 | -0.069 |
| 50 | 2 | 2.0 MHz | 9.9710e+06 | 5.20% | 59.21% | -0.028 | -0.189 |
| 60 | 2 | 2.0 MHz | 1.0794e+07 | 5.87% | 67.51% | -0.042 | -0.230 |
| 61 | 3 | 3.5 MHz | 1.4635e+08 | 6.05% | 68.18% | -0.033 | -0.215 |
| 70 | 3 | 3.5 MHz | 1.4624e+08 | 6.05% | 55.72% | 0.024 | -0.034 |
| 80 | 3 | 3.5 MHz | 1.4602e+08 | 6.33% | 66.44% | 0.041 | -0.211 |
| 90 | 3 | 3.5 MHz | 1.4603e+08 | 6.14% | 75.31% | 0.052 | -0.231 |
| 100 | 3 | 3.5 MHz | 1.4597e+08 | 6.26% | 82.93% | 0.045 | -0.241 |
| 110 | 3 | 3.5 MHz | 1.4588e+08 | 5.99% | 86.32% | 0.121 | -0.255 |
| 120 | 3 | 3.5 MHz | 1.4558e+08 | 6.09% | 89.49% | 0.115 | -0.263 |
| 130 | 3 | 3.5 MHz | 1.4573e+08 | 6.18% | 91.94% | 0.125 | -0.259 |
| 140 | 3 | 3.5 MHz | 1.4633e+08 | 5.79% | 91.41% | 0.197 | -0.254 |

