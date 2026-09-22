# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-25 21:54:14

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -2.50 mm, z = 11.00 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 3.60%
*   **Error L2 alpha0 final en ROI**: 91.00%
*   **SSIM Final c**: 0.344
*   **SSIM Final alpha0**: -0.278

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO004.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO004.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.1611e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 8.9128e+06 | 2.42% | 43.83% | 0.105 | -0.094 |
| 20 | 1 | 1.0 MHz | 8.3457e+06 | 2.37% | 41.11% | 0.039 | 0.095 |
| 30 | 1 | 1.0 MHz | 8.0365e+06 | 2.48% | 43.18% | -0.062 | 0.165 |
| 31 | 2 | 2.0 MHz | 8.2838e+06 | 2.51% | 43.70% | -0.079 | 0.158 |
| 40 | 2 | 2.0 MHz | 6.7104e+06 | 2.68% | 51.90% | -0.035 | 0.173 |
| 50 | 2 | 2.0 MHz | 6.1940e+06 | 2.83% | 46.98% | -0.019 | 0.108 |
| 60 | 2 | 2.0 MHz | 6.2598e+06 | 2.91% | 49.74% | -0.014 | -0.047 |
| 61 | 3 | 3.5 MHz | 2.3096e+08 | 2.87% | 50.28% | 0.011 | -0.052 |
| 70 | 3 | 3.5 MHz | 2.1733e+08 | 2.70% | 57.35% | 0.371 | 0.282 |
| 80 | 3 | 3.5 MHz | 2.0955e+08 | 2.85% | 60.33% | 0.368 | 0.235 |
| 90 | 3 | 3.5 MHz | 2.0184e+08 | 3.18% | 59.54% | 0.361 | -0.123 |
| 100 | 3 | 3.5 MHz | 1.7913e+08 | 3.31% | 73.36% | 0.360 | -0.352 |
| 110 | 3 | 3.5 MHz | 1.6974e+08 | 3.35% | 82.20% | 0.358 | -0.329 |
| 120 | 3 | 3.5 MHz | 1.6544e+08 | 3.42% | 86.35% | 0.362 | -0.308 |
| 130 | 3 | 3.5 MHz | 1.6174e+08 | 3.54% | 89.12% | 0.352 | -0.291 |
| 140 | 3 | 3.5 MHz | 1.5623e+08 | 3.60% | 91.00% | 0.344 | -0.278 |

