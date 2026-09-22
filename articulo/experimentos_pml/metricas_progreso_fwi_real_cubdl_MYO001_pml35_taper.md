# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-06 19:22:13

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 6.08%
*   **Error L2 alpha0 final en ROI**: 45.02%
*   **SSIM Final c**: -0.030
*   **SSIM Final alpha0**: 0.149

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_pml35_taper.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_pml35_taper.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.6550e+07 | 2.57% | 39.23% | -0.268 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.1029e+07 | 3.18% | 38.92% | 0.068 | 0.185 |
| 20 | 1 | 1.0 MHz | 1.0097e+07 | 3.53% | 43.24% | 0.149 | 0.033 |
| 30 | 1 | 1.0 MHz | 9.6608e+06 | 3.62% | 47.78% | 0.260 | -0.072 |
| 31 | 2 | 2.0 MHz | 1.0656e+07 | 3.85% | 48.50% | 0.245 | -0.074 |
| 40 | 2 | 2.0 MHz | 1.0641e+07 | 3.98% | 52.21% | 0.332 | -0.244 |
| 50 | 2 | 2.0 MHz | 1.0313e+07 | 5.66% | 46.66% | 0.136 | -0.114 |
| 60 | 2 | 2.0 MHz | 1.0754e+07 | 6.17% | 43.18% | -0.064 | 0.011 |
| 61 | 3 | 3.5 MHz | 1.4635e+08 | 6.27% | 43.61% | -0.061 | 0.008 |
| 70 | 3 | 3.5 MHz | 1.4556e+08 | 5.63% | 45.98% | -0.020 | -0.103 |
| 80 | 3 | 3.5 MHz | 1.4635e+08 | 4.46% | 43.97% | 0.112 | -0.043 |
| 90 | 3 | 3.5 MHz | 1.4646e+08 | 4.44% | 41.76% | 0.082 | 0.054 |
| 100 | 3 | 3.5 MHz | 1.4661e+08 | 4.50% | 41.71% | 0.137 | 0.064 |
| 110 | 3 | 3.5 MHz | 1.4645e+08 | 5.39% | 46.56% | 0.055 | -0.121 |
| 120 | 3 | 3.5 MHz | 1.4636e+08 | 6.03% | 43.16% | -0.044 | -0.039 |
| 130 | 3 | 3.5 MHz | 1.4737e+08 | 6.60% | 42.49% | -0.025 | 0.135 |
| 140 | 3 | 3.5 MHz | 1.4636e+08 | 6.08% | 45.02% | -0.030 | 0.149 |

