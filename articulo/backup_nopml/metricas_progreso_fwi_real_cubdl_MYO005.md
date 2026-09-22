# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-23 00:04:50

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 4.00 mm, z = 11.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 4.42%
*   **Error L2 alpha0 final en ROI**: 97.31%
*   **SSIM Final c**: -0.117
*   **SSIM Final alpha0**: -0.228

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO005.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO005.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.6821e+07 | 2.63% | 39.14% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.0239e+07 | 3.73% | 49.42% | -0.213 | -0.167 |
| 20 | 1 | 1.0 MHz | 9.8832e+06 | 3.48% | 52.78% | -0.187 | -0.178 |
| 30 | 1 | 1.0 MHz | 9.7874e+06 | 3.46% | 55.56% | -0.166 | -0.175 |
| 31 | 2 | 2.0 MHz | 9.2622e+06 | 3.49% | 55.83% | -0.169 | -0.177 |
| 40 | 2 | 2.0 MHz | 7.5720e+06 | 3.43% | 48.49% | -0.249 | -0.163 |
| 50 | 2 | 2.0 MHz | 7.1674e+06 | 3.60% | 54.45% | -0.260 | -0.262 |
| 60 | 2 | 2.0 MHz | 6.9167e+06 | 3.67% | 60.76% | -0.249 | -0.317 |
| 61 | 3 | 3.5 MHz | 2.9426e+08 | 3.72% | 61.30% | -0.248 | -0.318 |
| 70 | 3 | 3.5 MHz | 2.6815e+08 | 3.94% | 67.55% | -0.141 | -0.268 |
| 80 | 3 | 3.5 MHz | 2.4709e+08 | 3.90% | 79.34% | -0.096 | -0.234 |
| 90 | 3 | 3.5 MHz | 2.3393e+08 | 4.11% | 88.02% | -0.099 | -0.228 |
| 100 | 3 | 3.5 MHz | 2.2630e+08 | 4.14% | 93.05% | -0.086 | -0.223 |
| 110 | 3 | 3.5 MHz | 2.1893e+08 | 4.27% | 95.48% | -0.107 | -0.220 |
| 120 | 3 | 3.5 MHz | 2.1179e+08 | 4.29% | 96.30% | -0.106 | -0.230 |
| 130 | 3 | 3.5 MHz | 2.0530e+08 | 4.36% | 96.94% | -0.115 | -0.230 |
| 140 | 3 | 3.5 MHz | 2.0033e+08 | 4.42% | 97.31% | -0.117 | -0.228 |

