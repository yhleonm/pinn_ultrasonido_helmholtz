# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-06 18:25:36

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 4.91%
*   **Error L2 alpha0 final en ROI**: 94.89%
*   **SSIM Final c**: 0.193
*   **SSIM Final alpha0**: -0.245

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_pml35_sharp.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_pml35_sharp.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2894e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3720e+06 | 2.91% | 41.01% | -0.083 | 0.125 |
| 20 | 1 | 1.0 MHz | 8.8078e+06 | 2.88% | 42.28% | -0.077 | 0.188 |
| 30 | 1 | 1.0 MHz | 8.5995e+06 | 2.94% | 41.30% | -0.000 | 0.239 |
| 31 | 2 | 2.0 MHz | 9.4456e+06 | 3.00% | 41.56% | -0.025 | 0.227 |
| 40 | 2 | 2.0 MHz | 6.9123e+06 | 2.46% | 38.57% | 0.305 | 0.390 |
| 50 | 2 | 2.0 MHz | 6.7196e+06 | 2.82% | 37.02% | 0.254 | 0.278 |
| 60 | 2 | 2.0 MHz | 6.5347e+06 | 2.85% | 40.22% | 0.275 | 0.090 |
| 61 | 3 | 3.5 MHz | 1.3805e+08 | 2.84% | 40.84% | 0.277 | 0.071 |
| 70 | 3 | 3.5 MHz | 1.3699e+08 | 3.97% | 42.53% | 0.273 | 0.050 |
| 80 | 3 | 3.5 MHz | 1.3265e+08 | 4.26% | 59.86% | 0.233 | -0.176 |
| 90 | 3 | 3.5 MHz | 1.3013e+08 | 4.70% | 76.27% | 0.192 | -0.212 |
| 100 | 3 | 3.5 MHz | 1.2275e+08 | 4.67% | 83.63% | 0.188 | -0.227 |
| 110 | 3 | 3.5 MHz | 1.2206e+08 | 4.79% | 88.09% | 0.188 | -0.239 |
| 120 | 3 | 3.5 MHz | 1.1463e+08 | 4.69% | 90.25% | 0.200 | -0.248 |
| 130 | 3 | 3.5 MHz | 1.0861e+08 | 4.82% | 92.83% | 0.199 | -0.249 |
| 140 | 3 | 3.5 MHz | 1.0686e+08 | 4.91% | 94.89% | 0.193 | -0.245 |

