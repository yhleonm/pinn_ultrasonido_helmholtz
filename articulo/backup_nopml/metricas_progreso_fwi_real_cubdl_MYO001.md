# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-18 18:17:11

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 4.30%
*   **Error L2 alpha0 final en ROI**: 87.71%
*   **SSIM Final c**: -0.032
*   **SSIM Final alpha0**: -0.247

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 2.2791e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.0682e+07 | 5.15% | 38.52% | -0.234 | 0.078 |
| 20 | 1 | 1.0 MHz | 9.6314e+06 | 5.16% | 45.37% | -0.162 | 0.034 |
| 30 | 1 | 1.0 MHz | 8.6218e+06 | 5.26% | 56.11% | -0.175 | -0.124 |
| 31 | 2 | 2.0 MHz | 1.2592e+07 | 5.30% | 57.14% | -0.175 | -0.138 |
| 40 | 2 | 2.0 MHz | 1.0007e+07 | 6.43% | 75.31% | -0.175 | -0.252 |
| 50 | 2 | 2.0 MHz | 9.0608e+06 | 6.63% | 86.25% | -0.207 | -0.271 |
| 60 | 2 | 2.0 MHz | 7.8169e+06 | 6.63% | 90.16% | -0.210 | -0.243 |
| 61 | 3 | 3.5 MHz | 1.4555e+08 | 6.65% | 90.27% | -0.211 | -0.246 |
| 70 | 3 | 3.5 MHz | 1.4592e+08 | 4.39% | 96.34% | -0.031 | 0.044 |
| 80 | 3 | 3.5 MHz | 1.4293e+08 | 4.04% | 83.96% | -0.032 | 0.040 |
| 90 | 3 | 3.5 MHz | 1.3084e+08 | 4.20% | 71.28% | -0.058 | -0.091 |
| 100 | 3 | 3.5 MHz | 1.1974e+08 | 4.29% | 72.07% | -0.040 | -0.239 |
| 110 | 3 | 3.5 MHz | 1.1725e+08 | 4.16% | 79.13% | -0.035 | -0.255 |
| 120 | 3 | 3.5 MHz | 1.0336e+08 | 4.20% | 83.66% | -0.036 | -0.255 |
| 130 | 3 | 3.5 MHz | 9.7872e+07 | 4.27% | 86.04% | -0.033 | -0.244 |
| 140 | 3 | 3.5 MHz | 9.7598e+07 | 4.30% | 87.71% | -0.032 | -0.247 |

