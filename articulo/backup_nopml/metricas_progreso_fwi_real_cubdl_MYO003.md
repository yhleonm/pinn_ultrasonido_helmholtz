# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-22 19:45:44

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 16.00 mm, z = 18.00 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.59%
*   **Error L2 alpha0 final en ROI**: 92.07%
*   **SSIM Final c**: -0.163
*   **SSIM Final alpha0**: -0.246

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO003.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO003.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.7569e+07 | 3.00% | 43.42% | -0.287 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.4485e+07 | 5.59% | 41.88% | -0.263 | 0.342 |
| 20 | 1 | 1.0 MHz | 1.4066e+07 | 5.43% | 42.23% | -0.258 | 0.110 |
| 30 | 1 | 1.0 MHz | 1.3838e+07 | 5.84% | 52.20% | -0.195 | -0.178 |
| 31 | 2 | 2.0 MHz | 1.0815e+07 | 5.87% | 53.27% | -0.189 | -0.191 |
| 40 | 2 | 2.0 MHz | 1.0369e+07 | 4.14% | 47.13% | 0.149 | 0.346 |
| 50 | 2 | 2.0 MHz | 9.6765e+06 | 4.30% | 58.18% | 0.072 | 0.209 |
| 60 | 2 | 2.0 MHz | 9.0368e+06 | 4.11% | 66.67% | 0.095 | 0.105 |
| 61 | 3 | 3.5 MHz | 1.6890e+08 | 4.13% | 67.71% | 0.085 | 0.064 |
| 70 | 3 | 3.5 MHz | 1.6915e+08 | 4.45% | 66.21% | -0.051 | -0.093 |
| 80 | 3 | 3.5 MHz | 1.6571e+08 | 5.12% | 69.95% | -0.199 | -0.207 |
| 90 | 3 | 3.5 MHz | 1.6137e+08 | 5.25% | 72.04% | -0.183 | -0.207 |
| 100 | 3 | 3.5 MHz | 1.5463e+08 | 5.43% | 78.72% | -0.193 | -0.269 |
| 110 | 3 | 3.5 MHz | 1.5037e+08 | 5.48% | 83.02% | -0.180 | -0.289 |
| 120 | 3 | 3.5 MHz | 1.5082e+08 | 5.47% | 86.07% | -0.148 | -0.290 |
| 130 | 3 | 3.5 MHz | 1.4890e+08 | 5.64% | 90.25% | -0.165 | -0.276 |
| 140 | 3 | 3.5 MHz | 1.4331e+08 | 5.59% | 92.07% | -0.163 | -0.246 |

