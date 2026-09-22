# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-18 22:11:16

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 0.00 mm, z = 23.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 3.43%
*   **Error L2 alpha0 final en ROI**: 92.31%
*   **SSIM Final c**: -0.215
*   **SSIM Final alpha0**: -0.168

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO002.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO002.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 2.2529e+07 | 2.63% | 39.14% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.0618e+07 | 4.88% | 39.83% | -0.196 | 0.241 |
| 20 | 1 | 1.0 MHz | 9.3600e+06 | 4.57% | 43.50% | -0.092 | -0.009 |
| 30 | 1 | 1.0 MHz | 8.7037e+06 | 4.61% | 49.67% | 0.003 | -0.214 |
| 31 | 2 | 2.0 MHz | 1.5039e+07 | 4.60% | 50.37% | 0.019 | -0.234 |
| 40 | 2 | 2.0 MHz | 9.7138e+06 | 5.47% | 56.00% | -0.196 | -0.314 |
| 50 | 2 | 2.0 MHz | 8.4029e+06 | 5.49% | 76.68% | -0.232 | -0.318 |
| 60 | 2 | 2.0 MHz | 7.7685e+06 | 5.58% | 85.02% | -0.234 | -0.299 |
| 61 | 3 | 3.5 MHz | 1.7158e+08 | 5.61% | 85.71% | -0.234 | -0.297 |
| 70 | 3 | 3.5 MHz | 1.6622e+08 | 3.07% | 69.58% | -0.347 | -0.261 |
| 80 | 3 | 3.5 MHz | 1.4988e+08 | 3.03% | 65.30% | -0.296 | -0.310 |
| 90 | 3 | 3.5 MHz | 1.2989e+08 | 3.12% | 74.63% | -0.269 | -0.293 |
| 100 | 3 | 3.5 MHz | 1.1893e+08 | 3.26% | 84.67% | -0.276 | -0.270 |
| 110 | 3 | 3.5 MHz | 1.1865e+08 | 3.40% | 89.33% | -0.261 | -0.239 |
| 120 | 3 | 3.5 MHz | 1.1265e+08 | 3.41% | 91.84% | -0.258 | -0.212 |
| 130 | 3 | 3.5 MHz | 1.0629e+08 | 3.42% | 92.63% | -0.237 | -0.187 |
| 140 | 3 | 3.5 MHz | 1.0331e+08 | 3.43% | 92.31% | -0.215 | -0.168 |

