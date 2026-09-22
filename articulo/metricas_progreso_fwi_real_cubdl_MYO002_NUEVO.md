# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-08 20:16:35

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 0.00 mm, z = 23.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.32%
*   **Error L2 alpha0 final en ROI**: 95.56%
*   **SSIM Final c**: 0.011
*   **SSIM Final alpha0**: -0.238

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO002_NUEVO.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO002_NUEVO.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2864e+07 | 2.63% | 39.14% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.4505e+06 | 3.76% | 36.41% | -0.282 | 0.192 |
| 20 | 1 | 1.0 MHz | 9.1864e+06 | 3.76% | 36.88% | -0.240 | 0.241 |
| 30 | 1 | 1.0 MHz | 8.5811e+06 | 3.45% | 37.49% | -0.234 | 0.244 |
| 31 | 2 | 2.0 MHz | 1.0218e+07 | 3.50% | 37.81% | -0.219 | 0.238 |
| 40 | 2 | 2.0 MHz | 8.0888e+06 | 4.04% | 38.13% | -0.272 | 0.254 |
| 50 | 2 | 2.0 MHz | 6.8527e+06 | 4.10% | 47.41% | -0.228 | -0.090 |
| 60 | 2 | 2.0 MHz | 6.4001e+06 | 4.41% | 58.07% | -0.196 | -0.209 |
| 61 | 3 | 3.5 MHz | 1.6449e+08 | 4.47% | 58.88% | -0.192 | -0.215 |
| 70 | 3 | 3.5 MHz | 1.5795e+08 | 4.17% | 52.05% | -0.091 | 0.020 |
| 80 | 3 | 3.5 MHz | 1.4752e+08 | 4.92% | 61.97% | -0.059 | -0.108 |
| 90 | 3 | 3.5 MHz | 1.4210e+08 | 5.22% | 66.70% | -0.045 | -0.222 |
| 100 | 3 | 3.5 MHz | 1.3545e+08 | 5.35% | 77.56% | -0.026 | -0.243 |
| 110 | 3 | 3.5 MHz | 1.2780e+08 | 5.30% | 86.39% | -0.003 | -0.244 |
| 120 | 3 | 3.5 MHz | 1.2608e+08 | 5.44% | 91.26% | 0.001 | -0.241 |
| 130 | 3 | 3.5 MHz | 1.2106e+08 | 5.37% | 93.96% | 0.004 | -0.240 |
| 140 | 3 | 3.5 MHz | 1.1856e+08 | 5.32% | 95.56% | 0.011 | -0.238 |

