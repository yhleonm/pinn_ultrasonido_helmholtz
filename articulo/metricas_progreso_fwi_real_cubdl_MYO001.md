# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-25 19:12:56

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 4.67%
*   **Error L2 alpha0 final en ROI**: 95.03%
*   **SSIM Final c**: 0.129
*   **SSIM Final alpha0**: -0.244

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2831e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3736e+06 | 2.82% | 40.67% | -0.055 | 0.092 |
| 20 | 1 | 1.0 MHz | 8.8042e+06 | 2.84% | 41.22% | -0.071 | 0.144 |
| 30 | 1 | 1.0 MHz | 8.5912e+06 | 2.97% | 41.02% | -0.023 | 0.176 |
| 31 | 2 | 2.0 MHz | 9.0793e+06 | 2.97% | 41.05% | -0.021 | 0.177 |
| 40 | 2 | 2.0 MHz | 7.0350e+06 | 2.41% | 38.17% | 0.300 | 0.318 |
| 50 | 2 | 2.0 MHz | 6.6314e+06 | 2.70% | 37.79% | 0.257 | 0.216 |
| 60 | 2 | 2.0 MHz | 6.4654e+06 | 2.72% | 39.18% | 0.283 | 0.159 |
| 61 | 3 | 3.5 MHz | 1.3930e+08 | 2.68% | 39.72% | 0.285 | 0.141 |
| 70 | 3 | 3.5 MHz | 1.3909e+08 | 4.16% | 40.36% | 0.111 | 0.066 |
| 80 | 3 | 3.5 MHz | 1.2980e+08 | 4.24% | 49.66% | 0.117 | -0.115 |
| 90 | 3 | 3.5 MHz | 1.2282e+08 | 4.34% | 71.17% | 0.130 | -0.226 |
| 100 | 3 | 3.5 MHz | 1.1917e+08 | 4.39% | 83.44% | 0.150 | -0.248 |
| 110 | 3 | 3.5 MHz | 1.1162e+08 | 4.30% | 88.91% | 0.159 | -0.251 |
| 120 | 3 | 3.5 MHz | 1.1134e+08 | 4.43% | 92.04% | 0.148 | -0.249 |
| 130 | 3 | 3.5 MHz | 1.0784e+08 | 4.50% | 93.75% | 0.140 | -0.246 |
| 140 | 3 | 3.5 MHz | 1.0628e+08 | 4.67% | 95.03% | 0.129 | -0.244 |

