# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-18 19:19:30

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -6.50 mm, z = 20.00 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.37%
*   **Error L2 alpha0 final en ROI**: 94.05%
*   **SSIM Final c**: -0.215
*   **SSIM Final alpha0**: -0.214

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_TSH002.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_TSH002.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 2.6663e+05 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.6145e+05 | 4.60% | 38.66% | -0.187 | 0.086 |
| 20 | 1 | 1.0 MHz | 1.3640e+05 | 5.62% | 53.44% | -0.172 | -0.235 |
| 30 | 1 | 1.0 MHz | 1.2463e+05 | 5.92% | 65.37% | -0.144 | -0.233 |
| 31 | 2 | 2.0 MHz | 2.9232e+05 | 5.88% | 66.17% | -0.144 | -0.238 |
| 40 | 2 | 2.0 MHz | 2.4194e+05 | 4.69% | 64.09% | -0.132 | -0.163 |
| 50 | 2 | 2.0 MHz | 2.0812e+05 | 4.41% | 70.07% | -0.205 | -0.237 |
| 60 | 2 | 2.0 MHz | 1.9927e+05 | 4.46% | 76.52% | -0.218 | -0.280 |
| 61 | 3 | 3.5 MHz | 2.5322e+08 | 4.48% | 77.25% | -0.220 | -0.283 |
| 70 | 3 | 3.5 MHz | 2.3354e+08 | 4.82% | 89.63% | -0.208 | -0.116 |
| 80 | 3 | 3.5 MHz | 2.2774e+08 | 4.75% | 91.64% | -0.191 | -0.097 |
| 90 | 3 | 3.5 MHz | 2.1364e+08 | 4.68% | 91.46% | -0.199 | -0.123 |
| 100 | 3 | 3.5 MHz | 1.9753e+08 | 5.00% | 89.45% | -0.216 | -0.191 |
| 110 | 3 | 3.5 MHz | 1.9386e+08 | 4.98% | 89.03% | -0.223 | -0.222 |
| 120 | 3 | 3.5 MHz | 1.8848e+08 | 5.13% | 92.01% | -0.225 | -0.213 |
| 130 | 3 | 3.5 MHz | 1.8460e+08 | 5.28% | 94.67% | -0.216 | -0.192 |
| 140 | 3 | 3.5 MHz | 1.7405e+08 | 5.37% | 94.05% | -0.215 | -0.214 |

