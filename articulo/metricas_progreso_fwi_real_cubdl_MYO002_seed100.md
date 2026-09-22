# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-10 19:01:10

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 0.00 mm, z = 23.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.17%
*   **Error L2 alpha0 final en ROI**: 94.46%
*   **SSIM Final c**: 0.103
*   **SSIM Final alpha0**: -0.245

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO002_seed100.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO002_seed100.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2864e+07 | 2.63% | 39.14% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3761e+06 | 3.44% | 37.10% | -0.104 | 0.163 |
| 20 | 1 | 1.0 MHz | 8.8618e+06 | 3.19% | 42.66% | 0.067 | -0.030 |
| 30 | 1 | 1.0 MHz | 8.4008e+06 | 3.47% | 46.23% | 0.073 | -0.103 |
| 31 | 2 | 2.0 MHz | 9.9082e+06 | 3.54% | 46.29% | 0.063 | -0.104 |
| 40 | 2 | 2.0 MHz | 7.3741e+06 | 2.81% | 43.84% | 0.048 | -0.076 |
| 50 | 2 | 2.0 MHz | 6.7752e+06 | 3.40% | 39.61% | -0.080 | 0.120 |
| 60 | 2 | 2.0 MHz | 6.3447e+06 | 3.78% | 41.40% | -0.029 | 0.090 |
| 61 | 3 | 3.5 MHz | 1.5815e+08 | 3.85% | 41.83% | -0.022 | 0.076 |
| 70 | 3 | 3.5 MHz | 1.5648e+08 | 4.35% | 48.35% | -0.009 | -0.064 |
| 80 | 3 | 3.5 MHz | 1.4824e+08 | 4.77% | 53.48% | 0.030 | -0.153 |
| 90 | 3 | 3.5 MHz | 1.4497e+08 | 4.97% | 65.04% | 0.068 | -0.248 |
| 100 | 3 | 3.5 MHz | 1.3880e+08 | 5.11% | 76.28% | 0.067 | -0.270 |
| 110 | 3 | 3.5 MHz | 1.3363e+08 | 5.17% | 84.36% | 0.083 | -0.264 |
| 120 | 3 | 3.5 MHz | 1.2947e+08 | 5.24% | 89.51% | 0.088 | -0.256 |
| 130 | 3 | 3.5 MHz | 1.3044e+08 | 5.19% | 92.61% | 0.098 | -0.249 |
| 140 | 3 | 3.5 MHz | 1.2450e+08 | 5.17% | 94.46% | 0.103 | -0.245 |

