# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-30 13:12:25

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 4.77%
*   **Error L2 alpha0 final en ROI**: 94.45%
*   **SSIM Final c**: 0.147
*   **SSIM Final alpha0**: -0.247

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_quick_test.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_quick_test.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2837e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3778e+06 | 2.82% | 40.67% | -0.056 | 0.092 |
| 20 | 1 | 1.0 MHz | 8.8126e+06 | 2.85% | 41.22% | -0.072 | 0.142 |
| 30 | 1 | 1.0 MHz | 8.5978e+06 | 2.98% | 40.95% | -0.024 | 0.172 |
| 31 | 2 | 2.0 MHz | 9.0676e+06 | 2.97% | 40.98% | -0.020 | 0.173 |
| 40 | 2 | 2.0 MHz | 6.9895e+06 | 2.40% | 38.22% | 0.284 | 0.299 |
| 50 | 2 | 2.0 MHz | 6.6371e+06 | 2.64% | 37.36% | 0.284 | 0.217 |
| 60 | 2 | 2.0 MHz | 6.4751e+06 | 2.71% | 39.56% | 0.301 | 0.137 |
| 61 | 3 | 3.5 MHz | 1.3822e+08 | 2.69% | 40.01% | 0.300 | 0.125 |
| 70 | 3 | 3.5 MHz | 1.3494e+08 | 3.93% | 42.31% | 0.160 | -0.052 |
| 80 | 3 | 3.5 MHz | 1.2438e+08 | 4.33% | 62.39% | 0.165 | -0.258 |
| 90 | 3 | 3.5 MHz | 1.1962e+08 | 4.52% | 77.78% | 0.150 | -0.235 |
| 100 | 3 | 3.5 MHz | 1.1274e+08 | 4.58% | 84.20% | 0.150 | -0.243 |
| 110 | 3 | 3.5 MHz | 1.0911e+08 | 4.48% | 87.95% | 0.161 | -0.257 |
| 120 | 3 | 3.5 MHz | 1.0596e+08 | 4.54% | 91.09% | 0.159 | -0.254 |
| 130 | 3 | 3.5 MHz | 1.0439e+08 | 4.67% | 93.13% | 0.156 | -0.251 |
| 140 | 3 | 3.5 MHz | 1.0254e+08 | 4.77% | 94.45% | 0.147 | -0.247 |

