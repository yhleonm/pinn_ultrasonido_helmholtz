# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-10 21:39:18

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.40%
*   **Error L2 alpha0 final en ROI**: 64.69%
*   **SSIM Final c**: 0.024
*   **SSIM Final alpha0**: -0.319

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_corr.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_corr.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2837e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3282e+06 | 2.77% | 40.15% | -0.030 | 0.096 |
| 20 | 1 | 1.0 MHz | 8.8916e+06 | 2.71% | 39.80% | -0.087 | 0.122 |
| 30 | 1 | 1.0 MHz | 8.7179e+06 | 2.76% | 39.61% | -0.071 | 0.131 |
| 31 | 2 | 2.0 MHz | 9.3781e+06 | 2.76% | 39.65% | -0.067 | 0.131 |
| 40 | 2 | 2.0 MHz | 7.2578e+06 | 2.72% | 41.00% | 0.058 | 0.121 |
| 50 | 2 | 2.0 MHz | 6.8540e+06 | 3.19% | 41.19% | -0.015 | 0.034 |
| 60 | 2 | 2.0 MHz | 6.6422e+06 | 3.31% | 41.40% | 0.018 | 0.002 |
| 61 | 3 | 3.5 MHz | 1.3588e+08 | 3.33% | 41.46% | 0.019 | 0.000 |
| 70 | 3 | 3.5 MHz | 1.3124e+08 | 4.69% | 43.09% | -0.025 | -0.076 |
| 80 | 3 | 3.5 MHz | 1.2517e+08 | 5.39% | 49.58% | 0.011 | -0.225 |
| 90 | 3 | 3.5 MHz | 1.2065e+08 | 5.52% | 54.96% | 0.016 | -0.291 |
| 100 | 3 | 3.5 MHz | 1.1822e+08 | 5.48% | 58.85% | 0.021 | -0.312 |
| 110 | 3 | 3.5 MHz | 1.1671e+08 | 5.45% | 61.45% | 0.023 | -0.318 |
| 120 | 3 | 3.5 MHz | 1.1580e+08 | 5.42% | 63.09% | 0.024 | -0.320 |
| 130 | 3 | 3.5 MHz | 1.1523e+08 | 5.40% | 64.09% | 0.024 | -0.319 |
| 140 | 3 | 3.5 MHz | 1.1488e+08 | 5.40% | 64.69% | 0.024 | -0.319 |

