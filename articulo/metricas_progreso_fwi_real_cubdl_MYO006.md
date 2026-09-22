# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-25 22:53:10

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -16.50 mm, z = 10.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.80%
*   **Error L2 alpha0 final en ROI**: 92.00%
*   **SSIM Final c**: 0.025
*   **SSIM Final alpha0**: -0.325

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO006.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO006.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.5972e+07 | 3.02% | 43.64% | -0.288 | 0.040 |
| 10 | 1 | 1.0 MHz | 1.1198e+07 | 3.79% | 41.16% | -0.083 | 0.190 |
| 20 | 1 | 1.0 MHz | 1.0555e+07 | 4.06% | 42.30% | -0.077 | 0.138 |
| 30 | 1 | 1.0 MHz | 1.0335e+07 | 4.29% | 41.36% | -0.086 | 0.195 |
| 31 | 2 | 2.0 MHz | 1.1891e+07 | 4.31% | 41.38% | -0.092 | 0.195 |
| 40 | 2 | 2.0 MHz | 7.7999e+06 | 3.62% | 44.35% | -0.105 | 0.055 |
| 50 | 2 | 2.0 MHz | 7.6766e+06 | 3.84% | 45.30% | -0.084 | 0.010 |
| 60 | 2 | 2.0 MHz | 7.4435e+06 | 3.96% | 46.49% | -0.094 | -0.016 |
| 61 | 3 | 3.5 MHz | 2.4560e+08 | 3.98% | 46.61% | -0.085 | -0.019 |
| 70 | 3 | 3.5 MHz | 2.2644e+08 | 3.96% | 46.56% | 0.167 | -0.033 |
| 80 | 3 | 3.5 MHz | 2.0945e+08 | 4.86% | 51.14% | 0.097 | -0.064 |
| 90 | 3 | 3.5 MHz | 1.9947e+08 | 5.24% | 58.65% | 0.065 | -0.142 |
| 100 | 3 | 3.5 MHz | 1.9227e+08 | 5.57% | 68.39% | 0.045 | -0.242 |
| 110 | 3 | 3.5 MHz | 1.8761e+08 | 5.64% | 77.79% | 0.036 | -0.295 |
| 120 | 3 | 3.5 MHz | 1.8622e+08 | 5.70% | 84.96% | 0.029 | -0.319 |
| 130 | 3 | 3.5 MHz | 1.8238e+08 | 5.73% | 89.35% | 0.026 | -0.325 |
| 140 | 3 | 3.5 MHz | 1.8020e+08 | 5.80% | 92.00% | 0.025 | -0.325 |

