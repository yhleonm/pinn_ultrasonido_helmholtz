# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-07-07 14:38:43

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 0.00 mm, z = 23.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.53%
*   **Error L2 alpha0 final en ROI**: 95.73%
*   **SSIM Final c**: 0.113
*   **SSIM Final alpha0**: -0.247

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO002_pml35_sharp.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO002_pml35_sharp.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.0557e+07 | 2.63% | 39.14% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.9572e+06 | 4.46% | 37.73% | -0.248 | 0.215 |
| 20 | 1 | 1.0 MHz | 9.9102e+06 | 4.65% | 45.90% | -0.257 | 0.202 |
| 30 | 1 | 1.0 MHz | 9.7371e+06 | 4.83% | 55.78% | -0.218 | 0.211 |
| 31 | 2 | 2.0 MHz | 1.1880e+07 | 4.81% | 56.09% | -0.221 | 0.212 |
| 40 | 2 | 2.0 MHz | 1.0043e+07 | 4.59% | 80.88% | -0.286 | 0.276 |
| 50 | 2 | 2.0 MHz | 9.4400e+06 | 4.67% | 85.90% | -0.203 | 0.038 |
| 60 | 2 | 2.0 MHz | 8.8913e+06 | 4.70% | 94.33% | -0.137 | -0.120 |
| 61 | 3 | 3.5 MHz | 1.6429e+08 | 4.86% | 95.64% | -0.149 | -0.115 |
| 70 | 3 | 3.5 MHz | 1.5700e+08 | 4.87% | 76.17% | 0.008 | -0.210 |
| 80 | 3 | 3.5 MHz | 1.5265e+08 | 4.99% | 73.43% | 0.118 | -0.198 |
| 90 | 3 | 3.5 MHz | 1.4395e+08 | 5.30% | 76.64% | 0.093 | -0.166 |
| 100 | 3 | 3.5 MHz | 1.4062e+08 | 5.16% | 84.02% | 0.115 | -0.197 |
| 110 | 3 | 3.5 MHz | 1.3183e+08 | 5.34% | 88.50% | 0.117 | -0.224 |
| 120 | 3 | 3.5 MHz | 1.2702e+08 | 5.42% | 91.75% | 0.118 | -0.256 |
| 130 | 3 | 3.5 MHz | 1.2381e+08 | 5.48% | 94.18% | 0.118 | -0.250 |
| 140 | 3 | 3.5 MHz | 1.2226e+08 | 5.53% | 95.73% | 0.113 | -0.247 |

