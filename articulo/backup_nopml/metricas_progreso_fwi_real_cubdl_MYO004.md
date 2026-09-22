# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-22 20:42:52

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -2.50 mm, z = 11.00 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 6.14%
*   **Error L2 alpha0 final en ROI**: 94.63%
*   **SSIM Final c**: -0.211
*   **SSIM Final alpha0**: -0.260

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO004.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO004.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 2.4390e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.0458e+07 | 4.57% | 42.56% | -0.166 | 0.168 |
| 20 | 1 | 1.0 MHz | 7.8651e+06 | 5.45% | 51.42% | -0.136 | 0.136 |
| 30 | 1 | 1.0 MHz | 6.7058e+06 | 5.74% | 50.36% | -0.075 | -0.003 |
| 31 | 2 | 2.0 MHz | 1.2193e+07 | 5.73% | 49.00% | -0.077 | -0.007 |
| 40 | 2 | 2.0 MHz | 8.6611e+06 | 5.18% | 51.80% | -0.166 | -0.302 |
| 50 | 2 | 2.0 MHz | 7.2306e+06 | 5.39% | 67.77% | -0.177 | -0.335 |
| 60 | 2 | 2.0 MHz | 6.2247e+06 | 5.38% | 75.69% | -0.199 | -0.330 |
| 61 | 3 | 3.5 MHz | 2.6446e+08 | 5.35% | 76.07% | -0.204 | -0.332 |
| 70 | 3 | 3.5 MHz | 2.4938e+08 | 6.52% | 87.36% | -0.204 | -0.261 |
| 80 | 3 | 3.5 MHz | 2.2753e+08 | 6.32% | 93.45% | -0.202 | -0.209 |
| 90 | 3 | 3.5 MHz | 2.0912e+08 | 6.21% | 96.89% | -0.216 | -0.185 |
| 100 | 3 | 3.5 MHz | 1.9811e+08 | 6.14% | 93.09% | -0.205 | -0.206 |
| 110 | 3 | 3.5 MHz | 1.7576e+08 | 6.19% | 92.76% | -0.199 | -0.221 |
| 120 | 3 | 3.5 MHz | 1.7657e+08 | 6.11% | 92.87% | -0.210 | -0.249 |
| 130 | 3 | 3.5 MHz | 1.7126e+08 | 6.20% | 93.40% | -0.213 | -0.267 |
| 140 | 3 | 3.5 MHz | 1.7843e+08 | 6.14% | 94.63% | -0.211 | -0.260 |

