# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-22 22:24:02

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = -16.50 mm, z = 10.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.88%
*   **Error L2 alpha0 final en ROI**: 92.40%
*   **SSIM Final c**: -0.106
*   **SSIM Final alpha0**: -0.243

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO006.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO006.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 2.3662e+07 | 3.02% | 43.64% | -0.288 | 0.040 |
| 10 | 1 | 1.0 MHz | 1.4297e+07 | 4.77% | 42.53% | -0.276 | 0.089 |
| 20 | 1 | 1.0 MHz | 1.2691e+07 | 5.65% | 58.77% | -0.252 | -0.274 |
| 30 | 1 | 1.0 MHz | 1.1325e+07 | 6.25% | 74.55% | -0.271 | -0.333 |
| 31 | 2 | 2.0 MHz | 1.2483e+07 | 6.32% | 76.00% | -0.270 | -0.331 |
| 40 | 2 | 2.0 MHz | 1.1231e+07 | 7.52% | 86.21% | -0.281 | -0.401 |
| 50 | 2 | 2.0 MHz | 9.9194e+06 | 7.53% | 91.61% | -0.274 | -0.361 |
| 60 | 2 | 2.0 MHz | 8.7918e+06 | 7.73% | 98.18% | -0.259 | -0.262 |
| 61 | 3 | 3.5 MHz | 2.4555e+08 | 7.69% | 98.87% | -0.262 | -0.254 |
| 70 | 3 | 3.5 MHz | 2.4565e+08 | 8.11% | 97.00% | -0.209 | -0.219 |
| 80 | 3 | 3.5 MHz | 2.4309e+08 | 7.71% | 90.52% | -0.164 | -0.168 |
| 90 | 3 | 3.5 MHz | 2.4444e+08 | 6.20% | 78.60% | -0.053 | 0.056 |
| 100 | 3 | 3.5 MHz | 2.4296e+08 | 5.93% | 79.50% | -0.058 | 0.038 |
| 110 | 3 | 3.5 MHz | 2.4506e+08 | 6.09% | 87.19% | -0.101 | -0.134 |
| 120 | 3 | 3.5 MHz | 2.3039e+08 | 5.91% | 90.40% | -0.081 | -0.163 |
| 130 | 3 | 3.5 MHz | 2.2757e+08 | 5.97% | 92.41% | -0.078 | -0.186 |
| 140 | 3 | 3.5 MHz | 2.0941e+08 | 5.88% | 92.40% | -0.106 | -0.243 |

