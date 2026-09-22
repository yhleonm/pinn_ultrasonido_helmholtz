# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-25 21:02:42

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 16.00 mm, z = 18.00 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 6.50%
*   **Error L2 alpha0 final en ROI**: 96.37%
*   **SSIM Final c**: 0.040
*   **SSIM Final alpha0**: -0.303

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO003.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO003.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.4893e+07 | 3.00% | 43.42% | -0.287 | 0.038 |
| 10 | 1 | 1.0 MHz | 1.4322e+07 | 2.70% | 56.68% | 0.258 | -0.187 |
| 20 | 1 | 1.0 MHz | 1.4267e+07 | 3.08% | 62.68% | 0.156 | -0.148 |
| 30 | 1 | 1.0 MHz | 1.4117e+07 | 3.20% | 68.81% | 0.160 | -0.106 |
| 31 | 2 | 2.0 MHz | 9.2941e+06 | 3.22% | 69.27% | 0.153 | -0.095 |
| 40 | 2 | 2.0 MHz | 8.9102e+06 | 4.28% | 74.01% | 0.125 | -0.142 |
| 50 | 2 | 2.0 MHz | 8.7502e+06 | 5.22% | 69.61% | 0.099 | -0.147 |
| 60 | 2 | 2.0 MHz | 8.6651e+06 | 5.63% | 64.53% | 0.101 | -0.227 |
| 61 | 3 | 3.5 MHz | 1.6867e+08 | 5.69% | 64.24% | 0.096 | -0.227 |
| 70 | 3 | 3.5 MHz | 1.6782e+08 | 4.73% | 67.63% | 0.204 | -0.263 |
| 80 | 3 | 3.5 MHz | 1.6417e+08 | 5.42% | 75.15% | 0.135 | -0.202 |
| 90 | 3 | 3.5 MHz | 1.6631e+08 | 5.76% | 81.37% | 0.100 | -0.243 |
| 100 | 3 | 3.5 MHz | 1.6215e+08 | 5.80% | 80.27% | 0.032 | -0.177 |
| 110 | 3 | 3.5 MHz | 1.5858e+08 | 6.17% | 83.43% | 0.026 | -0.170 |
| 120 | 3 | 3.5 MHz | 1.5773e+08 | 6.28% | 90.46% | 0.020 | -0.263 |
| 130 | 3 | 3.5 MHz | 1.5611e+08 | 6.42% | 94.46% | 0.036 | -0.292 |
| 140 | 3 | 3.5 MHz | 1.5455e+08 | 6.50% | 96.37% | 0.040 | -0.303 |

