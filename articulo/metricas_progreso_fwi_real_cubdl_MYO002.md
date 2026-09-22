# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-25 20:03:46

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 0.00 mm, z = 23.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 4.29%
*   **Error L2 alpha0 final en ROI**: 94.01%
*   **SSIM Final c**: 0.040
*   **SSIM Final alpha0**: -0.248

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO002.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO002.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2861e+07 | 2.63% | 39.14% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.5122e+06 | 3.56% | 37.71% | -0.114 | 0.138 |
| 20 | 1 | 1.0 MHz | 8.8077e+06 | 3.26% | 41.81% | 0.000 | 0.006 |
| 30 | 1 | 1.0 MHz | 8.2114e+06 | 3.70% | 44.37% | 0.010 | -0.035 |
| 31 | 2 | 2.0 MHz | 9.8240e+06 | 3.75% | 44.56% | 0.009 | -0.035 |
| 40 | 2 | 2.0 MHz | 7.3884e+06 | 3.91% | 41.69% | -0.119 | 0.143 |
| 50 | 2 | 2.0 MHz | 6.7115e+06 | 3.91% | 42.30% | -0.058 | 0.148 |
| 60 | 2 | 2.0 MHz | 6.4476e+06 | 3.94% | 41.28% | -0.008 | 0.151 |
| 61 | 3 | 3.5 MHz | 1.6291e+08 | 4.00% | 41.30% | -0.009 | 0.152 |
| 70 | 3 | 3.5 MHz | 1.5662e+08 | 3.57% | 45.74% | 0.089 | 0.059 |
| 80 | 3 | 3.5 MHz | 1.5034e+08 | 4.21% | 45.95% | -0.012 | 0.041 |
| 90 | 3 | 3.5 MHz | 1.4708e+08 | 4.46% | 57.28% | -0.008 | -0.148 |
| 100 | 3 | 3.5 MHz | 1.4284e+08 | 4.58% | 70.40% | -0.010 | -0.247 |
| 110 | 3 | 3.5 MHz | 1.3831e+08 | 4.52% | 81.14% | 0.003 | -0.265 |
| 120 | 3 | 3.5 MHz | 1.3354e+08 | 4.47% | 87.61% | 0.005 | -0.260 |
| 130 | 3 | 3.5 MHz | 1.2842e+08 | 4.36% | 91.58% | 0.018 | -0.253 |
| 140 | 3 | 3.5 MHz | 1.2887e+08 | 4.29% | 94.01% | 0.040 | -0.248 |

