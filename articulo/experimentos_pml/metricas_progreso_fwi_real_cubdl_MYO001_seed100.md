# Integración FWI con CUBDL y Datos de Canal (Tarea C)
Fecha de Ejecución: 2026-06-27 16:56:03

## 1. Resumen de Características de la Integración
Se implementó un pipeline completo de inversión de forma de onda completa (FWI) preparado para la ingesta de datos reales CUBDL (RF canal). 

### Componentes Clave:
1.  **Beamforming Delay-and-Sum (DAS)**: Construye la imagen B-mode a partir del dataset de RF de canal y localiza automáticamente el centro del tumor mediante búsqueda del pico de envolvente (excluyendo el campo cercano).
2.  **Recorte de ROI Automático**: Define una ventana de 1.5 cm x 1.5 cm centrada en el tumor. La optimización mediante SIREN se limita exclusivamente a esta ROI, mientras que el resto del dominio permanece en los valores de fondo físico ($c_0=1540$ m/s, $\alpha_0=0.5$ dB/MHz/cm).
3.  **Remoción del Campo Incidente**: Resuelve de manera exacta la componente incidente $P_{\text{inc}}$ en un fondo homogéneo. La pérdida FWI se calcula sobre el campo disperso $P_{\text{scat}} = P_{\text{pred}} - P_{\text{inc}}$, permitiendo inversion a partir de mediciones puras de eco (sin onda directa).
4.  **Calibración Compleja de Amplitud/Fase**: Se calcula un factor $\beta$ óptimo a cada paso para compensar la ganancia y atenuación del canal físico sin requerir calibración absoluta de transductores.

## 2. Resultados en Mock CUBDL (Prueba Técnica)
*   **Localización del tumor detectado**: x = 2.00 mm, z = 18.50 mm (Teórico: 0.0 mm, 20.0 mm).
*   **Error L2 c final en ROI**: 5.29%
*   **Error L2 alpha0 final en ROI**: 93.13%
*   **SSIM Final c**: 0.195
*   **SSIM Final alpha0**: -0.241

El gráfico con los resultados de la reconstrucción y la localización de la ROI se encuentra en:
[resultados_fwi_real_cubdl_MYO001_seed100.png](file:///C:/Users/yelem/workspace/proyecto_pinn_ultrasonido_mama_fase_2/articulo/resultados_fwi_real_cubdl_MYO001_seed100.png)

## 3. Historial de Inversión por Etapas
| Época | Etapa | Frecuencia | Pérdida (MSE) | Err c (ROI) | Err alpha (ROI) | SSIM c | SSIM alpha0 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1 | 1.0 MHz | 1.2831e+07 | 2.63% | 39.23% | -0.246 | 0.038 |
| 10 | 1 | 1.0 MHz | 9.3241e+06 | 2.79% | 43.27% | -0.016 | -0.008 |
| 20 | 1 | 1.0 MHz | 8.8173e+06 | 2.81% | 44.11% | -0.088 | 0.029 |
| 30 | 1 | 1.0 MHz | 8.5574e+06 | 2.95% | 43.60% | -0.057 | 0.096 |
| 31 | 2 | 2.0 MHz | 8.8095e+06 | 2.96% | 43.53% | -0.059 | 0.100 |
| 40 | 2 | 2.0 MHz | 7.2033e+06 | 2.68% | 42.20% | 0.267 | 0.248 |
| 50 | 2 | 2.0 MHz | 6.7430e+06 | 3.03% | 40.82% | 0.280 | 0.045 |
| 60 | 2 | 2.0 MHz | 6.5984e+06 | 2.96% | 47.69% | 0.298 | -0.153 |
| 61 | 3 | 3.5 MHz | 1.3915e+08 | 2.97% | 48.31% | 0.301 | -0.163 |
| 70 | 3 | 3.5 MHz | 1.3649e+08 | 4.35% | 47.00% | 0.122 | -0.012 |
| 80 | 3 | 3.5 MHz | 1.3253e+08 | 5.15% | 62.36% | 0.183 | -0.172 |
| 90 | 3 | 3.5 MHz | 1.2891e+08 | 4.93% | 73.89% | 0.205 | -0.216 |
| 100 | 3 | 3.5 MHz | 1.2600e+08 | 5.12% | 79.11% | 0.230 | -0.255 |
| 110 | 3 | 3.5 MHz | 1.1800e+08 | 5.13% | 79.85% | 0.222 | -0.262 |
| 120 | 3 | 3.5 MHz | 1.1235e+08 | 5.20% | 85.79% | 0.209 | -0.258 |
| 130 | 3 | 3.5 MHz | 1.1302e+08 | 5.29% | 90.56% | 0.197 | -0.248 |
| 140 | 3 | 3.5 MHz | 1.0856e+08 | 5.29% | 93.13% | 0.195 | -0.241 |

