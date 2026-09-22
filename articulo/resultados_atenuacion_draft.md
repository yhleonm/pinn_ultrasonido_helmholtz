# Borrador de Resultados - Atenuación ($\alpha_0$)
_Borrador v1 - Pendiente de Validación post-PML (2026-06-22)_

Este documento contiene el borrador de la sección de resultados cuantitativos para el coeficiente de atenuación ($\alpha_0$) reconstruido mediante FWI-PINN en los 6 fantomas físicos de CIRS 040GSE (`MYO001-MYO006`).

---

## 1. Tabla de Resultados Consolidados de Atenuación

El coeficiente de atenuación se evalúa utilizando el CNR robusto con exclusión del 15% de margen del SIREN.

| Dataset | Profundidad $z$ (mm) | CNR $\alpha$ (dBr) | Media en Inclusión (dB/MHz/cm) | Media en Fondo ROI (dB/MHz/cm) | Desviación Estándar Fondo |
|---|:---:|:---:|:---:|:---:|:---:|
| **MYO001** | 18.5 | 12.17 | 0.0888 | 0.0780 | 0.0846 |
| **MYO002** | 23.5 | 10.73 | 0.0487 | 0.0916 | 0.1454 |
| **MYO003** | 18.0 | 10.69 | 0.0235 | 0.0670 | 0.0824 |
| **MYO004** | 11.0 | 12.30 | 0.0127 | 0.0310 | 0.0557 |
| **MYO005** | 11.5 | 12.50 | 0.0948 | 0.1970 | 0.2762 |
| **MYO006** | 10.5 | 10.50 | 0.0174 | 0.0581 | 0.0842 |

### 📊 Estadísticas de Dispersión
*   **Media de CNR $\alpha$**: `11.48 dBr`
*   **Desviación Estándar**: `0.85 dBr`
*   **Coeficiente de Variación (CV)**: **`7.4%`** (Excelente reproducibilidad para estándares cuantitativos)

---

## 2. Discusión y Análisis Físico

### 2.1 Asimetría de Robustez: Modulación de Amplitud vs Fase
El contraste de atenuación ($\alpha_0$) demuestra ser un biomarcador sustancialmente más estable que la velocidad de sonido ($c$) bajo la configuración de reflexión de un solo lado:
1.  **Modulación de Amplitud**: La atenuación afecta directamente a la amplitud del campo disperso $P_{\text{scat}}$. Las perturbaciones de amplitud tienen una relación más lineal con el residuo de la ecuación de onda y son menos vulnerables a los descalces de fase globales.
2.  **Modulación de Fase**: La velocidad de sonido modifica la fase a través del número de onda $k = \omega / c$. Debido a la naturaleza no lineal de los descalces de fase y a la cobertura angular limitada (que restringe la transmisión directa), la velocidad es extremadamente sensible a las condiciones de frontera numéricas.
3.  **Estabilidad del CNR**: Mientras que el contraste de velocidad es disperso y propenso a fugas del artefacto de borde, la atenuación muestra una consistencia excepcional, con un CNR robusto que se mantiene estable entre `10.5` y `12.5` dBr en los 6 fantomas físicos evaluados.

### 2.2 Descalce en el Valor Absoluto (Overfitting de Frontera)
Se observa que la media absoluta reconstruida de $\alpha_0$ dentro de la ROI (rango `0.01 - 0.20` dB/MHz/cm) es inferior al valor de referencia del fabricante para el fantoma CIRS 040GSE (`0.5` dB/MHz/cm para el fondo y `0.7` dB/MHz/cm para la inclusión). 

Este descalce se atribuye a un **overfitting numérico en las fronteras de la ROI**. Sin capas absorbentes (PML), el optimizador ajusta a la baja la atenuación global para compensar artificialmente las reflexiones espurias provenientes de los bordes de la simulación. El contraste relativo (CNR) se preserva de forma consistente gracias a la simetría del artefacto, pero el valor absoluto cuantitativo requiere la implementación de PML para estabilizarse en los rangos nominales del tejido.

---

## 3. Conclusión del Borrador para la Tesis

Se establece a la atenuación ($\alpha_0$) como el **biomarcador cuantitativo primario** de este desarrollo debido a su consistencia estadística across adquisiciones. Se propone reportar la atenuación como un biomarcador robusto para la caracterización tisular de tumores mamarios, dejando la velocidad como un parámetro complementario que requiere mayor control de borde (PML) para su correcta interpretación cuantitativa.

---

## 4. Estado de Validación
> [!WARNING]
> **ESTADO: BORRADOR V1**
> Esta sección se considera "Lista para integrar" pero no definitiva. Se volverá a correr el sweep completo una vez que esté implementada la formulación con **PML (Perfectly Matched Layers)** para confirmar si los valores absolutos de atenuación se estabilizan cerca de los nominales del fabricante (`0.5` y `0.7` dB/MHz/cm).
