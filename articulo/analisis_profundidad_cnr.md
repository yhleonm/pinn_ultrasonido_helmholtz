# Análisis de Dependencia de Profundidad (z) vs Contraste (CNR) - Definitivo
_Medical Insight 3D - Phase 2_

Este reporte analiza la relación entre la profundidad del tumor y el CNR robusto (margen 15%) tras realizar una **verificación física rigurosa y uniforme** sobre los 6 datasets del fantoma CIRS 040GSE (`MYO001-MYO006`).

---

## 1. Criterio de Verificación Física de Fuga de Borde

Para evitar clasificaciones inconsistentes, se aplicó un criterio objetivo y uniforme a las 6 reconstrucciones:
1.  **Exclusión de Borde**: Se recortó el 15% del margen de la ROI (interior de 10.5 x 10.5 mm).
2.  **Ubicación del Pico**: Se localizó el punto con el máximo valor de velocidad ($c$) dentro del interior.
3.  **Distancia Radial**: Se calculó la distancia física en milímetros desde este pico hasta el centro del tumor (0.0, 0.0) en la ROI.
    *   **LIMPIO**: Si la distancia es $\le 3.0$ mm (el pico máximo está dentro del radio físico del tumor).
    *   **FUGA**: Si la distancia es $> 3.0$ mm (el pico máximo está en la periferia, dominado por el artefacto de borde numérico).

---

## 2. Tabla Consolidada (Criterio Uniforme)

| Dataset | Profundidad z (mm) | CNR c (dBr) | Distancia del Pico al Centro (mm) | Clasificación Física | Comportamiento Reconstruido |
|---|---|---|---|---|---|
| **MYO001** | 18.5 | 1.51 | 5.30 mm | ❌ **FUGA** | Pico en el borde del interior |
| **MYO002** | 23.5 | 2.97 | 6.04 mm | ❌ **FUGA** | Pico en el borde del interior |
| **MYO003** | 18.0 | 6.26 | 5.43 mm | ❌ **FUGA** | Pico en el borde del interior |
| **MYO004** | 11.0 | 9.96 | 6.57 mm | ❌ **FUGA** | Pico en el borde del interior (Fuga severa) |
| **MYO005** | 11.5 | 3.53 | **3.00 mm** | **✅ LIMPIO** | **Pico central (dentro de la inclusión)** |
| **MYO006** | 10.5 | 9.20 | 5.32 mm | ❌ **FUGA** | Pico en el borde del interior |

---

## 3. Quiebre de la Hipótesis de Regresión Lineal

El análisis detallado e individual de las ROIs revela una conclusión crucial:

> [!CAUTION]
> **La correlación $r = -0.864$ era un espejismo numérico.**
> Los valores más altos de CNR (como **MYO004** con `9.96 dBr` y **MYO006** con `9.20 dBr`) no representan reconstrucciones físicas de alta calidad por su poca profundidad ($z \approx 11$ mm). Sus picos están ubicados a `6.57 mm` y `5.32 mm` del centro, respectivamente, indicando que el artefacto de borde del SIREN se filtró masivamente al interior del área evaluada, inflando de forma artificial el CNR.
> 
> Al contar con **solo 1 dataset verificado como limpio (MYO005)**, es estadísticamente imposible calcular una correlación de Pearson confiable ($N=1$).

---

## 4. Nueva Narrativa y Conclusiones para la Tesis

Esta verificación completa cambia por completo la perspectiva del artículo, aportando el rigor científico que un jurado de defensa o revisor de revista internacional exigirá:

1.  **Atenuación como Biomarcador Sólido:** La atenuación ($\alpha_0$) se mantiene reproducible y libre de dependencias espaciales complejas (CNR medio = 11.48 dBr, CV = 7.4%).
2.  **Velocidad Limitada por Artefactos Numéricos:** La dispersión en el CNR de velocidad ($c$) no es una degradación puramente física por cobertura de profundidad. Es la consecuencia de una **fuga no uniforme del artefacto de borde de SIREN** que contamina 5 de los 6 datasets evaluados.
3.  **Prioridad Absoluta de PML:** La prueba definitiva para el método FWI no es medir la dependencia con la profundidad bajo las condiciones actuales. Es **implementar capas de absorción PML** en la simulación directa y adjoint. 
    *   Si tras implementar PML la correlación con la profundidad aparece de forma limpia, se confirmará un límite físico de la geometría de reflexión.
    *   Si tras PML las velocidades se uniformizan, se confirmará que la dispersión era enteramente un artefacto numérico de la ROI.

*Reporte detallado de inspección visual:* [inspect_myo_comparison.md](file:///C:/Users/yelem/.gemini/antigravity-cli/brain/c09ad2ae-e67e-40d7-9db1-8b6956d3e523/inspect_myo_comparison.md)
