# Estado Final - Investigación de Correlación CNR-Profundidad (Cerrado)

_Medical Insight 3D - Fase 2 - Registro metodológico para la tesis_
_Cierre: 2026-06-22_

---

## Resumen ejecutivo

Se investigó si el contraste de velocidad de sonido (CNR_c) reconstruido por FWI-PINN
depende sistemáticamente de la profundidad de la inclusión, usando 6 adquisiciones del
mismo fantoma físico (CIRS 040GSE, datasets MYO001-MYO006). La hipótesis inicial de una
correlación fuerte (r=0.864) fue **retractada** tras verificación exhaustiva: se
determinó que es un artefacto del criterio de detección usado, no un hallazgo físico válido.

**Conclusión final, la única sostenible con la evidencia disponible:**
> El artefacto numérico de frontera en la parametrización SIREN contamina 5 de 6
> reconstrucciones de velocidad (83%), impidiendo cualquier análisis estadístico de
> dependencia con la profundidad hasta implementar condiciones de frontera absorbentes (PML).

---

## Cronología de la verificación (registro de trazabilidad)

| Paso | Criterio usado | Resultado | Estado |
|---|---|---|---|
| 1 | CNR robusto P90/P50 con margen 15% | CNR_c disperso (1.51 - 9.96 dBr), CNR_α consistente (10.50 - 12.50 dBr) | Observación inicial |
| 2 | Binario "pico cerca de borde" (sí/no) | 5/6 marcados "sí", MYO003 marcado "no" | **Criterio insuficiente - ver paso 4** |
| 3 | Correlación z vs CNR_c con los 6 puntos sin filtrar | r = -0.663 (global), r = -0.864 (excluyendo MYO005) | **Retractado - ver paso 4** |
| 4 | Distancia radial física del pico al centro del tumor (umbral: radio inclusión = 3.0 mm) | MYO003 reclasificado de "limpio" a "FUGA" (5.43 mm > 3.0 mm) | **Criterio definitivo** |

El paso 4 reveló que el criterio binario del paso 2 tuvo un **falso negativo** en MYO003,
que parecía limpio pero en realidad tenía fuga de borde tan severa como MYO004 y MYO006.

---

## Tabla final verificada - 6/6 datasets con criterio uniforme

| Dataset | z (mm) | CNR_c (dBr) | Distancia pico-centro (mm) | Diagnóstico |
|---|---|---|---|---|
| MYO001 | 18.5 | 1.51 | 5.30 | ❌ Fuga de borde |
| MYO002 | 23.5 | 2.97 | 6.04 | ❌ Fuga de borde |
| MYO003 | 18.0 | 6.26 | 5.43 | ❌ Fuga de borde |
| MYO004 | 11.0 | 9.96 | 6.57 | ❌ Fuga de borde (severa) |
| MYO005 | 11.5 | **3.53** | **3.00** | **✅ Limpio - único caso válido** |
| MYO006 | 10.5 | 9.20 | 5.32 | ❌ Fuga de borde |

**Criterio de clasificación:** distancia radial del pico de velocidad reconstruido al
centro DAS del tumor. Umbral = radio físico de la inclusión (3.0 mm, CIRS 040GSE).
Distancia > umbral -> el "pico" detectado no corresponde a la inclusión real, sino al
artefacto de frontera de la parametrización SIREN.

---

## Por qué no se puede calcular correlación con n=1

Con un único dataset confirmado libre de artefacto (MYO005), cualquier coeficiente de
correlación de Pearson calculado sería matemáticamente indefinido o estadísticamente vacío
(se requiere n>=3 como mínimo absoluto, n>=8 para potencia estadística razonable en estudios
exploratorios). Reportar un número de correlación en este escenario sería metodológicamente
incorrecto, independientemente de qué tan plausible sea la hipótesis física subyacente
(que sí tiene fundamento teórico en la pérdida de cobertura angular con la profundidad
en esquemas de reflexión de un solo lado - pero esa hipótesis requiere datos limpios para
ser probada, no datos contaminados).

---

## Lo que SÍ queda establecido con esta evidencia (afirmaciones sólidas para el artículo)

1. **El contraste físico de velocidad existe y es positivo en todos los casos** (rango 1.51 - 9.96 dBr, sin signos invertidos), confirmando que el parche de restricción de fase de β corrigió la ambigüedad de contraste invertido documentada en iteraciones previas.

2. **La atenuación α es un biomarcador sustancialmente más robusto que la velocidad** bajo esta geometría de adquisición: CV = 7.4% en α frente a una contaminación del 83% de los casos en c. Esto es consistente con el mecanismo físico ya identificado: α modula amplitud (relación directa, robusta a ambigüedad de fase), mientras c modula fase a través de k = ω/c (vulnerable a artefactos numéricos de frontera y a la ambigüedad de calibración de fase β).

3. **El artefacto de frontera SIREN es sistemático, no aleatorio** - aparece en 5/6 casos independientemente de la profundidad o posición lateral del tumor, lo cual señala una causa estructural (ausencia de condiciones de frontera absorbentes) más que ruido estocástico de entrenamiento.

4. **PML deja de ser una mejora opcional y se convierte en un requisito metodológico** para cualquier reporte cuantitativo futuro de CNR de velocidad. Sin PML, el 83% de las reconstrucciones de c no son confiables para análisis de correlación espacial o de dispersión, aunque sí lo son para demostración cualitativa de viabilidad (el contraste positivo en todos los casos sigue siendo evidencia de que el método localiza física real, solo que con magnitud contaminada).

---

## Redacción sugerida para la sección de Limitaciones del artículo

> "Si bien el contraste de velocidad reconstruido fue positivo en la totalidad de las
> adquisiciones evaluadas (n=6), la inspección de la posición radial del pico de
> reconstrucción reveló que en el 83% de los casos (5/6) el valor máximo se ubicaba a una
> distancia del centro de la inclusión mayor al radio físico conocido del fantoma,
> indicando contaminación por un artefacto numérico de frontera en la parametrización
> implícita (SIREN) del campo de velocidad. Este hallazgo impide, en el estado actual del
> método, establecer relaciones cuantitativas entre el contraste reconstruido y variables
> como la profundidad de la inclusión. La implementación de condiciones de frontera
> absorbentes (Perfectly Matched Layers) se identifica como el siguiente paso metodológico
> necesario para resolver esta limitación antes de extender el análisis cuantitativo de
> velocidad a cohortes más amplias."

---

## Próximos pasos (orden de prioridad actualizado)

1. **PML** - ya no es "deseable", es bloqueante para cualquier reporte cuantitativo de c más allá del actual.
2. Tras PML, **repetir el sweep completo MYO001-006** con el mismo criterio de distancia radial - si la proporción limpio/fuga mejora sustancialmente, valida que PML resuelve el problema raíz.
3. Solo entonces, **reintentar el análisis de correlación con profundidad** con un n verificado suficiente.
4. La atenuación α puede avanzar a la sección de Resultados del artículo **ya mismo**, sin esperar PML - su robustez está sólidamente establecida.
