"""
backend/app/services/clinical_interpretation.py
================================================
Traduce métricas cuantitativas de FWI (velocidad del sonido y atenuación)
en conclusiones en lenguaje clínico llano para el usuario (médico).
"""

from typing import Dict, Any

def interpret_fwi_metrics(dataset_id: str, metrics: Dict[str, Any], tumor_radius_mm: float = 3.0) -> Dict[str, Any]:
    """
    Toma métricas del sidecar y genera una interpretación clínica simplificada.
    """
    tumor_x = metrics.get("tumor_x_mm")
    tumor_z = metrics.get("tumor_z_mm")
    cnr_c = metrics.get("cnr_c_dbr")
    cnr_a = metrics.get("cnr_alpha_dbr")
    peak_dist = metrics.get("peak_distance_mm")
    
    # 1. Determinación de hallazgo principal (Finding)
    finding = "No se detecta inclusión o masa focalizada con contraste acústico significativo."
    
    def cnr_quality(val: float | None) -> str:
        if val is None:
            return "sin señal"
        if val > 12.0:
            return "excelente"
        if val > 6.0:
            return "significativo"
        if val > 1.5:
            return "moderado"
        return "bajo"

    if tumor_x is not None and tumor_z is not None:
        pos_str = f"a x = {tumor_x:.1f} mm, z = {tumor_z:.1f} mm"
        
        # Evaluar contraste en velocidad y atenuación
        has_speed_contrast = cnr_c is not None and cnr_c > 1.5
        has_atten_contrast = cnr_a is not None and cnr_a > 1.5
        is_velocity_unreliable = peak_dist is not None and peak_dist > 3.0
        
        if is_velocity_unreliable:
            if has_atten_contrast:
                finding = (
                    f"Masa focalizada detectada {pos_str}. Esta zona muestra atenuación elevada "
                    f"(contraste {cnr_quality(cnr_a)}, CNR α = {cnr_a:.2f} dBr); el canal de velocidad "
                    f"del sonido no pudo confirmarlo de forma independiente por una limitación técnica conocida "
                    f"(desplazamiento excesivo de {peak_dist:.2f} mm en el pico de velocidad)."
                )
            else:
                finding = (
                    f"Masa localizada {pos_str}. No se puede confirmar de forma fiable: el canal de velocidad "
                    f"del sonido presenta una desviación técnica ({peak_dist:.2f} mm) y no hay contraste de "
                    f"atenuación significativo que lo respalde."
                )
        else:
            if has_speed_contrast and has_atten_contrast:
                finding = (
                    f"Masa focalizada detectada {pos_str}. Presenta alta velocidad de sonido "
                    f"(contraste {cnr_quality(cnr_c)}, CNR c = {cnr_c:.2f} dBr) y alta atenuación "
                    f"(contraste {cnr_quality(cnr_a)}, CNR α = {cnr_a:.2f} dBr). "
                    f"Este perfil dual es compatible con una inclusión sólida de alta densidad acústica."
                )
            elif has_speed_contrast:
                finding = (
                    f"Masa focalizada detectada {pos_str}. Presenta contraste predominante en velocidad "
                    f"de sonido (contraste {cnr_quality(cnr_c)}, CNR c = {cnr_c:.2f} dBr) con atenuación normal. "
                    f"Compatible con una inclusión sólida de velocidad de sonido elevada pero baja absorción."
                )
            elif has_atten_contrast:
                finding = (
                    f"Masa focalizada detectada {pos_str}. Presenta contraste predominante en atenuación "
                    f"acústica (contraste {cnr_quality(cnr_a)}, CNR α = {cnr_a:.2f} dBr) con velocidad normal. "
                    f"Compatible con una zona altamente absorbente."
                )
            else:
                finding = f"Masa localizada {pos_str}, pero con bajo contraste acústico cuantitativo (CNR < 1.5 dBr)."

    # 2. Evaluación de fiabilidad de la reconstrucción (Reliability)
    is_reliable = True
    reliability_note = "Reconstrucción consistente. Excelente correspondencia entre los mapas físicos."
    
    if peak_dist is not None and peak_dist > 3.0:
        is_reliable = False
        reliability_note = (
            f"Advertencia de Artefacto: El pico de velocidad se encuentra a {peak_dist:.2f} mm del centro del tumor "
            f"(distancia > 3.0 mm esperada). Esto indica una fuga o artefacto numérico de borde en la "
            f"reconstrucción de velocidad del sonido (c). El mapa de atenuación (α) es físicamente más consistente."
        )
    elif cnr_c is not None and cnr_c < 1.0:
        is_reliable = False
        reliability_note = (
            "Advertencia de Baja Señal: Contraste de velocidad extremadamente bajo. "
            "Riesgo de falsos negativos en el mapa c."
        )
        
    # 3. Recomendación
    if not is_reliable:
        recommendation = "Utilizar el mapa de atenuación (α) como referencia primaria. Evitar diagnóstico cuantitativo basado únicamente en el mapa c."
    else:
        recommendation = "Correlacionar la posición de la masa con el ultrasonido convencional (B-mode) y evaluar la velocidad del sonido cuantitativa (m/s)."
        
    return {
        "finding": finding,
        "is_reliable": is_reliable,
        "reliability_note": reliability_note,
        "recommendation": recommendation
    }
