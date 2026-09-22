"""
backend/app/services/interpretation.py
=========================================
Traduce metricas fisicas (CNR en dBr, distancia en mm) a lenguaje llano
para un medico, con TRES restricciones de seguridad clinica no negociables:

1. NUNCA se combina un canal como evidencia positiva si su pico esta fuera
   del radio fisico esperado (artefacto de borde o limitacion de localizacion
   desconocida hasta ahora). En ese caso, el hallazgo se reporta SOLO desde
   el canal espacialmente consistente, explicitamente.

2. NUNCA se sugieren nombres de patologias (fibroadenoma, lesion maligna,
   quiste, etc.). Esta herramienta no tiene validacion contra biopsia,
   histologia, ni ningun estudio clinico -- solo fantomas fisicos y una
   adquisicion in vivo sin ground truth. Sugerir diagnosticos diferenciales
   es una afirmacion que ningun resultado de este proyecto respalda.

3. NUNCA se usa terminologia de otras modalidades clinicas ya establecidas
   (ej. "rigidez", que pertenece a elastografia) para nombrar cantidades
   fisicas distintas (velocidad de sonido compresional). Cada canal se
   nombra por lo que mide, no por una metafora que evoque otra tecnica.

Estas restricciones existen porque este es un proyecto de investigacion en
fase de validacion de metodo (n=6 fantomas, 1 adquisicion in vivo sin
ground truth), no una herramienta diagnostica validada clinicamente.
"""

from __future__ import annotations

from app.models.schemas import ClinicalInterpretation, ReliabilityLevel

TUMOR_RADIUS_MM = 3.0
CNR_USEFUL_THRESHOLD = 6.0
CNR_EXCELLENT_THRESHOLD = 12.0


def _velocity_is_reliable(peak_distance_mm: float | None) -> bool:
    """
    Criterio de fiabilidad del canal de velocidad: si el pico
    reconstruido esta dentro del radio fisico esperado del hallazgo.
    """
    if peak_distance_mm is None:
        return False
    return peak_distance_mm <= TUMOR_RADIUS_MM


def _attenuation_is_reliable(peak_distance_alpha_mm: float | None) -> bool:
    """
    Criterio de fiabilidad del canal de atenuacion: si el pico
    reconstruido esta dentro del radio fisico esperado del hallazgo.
    """
    if peak_distance_alpha_mm is None:
        return False
    return peak_distance_alpha_mm <= TUMOR_RADIUS_MM


def _location_phrase(tumor_x_mm: float, tumor_z_mm: float) -> str:
    side = "derecha" if tumor_x_mm >= 0 else "izquierda"
    return f"a {abs(tumor_x_mm):.0f} mm hacia la {side} del centro y {tumor_z_mm:.0f} mm de profundidad"


def build_interpretation(
    tumor_x_mm: float | None,
    tumor_z_mm: float | None,
    cnr_c_dbr: float | None,
    cnr_alpha_dbr: float | None,
    peak_distance_mm: float | None,
    peak_distance_alpha_mm: float | None,
) -> ClinicalInterpretation:
    """
    Punto unico de generacion de texto clinico. Cualquier cambio a como
    se describe un hallazgo o su nivel de confianza debe hacerse aqui,
    no en el frontend ni en ningun otro punto del backend.
    """
    if tumor_x_mm is None or tumor_z_mm is None:
        return ClinicalInterpretation(
            finding_summary="No se identifico una zona de interes localizable en esta adquisicion.",
            reliability_level=ReliabilityLevel.LOW,
            reliability_reason="No hay una localizacion espacial valida para evaluar.",
        )

    location = _location_phrase(tumor_x_mm, tumor_z_mm)
    velocity_ok = _velocity_is_reliable(peak_distance_mm)
    atten_ok = _attenuation_is_reliable(peak_distance_alpha_mm)

    if peak_distance_mm is None:
        vel_reason = "No se pudo verificar la posicion del hallazgo en el canal de velocidad."
    elif velocity_ok:
        vel_reason = (
            f"El punto de mayor contraste en velocidad de sonido esta dentro del "
            f"radio fisico esperado ({peak_distance_mm:.1f} mm de {TUMOR_RADIUS_MM:.0f} mm)."
        )
    else:
        vel_reason = (
            f"El punto de mayor contraste en velocidad de sonido aparece a "
            f"{peak_distance_mm:.1f} mm de la zona esperada -- mas lejos de lo "
            f"fisicamente razonable ({TUMOR_RADIUS_MM:.0f} mm). Esto corresponde a "
            f"un artefacto numerico de frontera o desfase local en el metodo. "
            f"Este canal se excluye como evidencia en este caso."
        )

    if peak_distance_alpha_mm is None:
        att_reason = "No se pudo verificar la posicion del hallazgo en el canal de atenuacion."
    elif atten_ok:
        att_reason = (
            f"El punto de mayor contraste en atenuacion acustica esta dentro del "
            f"radio fisico esperado ({peak_distance_alpha_mm:.1f} mm de {TUMOR_RADIUS_MM:.0f} mm)."
        )
    else:
        att_reason = (
            f"El punto de mayor contraste en atenuacion acustica aparece a "
            f"{peak_distance_alpha_mm:.1f} mm de la zona esperada -- mas lejos de lo "
            f"fisicamente razonable ({TUMOR_RADIUS_MM:.0f} mm). Esto indica una "
            f"limitacion de localizacion espacial en la inversion de atenuacion. "
            f"Este canal se excluye como evidencia en este caso."
        )

    has_speed_signal = cnr_c_dbr is not None and cnr_c_dbr > 1.5
    has_atten_signal = cnr_alpha_dbr is not None and cnr_alpha_dbr > 1.5

    if velocity_ok and atten_ok and has_speed_signal and has_atten_signal:
        level = ReliabilityLevel.HIGH
        finding = (
            f"Se identifica una zona con propiedades acusticas distintas al tejido "
            f"circundante, {location}. Tanto la velocidad de sonido como la "
            f"atenuacion acustica muestran contraste y localizacion consistentes."
        )
    elif atten_ok and has_atten_signal:
        level = ReliabilityLevel.MEDIUM
        finding = (
            f"Se identifica una zona con atenuacion acustica elevada {location}. "
            f"El canal de velocidad de sonido no pudo confirmar esta localizacion "
            f"de forma independiente por una limitacion tecnica conocida del metodo."
        )
    elif velocity_ok and has_speed_signal:
        level = ReliabilityLevel.MEDIUM
        finding = (
            f"Se identifica una zona con velocidad de sonido elevada {location}. "
            f"El canal de atenuacion acustica no pudo confirmar esta localizacion "
            f"de forma independiente por una limitacion tecnica conocida del metodo."
        )
    else:
        level = ReliabilityLevel.LOW
        finding = (
            f"Se marca una zona de referencia {location} a partir de la localizacion "
            f"inicial, pero ningun canal fisico confirma un contraste y localizacion "
            f"fiables en ella."
        )

    return ClinicalInterpretation(
        finding_summary=finding,
        reliability_level=level,
        reliability_reason=vel_reason if velocity_ok else att_reason,
        reliability_reason_velocity=vel_reason,
        reliability_reason_attenuation=att_reason,
        velocity_reliable=velocity_ok,
        attenuation_reliable=atten_ok,
    )
