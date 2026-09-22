import sys
from pathlib import Path

# Agregar backend al path
sys.path.append(str(Path(__file__).parent.parent / "medical_insight_webapp" / "backend"))

from app.services.interpretation import build_interpretation
import json

articulo_dir = Path(__file__).parent.parent / "articulo"

for ds in ["MYO005", "MYO006"]:
    metrics_path = articulo_dir / f"metrics_{ds}.json"
    with open(metrics_path, 'r') as f:
        metrics = json.load(f)
        
    interpretation = build_interpretation(
        tumor_x_mm=metrics["tumor_x_mm"],
        tumor_z_mm=metrics["tumor_z_mm"],
        cnr_c_dbr=metrics["cnr_c_dbr"],
        cnr_alpha_dbr=metrics["cnr_alpha_dbr"],
        peak_distance_mm=metrics["peak_distance_mm"],
        peak_distance_alpha_mm=metrics["peak_distance_alpha_mm"]
    )
    
    print(f"=== {ds} ===")
    print(f"  Reliability Level: {interpretation.reliability_level}")
    print(f"  Finding Summary: {interpretation.finding_summary}")
    print(f"  Reliability Reason: {interpretation.reliability_reason}")
    print(f"  Velocity Reason: {interpretation.reliability_reason_velocity}")
    print(f"  Attenuation Reason: {interpretation.reliability_reason_attenuation}")
    print()
