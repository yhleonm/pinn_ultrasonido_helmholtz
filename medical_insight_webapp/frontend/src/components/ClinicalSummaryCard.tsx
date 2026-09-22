import { ClinicalInterpretation } from '../types';

interface Props {
  interpretation: ClinicalInterpretation | null;
}

export default function ClinicalSummaryCard({ interpretation }: Props) {
  if (!interpretation) {
    return (
      <div className="clinical-summary-card">
        <p className="muted">No hay interpretación clínica disponible para este dataset.</p>
      </div>
    );
  }

  const {
    finding_summary,
    reliability_level,
    reliability_reason,
    reliability_reason_velocity,
    reliability_reason_attenuation,
  } = interpretation;

  let badgeClass = 'reliability-badge--low';
  let badgeText = '▲ Fiabilidad Baja';
  let isWarning = true;

  if (reliability_level === 'alta') {
    badgeClass = 'reliability-badge--high';
    badgeText = '● Fiabilidad Alta';
    isWarning = false;
  } else if (reliability_level === 'media') {
    badgeClass = 'reliability-badge--medium';
    badgeText = '◆ Fiabilidad Media';
    isWarning = true;
  }

  return (
    <div className="clinical-summary-card">
      <div className="clinical-card-header">
        <h3 className="clinical-card-title">Resumen Clínico (FWI)</h3>
        <span className={`reliability-badge ${badgeClass}`}>
          {badgeText}
        </span>
      </div>

      <div className="finding-section">
        <h4>Hallazgo Principal</h4>
        <div className={`finding-box ${isWarning ? 'finding-box--warning' : ''}`}>
          {finding_summary}
        </div>
      </div>

      <div className="recommendation-section">
        <h4>Evaluación de Fiabilidad</h4>
        <div className="recommendation-box">
          <p>{reliability_reason}</p>
          {(reliability_reason_velocity || reliability_reason_attenuation) && (
            <div style={{ marginTop: '0.75rem', fontSize: '0.82rem', opacity: 0.85 }}>
              {reliability_reason_velocity && (
                <div style={{ marginBottom: '0.4rem' }}>
                  <strong>Velocidad del sonido:</strong> {reliability_reason_velocity}
                </div>
              )}
              {reliability_reason_attenuation && (
                <div>
                  <strong>Atenuación acústica:</strong> {reliability_reason_attenuation}
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
