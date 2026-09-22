import type { DatasetSummary } from '../types';

interface Props {
  datasets: DatasetSummary[];
  selectedId: string | null;
  onSelect: (id: string) => void;
}

export default function DatasetList({ datasets, selectedId, onSelect }: Props) {
  if (datasets.length === 0) {
    return <p className="muted">Aún no hay reconstrucciones disponibles.</p>;
  }

  return (
    <ul className="dataset-list">
      {datasets.map((d) => (
        <li key={d.dataset_id}>
          <button
            className={`dataset-item ${selectedId === d.dataset_id ? 'is-selected' : ''}`}
            onClick={() => onSelect(d.dataset_id)}
          >
            <div className="dataset-item-top">
              <span className="dataset-id">{d.dataset_id}</span>
              <span className={`tissue-badge tissue-badge--${d.tissue_type}`}>
                {d.tissue_type === 'in_vivo' ? 'In vivo' : 'Fantoma'}
              </span>
            </div>
            {d.institution && <span className="dataset-meta">{d.institution}</span>}
            <div className="dataset-item-metrics">
              {d.cnr_alpha_dbr !== null && (
                <span className="metric-chip metric-chip--alpha">
                  α {d.cnr_alpha_dbr.toFixed(1)} dBr
                </span>
              )}
              {d.cnr_c_dbr !== null && (
                <span
                  className={`metric-chip ${d.is_reliable_velocity ? 'metric-chip--c' : 'metric-chip--warn'}`}
                  title={d.is_reliable_velocity ? undefined : 'Posible artefacto de borde en velocidad'}
                >
                  c {d.cnr_c_dbr.toFixed(1)} dBr {!d.is_reliable_velocity && '⚠'}
                </span>
              )}
            </div>
          </button>
        </li>
      ))}
    </ul>
  );
}
