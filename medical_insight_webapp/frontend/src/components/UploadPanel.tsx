import { useState, useEffect, useRef } from 'react';
import { uploadNewJob, fetchJobStatus } from '../api';
import type { JobStatusDetail } from '../types';

interface Props {
  onClose: () => void;
  onJobStarted: () => void;
}

export default function UploadPanel({ onClose, onJobStarted }: Props) {
  const [datasetId, setDatasetId] = useState('');
  const [file, setFile] = useState<File | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);
  const [status, setStatus] = useState<JobStatusDetail | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const pollRef = useRef<ReturnType<typeof setInterval> | null>(null);

  useEffect(() => {
    if (!jobId) return;
    pollRef.current = setInterval(() => {
      fetchJobStatus(jobId).then(setStatus).catch(() => {});
    }, 4000);
    return () => {
      if (pollRef.current) clearInterval(pollRef.current);
    };
  }, [jobId]);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (!file || !datasetId.trim()) return;
    setSubmitting(true);
    setError(null);
    try {
      const res = await uploadNewJob(datasetId.trim(), file);
      setJobId(res.job_id);
      onJobStarted();
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
    } finally {
      setSubmitting(false);
    }
  }

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal-panel" onClick={(e) => e.stopPropagation()}>
        <div className="modal-header">
          <h2>Nueva reconstrucción</h2>
          <button className="btn-icon" onClick={onClose} aria-label="Cerrar">✕</button>
        </div>

        {!jobId ? (
          <form onSubmit={handleSubmit} className="upload-form">
            <p className="upload-warning">
              Un cálculo completo toma entre 20 y 60 minutos. El trabajo corre
              en segundo plano; puedes cerrar esta ventana y volver más tarde
              para revisar el estado.
            </p>

            <label className="field-label" htmlFor="dataset-id">
              Identificador de la adquisición
            </label>
            <input
              id="dataset-id"
              type="text"
              placeholder="ej. MYO007"
              value={datasetId}
              onChange={(e) => setDatasetId(e.target.value)}
              required
            />

            <label className="field-label" htmlFor="rf-file">
              Archivo de canal RF (.h5, .hdf5, .mat)
            </label>
            <input
              id="rf-file"
              type="file"
              accept=".h5,.hdf5,.mat"
              onChange={(e) => setFile(e.target.files?.[0] ?? null)}
              required
            />

            {error && <p className="error-text">{error}</p>}

            <button className="btn-primary" type="submit" disabled={submitting}>
              {submitting ? 'Enviando...' : 'Iniciar reconstrucción'}
            </button>
          </form>
        ) : (
          <div className="job-status">
            <p><strong>Trabajo encolado:</strong> {jobId}</p>
            {status ? (
              <>
                <p><strong>Estado:</strong> {statusLabel(status.status)}</p>
                {status.queue_position > 0 && (
                  <p>Posición en cola: {status.queue_position}</p>
                )}
                {status.error && <p className="error-text">{status.error}</p>}
                <pre className="log-tail">{status.log_tail.join('\n')}</pre>
              </>
            ) : (
              <p className="muted">Consultando estado...</p>
            )}
            <button className="btn-secondary" onClick={onClose}>Cerrar y ver más tarde</button>
          </div>
        )}
      </div>
    </div>
  );
}

function statusLabel(status: string): string {
  switch (status) {
    case 'queued': return 'En cola';
    case 'running': return 'Calculando (esto toma tiempo)';
    case 'done': return 'Completado';
    case 'failed': return 'Falló';
    default: return status;
  }
}
