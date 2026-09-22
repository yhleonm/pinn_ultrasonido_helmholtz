import { useState, useEffect } from 'react';
import DatasetList from './components/DatasetList';
import ReconstructionViewer from './components/ReconstructionViewer';
import UploadPanel from './components/UploadPanel';
import type { DatasetSummary } from './types';
import { fetchDatasets } from './api';
import './App.css';

export default function App() {
  const [datasets, setDatasets] = useState<DatasetSummary[]>([]);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [showUpload, setShowUpload] = useState(false);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    fetchDatasets()
      .then((data) => {
        setDatasets(data);
        if (data.length > 0) setSelectedId(data[0].dataset_id);
      })
      .finally(() => setLoading(false));
  }, []);

  return (
    <div className="app-shell">
      <header className="app-header">
        <div className="header-brand">
          <span className="brand-mark" aria-hidden="true">λ</span>
          <div>
            <h1>Medical Insight 3D</h1>
            <p className="header-subtitle">Inversión física de onda completa · Fase 2</p>
          </div>
        </div>
        <button className="btn-primary" onClick={() => setShowUpload(true)}>
          Nueva reconstrucción
        </button>
      </header>

      <div className="app-body">
        <aside className="sidebar">
          <h2 className="sidebar-title">Adquisiciones</h2>
          {loading ? (
            <p className="muted">Cargando lista...</p>
          ) : (
            <DatasetList
              datasets={datasets}
              selectedId={selectedId}
              onSelect={setSelectedId}
            />
          )}
        </aside>

        <main className="main-panel">
          {selectedId ? (
            <ReconstructionViewer datasetId={selectedId} />
          ) : (
            <div className="empty-state">
              <p>Selecciona una adquisición para ver su reconstrucción.</p>
            </div>
          )}
        </main>
      </div>

      {showUpload && (
        <UploadPanel
          onClose={() => setShowUpload(false)}
          onJobStarted={() => setShowUpload(false)}
        />
      )}
    </div>
  );
}
