import { useState, useEffect } from 'react';
import ReconstructionCanvas from './ReconstructionCanvas';
import ConvergenceChart from './ConvergenceChart';
import ClinicalSummaryCard from './ClinicalSummaryCard';
import ErrorBoundary from './ErrorBoundary';
import type { ReconstructionMaps } from '../types';
import { fetchMaps } from '../api';

interface Props {
  datasetId: string;
}

export default function ReconstructionViewer({ datasetId }: Props) {
  const [maps, setMaps] = useState<ReconstructionMaps | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  // Estados para la vista clínica
  const [activeView, setActiveView] = useState<'clinical' | 'scientific'>('clinical');
  const [showOverlay, setShowOverlay] = useState(true);
  const [overlayScale, setOverlayScale] = useState<'velocity' | 'attenuation'>('velocity');
  const [overlayOpacity, setOverlayOpacity] = useState(0.45);
  const [isTechOpen, setIsTechOpen] = useState(false);
  const [kSigma, setKSigma] = useState(2.0);

  useEffect(() => {
    setLoading(true);
    setError(null);
    fetchMaps(datasetId)
      .then(setMaps)
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  }, [datasetId]);

  if (loading) return <p className="muted">Cargando reconstrucción de {datasetId}...</p>;
  if (error) return <p className="error-text">No se pudo cargar: {error}</p>;
  if (!maps) return null;

  const showWarningC = maps.peak_distance_mm !== null && maps.peak_distance_mm > maps.tumor_radius_mm;
  const showWarningAlpha = maps.peak_distance_alpha_mm !== null && maps.peak_distance_alpha_mm > maps.tumor_radius_mm;
  const isCurrentChannelReliable = maps.clinical_interpretation
    ? (overlayScale === 'velocity'
        ? maps.clinical_interpretation.velocity_reliable
        : maps.clinical_interpretation.attenuation_reliable)
    : true;
  const effectiveOpacity = isCurrentChannelReliable ? overlayOpacity : 0.05;

  // Renderizador de las métricas numéricas (CNR)
  const renderMetricCards = () => (
    <div className="viewer-metrics">
      {maps.cnr_c_dbr !== null && (
        <div className="metric-card">
          <span className="metric-label">CNR velocidad</span>
          <span className="metric-value">{maps.cnr_c_dbr.toFixed(2)} dBr</span>
        </div>
      )}
      {maps.cnr_alpha_dbr !== null && (
        <div className="metric-card metric-card--highlight">
          <span className="metric-label">CNR atenuación</span>
          <span className="metric-value">{maps.cnr_alpha_dbr.toFixed(2)} dBr</span>
        </div>
      )}
      {maps.peak_distance_mm !== null && (
        <div className="metric-card">
          <span className="metric-label">Dist. pico c</span>
          <span className="metric-value">{maps.peak_distance_mm.toFixed(1)} mm</span>
        </div>
      )}
      {maps.peak_distance_alpha_mm !== null && (
        <div className="metric-card">
          <span className="metric-label">Dist. pico α</span>
          <span className="metric-value">{maps.peak_distance_alpha_mm.toFixed(1)} mm</span>
        </div>
      )}
    </div>
  );

  const renderScientificMaps = () => (
    <div className="maps-grid">
      <ReconstructionCanvas
        data={maps.c_map}
        xMinMm={maps.x_min_mm} xMaxMm={maps.x_max_mm}
        zMinMm={maps.z_min_mm} zMaxMm={maps.z_max_mm}
        colorScale="velocity"
        tumorXMm={maps.tumor_x_mm} tumorZMm={maps.tumor_z_mm}
        tumorRadiusMm={maps.tumor_radius_mm}
        showWarningRing={showWarningC}
        label="Velocidad de sonido c(x,z)"
        unit="m/s"
        isReliable={!showWarningC}
        kSigma={kSigma}
      />
      <ReconstructionCanvas
        data={maps.alpha_map}
        xMinMm={maps.x_min_mm} xMaxMm={maps.x_max_mm}
        zMinMm={maps.z_min_mm} zMaxMm={maps.z_max_mm}
        colorScale="attenuation"
        tumorXMm={maps.tumor_x_mm} tumorZMm={maps.tumor_z_mm}
        tumorRadiusMm={maps.tumor_radius_mm}
        showWarningRing={showWarningAlpha}
        label="Atenuación α(x,z)"
        unit="dB/MHz/cm"
        isReliable={!showWarningAlpha}
        kSigma={kSigma}
      />
    </div>
  );

  return (
    <div className="viewer">
      <div className="viewer-header">
        <h2>{maps.dataset_id}</h2>
        
        {/* Selector de Pestañas de Vista */}
        <div className="view-tabs">
          <button
            className={`view-tab-btn ${activeView === 'clinical' ? 'is-active' : ''}`}
            onClick={() => setActiveView('clinical')}
          >
            📋 Vista Clínica
          </button>
          <button
            className={`view-tab-btn ${activeView === 'scientific' ? 'is-active' : ''}`}
            onClick={() => setActiveView('scientific')}
          >
            🔬 Vista Científica
          </button>
        </div>
      </div>

      {activeView === 'clinical' ? (
        /* VISTA CLÍNICA (PRINCIPAL) */
        <div className="clinical-view-container">
          <div className="clinical-viewer-layout">
            
            {/* Columna Izquierda: Canvas con B-Mode y Controles de Superposición */}
            <div>
              <ReconstructionCanvas
                bmodeData={maps.bmode_map}
                data={overlayScale === 'velocity' ? maps.c_map : maps.alpha_map}
                xMinMm={maps.x_min_mm} xMaxMm={maps.x_max_mm}
                zMinMm={maps.z_min_mm} zMaxMm={maps.z_max_mm}
                colorScale={overlayScale}
                showOverlay={showOverlay && !!maps.bmode_map}
                overlayOpacity={effectiveOpacity}
                tumorXMm={maps.tumor_x_mm} tumorZMm={maps.tumor_z_mm}
                tumorRadiusMm={maps.tumor_radius_mm}
                showWarningRing={overlayScale === 'velocity' ? showWarningC : showWarningAlpha}
                label={maps.bmode_map ? `Ultrasonido B-mode + Contraste (${overlayScale === 'velocity' ? 'Velocidad del sonido' : 'Atenuación'})${!isCurrentChannelReliable ? ' [Atenuado por baja fiabilidad]' : ''}` : 'Mapa de Contraste FWI'}
                findingLabel={maps.tumor_x_mm !== null ? 'Inclusión' : null}
                isReliable={isCurrentChannelReliable}
                kSigma={kSigma}
              />

              {maps.bmode_map && (
                <div className="overlay-controls">
                  <div className="control-row">
                    <label className="control-label">
                      <input
                        type="checkbox"
                        checked={showOverlay}
                        onChange={(e) => setShowOverlay(e.target.checked)}
                      />
                      Superponer Mapa FWI
                    </label>

                    {showOverlay && (
                      <div style={{ display: 'flex', gap: '0.75rem', fontSize: '0.8rem' }}>
                        <label style={{ cursor: 'pointer', display: 'flex', alignItems: 'center', gap: '4px' }}>
                          <input
                            type="radio"
                            name="overlay_scale"
                            checked={overlayScale === 'velocity'}
                            onChange={() => setOverlayScale('velocity')}
                          />
                          Velocidad del sonido
                        </label>
                        <label style={{ cursor: 'pointer', display: 'flex', alignItems: 'center', gap: '4px' }}>
                          <input
                            type="radio"
                            name="overlay_scale"
                            checked={overlayScale === 'attenuation'}
                            onChange={() => setOverlayScale('attenuation')}
                          />
                          Atenuación
                        </label>
                      </div>
                    )}
                  </div>

                  {showOverlay && (
                    <div className={`opacity-slider-container ${!isCurrentChannelReliable ? 'is-disabled' : ''}`}>
                      <span className="control-label" style={{ minWidth: '70px' }}>Transparencia:</span>
                      <input
                        type="range"
                        min="0.1"
                        max="0.8"
                        step="0.05"
                        value={isCurrentChannelReliable ? overlayOpacity : 0.05}
                        onChange={(e) => setOverlayOpacity(parseFloat(e.target.value))}
                        className="opacity-slider"
                        disabled={!isCurrentChannelReliable}
                        title={!isCurrentChannelReliable ? "Control desactivado: opacidad atenuada al 5% por baja fiabilidad" : ""}
                      />
                      <span className="opacity-val">{Math.round(effectiveOpacity * 100)}%</span>
                    </div>
                  )}

                  {showOverlay && isCurrentChannelReliable && (
                    <div className="k-slider-container" style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginTop: '0.5rem' }}>
                      <span className="control-label" style={{ minWidth: '70px' }}>Umbral (k-sigma):</span>
                      <input
                        type="range"
                        min="1.0"
                        max="3.5"
                        step="0.1"
                        value={kSigma}
                        onChange={(e) => setKSigma(parseFloat(e.target.value))}
                        className="k-slider"
                        style={{ flex: 1, accentColor: '#2ad4a8' }}
                      />
                      <span className="k-val" style={{ minWidth: '35px', textAlign: 'right', fontSize: '0.85rem' }}>{kSigma.toFixed(1)}σ</span>
                    </div>
                  )}
                </div>
              )}
            </div>

            {/* Columna Derecha: Tarjeta de Resumen Clínico */}
            <div>
              <ClinicalSummaryCard interpretation={maps.clinical_interpretation} />
            </div>
          </div>

          {/* Sección Colapsable para Detalles Técnicos al final */}
          <div className="technical-fold">
            <button
              className="technical-fold-trigger"
              onClick={() => setIsTechOpen(!isTechOpen)}
            >
              <span>{isTechOpen ? '▼ Ocultar Detalles Técnicos e Instrumentación' : '▶ Mostrar Detalles Técnicos e Instrumentación'}</span>
              <span style={{ fontSize: '0.75rem' }}>Ver CNR, heatmaps originales y convergencia</span>
            </button>
            {isTechOpen && (
              <div className="technical-fold-content">
                <ErrorBoundary>
                  {renderMetricCards()}
                </ErrorBoundary>
                <div style={{ marginTop: '0.5rem' }}>
                  <ErrorBoundary>
                    {renderScientificMaps()}
                  </ErrorBoundary>
                </div>
                <div style={{ marginTop: '0.5rem' }}>
                  <ErrorBoundary>
                    <ConvergenceChart datasetId={datasetId} />
                  </ErrorBoundary>
                </div>
              </div>
            )}
          </div>
        </div>
      ) : (
        /* VISTA CIENTÍFICA / TÉCNICA ORIGINAL */
        <div className="scientific-view-container">
          <div style={{ display: 'flex', justifyContent: 'flex-end', marginBottom: '1rem' }}>
            <ErrorBoundary>
              {renderMetricCards()}
            </ErrorBoundary>
          </div>
          {(showWarningC || showWarningAlpha) && maps.reliability_note && (
            <div className="reliability-banner" role="status">
              <strong>Nota de fiabilidad:</strong> {maps.reliability_note}
            </div>
          )}
          <ErrorBoundary>
            {renderScientificMaps()}
          </ErrorBoundary>
          <div style={{ marginTop: '1.5rem' }}>
            <ErrorBoundary>
              <ConvergenceChart datasetId={datasetId} />
            </ErrorBoundary>
          </div>
        </div>
      )}
    </div>
  );
}
