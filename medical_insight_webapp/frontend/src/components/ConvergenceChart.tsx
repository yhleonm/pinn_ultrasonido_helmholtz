import { useState, useEffect } from 'react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip,
  ReferenceLine, Legend,
} from 'recharts';
import type { ConvergenceHistory } from '../types';
import { fetchConvergence } from '../api';

interface Props {
  datasetId: string;
}

export default function ConvergenceChart({ datasetId }: Props) {
  const [history, setHistory] = useState<ConvergenceHistory | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetchConvergence(datasetId).then(setHistory).catch((e) => setError(e.message));
  }, [datasetId]);

  if (error) return null; // no bloquear la vista principal si no hay historial
  if (!history) return <p className="muted">Cargando curva de convergencia...</p>;

  const rows = history.epochs.map((ep, i) => ({
    epoch: ep,
    errC: history.err_c[i],
    errAlpha: history.err_alpha[i],
  }));

  return (
    <div className="convergence-panel" style={{ width: '100%', overflowX: 'auto' }}>
      <h3>Convergencia por época</h3>
      <LineChart width={580} height={220} data={rows} margin={{ top: 8, right: 35, bottom: 4, left: 15 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="var(--grid-line)" />
        <XAxis dataKey="epoch" tick={{ fontSize: 11 }} label={{ value: 'Época', position: 'insideBottom', offset: -2, fontSize: 11 }} />
        <YAxis
          yAxisId={0}
          tick={{ fontSize: 11 }}
          label={{ value: 'Error velocidad (%)', angle: -90, position: 'insideLeft', fontSize: 11, offset: 0 }}
          stroke="#3f7cd6"
        />
        <YAxis
          yAxisId={1}
          orientation="right"
          tick={{ fontSize: 11 }}
          label={{ value: 'Error atenuación (%)', angle: 90, position: 'insideRight', fontSize: 11, offset: 0 }}
          stroke="#c98a1f"
        />
        <Tooltip
          formatter={(value: any, name: string) => [`${(value || 0).toFixed(2)}%`, name === 'errC' ? 'Error velocidad' : 'Error atenuación']}
          labelFormatter={(ep) => `Época ${ep}`}
        />
        <Legend
          formatter={(value) => (value === 'errC' ? 'Error velocidad (c)' : 'Error atenuación (α)')}
        />
        {history.stage_boundaries.map((ep) => (
          <ReferenceLine key={ep} x={ep} yAxisId={0} stroke="var(--stage-line)" strokeDasharray="4 4" />
        ))}
        <Line yAxisId={0} type="monotone" dataKey="errC" stroke="#3f7cd6" strokeWidth={2} dot={false} />
        <Line yAxisId={1} type="monotone" dataKey="errAlpha" stroke="#c98a1f" strokeWidth={2} dot={false} />
      </LineChart>
    </div>
  );
}
