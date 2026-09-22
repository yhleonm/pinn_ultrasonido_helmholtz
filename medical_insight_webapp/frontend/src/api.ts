// frontend/src/api.ts
// Cliente delgado sobre la API FastAPI. Sin estado, sin caché --
// eso vive en los componentes que lo consumen.

import type {
  DatasetSummary, ReconstructionMaps, ConvergenceHistory,
  NewJobResponse, JobStatusDetail,
} from './types';

const BASE_URL = import.meta.env.VITE_API_URL ?? 'http://localhost:8000';

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(`${BASE_URL}${path}`);
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status}: ${body}`);
  }
  return res.json();
}

export function fetchDatasets(): Promise<DatasetSummary[]> {
  return getJSON('/api/reconstructions/');
}

export function fetchMaps(datasetId: string, resolution = 128): Promise<ReconstructionMaps> {
  return getJSON(`/api/reconstructions/${datasetId}/maps?resolution=${resolution}`);
}

export function fetchConvergence(datasetId: string): Promise<ConvergenceHistory> {
  return getJSON(`/api/reconstructions/${datasetId}/convergence`);
}

export async function uploadNewJob(datasetId: string, file: File): Promise<NewJobResponse> {
  const form = new FormData();
  form.append('dataset_id', datasetId);
  form.append('file', file);

  const res = await fetch(`${BASE_URL}/api/jobs/`, { method: 'POST', body: form });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status}: ${body}`);
  }
  return res.json();
}

export function fetchJobStatus(jobId: string): Promise<JobStatusDetail> {
  return getJSON(`/api/jobs/${jobId}`);
}
