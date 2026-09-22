// frontend/src/types.ts
// Espejo de los esquemas Pydantic en backend/app/models/schemas.py

export type JobStatus = 'queued' | 'running' | 'done' | 'failed';
export type TissueType = 'phantom' | 'in_vivo';

export interface DatasetSummary {
  dataset_id: string;
  tissue_type: TissueType;
  institution: string | null;
  status: JobStatus;
  tumor_x_mm: number | null;
  tumor_z_mm: number | null;
  cnr_c_dbr: number | null;
  cnr_alpha_dbr: number | null;
  is_reliable_velocity: boolean;
}

export type ReliabilityLevel = 'alta' | 'media' | 'baja';

export interface ClinicalInterpretation {
  finding_summary: string;
  reliability_level: ReliabilityLevel;
  reliability_reason: string;
  reliability_reason_velocity?: string;
  reliability_reason_attenuation?: string;
  velocity_reliable: boolean;
  attenuation_reliable: boolean;
}

export interface ReconstructionMaps {
  dataset_id: string;
  grid_shape: [number, number];
  x_min_mm: number;
  x_max_mm: number;
  z_min_mm: number;
  z_max_mm: number;
  c_map: number[][];
  alpha_map: number[][];
  bmode_map: number[][] | null;
  clinical_interpretation: ClinicalInterpretation | null;
  tumor_x_mm: number | null;
  tumor_z_mm: number | null;
  tumor_radius_mm: number;
  cnr_c_dbr: number | null;
  cnr_alpha_dbr: number | null;
  peak_distance_mm: number | null;
  peak_distance_alpha_mm: number | null;
  reliability_note: string | null;
}

export interface ConvergenceHistory {
  dataset_id: string;
  epochs: number[];
  loss: number[];
  err_c: number[];
  err_alpha: number[];
  stage_boundaries: number[];
}

export interface NewJobResponse {
  job_id: string;
  dataset_id: string;
  status: JobStatus;
  estimated_minutes: number;
}

export interface JobStatusDetail {
  job_id: string;
  dataset_id: string;
  status: JobStatus;
  queue_position: number;
  error: string | null;
  log_tail: string[];
}
