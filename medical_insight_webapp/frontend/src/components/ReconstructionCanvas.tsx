import { useRef, useEffect } from 'react';

interface Props {
  bmodeData?: number[][] | null;
  data: number[][] | null;
  xMinMm: number;
  xMaxMm: number;
  zMinMm: number;
  zMaxMm: number;
  colorScale?: 'velocity' | 'attenuation';
  showOverlay?: boolean;
  overlayOpacity?: number;
  tumorXMm: number | null;
  tumorZMm: number | null;
  tumorRadiusMm: number;
  roiContourMm?: number[][] | null;
  predContourMm?: number[][] | null;
  showWarningRing: boolean;
  label: string;
  unit?: string;
  findingLabel?: string | null;
  isReliable: boolean;
  kSigma?: number;
}

// Colormaps idénticos a Matplotlib:
// - velocity: 'RdBu_r' divergente de 11 paradas (azul profundo -> celeste -> blanco -> coral -> rojo intenso)
// - attenuation: 'hot' secuencial térmico (negro -> rojo vivo -> naranja/amarillo brillante -> blanco)
const RDBU_R_STOPS: [number, number, number][] = [
  [5, 48, 97],     // 0.0
  [33, 102, 172],  // 0.1
  [67, 147, 195],  // 0.2
  [146, 197, 222], // 0.3
  [209, 229, 240], // 0.4
  [247, 247, 247], // 0.5 (blanco neutro en ~1564 m/s)
  [253, 219, 199], // 0.6
  [244, 165, 130], // 0.7
  [214, 96, 77],   // 0.8
  [178, 24, 43],   // 0.9
  [103, 0, 31],    // 1.0
];

function colorFor(value: number, min: number, max: number, scale: 'velocity' | 'attenuation'): string {
  const t = Math.max(0, Math.min(1, (value - min) / (max - min || 1)));

  if (scale === 'velocity') {
    const pos = t * (RDBU_R_STOPS.length - 1);
    const idx = Math.min(RDBU_R_STOPS.length - 2, Math.floor(pos));
    const frac = pos - idx;
    const c0 = RDBU_R_STOPS[idx];
    const c1 = RDBU_R_STOPS[idx + 1];
    const r = Math.round(c0[0] + frac * (c1[0] - c0[0]));
    const g = Math.round(c0[1] + frac * (c1[1] - c0[1]));
    const b = Math.round(c0[2] + frac * (c1[2] - c0[2]));
    return `rgb(${r}, ${g}, ${b})`;
  }

  // Colormap 'hot' exacto de Matplotlib: negro -> rojo vivo -> amarillo brillante -> blanco
  if (t < 0.375) {
    const s = t / 0.375;
    return `rgb(${Math.round(10 + s * 245)}, 0, 0)`;
  } else if (t < 0.75) {
    const s = (t - 0.375) / 0.375;
    return `rgb(255, ${Math.round(s * 255)}, 0)`;
  } else {
    const s = (t - 0.75) / 0.25;
    return `rgb(255, 255, ${Math.round(s * 255)})`;
  }
}

function getTicks(minVal: number, maxVal: number): number[] {
  const span = maxVal - minVal;
  let step = 5;
  if (span <= 16) step = 2;
  else if (span <= 35) step = 5;
  else step = 10;

  const ticks: number[] = [];
  const start = Math.ceil(minVal / step) * step;
  for (let val = start; val <= maxVal; val += step) {
    ticks.push(Number(val.toFixed(1)));
  }
  return ticks;
}

/**
 * Aplica el algoritmo de segmentación clínica: exclusión del 15% de margen,
 * umbral estadístico mean + k * sigma, y Connected Component Analysis (CCA)
 * para seleccionar el componente conexo de mayor tamaño.
 * 
 * GATE DE SEGURIDAD CLÍNICA: Si el canal es no confiable (!isReliable),
 * el cálculo se aborta inmediatamente sin ejecutar ninguna operación de máscara.
 */
function computeCCAMask(data: number[][], isReliable: boolean, k: number): boolean[][] {
  const nz = data.length;
  const nx = data[0].length;
  const mask = Array.from({ length: nz }, () => Array(nx).fill(false));

  // Gate de fiabilidad: abortar de inmediato si el canal no es confiable
  if (!isReliable) {
    return mask;
  }

  // 1. Exclusión de margen del 15%
  const mz = Math.floor(0.15 * nz);
  const mx = Math.floor(0.15 * nx);

  if (nz <= 2 * mz || nx <= 2 * mx) {
    return mask;
  }

  // 2. Extraer valores del interior para estadísticas de fondo (background)
  const interiorVals: number[] = [];
  for (let i = mz; i < nz - mz; i++) {
    for (let j = mx; j < nx - mx; j++) {
      interiorVals.push(data[i][j]);
    }
  }

  if (interiorVals.length === 0) {
    return mask;
  }

  // Encontrar el percentil 90
  const sortedVals = [...interiorVals].sort((a, b) => a - b);
  const p90Idx = Math.floor(0.9 * sortedVals.length);
  const p90Val = sortedVals[p90Idx];

  // Identificar valores de fondo (menores al percentil 90)
  const bgVals = interiorVals.filter(v => v < p90Val);
  if (bgVals.length === 0) {
    return mask;
  }

  // Calcular la media y desviación estándar del fondo
  const meanBg = bgVals.reduce((sum, v) => sum + v, 0) / bgVals.length;
  const varianceBg = bgVals.reduce((sum, v) => sum + Math.pow(v - meanBg, 2), 0) / bgVals.length;
  const sigmaBg = Math.sqrt(varianceBg);

  const threshold = meanBg + k * sigmaBg;

  // 3. Crear máscara binaria cruda en el interior
  const rawMask = Array.from({ length: nz }, () => Array(nx).fill(false));
  for (let i = mz; i < nz - mz; i++) {
    for (let j = mx; j < nx - mx; j++) {
      if (data[i][j] > threshold) {
        rawMask[i][j] = true;
      }
    }
  }

  // 4. Connected Component Analysis (CCA) - BFS para buscar el mayor
  const visited = Array.from({ length: nz }, () => Array(nx).fill(false));
  let largestComponent: { r: number; c: number }[] = [];

  for (let i = mz; i < nz - mz; i++) {
    for (let j = mx; j < nx - mx; j++) {
      if (rawMask[i][j] && !visited[i][j]) {
        const component: { r: number; c: number }[] = [];
        const queue: { r: number; c: number }[] = [{ r: i, c: j }];
        visited[i][j] = true;

        let head = 0;
        while (head < queue.length) {
          const curr = queue[head++];
          component.push(curr);

          // Vecinos 4-conectados en el interior
          const neighbors = [
            { r: curr.r - 1, c: curr.c },
            { r: curr.r + 1, c: curr.c },
            { r: curr.r, c: curr.c - 1 },
            { r: curr.r, c: curr.c + 1 }
          ];

          for (const neighbor of neighbors) {
            const { r, c } = neighbor;
            if (r >= mz && r < nz - mz && c >= mx && c < nx - mx) {
              if (rawMask[r][c] && !visited[r][c]) {
                visited[r][c] = true;
                queue.push({ r, c });
              }
            }
          }
        }

        if (component.length > largestComponent.length) {
          largestComponent = component;
        }
      }
    }
  }

  // 5. Rellenar la máscara final
  for (const pixel of largestComponent) {
    mask[pixel.r][pixel.c] = true;
  }

  return mask;
}

export default function ReconstructionCanvas({
  bmodeData = null,
  data,
  xMinMm,
  xMaxMm,
  zMinMm,
  zMaxMm,
  colorScale = 'velocity',
  showOverlay = false,
  overlayOpacity = 0.45,
  tumorXMm,
  tumorZMm,
  tumorRadiusMm,
  roiContourMm = null,
  predContourMm = null,
  showWarningRing,
  label,
  unit = '',
  findingLabel = null,
  isReliable,
  kSigma = 2.0,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const nz = bmodeData ? bmodeData.length : data?.length || 0;
    const nx = bmodeData ? bmodeData[0].length : data?.[0].length || 0;
    if (nz === 0 || nx === 0) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const cw = canvas.width;
    const ch = canvas.height;

    // Configuración de padding para los ejes y reglas (padRight incrementado para colorbar)
    const padLeft = 45;
    const padBottom = 35;
    const padTop = 15;
    const padRight = 50;

    const plotW = cw - padLeft - padRight;
    const plotH = ch - padTop - padBottom;

    // Fondo del canvas completo (área de ejes)
    ctx.fillStyle = '#121b19'; // coincide con --bg-panel
    ctx.fillRect(0, 0, cw, ch);

    // Fondo negro del área de ploteo
    ctx.fillStyle = '#000000';
    ctx.fillRect(padLeft, padTop, plotW, plotH);

    // 1. Renderizar B-mode si está presente (rango clínico [-60, 0] dB en OASBUD o ventana de 45 dB)
    let bmodeMin = -60, bmodeMax = 0;
    if (bmodeData) {
      let actualMin = Infinity, actualMax = -Infinity;
      for (const row of bmodeData) for (const v of row) {
        if (v < actualMin) actualMin = v;
        if (v > actualMax) actualMax = v;
      }
      bmodeMax = actualMax;
      const dynRange = roiContourMm && roiContourMm.length > 0 ? 60.0 : 45.0;
      bmodeMin = Math.max(actualMin, actualMax - dynRange);
    }

    const cellW = plotW / nx;
    const cellH = plotH / nz;

    if (bmodeData) {
      for (let i = 0; i < nz; i++) {
        for (let j = 0; j < nx; j++) {
          const val = bmodeData[i][j];
          const tLin = Math.max(0, Math.min(1, (val - bmodeMin) / (bmodeMax - bmodeMin || 1)));
          const t = roiContourMm && roiContourMm.length > 0 ? tLin : Math.pow(tLin, 0.9);
          const gray = Math.round(t * 255);
          ctx.fillStyle = `rgb(${gray}, ${gray}, ${gray})`;
          ctx.fillRect(padLeft + j * cellW, padTop + i * cellH, cellW + 0.5, cellH + 0.5);
        }
      }
    }

    // 2. Establecer límites de escala idénticos a las figuras científicas de Matplotlib:
    //    - Velocidad c(x,z): [1535, 1595] m/s (con RdBu_r centrado en ~1565 m/s)
    //    - Atenuación α0(x,z): [0.45, 1.10] dB/MHz/cm (con hot colormap de alto contraste)
    const min = colorScale === 'velocity' ? 1535.0 : 0.45;
    const max = colorScale === 'velocity' ? 1595.0 : (roiContourMm && roiContourMm.length > 0 ? 1.10 : 0.90);

    if (showOverlay && data && nz > 0 && nx > 0) {
      // Calcular la máscara CCA con el gate de seguridad
      const mask = computeCCAMask(data, isReliable, kSigma);

      ctx.save();
      ctx.globalAlpha = overlayOpacity;
      for (let i = 0; i < nz; i++) {
        for (let j = 0; j < nx; j++) {
          if (mask[i][j]) {
            ctx.fillStyle = colorFor(data[i][j], min, max, colorScale);
            ctx.fillRect(padLeft + j * cellW, padTop + i * cellH, cellW + 0.5, cellH + 0.5);
          }
        }
      }
      ctx.restore();
    } else if (!bmodeData && data && nz > 0 && nx > 0) {
      // Si no hay B-mode, dibujar el heatmap directamente (Vista Científica)
      for (let i = 0; i < nz; i++) {
        for (let j = 0; j < nx; j++) {
          ctx.fillStyle = colorFor(data[i][j], min, max, colorScale);
          ctx.fillRect(padLeft + j * cellW, padTop + i * cellH, cellW + 0.5, cellH + 0.5);
        }
      }
    }

    // 3. Dibujar la barra de color (Colorbar) vertical en el margen derecho
    if (data) {
      const cbLeft = cw - padRight + 8; // posición X de la barra
      const cbWidth = 10;
      const cbHeight = plotH;

      // Dibujar el degradado de color vertical
      for (let y = 0; y < cbHeight; y++) {
        const t = 1 - (y / cbHeight); // 1 arriba (max), 0 abajo (min)
        const val = min + t * (max - min);
        ctx.fillStyle = colorFor(val, min, max, colorScale);
        ctx.fillRect(cbLeft, padTop + y, cbWidth, 1);
      }

      // Dibujar el contorno de la barra
      ctx.strokeStyle = '#24332f'; // --line
      ctx.lineWidth = 1;
      ctx.strokeRect(cbLeft, padTop, cbWidth, cbHeight);

      // Dibujar etiquetas de la barra de color
      ctx.fillStyle = '#7fa89c'; // --phosphor-dim
      ctx.font = '8px "IBM Plex Mono", monospace';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';

      // Etiqueta Max (arriba)
      ctx.fillText(`${max}`, cbLeft + cbWidth + 5, padTop);

      // Etiqueta Min (abajo)
      ctx.fillText(`${min}`, cbLeft + cbWidth + 5, padTop + cbHeight);

      // Etiqueta Media (medio)
      const midVal = (min + max) / 2;
      ctx.fillText(`${midVal.toFixed(colorScale === 'velocity' ? 0 : 2)}`, cbLeft + cbWidth + 5, padTop + cbHeight / 2);
    }

    // Dibujar líneas de los ejes
    ctx.strokeStyle = '#7fa89c'; // coincide con --phosphor-dim
    ctx.lineWidth = 1;

    // Eje X (inferior)
    ctx.beginPath();
    ctx.moveTo(padLeft, padTop + plotH);
    ctx.lineTo(padLeft + plotW, padTop + plotH);
    ctx.stroke();

    // Eje Z (izquierdo)
    ctx.beginPath();
    ctx.moveTo(padLeft, padTop);
    ctx.lineTo(padLeft, padTop + plotH);
    ctx.stroke();

    // Dibujar divisiones y etiquetas
    ctx.fillStyle = '#7fa89c';
    ctx.font = '10px "IBM Plex Mono", monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';

    // Ticks X
    const xTicks = getTicks(xMinMm, xMaxMm);
    xTicks.forEach(val => {
      const px = padLeft + ((val - xMinMm) / (xMaxMm - xMinMm)) * plotW;
      if (px >= padLeft && px <= padLeft + plotW) {
        ctx.beginPath();
        ctx.moveTo(px, padTop + plotH);
        ctx.lineTo(px, padTop + plotH + 4);
        ctx.stroke();
        ctx.fillText(`${val}`, px, padTop + plotH + 6);
      }
    });
    ctx.fillText('x (mm)', padLeft + plotW / 2, padTop + plotH + 20);

    // Ticks Z
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const zTicks = getTicks(zMinMm, zMaxMm);
    zTicks.forEach(val => {
      const pz = padTop + ((val - zMinMm) / (zMaxMm - zMinMm)) * plotH;
      if (pz >= padTop && pz <= padTop + plotH) {
        ctx.beginPath();
        ctx.moveTo(padLeft, pz);
        ctx.lineTo(padLeft - 4, pz);
        ctx.stroke();
        ctx.fillText(`${val}`, padLeft - 6, pz);
      }
    });

    // Título Z (rotado)
    ctx.save();
    ctx.translate(12, padTop + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    ctx.fillText('z (mm)', 0, 0);
    ctx.restore();

    // Overlay del tumor: polígono irregular real del radiólogo (OASBUD) + contorno punteado verde PINN
    if (roiContourMm && roiContourMm.length > 2) {
      ctx.save();
      ctx.beginPath();
      ctx.rect(padLeft, padTop, plotW, plotH);
      ctx.clip();

      // 1. Contorno continuo Cian: Ground Truth del Radiólogo (roi_mask)
      ctx.strokeStyle = showWarningRing ? '#e8542a' : '#00e5ff';
      ctx.lineWidth = 2.0;
      ctx.setLineDash(showWarningRing ? [5, 4] : []);
      ctx.beginPath();
      roiContourMm.forEach(([xPt, zPt], idx) => {
        const px = padLeft + ((xPt - xMinMm) / (xMaxMm - xMinMm)) * plotW;
        const pz = padTop + ((zPt - zMinMm) / (zMaxMm - zMinMm)) * plotH;
        if (idx === 0) ctx.moveTo(px, pz);
        else ctx.lineTo(px, pz);
      });
      ctx.closePath();
      ctx.stroke();

      // 2. Contorno punteado Verde: Detección Ciega PINN (Componente Conexa Principal)
      if (predContourMm && predContourMm.length > 2) {
        ctx.strokeStyle = '#22c55e';
        ctx.lineWidth = 1.8;
        ctx.setLineDash([5, 4]);
        ctx.beginPath();
        predContourMm.forEach(([xPt, zPt], idx) => {
          const px = padLeft + ((xPt - xMinMm) / (xMaxMm - xMinMm)) * plotW;
          const pz = padTop + ((zPt - zMinMm) / (zMaxMm - zMinMm)) * plotH;
          if (idx === 0) ctx.moveTo(px, pz);
          else ctx.lineTo(px, pz);
        });
        ctx.closePath();
        ctx.stroke();
      }
      ctx.setLineDash([]);

      // Leyenda compacta en la esquina superior izquierda del área de ploteo
      const legX = padLeft + 6;
      const legY = padTop + 6;
      const legH = predContourMm && predContourMm.length > 2 ? 28 : 16;
      ctx.fillStyle = 'rgba(10, 16, 15, 0.72)';
      ctx.fillRect(legX, legY, 118, legH);
      ctx.font = 'bold 8px "IBM Plex Mono", sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';

      // Línea cian en leyenda
      ctx.strokeStyle = '#00e5ff';
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      ctx.moveTo(legX + 5, legY + 8);
      ctx.lineTo(legX + 19, legY + 8);
      ctx.stroke();
      ctx.fillStyle = '#e0f7fa';
      ctx.fillText('ROI Radiólogo', legX + 23, legY + 8);

      // Línea verde punteada en leyenda
      if (predContourMm && predContourMm.length > 2) {
        ctx.strokeStyle = '#22c55e';
        ctx.lineWidth = 1.8;
        ctx.setLineDash([3, 2]);
        ctx.beginPath();
        ctx.moveTo(legX + 5, legY + 20);
        ctx.lineTo(legX + 19, legY + 20);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = '#bbf7d0';
        ctx.fillText('Detección PINN', legX + 23, legY + 20);
      }

      if (findingLabel && tumorXMm !== null && tumorZMm !== null) {
        const px = padLeft + ((tumorXMm - xMinMm) / (xMaxMm - xMinMm)) * plotW;
        const pz = padTop + ((tumorZMm - zMinMm) / (zMaxMm - zMinMm)) * plotH;
        ctx.fillStyle = '#00e5ff';
        ctx.font = 'bold 9px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        ctx.fillText(findingLabel, px, Math.max(padTop + 14, pz - 12));
      }
      ctx.restore();
    } else if (tumorXMm !== null && tumorZMm !== null) {
      const px = padLeft + ((tumorXMm - xMinMm) / (xMaxMm - xMinMm)) * plotW;
      const pz = padTop + ((tumorZMm - zMinMm) / (zMaxMm - zMinMm)) * plotH;
      const radiusPx = (tumorRadiusMm / (xMaxMm - xMinMm)) * plotW;

      ctx.strokeStyle = showWarningRing ? '#e8542a' : '#2ad4a8';
      ctx.lineWidth = 2;
      ctx.setLineDash(showWarningRing ? [5, 4] : []);
      ctx.beginPath();
      ctx.arc(px, pz, radiusPx, 0, Math.PI * 2);
      ctx.stroke();
      ctx.setLineDash([]);

      // Dibujar etiqueta del hallazgo
      if (findingLabel) {
        ctx.fillStyle = showWarningRing ? '#e8542a' : '#2ad4a8';
        ctx.font = 'bold 9px sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'bottom';
        ctx.fillText(findingLabel, px + radiusPx + 4, pz - 2);
      }
    }
  }, [bmodeData, data, xMinMm, xMaxMm, zMinMm, zMaxMm, colorScale, showOverlay, overlayOpacity, tumorXMm, tumorZMm, tumorRadiusMm, roiContourMm, predContourMm, showWarningRing, findingLabel, isReliable, kSigma]);

  return (
    <figure className="recon-canvas-figure">
      <figcaption className="recon-canvas-caption">
        <span>{label}</span>
        {unit && <span className="recon-canvas-unit">{unit}</span>}
      </figcaption>
      <canvas ref={canvasRef} width={360} height={360} className="recon-canvas" />
    </figure>
  );
}
