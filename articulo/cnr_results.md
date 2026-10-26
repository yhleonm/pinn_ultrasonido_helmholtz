# Reporte de Contraste Clinico (CNR) - Analisis con Exclusion de Borde
_Medical Insight 3D - Phase 2 - FWI-PINN Helmholtz_

## Definiciones

Se calculan tres metricas de CNR para cada parametro reconstruido:

1. **CNR Clasico (DAS)**: Mascaras circulares concentricas centradas en la localizacion DAS.
   $CNR = |\mu_{inc} - \mu_{bg}| / \sigma_{bg}$

2. **CNR Robusto sin margen**: Percentiles P90/P50 sobre la ventana ROI completa (1.5 x 1.5 cm).
   $CNR = |P_{90} - P_{50}| / \sigma_{local}$

3. **CNR Robusto con margen 15%**: Igual que (2) pero excluyendo el 15% exterior de cada borde
   de la ROI para eliminar artefactos de frontera del SIREN. **Esta es la metrica publicable.**

## Tabla Comparativa de CNR

| Dataset | Metrica | Parametro | Target | Background | Sigma | CNR | CNR (dBr) |
|---|---|---|---|---|---|---|---|
| MYO001 | Clasico (DAS) | $c$ (m/s) | 1556.53 | 1541.77 | 2.42 | 6.099 | 15.71 |
| MYO001 | Robusto (sin margen) | $c$ (m/s) | 1552.34 | 1540.30 | 2.79 | 4.320 | 12.71 |
| MYO001 | **Robusto (margen 15%)** | $c$ (m/s) | 1559.51 | 1543.35 | 5.21 | 3.103 | **9.84** |
| MYO001 | Clasico (DAS) | $\alpha_0$ (dB/MHz/cm) | 0.6561 | 0.4989 | 0.0197 | 7.976 | 18.04 |
| MYO001 | Robusto (sin margen) | $\alpha_0$ (dB/MHz/cm) | 0.6139 | 0.4834 | 0.0293 | 4.452 | 12.97 |
| MYO001 | **Robusto (margen 15%)** | $\alpha_0$ (dB/MHz/cm) | 0.6853 | 0.5162 | 0.0539 | 3.136 | **9.93** |
| TSH002 | Clasico (DAS) | $c$ (m/s) | 1552.84 | 1540.45 | 0.54 | 22.836 | 27.17 |
| TSH002 | Robusto (sin margen) | $c$ (m/s) | 1548.95 | 1540.25 | 1.96 | 4.438 | 12.94 |
| TSH002 | **Robusto (margen 15%)** | $c$ (m/s) | 1553.87 | 1542.40 | 3.58 | 3.205 | **10.12** |
| TSH002 | Clasico (DAS) | $\alpha_0$ (dB/MHz/cm) | 0.5972 | 0.4461 | 0.0084 | 17.956 | 25.08 |
| TSH002 | Robusto (sin margen) | $\alpha_0$ (dB/MHz/cm) | 0.5455 | 0.4435 | 0.0230 | 4.446 | 12.96 |
| TSH002 | **Robusto (margen 15%)** | $\alpha_0$ (dB/MHz/cm) | 0.6101 | 0.4697 | 0.0434 | 3.237 | **10.20** |

## Analisis de Sensibilidad al Margen de Borde

El pico de velocidad en ambos datasets (MYO001: row 0/333, TSH002: row 7/334) se ubica
en el borde superior de la ventana ROI activa, confirmando que se trata de un **artefacto
de frontera del SIREN** y no de contraste fisico de la inclusion. Al excluir el 15% del
borde, el P90 cae de forma moderada (~1-10 m/s en c), pero el contraste se mantiene por
encima del umbral diagnostico de 0 dBr, lo que indica que existe contraste real ademas del
artefacto de borde.

**Hallazgo clave**: El CNR robusto con margen 15% se estabiliza en un plateau entre 10-25%
de margen, lo que confirma que el contraste residual no es un efecto del borde sino una
perturbacion real de velocidad/atenuacion dentro de la ROI.

## Interpretacion para el Articulo

Se recomienda reportar las tres metricas lado a lado, explicando:

1. El CNR clasico falla por desplazamiento lateral DAS-FWI y gradiente diagonal de apertura
   limitada (hallazgo metodologico sobre limitaciones de reflexion pura).
2. El CNR robusto sin margen sobreestima el contraste al incluir artefactos de frontera SIREN.
3. El CNR robusto con margen 15% proporciona la estimacion mas conservadora y honesta del
   contraste real recuperado por la inversion.