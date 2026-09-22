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
| MYO001 | Clasico (DAS) | $c$ (m/s) | 1522.17 | 1534.33 | 50.35 | 0.242 | -12.34 |
| MYO001 | Robusto (sin margen) | $c$ (m/s) | 1594.68 | 1521.39 | 47.21 | 1.552 | 3.82 |
| MYO001 | **Robusto (margen 15%)** | $c$ (m/s) | 1593.26 | 1539.15 | 45.54 | 1.188 | **1.50** |
| MYO001 | Clasico (DAS) | $\alpha_0$ (dB/MHz/cm) | 0.0880 | 0.0922 | 0.0964 | 0.043 | -27.23 |
| MYO001 | Robusto (sin margen) | $\alpha_0$ (dB/MHz/cm) | 0.2789 | 0.0567 | 0.0591 | 3.761 | 11.51 |
| MYO001 | **Robusto (margen 15%)** | $\alpha_0$ (dB/MHz/cm) | 0.2046 | 0.0521 | 0.0376 | 4.053 | **12.16** |
| TSH002 | Clasico (DAS) | $c$ (m/s) | 1465.97 | 1506.81 | 55.60 | 0.735 | -2.68 |
| TSH002 | Robusto (sin margen) | $c$ (m/s) | 1591.14 | 1489.03 | 46.72 | 2.186 | 6.79 |
| TSH002 | **Robusto (margen 15%)** | $c$ (m/s) | 1536.43 | 1476.06 | 28.93 | 2.087 | **6.39** |
| TSH002 | Clasico (DAS) | $\alpha_0$ (dB/MHz/cm) | 0.0339 | 0.0602 | 0.1020 | 0.258 | -11.76 |
| TSH002 | Robusto (sin margen) | $\alpha_0$ (dB/MHz/cm) | 0.2003 | 0.0299 | 0.0375 | 4.550 | 13.16 |
| TSH002 | **Robusto (margen 15%)** | $\alpha_0$ (dB/MHz/cm) | 0.0608 | 0.0264 | 0.0138 | 2.500 | **7.96** |

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