# Sweep de consistencia - Fantoma CIRS 040GSE (MYO001-MYO006)
_Medical Insight 3D - Phase 2_

Margen de borde excluido: 15%

## Resultados por dataset
| Dataset | Disponible | CNR c (dBr) | CNR a (dBr) | Pico cerca borde (c) | Pico cerca borde (a) |
|---|---|---|---|---|---|
| MYO001 | si | 3.23 | 6.02 | si | si |
| MYO002 | si | 3.52 | 6.03 | si | si |
| MYO003 | si | 3.84 | 7.69 | si | si |
| MYO004 | si | 0.59 | 8.58 | si | si |
| MYO005 | si | 4.41 | 9.71 | si | si |
| MYO006 | si | 7.57 | 7.22 | si | si |

## Estadisticas de consistencia
- CNR c:  media = 3.86 dBr, std = 2.05 dBr
- CNR a:  media = 7.54 dBr, std = 1.32 dBr
- Veredicto c: DISPERSO - investigar causa
- Veredicto a: CONSISTENTE (std < 30% de la media)

## Interpretacion
Si ambos parametros son consistentes: el contraste fisico es reproducible
across multiples adquisiciones del mismo fantoma, y el artefacto de borde
(documentado en runs previos) puede tratarse como limitacion conocida sin
bloquear la publicacion del resultado principal.

Si alguno es disperso: revisar si los archivos 'pico cerca borde' coinciden
con los outliers - eso indicaria que el artefacto de borde no es uniforme
y SI requiere PML antes de reportar cualquier CNR.