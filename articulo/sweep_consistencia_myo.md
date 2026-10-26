# Sweep de consistencia - Fantoma CIRS 040GSE (MYO001-MYO006)
_Medical Insight 3D - Phase 2_

Margen de borde excluido: 15%

## Resultados por dataset
| Dataset | Disponible | CNR c (dBr) | CNR a (dBr) | Pico cerca borde (c) | Pico cerca borde (a) |
|---|---|---|---|---|---|
| MYO001 | si | 9.82 | 9.91 | no | no |
| MYO002 | si | 9.40 | 9.38 | no | no |
| MYO003 | si | 9.21 | 7.66 | no | no |
| MYO004 | si | 9.88 | 9.66 | no | no |
| MYO005 | si | 9.46 | 9.91 | no | no |
| MYO006 | si | 8.04 | 8.70 | no | no |

## Estadisticas de consistencia
- CNR c:  media = 9.30 dBr, std = 0.61 dBr
- CNR a:  media = 9.20 dBr, std = 0.80 dBr
- Veredicto c: CONSISTENTE (std < 30% de la media)
- Veredicto a: CONSISTENTE (std < 30% de la media)

## Interpretacion
Si ambos parametros son consistentes: el contraste fisico es reproducible
across multiples adquisiciones del mismo fantoma, y el artefacto de borde
(documentado en runs previos) puede tratarse como limitacion conocida sin
bloquear la publicacion del resultado principal.

Si alguno es disperso: revisar si los archivos 'pico cerca borde' coinciden
con los outliers - eso indicaria que el artefacto de borde no es uniforme
y SI requiere PML antes de reportar cualquier CNR.