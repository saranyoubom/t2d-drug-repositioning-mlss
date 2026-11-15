# Wider search: concordance of kept perturbation series with the INS-1 arms (rules fixed in PRESPEC_WP4b... B4/B5 before any value was read)

| series | arm | expected sign | status | genes | same sign | exact binomial P | Spearman rho (P) | strong INS-1 genes same sign |
|---|---|---|---|---|---|---|---|---|
| E-MTAB-2034 | DKK-1 | +1 | not analysable (K4): none: the ArrayExpress record holds raw reads only (ENA); no processed expression table;  | - | - | - | - | - |
| GSE167073 | DKK-1 | +1 | analysed | 17 | 11 | 0.332 | +0.184 (0.480) | 4 of 6 |
| GSE193888 | DAPT | -1 (mirror) | not analysable (K4): none by plain download: no supplementary file list (HTTP 404) and no NCBI-generated count | - | - | - | - | - |
| GSE249170 | DAPT | -1 (mirror) | not analysable (K4): none by plain download: no supplementary file list (HTTP 404) and no NCBI-generated count | - | - | - | - | - |
| GSE60803 | DKK-1 | +1 | analysed | 17 | 7 | 0.629 | -0.037 (0.889) | 3 of 6 |
| GSE99503 | DKK-1 | -1 (mirror) | analysed | 17 | 8 | 1.000 | +0.037 (0.889) | 1 of 6 |

## Pooled (REML + Hartung-Knapp on Fisher z)

| set | k | pooled rho | 95% CI | I2 | outcome |
|---|---|---|---|---|---|
| PRIMARY: all analysable series | 3 | +0.062 | [-0.216, +0.331] | 0% | inconclusive |
| Wnt series (DKK-1 arm) | 3 | +0.062 | [-0.216, +0.331] | 0% | inconclusive |
| Notch series (DAPT arm) | 0 | - | - | - | rule not applied (fewer than 3 analysable series) |
| mirror designs excluded | 2 | - | - | - | rule not applied (fewer than 3 analysable series) |
| GSE60803 excluded | 2 | - | - | - | rule not applied (fewer than 3 analysable series) |
