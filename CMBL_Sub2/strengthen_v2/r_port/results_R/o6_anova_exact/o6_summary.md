# Impact check O6: route 3 ANOVA artefact (constant genes), exact arithmetic against numpy's order

Written 2026-11-26 21:06 by o6_anova_exact/compare_o6.R. Dry runs only (DE counts and the frozen selection loop).

## Route 3 / development set (frozen_datasets.csv, 20 analysable per platform)

Gate (R numpy order = archived full run, status and DE count of every series): PASS (58 archived series, 58 R series, 0 differing, 0 of them known by exact key)

| platform | analysed (numpy) | analysed (exact) | screened (numpy) | screened (exact) | DE count changed | lost analysability | gained (incl. newly screened) |
|---|---|---|---|---|---|---|---|
| GPL1261 | 20 | 20 | 33 | 34 | 3 | 1 | 1 |
| GPL570 | 17 | 16 | 25 | 25 | 2 | 1 | 0 |

Series whose DE count changes (5):
- GPL1261 GSE49079 (#8): DE 262 -> 0; analysable (dry run) -> not analysable: 0 DE genes (< 34)
- GPL1261 GSE15892 (#53): DE 2881 -> 1874; analysable (dry run) -> analysable (dry run)
- GPL1261 GSE28477 (#82): DE 8834 -> 8619; analysable (dry run) -> analysable (dry run)
- GPL570 GSE53183 (#117): DE 10262 -> 10206; analysable (dry run) -> analysable (dry run)
- GPL570 GSE20677 (#122): DE 1320 -> 0; analysable (dry run) -> not analysable: 0 DE genes (< 34)

Series whose role changes (analysed in one run only, or screened in one run only) (3):
- GPL1261 GSE49079 (#8): numpy 'analysable (dry run)' (DE 262); exact 'not analysable: 0 DE genes (< 34)' (DE 0)
- GPL1261 GSE30164 (#135): numpy 'not screened' (DE NA); exact 'analysable (dry run)' (DE 650)
- GPL570 GSE20677 (#122): numpy 'analysable (dry run)' (DE 1320); exact 'not analysable: 0 DE genes (< 34)' (DE 0)

Analysed set identical under both rules: FALSE

## WP2 held-out arrays (frozen_heldout_arrays.csv as archived, 20 analysable per platform)

Gate (R numpy order = archived full run, status and DE count of every series): PASS (59 archived series, 59 R series, 1 differing, 1 of them known by exact key)

Differing series:
- GSE2634: archived 'not analysable: ValueError samples missing: ['GSM50690', 'GSM50691', 'GSM50692']' DE NA; R 'not analysable: 0 DE genes (< 34)' DE 0 | KNOWN: bug of the archived Python reader (pandas skiprows; audit_261126 O1); not analysable in both, no number moves

| platform | analysed (numpy) | analysed (exact) | screened (numpy) | screened (exact) | DE count changed | lost analysability | gained (incl. newly screened) |
|---|---|---|---|---|---|---|---|
| GPL1261 | 20 | 20 | 24 | 24 | 0 | 0 | 0 |
| GPL570 | 20 | 20 | 35 | 35 | 1 | 0 | 0 |

Series whose DE count changes (1):
- GPL570 GSE68424 (#206): DE 10332 -> 10326; analysable (dry run) -> analysable (dry run)

Series whose role changes (analysed in one run only, or screened in one run only) (0):
- none

Analysed set identical under both rules: TRUE

## WP2 held-out RNA-seq

Cached count files read: 27; series with at least one gene constant in all 12 samples: 0 (largest number 0). The exact rule changes only such genes, so the RNA-seq DE counts are unchanged where this is 0.

