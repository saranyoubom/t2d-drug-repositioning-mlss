# Example (synthetic data)

N = 12 samples, 3 groups (sizes [4, 4, 4]), 17 genes, denominators: none; residual df = 8; 2000 permutations (run test 2000), seed 2026.

## Warnings

- group means carry 98% of the pooled covariance of significant raw pairs; without any within-group coupling 97% would be expected (95% null interval 93% to 100%), so a share near 1 reflects the group effects of the genes and is not by itself evidence about within-group coupling
- T_r2 = 0.175 is above the within-group permutation null (P = 0.0005): within-group coupling exists; it does not by itself identify biology
- consensus count 35 is within the permutation null (mean 31.5, P = 0.13): weak evidence of no coupling, because the count test is insensitive to injected coupling in public benchmark data; read T_r2

## Group-mean share

Median share of the pooled covariance carried by group means (raw, pooled-significant pairs, n = 35): 0.98 (mean 0.89).
Expected without within-group coupling (2000 permutations, mean of the null medians): 0.97 (2.5th to 97.5th percentile 0.93 to 1.00); observed minus expected +0.013.

## Within-group coupling (log-scale residuals, raw scheme)

| statistic | observed | null mean | null 95th pct | P (one-sided) |
|---|---|---|---|---|
| T_r2 | 0.1752 | 0.109 | 0.128 | 0.0005 |
| T_p | 21 | 6.245 | 10 | 0.0005 |
| T_eig | 6.684 | 4.243 | 4.994 | 0.0005 |

## Consensus count

Observed 35 (same sign 35) of 136 eligible pairs; permutation null mean 31.5 (SD 2.7), P(null >= obs) = 0.130.

## Geometry

r_crit 0.576; expected share if group means carry the correlation: 0.61; observed share of significant raw pairs: all 0.26, genes with a group effect 0.35 (n genes 13).
