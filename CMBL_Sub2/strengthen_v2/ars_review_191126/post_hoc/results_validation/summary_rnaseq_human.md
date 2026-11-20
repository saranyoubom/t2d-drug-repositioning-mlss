# Post hoc checks of the benchmark: results (validation sets, rnaseq human)

Method: `POSTHOC_METHODS_191126.md` (SHA-256 dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f, self-recorded), 2000 permutations per matrix, seed 2026 + GSE number.

Gates: 118 matrices compared with the stored results; largest difference in the observed share 0.00e+00; 0 mismatches of the stored T_r2 P.

## Gene set A (19 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (19 matrices with a defined value): 0.949 / 0.982 / 0.990
- Expected without within-group coupling (mean of the null medians), quartiles: 0.986 / 0.995 / 0.999
- Excess (observed minus expected), quartiles: -0.046 / -0.013 / -0.001; observed above the 97.5th null percentile in 0 of 19
- Tr2 flagged (P < 0.05): original 15 of 19 (95% CI 0.57 to 0.91); after removal of PC1 12 of 19 (0.41 to 0.81)
- Paired: flagged both times 12, only before removal 3, only after removal 0, neither 4
- Share of the residual variance carried by PC1, quartiles: 0.504 / 0.613 / 0.713

## Gene set B (19 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (19 matrices with a defined value): 0.305 / 0.797 / 0.970
- Expected without within-group coupling (mean of the null medians), quartiles: 0.557 / 0.933 / 0.994
- Excess (observed minus expected), quartiles: -0.229 / -0.092 / -0.021; observed above the 97.5th null percentile in 1 of 19
- Tr2 flagged (P < 0.05): original 17 of 19 (95% CI 0.69 to 0.97); after removal of PC1 12 of 19 (0.41 to 0.81)
- Paired: flagged both times 11, only before removal 6, only after removal 1, neither 1
- Share of the residual variance carried by PC1, quartiles: 0.450 / 0.660 / 0.727

## Gene set A+B (38 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (38 matrices with a defined value): 0.799 / 0.951 / 0.989
- Expected without within-group coupling (mean of the null medians), quartiles: 0.926 / 0.989 / 0.998
- Excess (observed minus expected), quartiles: -0.100 / -0.040 / -0.002; observed above the 97.5th null percentile in 1 of 38
- Tr2 flagged (P < 0.05): original 32 of 38 (95% CI 0.70 to 0.93); after removal of PC1 24 of 38 (0.47 to 0.77)
- Paired: flagged both times 23, only before removal 9, only after removal 1, neither 5
- Share of the residual variance carried by PC1, quartiles: 0.460 / 0.613 / 0.739
