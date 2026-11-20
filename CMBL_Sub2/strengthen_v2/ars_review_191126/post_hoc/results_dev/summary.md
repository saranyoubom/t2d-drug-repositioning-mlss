# Post hoc checks of the benchmark: results (development series (37), stored matrices)

Method: `POSTHOC_METHODS_191126.md` (SHA-256 dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f, self-recorded), 2000 permutations per matrix, seed 2026 + GSE number.

Gates: 74 matrices compared with the stored results; largest difference in the observed share 0.00e+00; 0 mismatches of the stored T_r2 P.

## Gene set A (37 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (37 matrices with a defined value): 0.955 / 0.976 / 0.992
- Expected without within-group coupling (mean of the null medians), quartiles: 0.991 / 0.995 / 0.998
- Excess (observed minus expected), quartiles: -0.031 / -0.015 / -0.006; observed above the 97.5th null percentile in 0 of 37
- Tr2 flagged (P < 0.05): original 25 of 37 (95% CI 0.51 to 0.80); after removal of PC1 16 of 37 (0.29 to 0.59)
- Paired: flagged both times 15, only before removal 10, only after removal 1, neither 11
- Share of the residual variance carried by PC1, quartiles: 0.396 / 0.496 / 0.623

## Gene set B (37 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (37 matrices with a defined value): 0.424 / 0.691 / 0.864
- Expected without within-group coupling (mean of the null medians), quartiles: 0.511 / 0.769 / 0.937
- Excess (observed minus expected), quartiles: -0.106 / -0.057 / -0.003; observed above the 97.5th null percentile in 0 of 37
- Tr2 flagged (P < 0.05): original 29 of 37 (95% CI 0.63 to 0.89); after removal of PC1 18 of 37 (0.33 to 0.64)
- Paired: flagged both times 18, only before removal 11, only after removal 0, neither 8
- Share of the residual variance carried by PC1, quartiles: 0.392 / 0.442 / 0.544

## Gene set A+B (74 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (74 matrices with a defined value): 0.693 / 0.938 / 0.986
- Expected without within-group coupling (mean of the null medians), quartiles: 0.772 / 0.985 / 0.996
- Excess (observed minus expected), quartiles: -0.070 / -0.027 / -0.005; observed above the 97.5th null percentile in 0 of 74
- Tr2 flagged (P < 0.05): original 54 of 74 (95% CI 0.62 to 0.82); after removal of PC1 34 of 74 (0.35 to 0.57)
- Paired: flagged both times 33, only before removal 21, only after removal 1, neither 19
- Share of the residual variance carried by PC1, quartiles: 0.393 / 0.463 / 0.586
