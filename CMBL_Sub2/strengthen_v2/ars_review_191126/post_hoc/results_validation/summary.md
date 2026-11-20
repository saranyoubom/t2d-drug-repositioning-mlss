# Post hoc checks of the benchmark: results (validation sets (59), rebuilt matrices)

Method: `POSTHOC_METHODS_191126.md` (SHA-256 dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f, self-recorded), 2000 permutations per matrix, seed 2026 + GSE number.

Gates: 118 matrices compared with the stored results; largest difference in the observed share 0.00e+00; 0 mismatches of the stored T_r2 P.

## Gene set A (59 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (59 matrices with a defined value): 0.951 / 0.976 / 0.992
- Expected without within-group coupling (mean of the null medians), quartiles: 0.986 / 0.993 / 0.998
- Excess (observed minus expected), quartiles: -0.036 / -0.013 / -0.001; observed above the 97.5th null percentile in 0 of 59
- Tr2 flagged (P < 0.05): original 47 of 59 (95% CI 0.68 to 0.88); after removal of PC1 32 of 59 (0.42 to 0.66)
- Paired: flagged both times 32, only before removal 15, only after removal 0, neither 12
- Share of the residual variance carried by PC1, quartiles: 0.412 / 0.527 / 0.598

## Gene set B (59 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (59 matrices with a defined value): 0.405 / 0.845 / 0.952
- Expected without within-group coupling (mean of the null medians), quartiles: 0.714 / 0.920 / 0.984
- Excess (observed minus expected), quartiles: -0.147 / -0.073 / -0.007; observed above the 97.5th null percentile in 2 of 59
- Tr2 flagged (P < 0.05): original 48 of 59 (95% CI 0.70 to 0.89); after removal of PC1 32 of 59 (0.42 to 0.66)
- Paired: flagged both times 30, only before removal 18, only after removal 2, neither 9
- Share of the residual variance carried by PC1, quartiles: 0.430 / 0.532 / 0.727

## Gene set A+B (118 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (118 matrices with a defined value): 0.832 / 0.951 / 0.988
- Expected without within-group coupling (mean of the null medians), quartiles: 0.919 / 0.986 / 0.996
- Excess (observed minus expected), quartiles: -0.090 / -0.026 / -0.003; observed above the 97.5th null percentile in 2 of 118
- Tr2 flagged (P < 0.05): original 95 of 118 (95% CI 0.72 to 0.87); after removal of PC1 64 of 118 (0.45 to 0.63)
- Paired: flagged both times 62, only before removal 33, only after removal 2, neither 21
- Share of the residual variance carried by PC1, quartiles: 0.421 / 0.530 / 0.658
