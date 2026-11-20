# Post hoc checks of the benchmark: results (validation sets, arrays GPL570)

Method: `POSTHOC_METHODS_191126.md` (SHA-256 dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f, self-recorded), 2000 permutations per matrix, seed 2026 + GSE number.

Gates: 118 matrices compared with the stored results; largest difference in the observed share 0.00e+00; 0 mismatches of the stored T_r2 P.

## Gene set A (20 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (20 matrices with a defined value): 0.960 / 0.976 / 0.997
- Expected without within-group coupling (mean of the null medians), quartiles: 0.985 / 0.991 / 0.998
- Excess (observed minus expected), quartiles: -0.026 / -0.007 / -0.000; observed above the 97.5th null percentile in 0 of 20
- Tr2 flagged (P < 0.05): original 17 of 20 (95% CI 0.64 to 0.95); after removal of PC1 10 of 20 (0.30 to 0.70)
- Paired: flagged both times 10, only before removal 7, only after removal 0, neither 3
- Share of the residual variance carried by PC1, quartiles: 0.415 / 0.480 / 0.564

## Gene set B (20 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (20 matrices with a defined value): 0.335 / 0.821 / 0.932
- Expected without within-group coupling (mean of the null medians), quartiles: 0.650 / 0.888 / 0.966
- Excess (observed minus expected), quartiles: -0.142 / -0.067 / -0.007; observed above the 97.5th null percentile in 0 of 20
- Tr2 flagged (P < 0.05): original 16 of 20 (95% CI 0.58 to 0.92); after removal of PC1 9 of 20 (0.26 to 0.66)
- Paired: flagged both times 9, only before removal 7, only after removal 0, neither 4
- Share of the residual variance carried by PC1, quartiles: 0.447 / 0.530 / 0.718

## Gene set A+B (40 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (40 matrices with a defined value): 0.826 / 0.958 / 0.995
- Expected without within-group coupling (mean of the null medians), quartiles: 0.895 / 0.982 / 0.997
- Excess (observed minus expected), quartiles: -0.088 / -0.024 / -0.002; observed above the 97.5th null percentile in 0 of 40
- Tr2 flagged (P < 0.05): original 33 of 40 (95% CI 0.68 to 0.91); after removal of PC1 19 of 40 (0.33 to 0.63)
- Paired: flagged both times 19, only before removal 14, only after removal 0, neither 7
- Share of the residual variance carried by PC1, quartiles: 0.429 / 0.513 / 0.613
