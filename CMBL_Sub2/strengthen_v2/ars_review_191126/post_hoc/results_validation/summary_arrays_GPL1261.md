# Post hoc checks of the benchmark: results (validation sets, arrays GPL1261)

Method: `POSTHOC_METHODS_191126.md` (SHA-256 dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f, self-recorded), 2000 permutations per matrix, seed 2026 + GSE number.

Gates: 118 matrices compared with the stored results; largest difference in the observed share 0.00e+00; 0 mismatches of the stored T_r2 P.

## Gene set A (20 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (20 matrices with a defined value): 0.947 / 0.973 / 0.982
- Expected without within-group coupling (mean of the null medians), quartiles: 0.986 / 0.991 / 0.996
- Excess (observed minus expected), quartiles: -0.047 / -0.020 / -0.009; observed above the 97.5th null percentile in 0 of 20
- Tr2 flagged (P < 0.05): original 15 of 20 (95% CI 0.53 to 0.89); after removal of PC1 10 of 20 (0.30 to 0.70)
- Paired: flagged both times 10, only before removal 5, only after removal 0, neither 5
- Share of the residual variance carried by PC1, quartiles: 0.376 / 0.481 / 0.554

## Gene set B (20 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (20 matrices with a defined value): 0.637 / 0.896 / 0.946
- Expected without within-group coupling (mean of the null medians), quartiles: 0.822 / 0.928 / 0.966
- Excess (observed minus expected), quartiles: -0.105 / -0.047 / -0.006; observed above the 97.5th null percentile in 1 of 20
- Tr2 flagged (P < 0.05): original 15 of 20 (95% CI 0.53 to 0.89); after removal of PC1 11 of 20 (0.34 to 0.74)
- Paired: flagged both times 10, only before removal 5, only after removal 1, neither 4
- Share of the residual variance carried by PC1, quartiles: 0.380 / 0.514 / 0.661

## Gene set A+B (40 matrices)

- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (40 matrices with a defined value): 0.892 / 0.946 / 0.977
- Expected without within-group coupling (mean of the null medians), quartiles: 0.932 / 0.983 / 0.994
- Excess (observed minus expected), quartiles: -0.063 / -0.025 / -0.007; observed above the 97.5th null percentile in 1 of 40
- Tr2 flagged (P < 0.05): original 30 of 40 (95% CI 0.60 to 0.86); after removal of PC1 21 of 40 (0.37 to 0.67)
- Paired: flagged both times 20, only before removal 10, only after removal 1, neither 9
- Share of the residual variance carried by PC1, quartiles: 0.379 / 0.484 / 0.572
