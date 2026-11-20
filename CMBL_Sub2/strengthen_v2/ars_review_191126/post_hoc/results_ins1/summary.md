# Post hoc checks on the INS-1 matrix (addendum; 10000 permutations, seed 2026)

Methods: `POSTHOC_METHODS_191126.md` (SHA-256 dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f) and the addendum (SHA-256 02c7a9593fa5b3a71ac4af86a840dc81dd26780ce9d6d8c019c6939215e99f3d), both self-recorded.

Gates: observed share 0.971, T_r2 0.133 and P within 0.004 of 0.061 (manuscript): {'share_ok': True, 'tr2_ok': True, 'p_ok': True}.

## Check 1: share of the pooled covariance carried by group means
- Observed (median over 66 pooled-significant pairs): 0.971
- Expected without within-group coupling (mean of the null medians, 10000 permutations): 0.983 (2.5th to 97.5th percentile 0.964 to 0.998); observed minus expected -0.011; above the 97.5th percentile: False
- Pooled-significant pairs in the permuted matrices: mean 64.9 (observed 66); permutations without a significant pair: 0.000

## Check 2: T_r2 before and after removal of the first principal component
- Original: T_r2 0.1329, P = 0.0607
- After PC1 removal: T_r2 0.1431, P = 0.1622 (null mean 0.1303)
- Share of the residual variance carried by PC1: 0.356
