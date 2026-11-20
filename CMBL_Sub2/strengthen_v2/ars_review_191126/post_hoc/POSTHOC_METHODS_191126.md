# Post hoc checks of the benchmark: method fixed before any result (19 Nov 2026)

**Status: post hoc.** These two checks are NOT part of the pre-specified protocol (Additional file 5). They answer two questions raised by the simulated review of 19 Nov 2026 (`../ARS_review_JTM_manuscript_191126.md`, reviewer 1 weakness W1 and W2). This file was written and its SHA-256 was recorded by `post_hoc_checks_191126.py` before the script was run on any matrix; the hash is self-recorded (no third-party time stamp), so it shows that the file did not change after the run, not when it was made.

## Data
- Development series: the 37 series of the benchmark, the stored 12 x 17 matrices of gene set A and gene set B (`../../wp2_benchmark2/compact/before_lean/dev_sets.npz`), linear scale, 3 groups x 4 replicates in file order.
- Validation series (59: 40 held-out arrays, 19 RNA-seq): the same two sets rebuilt from the public GEO files by the benchmark's own loaders (`bench_data_io.py`, `run_rnaseq.py`). Only run if the user permits the download; the report says which sets were used.

## Code
`coexp_checks.py` (version 1.0) and the design of `bench2_core.py` (`Design([0]*4 + [1]*4 + [2]*4)`), read only. Permutations: 2000 per matrix (as `nperm_cp` and `nperm_cons` of `analyse_matrix`), seed = 2026 + GSE number (as `analyse_matrix`), alpha = 0.05. No new random seeds are tried.

## Check 1: the share of pooled covariance carried by group means, against its value without within-group coupling
- Observed: `group_mean_share(X, design)` (raw, Pearson on the linear scale), the median over the pairs with pooled |r| >= r_crit(N - 2) of cov(group means) / cov(all samples). This is the number that the manuscript reports (97% for INS-1; 98% and 85% for gene sets A and B).
- Expected without within-group coupling: the same statistic on each of 2000 matrices in which every gene's values are permuted independently within its group (group means kept, within-group coupling removed; the null of `consensus_null` and `coupling_test`). The set of pooled-significant pairs is recomputed in every permuted matrix, so the selection of significant pairs is part of the null.
- Reported per matrix: observed median share, mean of the null medians, 2.5th and 97.5th percentiles, excess = observed minus null mean, number of significant pairs (observed, null mean), share of permutations without a significant pair (their median is undefined and they are left out of the null summary).
- Across series (set A, set B, both): median, quartiles of the observed share, of the null mean and of the excess; the number of series whose observed share is above the 97.5th null percentile.

## Check 2: the within-group coupling test after removing one shared sample-level axis
- Residuals: log-scale values with group means removed (`design.Q @ log(X)`, as in `coupling_test`), a genes x samples matrix.
- PC1 removal: the first singular component of the residual matrix (largest singular value) is subtracted; T_r2' = mean squared off-diagonal correlation of the remaining rows. Exactly one component is removed; no other number of components is run.
- Null: the same permuted residual matrices as the original test (same seed, same permutations); PC1 is removed from every permuted matrix before T_r2' is computed, so the null is matched to the statistic. P' = (1 + #{null >= observed}) / (2000 + 1); a series is flagged when P' < 0.05.
- Reported per matrix: T_r2 and its P (original), T_r2' and P' (after removal), share of the residual variance carried by PC1.
- Across series (set A, set B): number flagged by the original test and after removal with Wilson 95% intervals, and the paired table (flagged both times, only before, only after, neither).

## Gates before any result is read
1. For the development matrices the script must reproduce the stored `share_between` of `dev_core_results.jsonl` (maximum absolute difference <= 1e-9) and the stored `T_r2_p` (exact equality), because it uses the same code and seeds. If a gate fails, nothing else is reported and the cause is investigated.
2. The type I error of the PC1-removed test is not re-derived: it is calibrated by construction (the same removal is applied to every permuted matrix).

## Rules for reading the results (fixed now)
- Everything listed above is reported, whichever direction it goes. No variant (two components, another threshold, another statistic) is added after seeing a result as part of this check; an extra analysis, if one is made, is labeled exploratory in the report.
- Check 1 is read through the excess (observed minus expected without coupling). If the median excess over series is small, the statement "group means carry 98%" is reworded as a property of the selected genes, not as evidence that the design dominates.
- Check 2 is read as a paired comparison. The text states both fractions side by side and does not name the removed axis "technical": PC1 removal also removes shared biology, so the comparison cannot separate the two.
- The manuscript, the abstract and Fig. 10 take their numbers from the report of this folder only.
