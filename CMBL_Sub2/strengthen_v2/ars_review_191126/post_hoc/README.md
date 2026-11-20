# Post hoc checks of the benchmark (19 to 20 November 2026)

Added in release v2.1.0. Two questions came out of a simulated review of the manuscript (written with the AI assistant described in the article, *Use of AI tools*; the review itself is not part of this repository):
1. Is "group means carry 98% (gene set A) and 85% (gene set B) of the pooled covariance" near its ceiling because the genes were selected for a group effect?
2. How much of the within-group coupling that T_r2 flags in public data survives the removal of one shared sample-level axis (the first principal component, PC1, of the residuals)?

**Status: post hoc.** Not part of the pre-specified protocol (`../../PROTOCOL_v2_prespecified.md`). The method was written before any result (`POSTHOC_METHODS_191126.md`, SHA-256 `dc88e898590bbf3c3b1a95ffd33c8ff6097d98d86387d972d24211ae52b18e2f`); the INS-1 matrix was added afterwards, after the development results were known (`POSTHOC_METHODS_ADDENDUM_191126.md`, SHA-256 `02c7a9593fa5b3a71ac4af86a840dc81dd26780ce9d6d8c019c6939215e99f3d`). The hashes are self-recorded, with no third-party time stamp (`FINGERPRINTS.txt`; `python -I verify_posthoc_fingerprints.py`). The article reports them as post hoc (Additional file 5, part 4.4; Additional file 6, sheets PostHoc_per_matrix and PostHoc_summary).

## Files
* `post_hoc_checks_191126.py` (the two checks; development mode on the stored matrices `../../wp2_benchmark2/compact/before_lean/dev_sets.npz`), `post_hoc_validation_191126.py` (the 59 validation sets, rebuilt from public GEO files with the benchmark's own loaders), `post_hoc_ins1_191126.py` (INS-1 matrix; the matrix is `../../../R_17genes/Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx`, or the environment variable `INS1_MATRIX`). The scripts import `coexp_checks.py` (version 1.0, unchanged) and `bench2_core.py` and write only into this folder.
* `results_dev/`, `results_ins1/`, `results_validation/`: `per_matrix.csv` (one row per matrix), `summary.json` and `summary.md`; the validation folder also has the summaries by arm.
* `FINGERPRINTS.txt`, `verify_posthoc_fingerprints.py`.
* `post_hoc_ins1_191126.py` differs from the file that produced the results in one line only: the default of the environment variable `INS1_MATRIX` was a path of the authors' computer and is now the path relative to the script.

## Gates (passed before any summary was written)
* Development series: the script reproduced the stored `share_between` of `dev_core_results.jsonl` (largest difference 0) and the stored `T_r2_p` (0 mismatches) for the 74 matrices.
* Validation sets: the rebuilt gene sets equalled the stored ones in all 118 matrices, the observed share equalled the stored share (largest difference 0) and the original T_r2 P equalled the stored value (0 mismatches). The 12 x 17 matrices were rebuilt from 61 public files (40 series matrices, 19 NCBI-generated count files, 2 platform annotation files; about 210 MB).
* INS-1: the script reproduced the share (0.971), T_r2 (0.133) and P (0.061) of the article.

## Results: 59 validation sets (40 held-out arrays and 19 human RNA-seq series; 2000 permutations per matrix, seed 2026 + GSE number)
| | Gene set A (17 genes drawn from those with a significant group effect) | Gene set B (17 random expressed genes) | Both |
|---|---|---|---|
| Observed share carried by group means, median (Q1 to Q3) | 0.976 (0.951 to 0.992) | 0.845 (0.405 to 0.952) | 0.951 (0.832 to 0.988) |
| Expected without within-group coupling, median of the null medians (Q1 to Q3) | 0.993 (0.986 to 0.998) | 0.920 (0.714 to 0.984) | 0.986 (0.919 to 0.996) |
| Observed minus expected, median (Q1 to Q3) | -0.013 (-0.036 to -0.001) | -0.073 (-0.147 to -0.007) | -0.026 (-0.090 to -0.003) |
| Matrices with the observed share above the 97.5th null percentile | 0 of 59 | 2 of 59 | 2 of 118 |
| T_r2 flagged (P < 0.05), original test | 47 of 59 (95% CI 0.68 to 0.88) | 48 of 59 (95% CI 0.70 to 0.89) | 95 of 118 (95% CI 0.72 to 0.87) |
| T_r2 flagged after removal of PC1 | 32 of 59 (0.42 to 0.66) | 32 of 59 (0.42 to 0.66) | 64 of 118 (0.45 to 0.63) |
| Paired: both / only before removal / only after removal / neither | 32 / 15 / 0 / 12 | 30 / 18 / 2 / 9 | 62 / 33 / 2 / 21 |
| Share of the residual variance in PC1, median (Q1 to Q3) | 0.527 (0.412 to 0.598) | 0.532 (0.430 to 0.727) | 0.530 (0.421 to 0.658) |

## Results: 37 development series (stored matrices; same settings)
| | Gene set A (17 genes drawn from those with a significant group effect) | Gene set B (17 random expressed genes) | Both |
|---|---|---|---|
| Observed share carried by group means, median (Q1 to Q3) | 0.976 (0.955 to 0.992) | 0.691 (0.424 to 0.864) | 0.938 (0.693 to 0.986) |
| Expected without within-group coupling, median of the null medians (Q1 to Q3) | 0.995 (0.991 to 0.998) | 0.769 (0.511 to 0.937) | 0.985 (0.772 to 0.996) |
| Observed minus expected, median (Q1 to Q3) | -0.015 (-0.031 to -0.006) | -0.057 (-0.106 to -0.003) | -0.027 (-0.070 to -0.005) |
| Matrices with the observed share above the 97.5th null percentile | 0 of 37 | 0 of 37 | 0 of 74 |
| T_r2 flagged (P < 0.05), original test | 25 of 37 (95% CI 0.51 to 0.80) | 29 of 37 (95% CI 0.63 to 0.89) | 54 of 74 (95% CI 0.62 to 0.82) |
| T_r2 flagged after removal of PC1 | 16 of 37 (0.29 to 0.59) | 18 of 37 (0.33 to 0.64) | 34 of 74 (0.35 to 0.57) |
| Paired: both / only before removal / only after removal / neither | 15 / 10 / 1 / 11 | 18 / 11 / 0 / 8 | 33 / 21 / 1 / 19 |
| Share of the residual variance in PC1, median (Q1 to Q3) | 0.496 (0.396 to 0.623) | 0.442 (0.392 to 0.544) | 0.463 (0.393 to 0.586) |

## Results: the INS-1 matrix (10000 permutations, seed 2026)
* Share carried by group means: observed 0.971 (66 pairs), expected without within-group coupling 0.983 (2.5th to 97.5th percentile 0.964 to 0.998); observed minus expected -0.011.
* T_r2: 0.1329, P = 0.0607 (original); after PC1 removal 0.1431, P = 0.1622; PC1 carries 0.356 of the residual variance.

## How to read it (fixed in the method note before the run)
1. The share of covariance carried by group means is near its ceiling for genes selected for a group effect even when within-group residuals are independent. The observed share was at or below its expectation in the median matrix and above the 97.5th percentile of its null in 2 of 192 matrices. "Group means carry 98%" is therefore a consequence of the selection of the genes and of the geometry of 3 groups of 4, not independent evidence that the design dominates a network.
2. Removing one shared axis lowered the number of validation sets flagged by T_r2 from 47 to 32 of 59 (gene set A) and from 48 to 32 of 59 (gene set B); about a third of the flagged sets lose the flag, and 54% of all validation sets are still flagged. The comparison does not say whether the removed axis is technical or biological (PC1 removal also removes shared biology).
3. The INS-1 matrix is typical in the first respect (97.1% against 98.3% expected). Its borderline coupling (P = 0.061) is not significant after PC1 removal (P = 0.16).

## Running it
* Development and INS-1 modes need only the files of this repository: `python -I post_hoc_checks_191126.py dev`, `python -I post_hoc_ins1_191126.py` (run in a copy: the scripts write `results_*` into their own folder).
* Re-running reproduces the tables (`summary.md`) and the `per_matrix.csv` files exactly; `summary.json` of the INS-1 run can differ from the stored file in the last binary digits of floating-point values (summation order; differences of about 1e-16), while the permutation P values are identical. A scratch copy of the repository layout gave this on 20 November 2026 (development mode 97 s, INS-1 run 3 minutes).
* The validation mode needs the 61 public files. `python -I post_hoc_validation_191126.py DATA_DIR --list` lists them without any network access; nothing is downloaded unless `--allow-download` is given. Downloaded files are data, never executed.
