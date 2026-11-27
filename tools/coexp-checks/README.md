# coexp-checks

Calibrated checks for pooled co-expression analyses of small designed experiments, for example 3 groups x 4 replicates and about 17 genes (RT-qPCR panels, small expression experiments). One Python module (NumPy, SciPy; pandas for the command line), with a command line and a Python API.

**Status: candidate checks.** Their calibration in public data is reported in the associated article. Whether using them improves the conclusions drawn from small experiments has not been tested, and passing them does not show that a network is biologically right: they detect some design artifacts, not all.

Version 1.1.1 (module `coexp_checks.py`; 1.1.0 plus a fix for genes with the same value in every sample, see `CHANGES.md`). The earlier module (version 1.0, 4 Nov 2026) is archived in the study repository; `CHANGES.md` lists what is different (version 1.1.0: a new function and two warnings, no existing statistic changed; version 1.1.1: a gene with the same value in every sample gets p = 1 in the group-effect test).

## Install and run
```
pip install .
coexp-checks examples/example_data.csv --groups 4,4,4 --nperm 2000 --out result.json --report report.md
```
or, without installing, `python src/coexp_checks.py DATA.csv --groups 4,4,4`. In Python:
```
import coexp_checks as cc
res = cc.run_checks(X, groups, genes=genes, den=None, nperm=2000)   # X: samples x genes, linear scale, strictly positive
print(cc.to_markdown(res))
```
`DATA.csv`: header = gene names, one row per sample (a non-numeric first column is read as sample ids), values on the linear scale (2^-dCq, CPM, intensity). `--groups 4,4,4` gives the sizes of consecutive groups, one label per sample, or `col:NAME` for a column of the file. Optional: `--denominators G1,G2,G3` (ratio schemes), `--runs` (run or plate label per gene), `--nperm`, `--seed`. The docstring at the top of `src/coexp_checks.py` documents every function.

`examples/make_example.py` writes a synthetic example (independent data with group effects and one injected within-group factor) and the report that the checks give for it; in that report the consensus count stays inside its null while the within-group coupling test flags the injected factor.

## What it reports
- **Group-mean share** (`group_mean_share`, `group_mean_share_expected`): the share of the pooled covariance of the significant pairs that is carried by the group means, and the share expected when within-group coupling is removed.
- **Permutation nulls** that keep each gene's group means (`consensus_null`, `adjusted_null`): the number of pairs significant under all normalization schemes and the treatment-adjusted counts, each against its null.
- **Within-group coupling** (`coupling_test`): the mean squared off-diagonal correlation of the within-group residuals (T_r2, primary), the number of residual pairs with P < 0.05 and the largest eigenvalue, against independent within-group residuals.
- **Run structure** (`run_structure_test`): whether within-group correlation is higher among genes measured in the same run or plate.
- **Geometry** (`geometry_check`), **simulation helpers** (`simulate`, `decouple`, `inject_factor`).

## How to read the results, and the limits
1. **A share of the pooled covariance near 1 is expected for genes selected for a group effect, even with no within-group coupling.** Read the excess over the expected share (`group_mean_share_expected`), not the share. In the study's benchmark (post hoc analyses) the median share for genes selected for a group effect was 0.976 against 0.995 expected without within-group coupling in the 37 development series, and 0.976 against 0.993 in the 59 validation sets; for random expressed genes it was 0.691 against 0.769 (development) and 0.845 against 0.920 (validation). The observed share was above the 97.5th percentile of its null in 2 of the 192 matrices.
2. **The count of pairs significant under all schemes is a weak check.** In the study's benchmark it detected injected within-group coupling far less often than T_r2 (28% against 91% at a within-group r of 0.7 in 6 genes), so a count inside its null does not show that coupling is absent.
3. **T_r2 did not meet the type I error requirement that was fixed in the study protocol.** In decoupled development data 276 of 3700 tests were significant (7.5%, requirement at most 6%); without the two series that contained genes with identical values in all samples it was 175 of 3500 (5.0%; a post hoc analysis). Genes without within-group variation get residual correlations that are numerical noise: `run_checks` warns, and they should be removed before T_r2 is read.
4. **A flagged coupling does not say where it comes from.** Within-group coupling was flagged in most public data sets of the benchmark; shared sample-level variation (batch, RNA quality, library size, normalization) and shared biology cannot be told apart by T_r2. After removal of the first principal component of the residuals (post hoc), T_r2 still flagged coupling in 32 of 59 validation sets for each gene set, against 47 and 48 before the removal, so about a third of the flagged sets lose the flag; the comparison does not say whether that component is technical or biological.
5. **Designs and genes.** Calibrated for one design (3 groups x 4 replicates) and 17 genes, on human and mouse arrays and human RNA-seq. Other designs, other gene numbers and other platforms are not calibrated. The permutation nulls assume that replicates are exchangeable within a group.
6. **Not shown:** that the checks improve conclusions, that they apply to toxicogenomic data, or that a pass means a network is right.

## Versions and tests
Tested only with Python 3.14.4, NumPy 2.4.4 and SciPy 1.17.1; `requires-python` is a declared minimum that was not tested. `tests/test_smoke.py` is a small self-contained suite of known-answer tests on simulated data (about a minute). The validation suite of the study (56 checks against an independent loop implementation and against the INS-1 reference numbers) is `test_coexp_checks.py` in the study repository and needs that repository's other files. The package was not built or installed in a clean environment.

## Benchmark kit
`benchmark_kit/` lists the 96 public designs of the benchmark (accessions, samples, genes of the two gene sets) so that others can rebuild the same matrices from GEO and test their own method on them. No expression data are included.

## Citation and licence
Please cite the associated article and this software (`CITATION.cff`). MIT licence (`LICENSE`).
