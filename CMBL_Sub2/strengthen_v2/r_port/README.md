# R port of the version-2 analyses (`r_port/`)

This folder holds an independent reimplementation, in R, of the version-2 analyses of this repository. Those analyses were first run in Python (`strengthen_v2/`, `benchmark_route3/`, `reanalysis_17genes/`). The port was written after all results were known, between 23 and 27 November 2026. It has two uses:

- It reproduces the archived Python results. The gates below check this.
- It computed the post hoc correction of the benchmark that the article reports. The scripts are in `o6_anova_exact/` and the results in `results_R/o6_anova_exact/`.

## Environment
- **R version.** R 4.3.1. The packages used include jsonlite, data.table, Matrix, openxlsx, readxl, curl and ragg. Their exact versions were not recorded; please keep the `sessionInfo()` of your own run.
- **Random numbers.** `rng/numpy_rng.R` reproduces numpy's generator (SeedSequence and PCG64) bit for bit, and `rng/pyrandom.R` reproduces Python's `random.Random`. Every permutation and simulation result therefore equals the Python one.
- **Number formatting.** `pycompat.R` writes numbers exactly as Python prints them.

## How to run
- `Rscript run_all_v2.R --list` lists the notes. Add `--dry-run` to check the scripts and inputs and print the steps without running them.
- `Rscript run_all_v2.R --source R` runs the pipeline from the data, in R:
  - downloads from GEO, EBI and NCBI, cached in `results_R/*_data/`;
  - the route-3 benchmark;
  - the WP2 data runs (part-1 gate, core analysis, held-out arrays, RNA-seq, calibration, design table);
  - the post hoc checks and the islet analyses;
  - every summary built from the R results.
  - It is estimated to take about 30 h on an 8 GB laptop with one worker; the full run has not been done (see below). The calibration and the design table can be run by hand with 2 workers (`calibrate_t2.R run 1 2` and `run 2 2`).
- `Rscript run_all_v2.R` (default `--source archive`) re-summarizes the archived Python results, so that `Rscript agreement_report.R` can compare the R results with the Python results.
- **Inputs that are not plain files of this repository:**
  - The saved openFDA responses of WP9 and WP12 are packed in `raw_responses.zip`. The first step, `unpack_openfda.R`, unpacks them into `results_R/openfda_raw/` and checks the WP9 folder digest against `RAW_SHA256.txt`.
  - The DrugCentral drug-target file of 25 Sep 2026 and the STRING 12.5 web export, which WP10 and WP12 read, are in `Database/` next to `CMBL_Sub2/`. They are third-party files under their own licenses (CC BY-SA 4.0 and CC BY 4.0); see `Database/README.md`. Without them the driver skips the steps that need them and names the missing file.
  - The Additional file 6 workbook of the article is not in the repository, so the AF6 steps are skipped in the same way.
  - A few steps read files that are not in this repository: the qPCR instrument exports (WP11), Additional files 3 and 4 of the article (WP5a), and the authors' original INS-1 workbook (optional cross-checks). The header of each script names these files.
  - The live queries (dataset selection, database retrievals) are not re-run: their results are frozen in the repository.
- **Tested from a clone of the repository (27 November 2026).**
  - The dry run of both modes passed.
  - WP9 ran: the responses were unpacked with the digest checked, then the analysis, the fallback, and 120 checks with 0 failures.
  - WP10 ran: STRING analysis, DrugCentral sensitivity, web-export check, and 54 checks with 0 failures.
  - WP12 ran. The checked WP9, WP10 and WP12 outputs are byte-identical to the archived ones.
  - The full run has not been done.

## Gates (results of 27 November 2026)
- **Gate G0** (`rng/test_numpy_rng.R`): pass, 0 failing checks.
- **Gate G0b** (`rng/test_pyrandom.R`): pass.
- **Number formatting** (`rng/test_pyformat.R`): 331 values × 5 formats, 0 differences.
- **Coexpression checks** (`test_coexp_checks.R`): 52 checks, 0 failed.
- **Route-3 validation on the INS-1 data** (`benchmark_route3/validate_pipeline.R`): 16 of 16.
- **Agreement report** (`results_R/agreement_report.md`, `agreement_lines.txt`).
  - 249 items agree and 25 differ, each difference with its reason in `agreement_known_reasons.tsv`; 0 missing, 0 errors. In addition, 27 files were compared byte for byte, and all are identical.
  - Reasons for the 25 differences (21 of them are floating-point rounding):
    - 10 items: in two development series (GSE20677, GSE28477) some genes have the same value in all 12 samples. Statistics computed on such genes are rounding noise in both languages. In addition, the archived Python count of GSE20677 set A includes one pair of two such genes (7 pairs; R gives the exact count, 6).
    - 6 items: rounding noise printed in check lines.
    - 3 items: a numpy result one unit in the last place from the correctly rounded value (the float32 logarithm in GSE84133, 2 items; a null draw that ties the observed value of an empirical P, 1 item).
    - 2 items: the WP1 comparison read the input copy in `R_17genes`, whose values are rounded to 16 digits; on the original matrix R gives the archived values.
    - 2 items: live GEO answers (a download that failed during the R selection run; HTTP status 403 instead of 404).
    - 2 items: an error of the archived Python reader (GSE2634, held-out arrays): the series is not analyzable in both, with a different reason.
- **Degenerate inputs** (`test_degenerate_inputs.R`, summary in `results_R/test_degenerate_inputs.txt`): 12 PASS, 2 KNOWN, 2 NOISE, 0 FAIL.

## Floating-point behavior reproduced on purpose, and the switches that give the exact rule
The archived Python analyses contain three floating-point artifacts. The port reproduces them by default, so that it gives the archived numbers. Each has a switch that applies the rule exactly instead.

| Where | Archived behavior (default) | Switch for the exact rule |
|---|---|---|
| Gene selection by analysis of variance (route 3 and WP2) | The grand mean is summed in numpy's order. A gene with the same value in all 12 samples can then get P = 0. | `ANOVA_EXACT=1`: such genes get P = 1 |
| Metacell detection (islet analysis) | The detected share is summed as scipy sums it. A gene detected in exactly half of the metacells falls below 50%. | `MC_EXACT=1`: count / n |
| Group-effect analysis of variance of the coexpression checks | As in the archived module. A gene with the same value in every sample gets P = 0, or every q becomes NaN. | `CC_CONSTANT_FIX=1`: the rule of coexp-checks 1.1.1 (`tools/coexp-checks`) |

## Post hoc correction of the benchmark (`o6_anova_exact/`)
- **When and how it was found.** While the benchmark was being reimplemented, after all results were known, the gene selection was found to depend on the floating-point artifact in the first row of the table above. This affected which series could be analyzed and which genes were drawn.
- **The pre-specified result.** It stays as run: type I error of T_r2 7.5% (276 of 3700), above the 6% limit, so it FAILS.
- **The corrected run.** With `ANOVA_EXACT=1`:
  - the development set has 36 series: GSE20677 and GSE49079 drop out, GSE30164 enters;
  - T_r2 type I error is 172 of 3600 = 4.8%, which passes;
  - power at r = 0.7 with 6 coupled genes is 656 of 720 = 91%.
  - The verdict file is `results_R/o6_anova_exact/wp2_exact/calibration_verdict.txt`, and the summary is `o6_calibration_summary.md`.
- **Run order.** `run_o6_route3.sh`, then `run_o6_heldout.sh`, then `run_o6_calib.sh`.
- **In the article.** See Additional file 5, part 4.1 and part 8.

## Data
- Downloaded data are cached in `results_R/wp2_data/` and `results_R/wp4_data/`. They are not in the repository and are read as numbers only.
- For a time, GEO answered HTTP 403 for GSE31851. When its download fails, the archived matrices (`oracle/dev_sets_archive.json`) stand in, and the log says so.
