# R/Python Analysis Scripts — Revision 3

**Repository:** [saranyoubom/t2d-drug-repositioning-mlss](https://github.com/saranyoubom/t2d-drug-repositioning-mlss)
**Subfolder:** `R_scripts_Rev3/`
**Generated:** 2026-09-15

---

## Overview

This subfolder contains the analysis scripts, outputs, and regenerated figure assets added or changed in Revision 3 for:

> **Notch and Wnt inhibitions reveal beta cell dysregulation and predict synergistic diabetes drug candidates**
> *npj Systems Biology and Applications* — Revision 3 (Submission ID `c1ae7071-d445-422d-baf3-85229ddb6cc1`)

Revision 3 responded to a new third reviewer (Reviewer 3) in addition to Reviewer 2's outstanding comments from Revision 1. It does **not** supersede `R_scripts_Rev1/` — this folder contains only what is new or changed this round; unmodified figure scripts (Fig. 1–3, 5, 6c–e, S1) are unchanged from `R_scripts_Rev1/` and are not duplicated here. Revision 2 did not produce a separate deposited subfolder; changes described below as "Revision 3" are therefore changes since the Revision 1 deposit.

---

## Revision 3 Changes Summary

| Change | Scripts Affected | Reviewer Comment |
|--------|-------------------|-------------------|
| Treatment-adjusted correlation reanalysis (regress out treatment group, correlate residuals, BH-corrected across all consensus pairs) — new analysis, not in Rev1 | `Fig4_treatment_adjusted_reanalysis.py` | R3C2 |
| Normalization-gene sensitivity extended to Tcf7 (Wnt pathway) and Cacna1c (calcium pathway) as alternative normalizers from previously-unused functional categories | `R2C1_normalization_sensitivity.py` | R2C1 |
| **Data bug fix**: `Fig_4d_consensus_pairs.csv` incorrectly retained a `Wnt9a_Ryr1` row (Ryr1 is deliberately excluded from the 20-gene consensus panel everywhere else in the study, per Rev1's own exclusion). Corrected file has 27 data rows, not 28 — same filename as the Rev1 version, content corrected. | `Fig_4d_consensus_pairs.csv` (this folder) | Author-identified (R3C2 reconciliation) |
| Fig. 6b hub network regenerated from the corrected 27-pair consensus data; `print()` call removed so the script runs headless via `Rscript` (was previously interactive-device-only) | `Fig_6b_consensus_hub_network_revised.r` | Author-identified |
| Fig. 6f Sankey diagram: label/axis font sizes increased for legibility (stratum label 2.8pt, title 13pt, axis 10pt), and gene-symbol labels (Hub Genes axis only) italicized via `ggalluvial`'s exposed `stat_stratum` axis index | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Author-identified |
| LinRegPCR/qPCR amplicon-efficiency reanalysis attempted via the real published `sliwin()` algorithm (`qpcR` package, Ritz & Spiess 2008) on raw CFX well data, bypassing `qpcR`'s broken `rgl` install dependency by sourcing its R files directly. Included as evidence, not as a source of reported efficiency numbers — see script header and manuscript Methods/Response letter for why. | `R2C3_load_qpcR_functions.R`, `R2C3_test_sliwin_real.R` | R2C3 |

---

## Folder Contents

### Reanalysis Scripts (Python)

| File | Description |
|------|-------------|
| `Fig4_treatment_adjusted_reanalysis.py` | Regresses out treatment-group mean from each gene's expression, correlates residuals pairwise across the 20-gene consensus panel, applies Benjamini-Hochberg correction across all tests. Reports how many of the 27 pooled-consensus pairs survive treatment-adjustment. Output: `Fig4_treatment_adjusted_reanalysis_output.txt`, `SuppTable3_treatment_adjusted_data.csv` (Supplementary Table 3) |
| `R2C1_normalization_sensitivity.py` | Recomputes the four-normalization-method consensus-pair overlap using Tcf7 and Cacna1c as additional normalizer genes (Wnt and calcium-signalling categories, respectively — previously untested). Reports evaluable-pair retention (excludes structurally non-evaluable self-pairs from the denominator). Output: `R2C1_normalization_sensitivity_output.txt` |

### Evidence Scripts (R) — LinRegPCR / qPCR efficiency (R2C3)

| File | Description |
|------|-------------|
| `R2C3_load_qpcR_functions.R` | Sources the `qpcR` package's R functions directly from its installed `R/` source files and `sysdata.rda`, bypassing the package's normal `library(qpcR)` load path (blocked locally by a broken `rgl` Depends). Provides a working `pcrfit()`/`sliwin()` without a full package install. |
| `R2C3_test_sliwin_real.R` | Runs `sliwin()` (Ritz & Spiess 2008 window-of-linearity method) on real raw CFX well fluorescence data for a sample of genes/wells. Documents non-convergence / implausible fitted windows on CFX Maestro's baseline-subtracted export format — used in the manuscript response as evidence that the data format, not algorithm choice, is the blocker for a full per-amplicon efficiency reanalysis. Requires a local folder of raw CFX `.xlsx` exports (see script header; not included, provided separately per journal raw-data policy). |
| `R2C3_sliwin_test_output.txt` | Captured console output of `R2C3_test_sliwin_real.R` — 9 of 11 tested wells fail to converge (`SLIWIN FAILED: NA/NaN argument`); the 2 that do converge return implausible efficiencies (~0.1%). This is the evidence referenced above. |

### Corrected Data / Regenerated Figure Scripts (R)

| File | Figure Panel | Description |
|------|-------------|-------------|
| `Fig_4d_consensus_pairs.csv` | Fig. 4d, 6b | Corrected consensus gene-pair list — 27 rows (Ryr1-containing row removed; see Changes Summary above). Same filename as `R_scripts_Rev1/Fig_4d_Venn_diagram.r`'s output; content corrected. |
| `Fig_6b_consensus_hub_network_revised.r` | Fig. 6b | Regenerated from the corrected 27-pair data; hub-gene identities and degree values confirmed unchanged from the Rev1 (28-pair, buggy) version except for the removal of Ryr1 itself and one lost edge for Wnt9a. Outputs `Fig_6b_ConsensusHubNetwork_NPG.png/.svg`, `Fig_6b_Gene_Centrality_Summary.csv` |
| `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Fig. 6f | Six-layer mechanistic Sankey — font-size and gene-italics fixes only (see Changes Summary); underlying data unchanged from `R_scripts_Rev1/Fig_6f_regulatory_mechanistic_sankey.r`. Outputs `Fig_6f_Regulatory_Mechanistic_Sankey.png/.svg` |

---

## Required R / Python Packages

R scripts use the same dependencies as `R_scripts_Rev1/` (see that folder's README), plus:

```r
# R2C3 evidence scripts only — qpcR loaded via direct source(), not install.packages()
# Standard qpcR runtime deps (already required by a normal qpcR install):
install.packages(c("minpack.lm", "robustbase", "lattice", "MASS"))
```

Python scripts require:

```
pandas
numpy
scipy
statsmodels
```

---

## Execution Order

These scripts are independent of each other and of `R_scripts_Rev1/`'s execution order; each can be run standalone provided its input files are present.

```
1. Fig_4d_consensus_pairs.csv (this folder's corrected version) must exist before:
2. Fig_6b_consensus_hub_network_revised.r
3. Fig4_treatment_adjusted_reanalysis.py   (independent)
4. R2C1_normalization_sensitivity.py       (independent)
5. Fig_6f_regulatory_mechanistic_sankey_revised.r  (independent; data unchanged from Rev1)
6. R2C3_load_qpcR_functions.R, then R2C3_test_sliwin_real.R  (evidence only, not part of the figure pipeline)
```

---

## Input Data Files Required

All scripts here use relative paths and expect to be run with the working directory set to a folder that also contains the relevant `R_scripts_Rev1/` pipeline outputs below (this folder is an incremental overlay on Rev1, not a standalone pipeline — see Overview).

| File | Used by | Source |
|------|---------|--------|
| `Fig_4d_consensus_pairs.csv` (corrected, 27 rows) | `Fig_6b_consensus_hub_network_revised.r` | Included in this folder |
| `Fig_4e_Hes1_significant_correlations_only.csv` | `Fig_6b_consensus_hub_network_revised.r` | Output of `R_scripts_Rev1/Fig_4e_Hes1-normalized_correlation_matrix.R` |
| `Fig_6b_Gene_Centrality_Summary.csv` | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Included in this folder (output of `Fig_6b_consensus_hub_network_revised.r` above) |
| `Fig_6c_MLSS_Top50.csv` | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Output of the current MLSS scoring script (`Figures/Fig_6c_MLSS_v4.0.r` in the working tree; **note**: this postdates and is not identical to `R_scripts_Rev1/Fig_6c_MLSS_drug_combinations.r`'s `Fig_6c_mlss_top15.csv` output — the MLSS script itself was not part of this round's changes and has not been re-deposited; flagged here for awareness, not fixed in this round) |
| `FDA_Approved.csv`, `EMA_Approved.csv`, `PMDA_Approved.csv`, `FDA-EMA-PMDA_Approved.csv` | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Same files as `R_scripts_Rev1/README.md`'s Input Data Files table |
| `Fig_4c_Hes1_normalized_data.csv`, `Fig_4c_Kcnj11_normalized_data.csv`, `Fig_4c_Ins1_normalized_data.csv`, `Fig_4c_Raw_expression_data.csv` | `Fig4_treatment_adjusted_reanalysis.py`, `R2C1_normalization_sensitivity.py` (raw only) | Output of `R_scripts_Rev1/Fig_4c_normalization_methods_revised.r` |
| Raw CFX Maestro well-level fluorescence exports (baseline-subtracted, `.xlsx`) | `R2C3_test_sliwin_real.R` | Provided separately per journal raw-data policy, not included in this deposit |

---

## Citation

If you use these scripts, please cite:

> [Manuscript citation to be added upon acceptance]

---

*Revision 3 scripts compiled: 2026-09-15*
