# R/Python Analysis Scripts — Revision 3

**Repository:** [saranyoubom/t2d-drug-repositioning-mlss](https://github.com/saranyoubom/t2d-drug-repositioning-mlss) (`Rev3` branch)
**Generated:** 2026-09-15, corrected 2026-09-21

---

## Overview

This branch contains the analysis scripts, outputs, and regenerated figure assets added or changed in Revision 3 for:

> **Notch and Wnt inhibitions reveal beta cell dysregulation and prioritize candidate diabetes drug combinations**
> *npj Systems Biology and Applications* — Revision 3 (Submission ID `c1ae7071-d445-422d-baf3-85229ddb6cc1`)

Revision 3 responded to a new third reviewer (Reviewer 3) in addition to Reviewer 2's three outstanding comments from earlier rounds. The files listed below are only what is new or changed this round; unmodified figure scripts (Fig. 1-3, 5, 6c-e, S1) are unchanged from the main branch and are not duplicated here.

---

## Revision 3 Changes Summary

| Change | Scripts Affected | Reviewer Comment |
|--------|-------------------|-------------------|
| Treatment-adjusted correlation reanalysis (regress out treatment group, correlate residuals, BH-corrected across all consensus pairs) — new analysis | `Fig4_treatment_adjusted_reanalysis.R` | R3C2 |
| Normalization-gene sensitivity extended to Tcf7 (Wnt pathway) and Cacna1c (calcium pathway) as alternative normalizers from previously-unused functional categories | `R2C1_normalization_sensitivity.py` | R2C1 |
| **Data bug fix**: `Fig_4d_consensus_pairs.csv` incorrectly retained a `Wnt9a_Ryr1` row (Ryr1 is deliberately excluded from the 20-gene consensus panel everywhere else in the study). Corrected file has 27 data rows, not 28. | `Fig_4d_consensus_pairs.csv` | Author-identified, found while addressing R3C2's edge-count reconciliation request |
| Fig. 6b hub network regenerated from the corrected 27-pair consensus data; `print()` call removed so the script runs headless via `Rscript` | `Fig_6b_consensus_hub_network_revised.r` | Author-identified |
| Fig. 6f Sankey diagram: label/axis font sizes increased for legibility, and gene-symbol labels (Hub Genes axis only) italicized | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Author-identified |
| **LinRegPCR/qPCR amplification-efficiency re-derivation — succeeded.** Applying the published `sliwin()` algorithm (`qpcR` package, Ritz & Spiess 2008) to every well of the actual experimental runs (384 wells, 22 genes, 5 qPCR runs) converged for 379 wells (98.7%, mean R² = 0.992). Per-gene mean efficiency is genuinely low and non-uniform (52.9-74.6%), now reported in the manuscript's Supplementary Table 1 and Supplementary Data 3. Retrospective correction of the already-published fold-change values was not possible, because the well-to-replicate correspondence needed to apply it was lost after the original analysis (`R2C3_match_well_to_replicate.py` documents this: only 16 of 63 gene x treatment groups could be resolved). | `R2C3_extract_all_wells.py`, `R2C3_run_sliwin_all.R`, `R2C3_load_qpcR_functions.R`, `R2C3_extract_cq_all.py`, `R2C3_match_well_to_replicate.py` | R2C3 |

**Note on the 2026-09-21 correction**: files first pushed to this branch on 2026-09-15/17 mistakenly included an early, small-scale feasibility test (`R2C3_test_sliwin_real.R`, 11 wells, mostly non-convergent) that predated and was superseded by the comprehensive 384-well re-analysis above. That file has been removed and replaced with the actual pipeline that produced the manuscript's reported numbers. `Fig4_treatment_adjusted_reanalysis.py` has likewise been replaced with an R port (`Fig4_treatment_adjusted_reanalysis.R`, matching the filename cited in the manuscript's Supplementary Table 3); both versions produce numerically identical output.

---

## Folder Contents

### Reanalysis Scripts

| File | Description |
|------|-------------|
| `Fig4_treatment_adjusted_reanalysis.R` | Regresses out treatment-group mean from each gene's expression, correlates residuals pairwise across the 20-gene consensus panel, applies Benjamini-Hochberg correction across all tests. Reports how many of the 27 pooled-consensus pairs survive treatment-adjustment. Output: `Fig4_treatment_adjusted_reanalysis_output.txt`, `SuppTable3_treatment_adjusted_data.csv` (Supplementary Table 3) |
| `R2C1_normalization_sensitivity.py` | Recomputes the four-normalization-method consensus-pair overlap using Tcf7 and Cacna1c as additional normalizer genes (Wnt and calcium-signalling categories, respectively — previously untested). Reports evaluable-pair retention (excludes structurally non-evaluable self-pairs from the denominator). Output: `R2C1_normalization_sensitivity_output.txt` |

### Amplification-Efficiency Pipeline (R2C3)

| File | Description |
|------|-------------|
| `R2C3_extract_all_wells.py` | Extracts per-cycle fluorescence for every well of the actual RT-qPCR runs from the raw Bio-Rad CFX Maestro plate exports into one wide CSV (`all_wells_wide.csv`). Uses `sybr_well_gene_mapping.csv` (included) for well-to-gene identity on the three runs whose export format doesn't otherwise label it. |
| `R2C3_load_qpcR_functions.R` | Sources the `qpcR` package's R functions directly from a local qpcR source tree, bypassing the package's normal `library(qpcR)` load path where it's blocked by a broken `rgl` Depends. Only needed if `library(qpcR)` doesn't work in your environment; see the script header. |
| `R2C3_run_sliwin_all.R` | Runs `sliwin()` (Ritz & Spiess 2008 sliding-window log-linear regression) on every well in `all_wells_wide.csv`. Output: `sliwin_all_results.csv` (per-well efficiency, R², fit window), `R2C3_run_sliwin_all_output.txt` (per-gene convergence summary: 379/384 wells converged, 52.9-74.6% mean efficiency range). |
| `R2C3_per_gene_efficiency_summary.csv` | Per-gene mean efficiency, SD, range, and converged/total well counts, aggregated from `sliwin_all_results.csv` — the data underlying the manuscript's Supplementary Table 1 efficiency column and Supplementary Data 3. |
| `R2C3_extract_cq_all.py` | Extracts per-well Sample/Target/Cq records from the raw CFX "Quantification Cq Results" exports (`all_cq.csv`), for use by `R2C3_match_well_to_replicate.py`. |
| `R2C3_match_well_to_replicate.py` | Attempts to reconstruct which raw well corresponds to which technical-replicate label in the already-published fold-change table, by matching implied dCt values against the published values. Evidence for a disclosed manuscript limitation: only 16 of 63 gene x treatment groups could be resolved (`R2C3_match_well_to_replicate_output.txt`), confirming the well-to-replicate correspondence needed for retrospective efficiency correction cannot be reconstructed for most genes. |

### Corrected Data / Regenerated Figure Scripts (R)

| File | Figure Panel | Description |
|------|-------------|-------------|
| `Fig_4d_consensus_pairs.csv` | Fig. 4d, 6b | Corrected consensus gene-pair list — 27 rows (Ryr1-containing row removed; see Changes Summary above). |
| `Fig_6b_consensus_hub_network_revised.r` | Fig. 6b | Regenerated from the corrected 27-pair data; hub-gene identities and degree values confirmed unchanged except for the removal of Ryr1 itself and one lost edge for Wnt9a. Outputs `Fig_6b_ConsensusHubNetwork_NPG.png/.svg`, `Fig_6b_Gene_Centrality_Summary.csv` |
| `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Fig. 6f | Six-layer mechanistic Sankey — font-size and gene-italics fixes only; underlying data unchanged. Outputs `Fig_6f_Regulatory_Mechanistic_Sankey.png/.svg` |

---

## Required R / Python Packages

R scripts use the same dependencies as the main branch's figure scripts, plus:

```r
install.packages(c("openxlsx", "minpack.lm", "robustbase", "MASS"))
# qpcR itself if `library(qpcR)` works in your environment; otherwise see
# R2C3_load_qpcR_functions.R's header for the source()-based workaround.
```

Python scripts require:

```
pandas
numpy
scipy
statsmodels
openpyxl
```

---

## Execution Order

```
Treatment-adjustment / normalization-sensitivity reanalyses (independent of each other and of the R2C3 pipeline below):
1. Fig4_treatment_adjusted_reanalysis.R
2. R2C1_normalization_sensitivity.py

Consensus-pairs correction:
3. Fig_4d_consensus_pairs.csv (this folder's corrected version) must exist before:
4. Fig_6b_consensus_hub_network_revised.r
5. Fig_6f_regulatory_mechanistic_sankey_revised.r

Amplification-efficiency pipeline (each step's output feeds the next):
6. R2C3_extract_all_wells.py            -> all_wells_wide.csv
7. R2C3_load_qpcR_functions.R (if needed), then R2C3_run_sliwin_all.R
                                         -> sliwin_all_results.csv
8. (derive R2C3_per_gene_efficiency_summary.csv from sliwin_all_results.csv)
9. R2C3_extract_cq_all.py               -> all_cq.csv
10. R2C3_match_well_to_replicate.py     -> resolved_replicate_mapping.csv
```

---

## Input Data Files Required

| File | Used by | Source |
|------|---------|--------|
| `Fig_4d_consensus_pairs.csv` (corrected, 27 rows) | `Fig_6b_consensus_hub_network_revised.r` | Included in this deposit |
| `Fig_4e_Hes1_significant_correlations_only.csv` | `Fig_6b_consensus_hub_network_revised.r` | Output of the main branch's `Fig_4e_Hes1-normalized_correlation_matrix.R` |
| `Fig_6b_Gene_Centrality_Summary.csv` | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Included in this deposit (output of `Fig_6b_consensus_hub_network_revised.r` above) |
| `Fig_6c_MLSS_Top50.csv` | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Output of the current MLSS scoring script; not part of this round's changes and not re-deposited here |
| `FDA_Approved.csv`, `EMA_Approved.csv`, `PMDA_Approved.csv`, `FDA-EMA-PMDA_Approved.csv` | `Fig_6f_regulatory_mechanistic_sankey_revised.r` | Same files as the main branch's README Input Data Files table |
| `Fig_4c_Hes1_normalized_data.csv`, `Fig_4c_Kcnj11_normalized_data.csv`, `Fig_4c_Ins1_normalized_data.csv`, `Fig_4c_Raw_expression_data.csv` | `Fig4_treatment_adjusted_reanalysis.R`, `R2C1_normalization_sensitivity.py` (raw only) | Output of the main branch's `Fig_4c_normalization_methods_revised.r` |
| `sybr_well_gene_mapping.csv` | `R2C3_extract_all_wells.py` | Included in this deposit |
| Raw CFX Maestro well-level fluorescence exports (`.xlsx`) | `R2C3_extract_all_wells.py` | Raw instrument exports, not included per journal raw-data policy; provided separately on request |
| Raw CFX Maestro "Quantification Cq Results" exports (`.xlsx`) | `R2C3_extract_cq_all.py` | Same as above |
| Published per-replicate fold-change table (reproduces Supplementary Data 1) | `R2C3_match_well_to_replicate.py` | Not included; see manuscript Supplementary Data 1 |

---

## Citation

If you use these scripts, please cite:

> [Manuscript citation to be added upon acceptance]

---

*Revision 3 scripts compiled: 2026-09-15; corrected 2026-09-21 (see "Note on the 2026-09-21 correction" above).*
