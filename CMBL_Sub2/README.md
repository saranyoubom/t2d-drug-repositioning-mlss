# CMBL_Sub2: 17-gene analysis (November 2026)

Scripts, inputs and results behind the revised manuscript, which analyses **17 genes** (Hey1, Ryr1, Wnt2 and Wnt2b were measured but are excluded). The top-level scripts of this repository (npj SBA submission, 20 genes) and `BiologyDirect_Sub1/` are older versions and are left unchanged.

## Folders

| Folder | Content |
|---|---|
| `R_17genes/` | R scripts for Fig. 3f, 4b-4e, 5a-5d and 6b-6f with their inputs (relative-expression matrix `.xlsx`, drug-target and approval tables, STRING export) and CSV outputs. `data_export/` holds the Fig. 5 exports. Figures (PNG/SVG) are not included; the scripts regenerate them. |
| `R_panels_authors/` | R scripts for Fig. 1c, 1d, 2b-2f and 3b-3f (italic *p* labels, 7-gene Fig. 1d with key), with inputs and statistics CSV. |
| `reanalysis_17genes/` | Python re-implementation of the pooled co-expression pipeline (`pipeline.py`: four normalization schemes, group-preserving permutation null with seed 2026, treatment-adjusted partial correlations with Benjamini-Hochberg q), table counts, and the Gapdh-reference work for the 26 August 2020 DAPT group (`gapdh_exact.py`, `gapdh_qc_amplification.py`, `gapdh_qc_replicates.py`). `README_reanalysis.md` explains each result. |
| `primer_check/` | Primer evidence and the group-level efficiency sensitivity tables. |
| `additional_files/` | Scripts that build Additional files 2 and 4 from the outputs above. |
| `strengthen_v2/` | Version-2 analyses of November 2026: pre-specified protocol with SHA-256 fingerprints, calibrated within-group coupling checks, benchmark in public expression data, cross-check in public islet data, MLSS recomputation, pharmacovigilance annotation with openFDA, database-version sensitivity, plate provenance, database annotations, Gapdh reference sensitivity, recomputed Table S8, and, from release v2.1.0, the post hoc checks of the benchmark (`ars_review_191126/post_hoc/`) and the code of Fig. 10 (`fig10_201126/`). See its own `README.md`. |
| `benchmark_route3/` | First part of the benchmark (route 3) that the version-2 benchmark builds on, with its pre-specification and fingerprints. |

## Before you run anything

* Many scripts contain **absolute Windows paths and `setwd()` calls** from the authors' machine. Edit them first; some scripts write into their working directory and can overwrite files there.
* R 4.4.0 was used for the figures (patchwork >= 1.3.2 is needed with ggplot2 4.x). The Python scripts use numpy, pandas, scipy, statsmodels and openpyxl; the imports are at the top of each file.
* Raw CFX Maestro exports and the authors' calculation sheets are **not** included. Per-well efficiency and melt-curve data are in Additional file 2 of the manuscript.
* `drug.target.interaction.csv`, `FDA_Approved.csv`, `EMA_Approved.csv`, `PMDA_Approved.csv`, `FDA-EMA-PMDA_Approved.csv` and `string_interactions.csv` are derived from DrugCentral and STRING. They are included for reproducibility under the licenses of those sources; please cite the original databases.
* `reanalysis_17genes/insilico.py` queries public sequence services; replace the `YOUR_EMAIL_HERE` placeholder in its User-Agent string with your own address before running it.
* The scripts of `strengthen_v2/` and `benchmark_route3/` are the exception to the first point: they find their inputs relative to their own folder or through environment variables (see `strengthen_v2/README.md`).
* The Fig. 5 exports in `R_17genes/data_export/` and the script `Fig_5abcd_Gene_Interaction_Network_Analysis.r` are the version as first reported: gene names were matched to STRING literally, so the matrix gene `Glut2` (called `Slc2a2` in STRING) found no partners and 34 of 120 pairs were in STRING. The corrected classification, with the alias, is in `strengthen_v2/wp10_versions/fig5_corrected/alias/` (38 of 120 pairs in STRING; classes 5 / 16 / 33 / 66); the revised manuscript uses the corrected numbers.
* The numbers reported in the manuscript and its Additional files are authoritative; this folder documents how they were produced.
* Release v2.1.0 (November 2026) adds the post hoc checks of the benchmark and the code of Fig. 10 to `strengthen_v2/`, and the calibrated checks as a small Python package in `../tools/coexp-checks/`.
