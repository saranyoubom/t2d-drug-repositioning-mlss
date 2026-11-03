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

## Before you run anything

* Many scripts contain **absolute Windows paths and `setwd()` calls** from the authors' machine. Edit them first; some scripts write into their working directory and can overwrite files there.
* R 4.4.0 was used for the figures (patchwork >= 1.3.2 is needed with ggplot2 4.x). The Python scripts use numpy, pandas, scipy, statsmodels and openpyxl; the imports are at the top of each file.
* Raw CFX Maestro exports and the authors' calculation sheets are **not** included. Per-well efficiency and melt-curve data are in Additional file 2 of the manuscript.
* `drug.target.interaction.csv`, `FDA_Approved.csv`, `EMA_Approved.csv`, `PMDA_Approved.csv`, `FDA-EMA-PMDA_Approved.csv` and `string_interactions.csv` are derived from DrugCentral and STRING. They are included for reproducibility under the licenses of those sources; please cite the original databases.
* `reanalysis_17genes/insilico.py` queries public sequence services; replace the `YOUR_EMAIL_HERE` placeholder in its User-Agent string with your own address before running it.
* The numbers reported in the manuscript and its Additional files are authoritative; this folder documents how they were produced.
