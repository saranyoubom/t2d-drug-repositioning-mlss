## effect_sizes.R - WP5a in R: effect sizes with 95% intervals, across-gene Benjamini-Hochberg q and the 26 August Gapdh sensitivity for the INS-1 RT-qPCR panel,
## plus ratios for DNA content, 2-NBDG and C-peptide per ug DNA where the data support them (R port of wp5_reporting/effect_sizes.py; nothing random is used).
## Inputs (read only): Submission_package_JTM/Additional file 3.xlsx and 4.xlsx, Analysis_outputs/reanalysis_17genes/gapdh_exact.json, SBA/Figures/Fig_4b_...xlsx (optional cross-check).
## Outputs: <out>/effect_sizes.csv, effect_sizes.md, gate_reproduction.json. Usage: Rscript effect_sizes.R   (ES_OUT = output folder, ES_LABEL = script name printed in effect_sizes.md, CMBL_ROOT = project root)
ES_HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
suppressMessages(library(jsonlite))
source(file.path(ES_HERE, "pycompat.R"))
for (p in 1:4) source(file.path(ES_HERE, sprintf("effect_sizes_part%d.R", p)))
es_main(Sys.getenv("ES_OUT", file.path(ES_HERE, "results_R", "wp5")))
