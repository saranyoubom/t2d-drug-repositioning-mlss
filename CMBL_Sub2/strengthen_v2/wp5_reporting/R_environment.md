# R environment record (6 Nov 2026; R was not started and nothing was installed to write this)

Sources: the Methods and Declarations of `Manuscript/Manuscript_JTM_working.docx` (read only) and the `DESCRIPTION` files of the R libraries on this machine. The v2 analyses (WP1 to WP5) are in Python; R is used only by the authors' scripts (`../../R_17genes/`, `../../R_panels_authors/`) and by nothing in `run_all_v2.py`.

## What the manuscript says was used
| use | R version | packages stated in the Methods |
|---|---|---|
| amplification-efficiency re-estimation | not stated for this step | qpcR 1.4-1 (`sliwin()`) |
| networks, correlations and figures of the original analysis | 4.3.1 | igraph 2.1.4, ggraph 2.2.2, factoextra 1.0.7, ggcorrplot 0.1.4.1, ggpubr 0.6.0 |
| treatment-adjusted correlation, permutation null, group statistics, 17-gene reanalysis and the figures derived from it | 4.4.0 | not stated separately |
The Declarations (Availability of data and materials, software details) list: R (4.3.1 and 4.4.0); Python 3 (NumPy, pandas, SciPy, statsmodels) for cross-checks and the islet analyses; R packages qpcR, igraph, ggraph, ggplot2, factoextra, ggcorrplot, ggpubr.

## What is installed here (read from `DESCRIPTION`; the build field names the R version that compiled the package)
| package | R 4.3 user library | R 4.4 user library | R-4.3.1 system library | manuscript |
|---|---|---|---|---|
| igraph | 2.1.4 (built under R 4.3.3) | 2.3.0 (4.4.3) | - | 2.1.4 |
| ggraph | 2.2.2 (4.3.1) | 2.2.2 (4.4.3) | - | 2.2.2 |
| ggplot2 | 4.0.2 (4.3.1) | 4.0.3 (4.4.3) | 3.4.2 | not stated |
| factoextra | 1.0.7 (4.3.3) | 2.0.0 (4.4.3) | - | 1.0.7 |
| ggcorrplot | 0.1.4.1 (4.3.1) | 0.1.4.1 (4.4.3) | - | 0.1.4.1 |
| ggpubr | 0.6.3 (4.3.1) | 0.6.3 (4.4.3) | 0.6.0 | 0.6.0 |
| qpcR | not found in any library on this machine | not found | not found | 1.4-1 |
The project library `../../R_panels_authors/Rlib/` holds only patchwork. System R installations: R-4.2.2, R-4.3.1 and R-4.4.0 under `C:\Program Files\R`; R 4.4.0 is the interpreter named in the plan for the authors' scripts (`Rscript.exe`, `R_LIBS` pointing at `Rlib`).

## Points for the authors (not resolved here)
1. The 17-gene reanalysis is stated to have run under R 4.4.0, but the only 4.4 packages found here are igraph 2.3.0 and factoextra 2.0.0, which differ from the versions printed in the Methods (2.1.4 and 1.0.7, found in the 4.3 library). Say which versions produced the 17-gene outputs (`sessionInfo()` of the run) before the Methods or the repository README state versions.
2. qpcR 1.4-1 is cited for the efficiency re-estimation but is not installed in any library here; the scripts and outputs of that step (`../../primer_check/`) should be accompanied by the package version used.
3. A `sessionInfo()` printout from each of the two R versions should be saved next to the scripts when they are run on the authors' machine; none exists in the project.
Python packages: `requirements.txt` (exact versions read from the installed packages).
