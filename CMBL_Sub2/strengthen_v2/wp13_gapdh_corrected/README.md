# WP13: the corrected version of the 26 August Gapdh reference (9 Nov 2026)

> Release copy: working document of the analysis. Sections that only listed points for the authors or proposed wording for the manuscript are omitted and marked; the reading key is in `../README.md`.

**Strength of this evidence.** A descriptive recomputation, post hoc, no new rule: the three versions were defined on 2 Nov 2026 in `reanalysis_17genes/gapdh_exact.py` and are used unchanged. The reported analysis stays primary until the authors decide (section 4). Request of the user (9 Nov 2026): the reason why Gapdh wells A10 and B10 were replaced is not known; "can you try to calculate the corrected version". Independent check: `check_wp13_numbers.py` (`check_wp13_last_run.txt`).

## 1. The answer in plain words
* **What was done in 2020.** In the DAPT group of the 26 August plate, Gapdh wells A10 (Cq 24.99) and B10 (24.13) are 1.3 to 2.4 cycles later than C10 (22.57) and D10 (22.80). The calculation sheet does not use A10 and B10; it enters the values of C10 and D10 twice (22.5685, 22.7983, 22.5685, 22.7983), and the four DAPT replicates of each of the five genes of that run are divided by these four numbers. The reason is not recorded anywhere in the project.
* **The data decide the question better than memory can.** The plateau fluorescence of A10 and B10 is 513 and 725 RFU, which is 38% and 54% of the median of the other ten Gapdh wells (all between 1279 to 1507 RFU). A smaller amount of template delays the Cq but does not lower the plateau (the control wells, one cycle later than the treated ones, have a normal plateau), so these two reactions failed or were inhibited. The Gapdh replicate SD in the DAPT group is 1.14 cycles with all four wells, against 0.10 (CTRL) and 0.15 (DKK-1). Both melt curves still show one peak at 82.5 degrees C (the peaks are lower: 169 and 245 against 453 to 544; earlier QC `reanalysis_17genes/README_reanalysis.md`, not re-derived here). **Excluding A10 and B10 is therefore well supported; duplicating C10 and D10 to fill the empty rows is the part that is not correct.**
* **The corrected version is "C0": exclude A10 and B10 without duplicating anything**, i.e. every DAPT replicate of the five genes is divided by the mean of the two valid wells (22.6834). "C1" (all four wells, mean 23.6215) is the worst case in which the failed wells are kept.

## 2. What changes (`results/gapdh_versions_summary.csv`, `gapdh_versions.json`)
Gate: the version "reported" reproduces the reported numbers (21 consensus pairs; pooled 66 / 80 / 84 / 120; adjusted 0 / 51 / 11 / 15; hubs Glp1r 6 and Tcf7l2, Cacna1d, Ptbp1 5; Fig. 5 with the Glut2 alias 5 / 16 / 33 / 66).
| | reported (duplicated C10, D10) | C0 (corrected: A10, B10 excluded) | C1 (all four wells) |
|---|---|---|---|
| DAPT reference, Gapdh Cq | 22.5685, 22.7983, 22.5685, 22.7983 | 22.6834 for all four | 23.6215 for all four |
| DAPT fold change of Tcf7l2, Kcnj11, Cacna1d, Cacna1c, Ptbp1 | 0.78, 0.40, 0.76, 0.71, 0.62 | 0.78, 0.40, 0.76, 0.72 and 0.62 (no change) | 1.49, 0.76, 1.46, 1.37 and 1.18 (DIRECTION REVERSES for four genes) |
| Consensus pairs of the 4 normalization schemes | 21 | 24 consensus pairs; 20 of the 21 kept, lost Itpr1-Tcf7l2; new Itpr1-Ryr3, Itpr1-Wnt5a, Ryr2-Tcf7l2, Ryr2-Wnt5b | 21, but only 11 of the 21 kept (10 lost, 10 new) |
| Hubs (degree 5 or more) | Glp1r 6; Tcf7l2, Cacna1d, Ptbp1 5 | the same four plus Wnt5a 5 | Glp1r 5 and Cacna1d 5 only |
| Pooled significant pairs, raw / Ins1 / Kcnj11 / Hes1 | 66 / 80 / 84 / 120 | 63 / 80 / 86 / 120 | 62 / 49 / 44 / 120 |
| Adjusted q < 0.05, same order | 0/51/11/15 | 0/47/36/12 | 2/41/12/12 |
| r(Cacna1d, Glp1r), Hes1-normalized | 0.960 | 0.976 | 0.974 |
| Fig. 5 with the alias: Known / Novel / STRING Only / Not analytically consistent | 5 / 16 / 33 / 66 (38 in STRING) | 6 / 18 / 32 / 64 (38 in STRING) | 5 / 16 / 33 / 66 (38 in STRING; different pairs) |
| MLSS ranking (WP3 scenarios S01) | reported | Spearman 1.000000, top 15 kept 15, largest score change 0.0152 | Spearman 1.000000, top 15 kept 15, largest score change 0.0141 |
**Reading.** Under the corrected version (C0) the biology of the five genes does not change: the DAPT fold changes are the same to the second decimal, and the MLSS ranking is identical in order. What changes is the composition of the network at the edges: 4 new pairs and 1 lost, Wnt5a joins the hubs, and the Fig. 5 classes "known" and "database-unregistered" hold 6 and 18 pairs instead of 5 and 16 (three more consensus pairs in all). Only if the failed wells were counted (C1) would the DAPT effect on four genes reverse and half the consensus pairs change; the plateau data argue against doing that. The statement "the network is not stable to the reference" of the manuscript stands, and C0 shows that it is stable to the correct handling of A10 and B10 within a few pairs.

## 3. Wording that the data allow (omitted)

Omitted from this copy: this section proposed wording for the article; the wording used is that of the article.

## 4. Decision for the authors (not taken here)
* **Option A (recommended): keep the reported analysis as the primary one and report C0 and C1 as sensitivity analyses.** Cost: two or three sentences, one more part in Additional file 5 and one sheet in Additional file 6; all figures and tables stay.
* **Option B: make C0 the primary analysis.** Every number that depends on the INS-1 matrix changes slightly (Figs 2 to 6, S-figures, Tables, Additional file 4, the benchmark of WP1 to WP4 that uses the 21 consensus pairs, the islet cross-check with its evaluable pairs): it would repeat most of the work of plan v2. The conclusions would not change, so this is not recommended for this submission.

## 5. Limits
Post hoc; one preparation; the pairing of target and Gapdh replicates of the sheet is by rank for this run (WP1, WP11), which C0 sidesteps for the DAPT group (a constant reference). The versions change only the 26 August DAPT group; the other runs, the other groups and the 24 June thresholds are as reported. The melt-peak heights are quoted from the earlier QC, not re-derived. The docstring of the earlier script `reanalysis_17genes/gapdh_exact.py` gives the C1 reference as 23.6225; the value that the script computes, and that is used here, is the mean of the four raw wells, 23.6215 (`results/gapdh_versions.json`, key `ref`).

## 6. Files
`wp13_gapdh_versions.py` (reads `reanalysis_17genes/pipeline.py`; no file of that folder is changed), `check_wp13_numbers.py`, `results/gapdh_versions.json`, `results/gapdh_versions_summary.csv`. Use `python -I`.
