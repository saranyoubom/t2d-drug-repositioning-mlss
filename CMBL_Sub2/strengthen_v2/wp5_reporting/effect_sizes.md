# Effect sizes with 95% intervals, across-gene BH q and Gapdh sensitivity (WP5a)

Produced by `effect_sizes.py` from Additional files 3 and 4 and `gapdh_exact.json`. Every interval below is **within-experiment, one preparation, wells are not independent biological replicates**.

## How to read the tables

* **Values used.** For each gene, the 12 relative expression values (2^-dCt, one per well-level replicate; Additional file 4, sheet `Fig_4c_Raw_expression_data`, full precision) divided by the arithmetic mean of the four CTRL values of that gene. These are the replicate fold changes R1-R4 of Additional file 3 (reproduced to rounding, see README gate). Rows 1-4 are CTRL, 5-8 DAPT, 9-12 DKK-1. DAPT values of the five 26 August genes are as reported unless a variant is named.
* **log2 FC [95% CI].** Each replicate fold change is log2-transformed. Effect = mean log2 (treated) minus mean log2 (CTRL), that is log2 of the ratio of geometric means. Interval = Welch t-interval (n = 4 vs 4; Welch-Satterthwaite df, range in the summary below), so no equal-variance assumption is made. Not adjusted for the three pairwise contrasts per gene or for the 17 genes. The same values rounded to the 3 decimals printed in Additional file 3 move the interval limits by at most 0.0016 log2 units.
* **Mean +/- SEM.** As printed in Additional file 3 (fold change of the arithmetic mean relative to the CTRL arithmetic mean; SEM of the four replicate fold changes). The manuscript fold change is the ratio of arithmetic means, so it differs slightly from 2^(log2 FC).
* **p (Bonf.).** The manuscript's Bonferroni-corrected pairwise p (Additional file 3): pooled-variance t-tests after one-way ANOVA, or Wilcoxon rank-sum tests after Kruskal-Wallis for Wnt5a and Lef1 (marked KW; the smallest attainable adjusted p is 0.086).
* **BH q.** Benjamini-Hochberg adjustment over the 17 raw pairwise p of the same contrast (CTRL vs DAPT, or CTRL vs DKK-1). The raw p are those of the manuscript's tests (pooled-variance t-test, or Wilcoxon rank-sum for the two KW genes) recomputed from the 12 full-precision values. The p printed in Additional file 3 were computed from its 3-decimal replicate values (all 102 of 102 printed pairwise p are reproduced from those values); the full-precision p differ from them by a median of 0.2% and at most 2% (relative), and q changes by at most 0.0025 with no change in whether q < 0.05. q is reported beside the manuscript's tests and does not replace them; it answers a different question (17 genes) from the Bonferroni p (3 contrasts per gene). q < 0.05 is marked in bold.
* **Excludes 0 (pooled).** Sensitivity: whether the interval excludes 0 when the within-group variance is pooled over the three groups (df = 9), the model behind the manuscript's pairwise tests.
* **†** The DAPT value depends on the Gapdh wells of the 26 August run (Tcf7l2, Kcnj11, Cacna1d, Cacna1c, Ptbp1); see the section on Gapdh variants.

## Summary of the 17-gene panel

Welch-Satterthwaite df range 3.1 to 6.0; median width of the 95% interval 0.62 log2 units (a 1.5-fold range). Counts below are out of 17 genes per contrast.

| Contrast | Manuscript Bonferroni p < 0.05 | BH q < 0.05 | Welch interval excludes 0 | Pooled interval excludes 0 | Genes where Bonferroni and BH disagree | Genes where Welch interval and Bonferroni disagree |
|---|---|---|---|---|---|---|
| DAPT vs CTRL | 11 | 14 | 13 | 13 | *Wnt5a* (q only), *Ryr3* (q only), *Ins1* (q only) | *Wnt5a* (interval only), *Ins1* (interval only) |
| DKK-1 vs CTRL | 6 | 6 | 7 | 7 | none | *Wnt5a* (interval only) |

Bonferroni p < 0.05 is the manuscript's own criterion. The intervals are unadjusted, so a few genes whose Bonferroni p is just above 0.05 have intervals that exclude 0.

## RT-qPCR: DAPT vs CTRL (17 genes)

Interval label: within-experiment, one preparation, wells are not independent biological replicates.

| Gene | Test | CTRL mean +/- SEM | DAPT mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.), manuscript | BH q | Excludes 0 (Welch / pooled) |
|---|---|---|---|---|---|---|---|---|
| *Wnt5a* | KW | 1 +/- 0.056 | 0.552 +/- 0.032 | -0.86 [-1.14, -0.58] | 0.55 [0.46, 0.67] | 0.086 | **0.037** | yes / yes |
| *Wnt5b* | ANOVA | 1 +/- 0.054 | 0.696 +/- 0.07 | -0.54 [-1.00, -0.08] | 0.69 [0.50, 0.95] | 0.023 | **0.012** | yes / yes |
| *Wnt9a* | ANOVA | 1 +/- 0.027 | 1.93 +/- 0.11 | 0.94 [0.68, 1.20] | 1.92 [1.60, 2.30] | 3.44e-05 | **4.88e-05** | yes / yes |
| *Lef1* | KW | 1 +/- 0.097 | 1.171 +/- 0.036 | 0.24 [-0.16, 0.65] | 1.19 [0.89, 1.57] | 0.600 | 0.200 | no / no |
| *Tcf7* | ANOVA | 1 +/- 0.041 | 1.296 +/- 0.145 | 0.35 [-0.15, 0.85] | 1.27 [0.90, 1.80] | 0.362 | 0.128 | no / no |
| *Tcf7l2*† | ANOVA | 1 +/- 0.017 | 0.778 +/- 0.038 | -0.37 [-0.59, -0.15] | 0.78 [0.67, 0.90] | 0.002 | **0.001** | yes / yes |
| *Hes1* | ANOVA | 1 +/- 0.043 | 0.136 +/- 0.003 | -2.87 [-3.06, -2.68] | 0.14 [0.12, 0.16] | 3.03e-08 | **1.69e-07** | yes / yes |
| *Itpr1* | ANOVA | 1 +/- 0.038 | 0.554 +/- 0.098 | -0.91 [-1.65, -0.17] | 0.53 [0.32, 0.89] | 0.007 | **0.005** | yes / yes |
| *Ryr2* | ANOVA | 1 +/- 0.128 | 0.531 +/- 0.03 | -0.88 [-1.44, -0.33] | 0.54 [0.37, 0.79] | 0.010 | **0.007** | yes / yes |
| *Ryr3* | ANOVA | 1 +/- 0.15 | 0.612 +/- 0.027 | -0.66 [-1.35, 0.02] | 0.63 [0.39, 1.02] | 0.062 | **0.029** | no / yes |
| *Cacna1c*† | ANOVA | 1 +/- 0.061 | 0.712 +/- 0.066 | -0.50 [-0.93, -0.07] | 0.71 [0.52, 0.95] | 0.013 | **0.007** | yes / yes |
| *Cacna1d*† | ANOVA | 1 +/- 0.038 | 0.763 +/- 0.034 | -0.39 [-0.60, -0.18] | 0.76 [0.66, 0.88] | 0.003 | **0.002** | yes / yes |
| *Kcnj11*† | ANOVA | 1 +/- 0.038 | 0.396 +/- 0.031 | -1.35 [-1.70, -1.00] | 0.39 [0.31, 0.50] | 1.03e-06 | **2.91e-06** | yes / yes |
| *Ptbp1*† | ANOVA | 1 +/- 0.007 | 0.616 +/- 0.027 | -0.70 [-0.90, -0.51] | 0.61 [0.54, 0.70] | 9.17e-06 | **1.72e-05** | yes / yes |
| *Glut2* | ANOVA | 1 +/- 0.044 | 0.848 +/- 0.052 | -0.24 [-0.52, 0.03] | 0.85 [0.70, 1.02] | 0.205 | 0.077 | no / no |
| *Ins1* | ANOVA | 1 +/- 0.056 | 1.276 +/- 0.056 | 0.36 [0.10, 0.61] | 1.28 [1.07, 1.53] | 0.113 | **0.046** | yes / no |
| *Glp1r* | ANOVA | 1 +/- 0.038 | 0.577 +/- 0.038 | -0.80 [-1.08, -0.52] | 0.57 [0.47, 0.70] | 0.0003 | **3.49e-04** | yes / yes |

## RT-qPCR: DKK-1 vs CTRL (17 genes)

Interval label: within-experiment, one preparation, wells are not independent biological replicates.

| Gene | Test | CTRL mean +/- SEM | DKK-1 mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.), manuscript | BH q | Excludes 0 (Welch / pooled) |
|---|---|---|---|---|---|---|---|---|
| *Wnt5a* | KW | 1 +/- 0.056 | 0.611 +/- 0.03 | -0.71 [-0.97, -0.45] | 0.61 [0.51, 0.73] | 0.086 | 0.069 | yes / yes |
| *Wnt5b* | ANOVA | 1 +/- 0.054 | 0.925 +/- 0.064 | -0.12 [-0.45, 0.22] | 0.92 [0.73, 1.16] | 1.000 | 0.653 | no / no |
| *Wnt9a* | ANOVA | 1 +/- 0.027 | 1.079 +/- 0.067 | 0.10 [-0.17, 0.37] | 1.07 [0.89, 1.29] | 1.000 | 0.681 | no / no |
| *Lef1* | KW | 1 +/- 0.097 | 1.006 +/- 0.144 | -0.01 [-0.58, 0.56] | 0.99 [0.67, 1.48] | 1.000 | 0.897 | no / no |
| *Tcf7* | ANOVA | 1 +/- 0.041 | 0.689 +/- 0.149 | -0.67 [-1.85, 0.52] | 0.63 [0.28, 1.43] | 0.314 | 0.197 | no / no |
| *Tcf7l2* | ANOVA | 1 +/- 0.017 | 1.003 +/- 0.031 | 0.00 [-0.13, 0.14] | 1.00 [0.91, 1.10] | 1.000 | 0.947 | no / no |
| *Hes1* | ANOVA | 1 +/- 0.043 | 0.517 +/- 0.031 | -0.95 [-1.24, -0.67] | 0.52 [0.42, 0.63] | 4.68e-06 | **2.61e-05** | yes / yes |
| *Itpr1* | ANOVA | 1 +/- 0.038 | 0.838 +/- 0.075 | -0.27 [-0.64, 0.10] | 0.83 [0.64, 1.07] | 0.475 | 0.269 | no / no |
| *Ryr2* | ANOVA | 1 +/- 0.128 | 1.038 +/- 0.064 | 0.08 [-0.47, 0.63] | 1.06 [0.72, 1.55] | 1.000 | 0.921 | no / no |
| *Ryr3* | ANOVA | 1 +/- 0.15 | 1.01 +/- 0.075 | 0.05 [-0.62, 0.72] | 1.04 [0.65, 1.65] | 1.000 | 0.947 | no / no |
| *Cacna1c* | ANOVA | 1 +/- 0.061 | 0.823 +/- 0.026 | -0.28 [-0.55, 0.00] | 0.83 [0.68, 1.00] | 0.135 | 0.097 | no / no |
| *Cacna1d* | ANOVA | 1 +/- 0.038 | 0.732 +/- 0.031 | -0.45 [-0.66, -0.25] | 0.73 [0.63, 0.84] | 0.001 | **0.002** | yes / yes |
| *Kcnj11* | ANOVA | 1 +/- 0.038 | 0.659 +/- 0.028 | -0.60 [-0.80, -0.40] | 0.66 [0.57, 0.76] | 0.0001 | **2.21e-04** | yes / yes |
| *Ptbp1* | ANOVA | 1 +/- 0.007 | 0.66 +/- 0.037 | -0.60 [-0.85, -0.36] | 0.66 [0.55, 0.78] | 2.54e-05 | **7.21e-05** | yes / yes |
| *Glut2* | ANOVA | 1 +/- 0.044 | 0.635 +/- 0.059 | -0.67 [-1.06, -0.28] | 0.63 [0.48, 0.82] | 0.002 | **0.003** | yes / yes |
| *Ins1* | ANOVA | 1 +/- 0.056 | 0.991 +/- 0.114 | -0.03 [-0.53, 0.46] | 0.98 [0.69, 1.37] | 1.000 | 0.947 | no / no |
| *Glp1r* | ANOVA | 1 +/- 0.038 | 0.746 +/- 0.058 | -0.43 [-0.77, -0.10] | 0.74 [0.59, 0.94] | 0.010 | **0.010** | yes / yes |

## The five genes whose DAPT values depend on the 26 August Gapdh wells

DAPT vs CTRL; CTRL and DKK-1 values are identical in all three variants. Reported = Gapdh values as analyzed in the authors' calculation sheet (C10 and D10 entered twice, A10 and B10 set aside). C0 = reference is the mean of the two valid wells (Cq 22.68) for all four DAPT replicates. C1 = reference is the mean of all four wells (Cq 23.62), failed wells included (worst case). Same interval method as above, label: within-experiment, one preparation, wells are not independent biological replicates. BH q is recomputed over the 17 genes of the DAPT contrast with the five p replaced by the variant values.

| Gene | Variant | Mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.) | BH q | Sign of log2 FC |
|---|---|---|---|---|---|---|---|
| *Tcf7l2* | reported | 0.778 +/- 0.038 | -0.37 [-0.59, -0.15] | 0.78 [0.67, 0.90] | 0.002 | **0.001** | down |
| *Kcnj11* | reported | 0.396 +/- 0.031 | -1.35 [-1.70, -1.00] | 0.39 [0.31, 0.50] | 1.03e-06 | **2.91e-06** | down |
| *Cacna1d* | reported | 0.763 +/- 0.034 | -0.39 [-0.60, -0.18] | 0.76 [0.66, 0.88] | 0.003 | **0.002** | down |
| *Cacna1c* | reported | 0.712 +/- 0.066 | -0.50 [-0.93, -0.07] | 0.71 [0.52, 0.95] | 0.013 | **0.007** | down |
| *Ptbp1* | reported | 0.616 +/- 0.027 | -0.70 [-0.90, -0.51] | 0.61 [0.54, 0.70] | 9.09e-06 | **1.72e-05** | down |
| *Tcf7l2* | C0 (two valid wells) | 0.779 +/- 0.045 | -0.37 [-0.61, -0.13] | 0.78 [0.66, 0.92] | 0.003 | **0.003** | down |
| *Kcnj11* | C0 (two valid wells) | 0.398 +/- 0.038 | -1.35 [-1.80, -0.89] | 0.39 [0.29, 0.54] | 2.09e-06 | **5.91e-06** | down |
| *Cacna1d* | C0 (two valid wells) | 0.762 +/- 0.017 | -0.39 [-0.56, -0.22] | 0.76 [0.68, 0.86] | 0.001 | **9.46e-04** | down |
| *Cacna1c* | C0 (two valid wells) | 0.715 +/- 0.074 | -0.50 [-1.01, 0.01] | 0.71 [0.50, 1.00] | 0.020 | **0.012** | down |
| *Ptbp1* | C0 (two valid wells) | 0.616 +/- 0.028 | -0.70 [-0.92, -0.49] | 0.61 [0.53, 0.71] | 1.06e-05 | **2.01e-05** | down |
| *Tcf7l2* | C1 (all four wells) | 1.493 +/- 0.086 | 0.57 [0.33, 0.81] | 1.49 [1.26, 1.76] | 3.26e-04 | **3.70e-04** | up |
| *Kcnj11* | C1 (all four wells) | 0.763 +/- 0.074 | -0.41 [-0.86, 0.04] | 0.75 [0.55, 1.03] | 0.027 | **0.015** | down |
| *Cacna1d* | C1 (all four wells) | 1.459 +/- 0.033 | 0.55 [0.38, 0.71] | 1.46 [1.30, 1.64] | 1.69e-05 | **4.78e-05** | up |
| *Cacna1c* | C1 (all four wells) | 1.371 +/- 0.142 | 0.44 [-0.07, 0.94] | 1.35 [0.95, 1.93] | 0.054 | **0.028** | up |
| *Ptbp1* | C1 (all four wells) | 1.180 +/- 0.054 | 0.23 [0.02, 0.45] | 1.18 [1.01, 1.37] | 0.025 | **0.015** | up |

The "reported" rows are the main-table values; all p and q in this table are recomputed from the full-precision values (pooled-variance t-tests, df = 9), so their third digit can differ from the printed p of Additional file 3. Fold changes and Bonferroni p of all three variants match `gapdh_exact.json` (max absolute difference of fold changes 2.7e-15, max relative difference of p 9.7e-14). With C1 the sign of the log2 FC reverses for *Tcf7l2*, *Cacna1d*, *Cacna1c*, *Ptbp1*. The 95% interval includes 0 under C1 for *Kcnj11*, *Cacna1c*, under C0 for *Cacna1c* and under the reported values for none.

## DNA content, 2-NBDG and C-peptide per ug DNA

Interval label: within-experiment, one preparation, wells are not independent biological replicates. Intervals are Welch t-intervals on the log2 scale of the replicate values, as above.

| Readout | Contrast | Variant | n (CTRL / treated) | Mean ratio (arith.) | log2 ratio [95% CI] | 2^log2 ratio [95% CI] | p (Bonf.), manuscript |
|---|---|---|---|---|---|---|---|
| Total DNA (ug), wells at 2.8, 5.5, 22, 44 mM glucose | DAPT vs CTRL | as reported | 4 / 4 | 2.012 | 1.12 [0.02, 2.23] | 2.18 [1.02, 4.68] | 0.013 |
| Total DNA (ug), all five wells incl. basal (sensitivity) | DAPT vs CTRL | sensitivity: five wells | 5 / 5 | 2.274 | 1.35 [0.37, 2.33] | 2.55 [1.29, 5.04] |  |
| Total DNA (ug), wells at 2.8, 5.5, 22, 44 mM glucose | DKK-1 vs CTRL | as reported | 4 / 4 | 2.166 | 1.25 [0.11, 2.38] | 2.37 [1.08, 5.21] | 0.005 |
| Total DNA (ug), all five wells incl. basal (sensitivity) | DKK-1 vs CTRL | sensitivity: five wells | 5 / 5 | 2.566 | 1.53 [0.54, 2.52] | 2.89 [1.46, 5.72] |  |
| 2-NBDG fluorescence, scenario S1 | DAPT vs CTRL | published (Fig. 2b) | 3 / 4 | 0.946 | -0.13 [-1.14, 0.88] | 0.91 [0.45, 1.84] | 1.000 |
| 2-NBDG fluorescence, scenario S1 | DKK-1 vs CTRL | published (Fig. 2b) | 3 / 4 | 1.288 | 0.21 [-1.11, 1.54] | 1.16 [0.46, 2.90] | 1.000 |
| C-peptide per ug group-mean DNA (workbook units), 5.5 mM glucose | DAPT vs CTRL | as reported | 4 / 4 | 0.797 | -0.28 [-0.92, 0.36] | 0.82 [0.53, 1.29] | 0.077 |
| C-peptide per ug group-mean DNA (workbook units), 5.5 mM glucose | DAPT vs CTRL | with DNA-denominator uncertainty | 4 / 4 | 0.797 | -0.28 [-1.29, 0.73] | 0.82 [0.41, 1.66] |  |
| C-peptide per ug group-mean DNA (workbook units), 5.5 mM glucose | DAPT vs CTRL | sensitivity: CTRL artefact well left out | 3 / 4 | 0.948 | -0.08 [-0.12, -0.03] | 0.95 [0.92, 0.98] |  |
| C-peptide per ug group-mean DNA (workbook units), 5.5 mM glucose | DKK-1 vs CTRL | as reported | 4 / 4 | 0.617 | -0.65 [-1.29, -0.01] | 0.64 [0.41, 1.00] | 0.085 |
| C-peptide per ug group-mean DNA (workbook units), 5.5 mM glucose | DKK-1 vs CTRL | with DNA-denominator uncertainty | 4 / 4 | 0.617 | -0.65 [-1.66, 0.36] | 0.64 [0.32, 1.28] |  |
| C-peptide per ug group-mean DNA (workbook units), 5.5 mM glucose | DKK-1 vs CTRL | sensitivity: CTRL artefact well left out | 3 / 4 | 0.733 | -0.45 [-0.49, -0.41] | 0.73 [0.71, 0.75] |  |
| ELISA signal dOD, 5.5 mM glucose | DAPT vs CTRL | as reported | 4 / 4 | 1.813 | 0.90 [0.26, 1.55] | 1.87 [1.20, 2.92] | 0.077 |
| ELISA signal dOD, 5.5 mM glucose | DKK-1 vs CTRL | as reported | 4 / 4 | 1.582 | 0.71 [0.06, 1.35] | 1.63 [1.05, 2.55] | 0.085 |

Notes on this table.

* **DNA content.** The four values per group are wells that differed in their 20-minute glucose exposure (Methods), not four replicates of one condition, and the CTRL wells range from 0.091 to 0.239 ug. The lower interval limit is 0.02 (DAPT) and 0.11 (DKK-1) log2 units above 0, that is, barely above no change. Adding the fifth (basal) well, which the authors left out for an unrecorded reason, gives 1.35 and 1.53 log2 units with lower limits 0.37 and 0.54 (post hoc sensitivity, not a replacement).
* **2-NBDG: not supported.** Across the 12 scenarios of Table S6, 0 of 24 intervals exclude 0, the log2 ratio of DAPT runs from -1.03 to 0.28 and that of DKK-1 from -0.69 to 0.70, so even the direction depends on the analysis choice. The published analysis (S1) uses the CTRL values 195, 331 and 331 (the value 331 appears twice because the analysis block of the workbook holds 331 for CTRL replicate 2, where the reading block holds 184; the reason is not recorded). Only S1 is shown in the table above; all twelve are in the scenario table below and in the CSV.
* **C-peptide at 5.5 mM glucose.** The raw signal (dOD) is higher than CTRL in both treated groups (mean ratio 1.81 and 1.58) whereas the value per ug of group-mean DNA is lower (0.80 and 0.62): the sign of the difference is set by the DNA denominator, which is 2.3 and 2.6 times larger in the treated groups (Table S7). Three of the four CTRL wells read 0.047 to 0.048, only 0.003 to 0.004 above the lowest standard (0.044), and the fourth (0.083) carries the optical artefact noted in Table S7. The "wells only" interval treats the DNA denominator as an exact constant; the row "with DNA-denominator uncertainty" adds the variance of the five DNA values per group (one well per glucose condition) and is much wider. The narrow intervals of the row "CTRL artefact well left out" come from near-identical readings (the readings are 0.001 apart in dOD) and describe the spread of the readings, not the uncertainty of the biology.

### 2-NBDG: range over the twelve analysis scenarios of Table S6

| Scenario | CTRL values used | DAPT mean ratio | DAPT log2 ratio [95% CI] | DKK-1 mean ratio | DKK-1 log2 ratio [95% CI] |
|---|---|---|---|---|---|
| S1 | 195, 331, 331 | 0.95 | -0.13 [-1.14, 0.88] | 1.29 | 0.21 [-1.11, 1.54] |
| S2 | 195, 184, 331 | 1.14 | 0.15 [-0.88, 1.19] | 1.55 | 0.50 [-0.84, 1.83] |
| S3 | 195, 331, 331, 3382 | 0.26 | -1.03 [-3.83, 1.76] | 0.35 | -0.69 [-3.44, 2.07] |
| S4 | 195, 184, 331, 3382 | 0.26 | -0.82 [-3.85, 2.21] | 0.36 | -0.48 [-3.45, 2.50] |
| S5 | 121, 90, 209 | 1.22 | 0.17 [-1.30, 1.64] | 1.92 | 0.70 [-1.06, 2.46] |
| S6 | 121, 90, 209, 3276 | 0.19 | -0.99 [-4.57, 2.59] | 0.29 | -0.46 [-4.01, 3.09] |
| S7 | 96, 85, 232, 3283 | 0.19 | -0.97 [-4.65, 2.71] | 0.29 | -0.49 [-4.17, 3.19] |
| S8 | 189, 161, 277, 3435 | 0.26 | -0.74 [-3.91, 2.43] | 0.34 | -0.41 [-3.52, 2.70] |
| S9 | 189, 161, 277 | 1.25 | 0.28 [-0.66, 1.22] | 1.67 | 0.61 [-0.62, 1.83] |
| S10 | 397, 701, 1020, 1746 | 1.15 | 0.08 [-1.74, 1.91] | 1.64 | 0.54 [-1.41, 2.50] |
| S11 | 453, 670, 1054, 1711 | 1.02 | -0.03 [-1.65, 1.59] | 1.43 | 0.41 [-1.25, 2.08] |
| S12 | 379, 576, 932, 1605 | 1.03 | -0.04 [-1.82, 1.74] | 1.48 | 0.46 [-1.33, 2.24] |

### C-peptide: which glucose conditions the data support

A comparison is estimated only where no well of any of the three groups lies below the lowest standard (dOD < 0.044). Wells below the lowest standard per group (CTRL / DAPT / DKK-1):

| Glucose (mM) | Below lowest standard | Estimated here |
|---|---|---|
| 0.0 | 0 / 4 / 2 | no: at least one group has wells below the lowest standard |
| 2.8 | 0 / 4 / 4 | no: at least one group has wells below the lowest standard |
| 5.5 | 0 / 0 / 0 | yes |
| 22.0 | 0 / 4 / 0 | no: at least one group has wells below the lowest standard |
| 44.0 | 3 / 1 / 0 | no: at least one group has wells below the lowest standard |
