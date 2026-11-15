# WP5a: effect sizes with intervals, across-gene BH q and MIQE checklist (4 Nov 2026)

Protocol section 6 (reporting assets). Nothing here changes the manuscript, the JTM package, the figures or the authors' folders; no commit or push. All intervals are labelled **within-experiment, one preparation, wells are not independent biological replicates**.

## Files

| File | Content |
|---|---|
| `effect_sizes.py` | Computes everything below from `Additional file 3.xlsx`, `Additional file 4.xlsx` and `reanalysis_17genes/gapdh_exact.json` (read only). `python -I effect_sizes.py`; the project root is found from the script location, `CMBL_ROOT` overrides it. Runs in under a minute. |
| `effect_sizes.csv` | 80 rows, 37 columns: 34 gene contrasts (17 genes x DAPT and DKK-1 vs CTRL, as reported), 10 rows for the Gapdh variants C0 and C1 of the five 26 August genes (DAPT), 4 DNA rows (two of them a five-well sensitivity), 24 2-NBDG rows (12 scenarios x 2), 8 C-peptide rows (5.5 mM only). |
| `effect_sizes.md` | The same as readable tables, with the "how to read" notes and a summary of the 17-gene panel. |
| `gate_reproduction.json` | Every number of the gate below. |
| `check_effect_sizes.py` | Independent re-computation (scipy, statsmodels OLS, hand-written BH) of the 34 intervals and the q values; differences below 5e-6. |
| `miqe_checklist.md` | MIQE essential items from the Methods text, Additional file 1 and the named files; ends with the block of items the authors must supply. |

## Gate: reported means reproduced before anything was added (passed)

The 12 full-precision relative expression values per gene (Additional file 4, sheet `Fig_4c_Raw_expression_data`; identical to the SBA `Fig_4b` matrix, largest relative difference 3.8e-15) divided by the CTRL arithmetic mean reproduce Additional file 3:

| Quantity | Largest absolute difference | Count |
|---|---|---|
| Replicate fold changes R1-R4 | 0.00049 | 204 values |
| Group means | 0.00075 (*Ptbp1* DKK-1: 0.66075 against 0.66 printed) | 51 means |
| Group SEMs | 0.00055 | 51 SEMs |
| CTRL mean | 2e-16 from 1 | 17 genes |

Additional file 3 prints 3-decimal replicate values and its means, SEMs and p were computed from those rounded values: from them the means and SEMs agree to 0.0005, and **all 102 printed pairwise p (raw and Bonferroni, three contrasts) and all 17 overall p are reproduced** at the printed precision. From the full-precision values 72 of the 102 and 16 of the 17 agree at the printed precision; the full-precision raw p differ from those computed from the printed replicates by at most 2.0% (median 0.2%, relative). Using the printed replicates instead of the full-precision values would move the interval limits by at most 0.0016 log2 units and q by at most 0.0025, with no change in whether q < 0.05 and no change in any significance call. The tests were identified as pooled-variance t-tests with the mean square over all three groups (df = 9, Bonferroni x 3) and Wilcoxon rank-sum tests for *Wnt5a* and *Lef1*.

Other reproductions: the Gapdh variants (reported, C0, C1) give the fold changes and Bonferroni p of `gapdh_exact.json` to 2.7e-15 (absolute) and 9.7e-14 (relative) and the fold changes of Additional file 1 Table S1c at two decimals; DNA means, SEMs and Bonferroni p match `Fig_1c_stats` to 4.4e-5 and 3.1e-7; the 12 2-NBDG scenarios match Table S6 (means and SEMs) to 4.4e-15; the C-peptide 5.5 mM means and SEMs match `Fig_3b_stats` to 5e-5.

## What the numbers are

* **Values used.** Per gene, the four replicate fold changes of each group (rows 1-4 CTRL, 5-8 DAPT, 9-12 DKK-1), log2-transformed. DAPT values of *Tcf7l2*, *Kcnj11*, *Cacna1d*, *Cacna1c* and *Ptbp1* are the reported ones except in the variant rows.
* **Effect and interval.** Difference of the group means of the log2 values (log2 of the ratio of geometric means) with a Welch t-interval (n = 4 vs 4, df 3.1 to 6.0), not adjusted for the three contrasts per gene or the 17 genes. The manuscript's fold change (ratio of arithmetic means) is listed beside it. Median interval width 0.62 log2 units (a 1.5-fold range). With the three-group pooled variance (df = 9) the verdict "interval excludes 0" differs for 2 of 34 contrasts (*Ryr3* and *Ins1*, DAPT); 20 of 34 intervals exclude 0 either way.
* **BH q.** Over the 17 raw pairwise p of each contrast, beside and not instead of the manuscript's Bonferroni p. Bonferroni p < 0.05: 11 genes (DAPT) and 6 (DKK-1); BH q < 0.05: 14 and 6. The three extra DAPT genes are *Wnt5a* (the rank-based test cannot reach a Bonferroni p below 0.086), *Ryr3* and *Ins1*. The two corrections answer different questions (3 contrasts against 17 genes), so this is a difference in family, not in evidence.

## The three intervals most worth knowing

1. ***Hes1* under DAPT: log2 FC -2.87 [-3.06, -2.68]**, that is 0.14-fold [0.12, 0.16] (manuscript 0.136 +/- 0.003; BH q 1.7e-07). The largest and tightest effect and the only Notch readout; it does not depend on the 26 August run. Under DKK-1 it is -0.95 [-1.24, -0.67] (0.52-fold).
2. ***Kcnj11* under DAPT depends on the Gapdh wells: reported -1.35 [-1.70, -1.00] (0.39-fold), C0 -1.35 [-1.80, -0.89], C1 -0.41 [-0.86, +0.04] (0.75-fold).** Under C1 the sign of the DAPT effect reverses for *Tcf7l2* (+0.57 [0.33, 0.81]), *Cacna1d* (+0.55 [0.38, 0.71]), *Cacna1c* (+0.44 [-0.07, 0.94]) and *Ptbp1* (+0.23 [0.02, 0.45]). The five genes have BH q < 0.05 under all three variants, but with opposite signs under C1 for four of them. C0 moves every interval limit by at most 0.11 log2 units, but the *Cacna1c* interval now reaches 0 (-0.50 [-1.01, +0.01]).
3. **DNA content: DAPT +1.12 [0.02, 2.23] and DKK-1 +1.25 [0.11, 2.38] log2 units**, that is 2.2-fold [1.0, 4.7] and 2.4-fold [1.1, 5.2]. The direction is clear, the size is not, and the lower limits are barely above 0 because the four CTRL wells (which differ in glucose exposure) range from 0.091 to 0.239 ug. Adding the omitted basal well (post hoc) gives +1.35 and +1.53.

Also useful: *Wnt9a* is the only clear increase with DAPT, +0.94 [0.68, 1.20] (1.9-fold [1.6, 2.3]); it does not depend on the 26 August run. 2-NBDG and C-peptide are not supported (next section).

## Where the data do not support an interval

* **2-NBDG.** All 24 intervals (12 scenarios x 2 contrasts) include 0 and the DAPT log2 ratio runs from -1.03 to +0.28 (DKK-1 -0.69 to +0.70) depending on the analysis choice. The published CTRL set (scenario S1) is 195, 331, 331: the value 331 appears twice because the analysis block holds 331 for CTRL replicate 2 where the reading block holds 184 (disclosed in the Methods, reason not recorded).
* **C-peptide.** Estimated only at 5.5 mM glucose, the only condition in which no well of any group lies below the lowest standard (22, 44, 0 and 2.8 mM have wells below it in at least one group). There the raw signal is higher than CTRL (ratio 1.81 DAPT, 1.58 DKK-1) while the value per ug group-mean DNA is lower (0.80, 0.62): the sign is set by the DNA denominator, which is 2.3 and 2.6 times larger in the treated groups. The wells-only intervals treat that denominator as a constant; the extra rows with its variance (five DNA values per group, one per glucose condition) are about 1.6 times as wide and include 0 for both groups.

## Unusual in the data files

1. **Rounded inputs to the printed statistics** (see the gate): harmless, but a reader recomputing from the full-precision sheet gets p that differ in the third digit.
2. **Normalization wording.** The Methods say the delta-delta Ct is relative to "the mean dCt of all CTRL replicates". Additional file 3 shows a CTRL mean of exactly 1.000 for all genes, which is the arithmetic mean of 2^-dCt (as the Table S1c remark says). Normalizing to the mean dCt would give a CTRL arithmetic mean of up to 1.036 (median 1.003). The sentence needs correcting (listed in the MIQE block).
3. **Near-identical ELISA wells in the treated groups.** At 5.5 mM glucose the within-group coefficient of variation of dOD is 2.0% (DAPT) and 1.6% (DKK-1) against 32% in CTRL (one artefact well), and the four DAPT wells have the same OD450 (0.139). The Methods describe them as four separate culture wells; the authors should confirm that they are not repeat readings of one supernatant.
4. **`Fig_3b_stats` title** says "ng/g DNA" while its table header says "ng/ug DNA" (the values are per ug of group-mean DNA, workbook units).
5. **DNA "replicates".** The four DNA values per group are wells that differ in their glucose exposure, not replicates of one condition (the Methods say so); the interval treats them as exchangeable.

## Choices and deviations

* The brief did not fix the interval type. Welch is primary because it makes no equal-variance assumption across groups with n = 4; the pooled-variance interval is in the CSV (`ci95_lo_pooled3`, `ci95_hi_pooled3`).
* Rows labelled "sensitivity" or "with DNA-denominator uncertainty" (five DNA wells, CTRL artefact well left out, denominator variance) are post hoc and are not replacements for the authors' analyses.
* No deviation from `PROTOCOL_v2_prespecified.md` (section 8 untouched; no fingerprint change). The extras beyond protocol section 6 (pooled interval, DNA, 2-NBDG and C-peptide rows) come from the WP5a brief.
* Software: Python 3.14.4, numpy 2.4.4, scipy 1.17.1, statsmodels 0.14.6, openpyxl 3.1.5 (pandas only for the optional cross-check with the SBA matrix). Nothing random is used.

## For later work packages

`effect_sizes.csv` can feed a new sheet of Additional file 4 and a column of log2 FC with interval in Additional file 3 (WP6a), and the MIQE block can become a supplementary table (WP6b); the Methods wording in note 1 of `miqe_checklist.md` needs an author decision.
