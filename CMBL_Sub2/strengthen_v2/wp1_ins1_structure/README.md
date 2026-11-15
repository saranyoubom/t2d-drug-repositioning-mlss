# WP1: INS-1 run structure and calibrated checks (post hoc), 4 Nov 2026

**Strength of every result below.** These are post hoc analyses of one cell preparation (3 groups x 4 wells, n = 12) done after the manuscript analysis was fixed. They describe this data set; none is confirmatory, none is a pre-specified hypothesis test about biology, and "P" values are descriptive. The protocol is `../PROTOCOL_v2_prespecified.md` (section 2, fingerprint unchanged at the start of this work package); anything added beyond the protocol is marked ADD-ON.

## 1. What is here

| File | Content |
|---|---|
| `../coexp_checks.py` | reusable module and command line (API in its docstring; section 9 below) |
| `../test_coexp_checks.py`, `test_output.txt` | regression and known-answer tests, 56 checks; the output is from the last full run (afterwards only the speed note in the module docstring was edited, no code line) |
| `build_gene_run_map.py` -> `gene_run_map.csv`, `gene_run_evidence.csv`, `ins1_sheet_cq.csv` | item 2.1: verified gene-to-run map, every gene x candidate comparison, the matched sheet rows (target and reference Cq per replicate) |
| `check_sheet_vs_exports.py` -> `sheet_vs_exports.csv` | provenance check of the sheet Cq values against the CFX exports in the folders (read only) |
| `wp1_ins1_analysis.py` -> `results_ins1.json`, `tbl_2_*.csv`, `analysis_log.txt` | items 2.2 to 2.8 (about 4 minutes; seed 2026) |

Run order, all with `python -I` from this folder: `build_gene_run_map.py`, `check_sheet_vs_exports.py`, `wp1_ins1_analysis.py`; tests: `python -I ../test_coexp_checks.py`. Sources are only read: the three calculation workbooks in `Undergrad_Glucose uptake signaling/RT-PCR/` and `SBA/Figures/Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx`.

## 2. Gates

**Regression test: passed (56 checks, 0 failed, `test_output.txt`).** From the 17-gene matrix the module reproduces every number of `benchmark_route3/validate_pipeline.py`: 21 consensus pairs (any sign and same sign; the same 21 pairs as `pipeline_17genes.json`); pooled significant pairs 66/80/84/120; adjusted q < 0.05 0/51/11/15 (q matrices equal the statsmodels result to 1e-12); 13 consensus pairs with q < 0.05 in at least one scheme and 0 in all four; 91 eligible pairs; degrees Glp1r 6, Tcf7l2/Cacna1d/Ptbp1 5, Wnt5a/Cacna1c 4; null mean 22.10 (SD 3.03; manuscript 22.11, 3.06) and P(null >= 21) = 0.682 (manuscript 0.688) with 10000 permutations; covariance share 0.971. Also checked: vectorized primitives against explicit loops, statsmodels, scipy Spearman and Monte Carlo geometry; permutation engine properties; synthetic known answers (section 9).

**Gene-to-run map: passed (21 of 21 matrix columns resolved; none unresolved).** Section 3.

## 3. Item 2.1: verified gene-to-run map (`gene_run_map.csv`)

Method: every matrix column (12 values, sheet row order within group, no well-to-replicate mapping reconstructed) was compared with every candidate series of the three workbooks (53 blocks; 74 series versions because cached and recomputed values are both compared; 1554 comparisons in `gene_run_evidence.csv`), by value and not by label. Each gene has an exact match, and the matching blocks always carry the same gene name.

| Run (calculation sheet) | Genes of the 17 analyzed | Other matrix columns | Gapdh vectors on the sheet |
|---|---|---|---|
| 24 June 2020: sheet "24 June 2020"; the matrix equals cached column M of "Mix_3+24 June 2020" | Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7 (vector refA); Hes1, Glut2, Ins1, Glp1r (vector refC): 9 genes | Wnt2 (refA), Wnt2b (refB, own DAPT values), Hey1 (refC) | three |
| 26 August 2020: "Quantification Summary", sheet "0", column M | Tcf7l2, Cacna1c, Cacna1d, Kcnj11, Ptbp1: 5 genes | none | one |
| 9 November 2020: sheet "Box Plot", column N | Itpr1, Ryr2, Ryr3: 3 genes | Ryr1 | one |
| 3 June 2020 plates Gr1, Gr2 | none | none | not applicable |

* **Largest difference: 0.0.** The cached cells equal the matrix bit for bit for all 21 columns (largest absolute difference 0.0); recomputing 2^-(target Cq - reference Cq) from the typed Cq gives a relative difference of at most 1.4e-15 (floating point). No 3 June block comes close to any gene: for every gene the best 3 June block still differs by at least 99% in some replicate (`best_3June_plate_max_rel_diff` >= 0.993).
* **Table S1c lists Glp1r, Glut2, Hes1 and Ins1 under the 3 June plates only, and says which run supplied genes run more than once is not recorded. For the matrix it can now be stated:** all of them, and Lef1, Tcf7, Wnt5a, Wnt5b, Wnt9a, come from the "24 June 2020" calculation sheet. The nine June genes of the analysis therefore form one run of nine, not the runs of 5 + 4 assumed in the planning reconnaissance.
* **Formulas** (checked, not assumed): relative expression = 2^-(target Cq - reference Cq); 24 June: `M = 2^(-F)`, `F = D - E`, with D and E typed with two decimals; 26 Aug: `M = 2^(-K)`, `K = G - J`, J = the Gapdh rows (DAPT rows J30:J33 hold C10 and D10 twice, as already documented); 9 Nov: `N = 2^(-G)`, `G = E - F`. The summary columns S to U of the "Mix" sheet hold hand-typed values in some cells (for example U4, U10, T20, T148, T180) and are not used by the matrix. The 9 Nov sheet also holds copies of 3 June blocks and a second Ryr3 block with one Gapdh value per group; neither is used.
* **Gapdh reference vectors** (identical 12-value vectors, one value per replicate): 24 June refA (Gapdh Cq 18.9 to 20.1; used by Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7, Wnt2), refB (Wnt2b), refC (18.5 to 19.5; Hes1, Glut2, Ins1, Glp1r, Hey1); 26 Aug (22.5 to 23.8); 9 Nov (23.9 to 25.4). Genes sharing a vector share one Gapdh value per replicate.
* **Row order.** On the 26 Aug sheet the four rows of every target and group are in ascending Cq order (15 of 15 target blocks; the reference rows too in CTRL and DKK-1, not in DAPT, where C10 and D10 are entered twice). Replicate k of this run is therefore the k-th smallest Cq, not a well. The 24 June and 9 Nov rows are not sorted (0 or 1 of the 3 groups ascending per gene, as chance gives; `target_cq_ascending_groups` in the map).

## 4. Results (items 2.2 to 2.8; all post hoc)

### 2.2 Run-structure test (`tbl_2_2_run_structure.csv`)
D = mean within-group r over same-run pairs minus mean over between-run pairs; Pearson r on log-scale residuals after removing group means, raw scheme; run labels permuted across the 17 genes preserving run sizes (20000 permutations, seed 2026); two-sided P.

| Grouping | Same / between pairs | Mean r same | Mean r between | D | P | D for abs(r), P |
|---|---|---|---|---|---|---|
| Run = sheet (pre-specified, primary) | 49 / 87 | +0.149 | -0.019 | +0.168 | 0.020 | +0.102, 0.012 |
| Run = Gapdh vector (ADD-ON, secondary) | 29 / 107 | +0.385 | -0.052 | +0.436 | 0.00005 (the smallest possible with 20000 permutations) | +0.179, 0.00055 |

Per run (primary grouping; mean r inside the run, mean r between that run and the rest, D, P):

| Run | Genes | Pairs inside | Inside | Cross | D | P |
|---|---|---|---|---|---|---|
| 24 June | 9 | 36 | +0.014 | -0.032 | +0.047 | 0.55 |
| 26 August | 5 | 10 | +0.537 | -0.030 | +0.566 | 0.0006 |
| 9 November | 3 | 3 | +0.469 | +0.020 | +0.450 | 0.035 |

**Reading.** With the pre-specified definition (run = sheet) within-group correlation is higher among same-run pairs, but the excess is carried by the 26 August run (10 pairs) and the 9 November run (3 pairs); inside the 24 June run (36 pairs) there is none. The stronger result with run = Gapdh vector (ADD-ON) separates the 24 June genes by reference vector: refA pairs +0.333 (10 pairs, P = 0.006), refC pairs +0.176 (6 pairs, P = 0.055), whereas the 20 pairs between a refA gene and a refC gene, on the same sheet, average -0.194. This second grouping was not pre-specified and was chosen after the map showed the two vectors.

### 2.3 Consensus pairs and hubs by run (`tbl_2_3_*.csv`)
Of the 91 eligible pairs (14 non-denominator genes), 30 are same-run and 61 between-run. Consensus pairs: 9 of 30 same-run (30%) and 12 of 61 between-run (20%); Fisher odds ratio 1.75, P = 0.30 (one-sided 0.20); gene-label permutation (20000): observed 9 against a null mean of 6.9 (SD 1.9), P(null >= 9) = 0.19. **Consensus pairs are not enriched among same-run pairs.**
Hubs (degree >= 5, a post hoc rule): Glp1r (24 June, degree 6) and Tcf7l2, Cacna1d, Ptbp1 (26 Aug, 5 each). The 26 Aug run has 4 eligible genes, 3 of them hubs (1.1 expected; hypergeometric P = 0.041); mean degree 4.75 against 2.29 (24 June, 7 eligible genes) and 2.33 (9 Nov, 3 genes); label-permutation P for the 26 Aug mean degree 0.016; Kruskal-Wallis H = 5.0 for degree by run, label-permutation P = 0.023 (asymptotic 0.081). **15 of the 21 consensus pairs involve a 26 Aug gene** (10.6 expected; hypergeometric P = 0.026): 9 with a 24 June gene, 4 inside the run, 2 with a 9 Nov gene. These tests are descriptive and uncorrected.

### 2.4 Shared-reference contribution (`tbl_2_4_*.csv`)
Within-group r of same-run pairs (mean over pairs), on target Cq alone (no reference) and on dCq (the analysis scale); "expected" is the within-group r that two series sharing the observed reference residuals would show if their targets were independent (analytic).

| Pairs | n | r on target Cq | r on dCq | expected from shared reference alone |
|---|---|---|---|---|
| 24 June, same vector refA | 10 | +0.123 | +0.333 | +0.413 |
| 24 June, same vector refC | 6 | -0.164 | +0.176 | +0.297 |
| 24 June, refA x refC | 20 | +0.098 | -0.194 | +0.049 |
| 26 August | 10 | +0.777 | +0.536 | +0.343 |
| 9 November | 3 | +0.192 | +0.469 | +0.347 |
| All same-run | 49 | +0.215 | +0.149 | +0.232 |
| All between-run | 87 | -0.090 | -0.019 | +0.036 |

* In the 24 June run the reference raises the correlation of same-vector pairs (target Cq alone is near zero or negative), and the observed dCq correlations (+0.333, +0.176) are below what a shared reference alone would give (+0.413, +0.297). **The 24 June same-vector coupling is compatible with the shared Gapdh vector; it needs no co-regulation.** The 9 Nov run (3 pairs) is above the shared-reference value (+0.469 against +0.347), which three pairs cannot resolve.
* In the 26 August run the opposite holds: target Cq alone is already +0.777 (all 10 pairs have abs(r) >= 0.632) because the rows are in ascending Cq order, and the sorted reference partly cancels it (+0.536). **ADD-ON layout control:** independent genes simulated with the observed group means and SDs and the same ascending row order give a mean dCq correlation of +0.520 (95% interval 0.305 to 0.726; observed +0.536, P(simulated >= observed) = 0.45) and +0.741 on target Cq (observed +0.777). With random row order the same simulation gives +0.343 (0.032 to 0.632; P = 0.11). **The 26 Aug same-run coupling is reproduced by the row layout alone. Because each gene's four values per group are only a sorted list, the sheet carries no information on how replicates pair across genes for this run.**

Sensitivity versions of the pooled analysis (reference replaced as specified; 10000 permutations):

| Version | Consensus pairs (kept / lost / new vs reported) | Null mean, P(null >= obs) | Adjusted q < 0.05 (raw / Ins1 / Kcnj11 / Hes1) | Hubs (degree >= 5) | T_r2 (P) | D (P) |
|---|---|---|---|---|---|---|
| Reported (primary) | 21 | 22.1, 0.68 | 0 / 51 / 11 / 15 | Glp1r 6; Tcf7l2, Cacna1d, Ptbp1 5 | 0.133 (0.061) | +0.168 (0.020) |
| (a) reference = run-and-group mean | 27 (18 / 3 / 9) | 27.8, 0.66 | 2 / 17 / 27 / 18 | Glp1r 8; Tcf7l2 6; Ptbp1 6; Wnt5a, Ryr2, Cacna1c, Cacna1d 5 | 0.152 (0.001) | +0.305 (0.0015) |
| (b) reference = run-wide mean (no normalization) | 19 (8 / 13 / 11) | 19.4, 0.68 | 2 / 14 / 9 / 23 | Tcf7l2 6; Cacna1c 5 | 0.152 (0.001) | +0.305 (0.0015) |

Versions (a) and (b) have identical within-group residuals (both remove all per-replicate reference variation; a per-gene constant does not change r), hence the same T_r2 and D. **The consensus network is not stable to removing the per-replicate reference: only 8 to 18 of the 21 pairs remain and the hub list changes** (Glp1r's pairs with Cacna1d, Ptbp1 and Tcf7l2 are lost in (b)). Removing the reference raises, not lowers, the within-group coupling (T_r2 P = 0.001; two raw-scheme pairs reach q < 0.05, both inside the 26 Aug run: Cacna1c-Kcnj11 and Kcnj11-Ptbp1), again through the sorted rows: in these versions D is +0.853 for 26 Aug (P = 0.0004), +0.169 for 24 June (P = 0.039) and +0.271 for 9 Nov (P = 0.20).

### 2.5 Calibrated adjusted counts (`tbl_2_5_adjusted_calibrated.csv`, `tbl_2_5_q_per_pair.csv`)
Observed number of pairs with BH q < 0.05 after removing group means, against its distribution when within-group residuals are independent across genes (2000 permutations, seed 2026).

| Scheme | Observed | Null mean (SD, 95th pct) | Excess | P(null >= obs) | Status |
|---|---|---|---|---|---|
| raw | 0 | 0.07 (0.28, 1) | -0.07 | 1.00 | calibrated |
| Ins1-normalized | 51 | 29.3 (14.4, 54) | +21.7 | 0.079 | artefact-prone |
| Kcnj11-normalized | 11 | 11.3 (9.2, 30) | -0.3 | 0.44 | artefact-prone |
| Hes1-normalized | 15 | 6.5 (5.2, 17) | +8.5 | 0.078 | artefact-prone |

**The ratio-scheme counts (51, 11, 15) cannot be read as findings:** every ratio shares its denominator, so independent genes already give 29, 11 and 6.5 such pairs on average, and none of the three observed counts is above its own null at P < 0.05. In the calibrated raw scheme no pair reaches q < 0.05 (smallest q 0.114, Cacna1c-Kcnj11; per-pair q in `tbl_2_5_q_per_pair.csv`, ratio schemes marked as artefact-prone).

### 2.6 Within-group coupling and the run-block null (`tbl_2_6_coupling.csv`, 10000 permutations)

| Statistic | Observed | Independent-gene null: mean (95th pct), P | Run-block null: mean (95th pct), P |
|---|---|---|---|
| T_r2, all pairs | 0.133 | 0.114 (0.134), 0.061 | 0.132 (0.151), 0.42 |
| T_p (pairs with p < 0.05 at df 8) | 13 | 7.3 (12), 0.031 | 10.8 (15), 0.20 |
| T_eig | 5.06 | 4.39 (5.20), 0.081 | 5.01 (5.85), 0.43 |
| T_r2, same-run pairs (49) | 0.179 | 0.130 (0.166), 0.020 | invariant by construction |
| T_r2, between-run pairs (87) | 0.107 | 0.105 (0.128), 0.42 | 0.105 (0.135), 0.42 |
| Consensus count | 21 | 22.1 (27), P(null >= 21) = 0.68 | 23.6 (29), P = 0.85 |

Overall within-group coupling is borderline (T_r2 P = 0.061, T_p P = 0.031, T_eig P = 0.081). It is confined to same-run pairs, and the between-run pairs, which are the negative control, sit at the null. Once coupling inside runs is allowed (run-block null) nothing is left over: T_r2 is at its null mean and the observed consensus count (21; 12 of them between-run) is below the run-block null (mean 23.6, of which the between-run part is 14.6). By protocol section 7 this is read as compatible with shared reference, layout, cDNA or plate effects and not as co-regulation.

### 2.7 Geometric expectation (`tbl_2_7_geometry.csv`)
With 3 groups the centered group-mean vector of a gene lies in a plane; if pooled r is carried by the group means with uniform orientation, P(abs(r) >= 0.576) = (2/pi) arccos(0.576) = 0.609. Observed share of pooled-significant raw pairs: INS-1, pairs among the 15 genes with an ANOVA group effect (log scale, BH q < 0.05; the linear-scale ANOVA picks the same 15): 61 of 105 = **0.581**; all 17 genes: 66 of 136 = 0.485. Benchmark part 1, gene set A (17 DE genes, 37 series): median 0.794 (IQR 0.610 to 0.904; 28 of 37 above 0.609). INS-1 is near the anchor and below the benchmark median; the benchmark exceeds the anchor because group-mean vectors of DE genes are not uniformly oriented. An anchor, not a test.

### 2.8 Scale sensitivity (`tbl_2_8_scale_sensitivity.csv`, 10000 permutations)

| Scale | Pooled significant (raw / Ins1 / Kcnj11 / Hes1) | Consensus | Null mean (SD), P | Covariance share (median) | Adjusted q < 0.05 | Hubs (degree >= 5) |
|---|---|---|---|---|---|---|
| linear Pearson (manuscript) | 66 / 80 / 84 / 120 | 21 | 22.1 (3.0), 0.68 | 0.971 | 0 / 51 / 11 / 15 | Glp1r 6; Tcf7l2, Cacna1d, Ptbp1 5 |
| log-scale Pearson | 63 / 77 / 86 / 120 | 24 | 22.1 (2.9), 0.30 | 0.976 | 0 / 45 / 2 / 8 | Tcf7l2 8; Ptbp1 6; Wnt5a 5; Glp1r 5 |
| Spearman | 62 / 77 / 75 / 120 | 18 | 16.4 (3.0), 0.34 | 0.928 | 0 / 31 / 1 / 3 | none |

The statements "consensus count within the null" and "group means carry about 93 to 98% of the pooled covariance" hold on all three scales; the hub list does not (it differs on every scale).

## 5. What weakens the planning reconnaissance, and what does not

* **Weakened.** The reconnaissance (provisional map from Table S1c: 29 same-run pairs, mean +0.41 against -0.07 between-run, MWU P = 0.0007) suggested a general same-run effect. The provisional groups turn out to be the Gapdh-vector groups (its "3 June" genes are refC, its "24 June" genes refA), so its partition reproduces under the protocol statistic as D = +0.436 (listed as "provisional map" in `tbl_2_2_run_structure.csv`). With the verified runs (run = sheet) the effect is smaller (D = +0.168, P = 0.020), absent inside the 24 June run, and located in the 26 August and 9 November runs.
* **Not weakened.** Consensus pairs are not enriched among same-run pairs (reconnaissance Fisher P = 0.23; here 0.30). Ratio-scheme adjusted counts are not above their own null at P < 0.05 (reconnaissance null means 29/12/6, P 0.087/0.465/0.074; here 29/11/7, P 0.079/0.44/0.078). Hubs concentrate in the 26 Aug run, whose rows are in ascending Cq order; 15 of 21 consensus pairs involve one of its genes.
* **New.** For the 26 Aug run the layout alone reproduces the same-run correlation; removing the per-replicate reference makes the within-group coupling stronger in that run, not weaker; the consensus network and its hubs change when the reference is removed.

## 6. What the calculation sheets did and did not allow

Allowed: an exact value match for every matrix column; target and reference Cq per replicate for all 17 genes (so item 2.4 could be run in full); verification of every formula.
Not allowed or unknown:
* **Well-to-replicate pairing.** Not used and not reconstructed. The 26 Aug rows are in ascending Cq order within group, so replicate k is a rank there. For the other runs the row order is the only information; replicate k across runs has no verified physical meaning.
* **Raw source of the June values.** The sheet's 24 June Cq values (typed with two decimals; the same values are typed in `SeniorProject-2019/2020/Jun 24/2020-Jun-29_RT-PCR_Cal.xlsx`) were not found in any retained export (`sheet_vs_exports.csv`). Against the 2025 re-export "2020-Jun-24_Round2" (SBA/Rev2_250626/Export_qPCR), 20 of 21 gene x group blocks (Lef1, Tcf7, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a) differ from the sheet by +2.5 to +2.9 cycles (median +2.6); only the Wnt2b DAPT block (an excluded gene) is equal. That export contains no Hes1, Hey1, Ins1, Glp1r or Glut2; against the union of the two 3 June re-exports (Gr1, Gr2: two wells per group each) none of the 36 gene x group blocks of the sheet has an equal value (median difference +1.8 to +3.6 cycles per gene). Native `.pcrd` files of Round 2 exist in that senior-project folder (two copies); they were not opened. **This does not by itself show that a raw file is missing:** the 2025 re-exports also differ from the sheets elsewhere (26 Aug: 9 of 15 blocks differ by up to 0.8 cycles from the sheet, although the export kept in the RT-PCR folder equals the sheet in 60 of 60 values; 9 Nov: re-export and folder export both equal the sheet, 44 of 44), so analysis settings can move Cq. The June offset is larger and systematic, and the files do not say whether it is a setting or a different run. Table S1c quotes the Gapdh range of the Round 2 re-export (16.7 to 17.5) for 24 June, whereas the sheet that fed the matrix has 18.5 to 20.1 for the genes of the analysis.
* **Why two Gapdh vectors exist on the 24 June sheet** (refA and refC differ by 0.2 to 0.8 cycles per replicate) is not recorded.
* **Cross-run independence.** Whether the runs used different cDNA aliquots is not recorded. Between-run pairs are used as the negative control by design, but their independence cannot be verified from the sheets.

## 7. Limits

One preparation, 12 observations, 17 genes; all P values are uncorrected across the many descriptive tests listed; the hub rule is post hoc; simulations of the layout use the observed group means and SDs of four values per group.

## 8. Additions and deviations (also logged in the protocol, section 8)

* ADD-ON, not in the protocol: run = Gapdh vector as a second grouping (2.2); the analytic shared-reference expectation and the 26 Aug layout simulation (2.4); consensus pairs involving a run, degree by run tests (2.3); the provenance check against CFX exports (section 6).
* The protocol does not state the number of permutations for 2.6, 2.8 and the sensitivity versions: 10000 were used. Items 2.2 and 2.3 use 20000 and 2.5 uses 2000, as written. The ANOVA of 2.7 is on the log scale (linear scale gives the same 15 genes).
* No protocol rule, threshold or definition was changed after a result was seen; the verified map replaced the provisional one as the protocol required. Edited outside this folder: the module, its test, section 8 of the protocol and the fingerprint file (appended), and `PREP_PLAN_CMBL.md` (one line appended). No commit or push.

## 9. Using `coexp_checks.py` from other work packages

```python
import sys; sys.path.insert(0, r'.../Analysis_outputs/strengthen_v2')      # python -I removes the script folder
import coexp_checks as cc
res = cc.run_checks(X, groups, genes=genes, den=['Ins1','Kcnj11','Hes1'], runs=None, nperm=2000)   # X samples x genes, linear scale
cc.coupling_test(X, groups, nperm=500)        # T_r2, T_p, T_eig: observed, null, p_ge_plus1
cc.adjusted_null(X, groups, den, nperm=500)   # calibrated BH counts per scheme
Xd = cc.decouple(X, groups, rng); Xi = cc.inject_factor(Xd, groups, range(6), 0.7, rng)    # protocol 3.2
```
* Observed values are identical to `bench_pipeline.py`. The random streams differ (vectorized permutations), so null means and P agree with the part 1 loop within Monte Carlo error but not draw for draw (10000 vectorized permutations: 22.10 and 0.682; 300 loop permutations: 21.90). A gate that needs exact reproduction of part 1 null means must call `bench_pipeline.analyse` for that comparison.
* Speed at G = 17, N = 12 (3 x 4), one process, other jobs running on the machine: T_r2 + T_p 5100 permutations/s, with T_eig 3200/s (`eig=False` skips it); consensus null 2500/s (the loop does 48/s: 53x); adjusted counts over 4 schemes 1800/s; Spearman consensus 420/s; run-label permutations 5600/s. Memory is bounded by `MAX_ELEMS` chunking and results do not depend on it.
* Use `p_ge_plus1` for test P values (`p_ge` is the plain fraction, as in the manuscript's P = 0.69). Pass a prebuilt `cc.Design(groups)` to avoid rebuilding it in a loop.
* Synthetic checks in the test file (250 independent 3 x 4 data sets, 199 permutations): T_r2 type I error 0.052 (Wilson 0.031 to 0.087), P approximately uniform (KS P = 0.93); T_p 0.028 and the consensus count 0.052 with mean P 0.56 and 0.60 (discrete and conservative); injected factor, r = 0.7 in 6 genes: detected by T_r2 in 92% (T_p 95%, T_eig 92%) and by the consensus count in 20%; r = 0.9 in 10 genes: 100% for all; r = 0.5 in 3 genes: about 5% (no power). With noisy shared denominators and independent genes the ratio-scheme adjusted count is 45 of 120 pairs against 0.8 for the raw scheme, and the calibrated null reproduces it (share of P <= 0.05: 0.04). These are simulations on synthetic data, not the protocol 3.2 calibration on the development series (WP2).

## 10. For the authors

1. Which plate or round produced the second block of the "24 June 2020" sheet (Hes1, Hey1, Ins1, Glp1r, Glut2, Gapdh vector refC), which has no counterpart in the retained Round 2 export, and whether the Cq values of the first block (Lef1 to Wnt9a, refA), which are 2.6 cycles higher than in the 2025 Round 2 re-export, come from other analysis settings of that run or from another run; an original export or the `.pcrd` re-analyzed with the 2020 settings would settle it. Table S1c should then name the run that supplied Glp1r, Glut2, Hes1 and Ins1 (the matrix uses the 24 June sheet, not the 3 June plates).
2. Whether the ascending-Cq row order of the 26 August summary is a CFX export setting or a manual sort, and whether the plate layout file allows a documented pairing of target and Gapdh wells. Nothing here depends on it, and no pairing was reconstructed.
3. Why the 24 June sheet carries two different Gapdh vectors for the same samples.
