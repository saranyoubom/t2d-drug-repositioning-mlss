# WP4: public islet cross-check, version 2 (6 Nov 2026)

Protocol v2 section 4 (`../PROTOCOL_v2_prespecified.md`); brief `../briefs/WP4_islets_v2_brief.md`; choices that the protocol leaves open were fixed before the new analyses in `PRESPEC_WP4_choices.md` (fingerprints in `../PROTOCOL_v2_SHA256.txt`). Executed in the main session. Every table below is generated from the result files by `make_tables.py` (`WP4_summary_tables.md`, same numbering). Nothing here changes the manuscript, the JTM package, Fig. 1-6 or the earlier islet folders (`../../human_islet_GSE81608/`, `../../islet_GSE84133/`); no commit or push.

**What kind of evidence this is.** The cross-check asks whether the 21 consensus pairs of the INS-1 experiment are co-expressed from cell to cell in untreated human or mouse beta-cells. It cannot test the treatment response. The earlier GSE81608 and GSE84133 results stay as reported (gate passed); everything else is a pre-specified sensitivity layer on public data, except the add-ons that are labelled POST HOC (the gene-graph null, section 5, and the items of section 11).

## 1. Summary
1. **Gate (4.1): PASS.** Both earlier analyses were re-run on the cached data and every stored file and all 49 quoted numbers were reproduced exactly; the new module `wp4_core.py` reproduces them too (section 2).
2. **Primary endpoint (4.4).** Three further human datasets were obtained with author-provided beta-cell labels (GSE86469, GSE81547, E-MTAB-5061; GSE85241 has no per-cell labels and was skipped). With the 5 human datasets, the pooled excess co-expression of the consensus pairs over expression-matched random pairs, depth-adjusted, is **E = +0.013 (95% CI -0.052 to +0.078; REML with Hartung-Knapp, k = 5)**. The interval lies inside the pre-set margin of +/-0.10, so the fixed rule gives: **no meaningful excess co-expression**. This is a statement about the 11 of 21 consensus pairs that can be evaluated; the other 10 (all 7 Wnt pairs, i.e. every pair that contains WNT5A, WNT5B or WNT9A, and the 3 RYR3 pairs) cannot be evaluated in any human dataset because a gene is detected in fewer than 10% of beta-cells (LEF1 is in no consensus pair).
3. **How firm.** Heterogeneity is high (I2 71%); leaving out one dataset gives "inconclusive" (an interval beyond +/-0.10) for 3 of the 5 leave-one-out runs; the non-diabetic-only and not-depth-adjusted pooled estimates, and a post hoc null that keeps the shared genes of the consensus pairs (+0.012 [-0.050, +0.074]), give the same verdict as the primary (section 5).
4. **Depth adjustment (4.2).** Moves the pooled E from -0.001 to +0.013; per dataset E rises in 5 of 6 datasets (by 0.007 to 0.047) and falls in GSE84133 human (+0.262 to +0.192, by 0.070) (section 4).
5. **Metacells (4.5).** Aggregating random groups of 20 cells of a donor (donors with at least 40 beta-cells) makes 1 to 4 more consensus pairs evaluable (GSE84133 human 3 to 7 at size 20 and 11 at size 50; GSE81608 11 to 15; GSE86469 11 to 14; GSE81547 9 to 11; E-MTAB-5061 11 to 12; mouse 5 to 7), but the metacell correlations are small and not above the matched null: in the two datasets with many metacells the median consensus rho is +0.03 (GSE84133 human, 124 metacells: +0.032 against +0.239 at the cell level; E +0.016, range over 20 partitions -0.060 to +0.055); in the others 5 to 13 metacells from 2 to 4 donors make the correlations too imprecise to read (section 7).
6. **Perturbation concordance (4.6).** 43 candidate series from the two fixed GEO queries; one qualified (GSE60803, beta-catenin shRNA in human beta-cell-derived cells; no Notch-inhibition series qualified). Its sign agreement with the INS-1 DKK-1 arm is 7 of 17 genes (exact binomial P 0.63, Spearman rho -0.04): no concordance, but the design is confounded (controls are not expanded cells) and the model systems differ (section 8).

## 2. Gate (deliverable A, protocol 4.1): PASS
`gate_output.txt` holds the full text. Stage 1: `../../islet_GSE84133/islet_check_baron.py` was run unchanged on the cached parsed data, and `../../human_islet_GSE81608/human_islet_check.py` through `gate_human_islet_check_lean.py`, a copy in which only the loading block differs (the original reads the whole 39851 x 1600 RPKM table as a dense float64 frame, about 510 MB, three to four copies at peak, too much for this machine; the lean loader keeps the 503 beta-cell columns in chunks of 2500 rows, every statistic line is verbatim). All 21 output files are identical to the stored results (max absolute difference 0) and all 49 numbers quoted in the earlier READMEs and in the brief are reproduced (for example GSE81608 11 evaluable pairs, median rho -0.021, empirical P 0.985; GSE84133 human 3 pairs, 0.24, P 0.071; mouse 5 pairs, -0.04, P 0.805; non-diabetic donors only, 5% rule: 7 pairs, 0.28, P 0.09). Stage 2: the new module `wp4_core.py` in unadjusted mode reproduces rho, p, q of all 136 pairs, the detection table, the pair-level nulls and the summaries of GSE81608, GSE84133 human, mouse and the two non-diabetic variants (tolerance 1e-9, empirical P exactly).

## 3. Datasets obtained and not obtained (deliverable C, 4.3; `attempts.csv`, `frozen_islet_datasets.csv`, `frozen_beta_cells.csv`)
The dataset table and the beta-cell lists were frozen from metadata only (GEO series matrices, ArrayExpress SDRF) on 6 Nov at 00:17, before any expression value of these datasets was read (sha256 in the fingerprint file). The label is the cell-type annotation written by the data authors in the public record; no label was derived from expression.

| order | dataset | source and size | data | author label | beta-cells | donors | outcome |
|---|---|---|---|---|---|---|---|
| 1 | GSE86469 (Lawlor 2017) | GEO, RSEM expected counts, 12.7 MB | counts | "cell type" = Beta | 264 | 8 (3 with T2D, 96 cells) | obtained |
| 2 | GSE81547 (Enge 2017) | GEO, RAW.tar of 2544 per-cell files, 176 MB | htseq counts | "inferred_cell_type" = beta | 348 | 8 (age and sex define the donor; disease not stated, healthy donors per the series summary) | obtained |
| 3 | E-MTAB-5061 (Segerstolpe 2016) | EBI, processed.1.zip 159 MB; the 723 MB text was streamed from the zip, never written out | Smart-seq2 counts | SDRF "inferred cell type" = beta cell | 270 | 10 (4 with T2D, 99 cells) | obtained |
| 4 | GSE85241 (Muraro 2016) | GEO, counts csv.gz 17 MB | - | none (32 library-level GEO samples, no per-cell label) | 0 | - | skipped before download: the protocol requires author-provided labels and we do not assign cell types |

Gene labels: GSE86469 uses Ensembl ids (key `ensembl_ids_17genes.csv`, from the NCBI Datasets gene report); the other two use symbols. E-MTAB-5061 passed the cell-identity check (INS read count 51329 in beta-cells against at most 382 in every single-hormone type; GCG 95994 in alpha-cells against at most 760; `logs/identity_check_E-MTAB-5061.json`). Publications attached to the records (`dataset_citations.csv`, from PubMed via the records; to be verified in the revision before citing): Xin 2016 (GSE81608), Baron 2016 (GSE84133), Lawlor 2017 (GSE86469), Enge 2017 (GSE81547), Segerstolpe 2016 (E-MTAB-5061), Lenz 2014 (GSE60803).

## 4. Depth adjustment (deliverable B, 4.2)
Per gene, the donor-centred log expression was regressed on the donor-centred depth and its square (log10 library size for UMI and read counts, log10 number of detected genes for RPKM), the residuals were ranked and correlated (partial Spearman, df = cells - donors - 4); the matched null (same detection deciles, 2000 draws, seeds 2026 and 2027) was rebuilt on the adjusted ranks with the same draws. Table 1 gives both analyses side by side for each dataset (`results/`; `WP4_summary_tables.md`). Key numbers (median consensus rho; E):

| dataset | not adjusted | depth-adjusted |
|---|---|---|
| GSE81608 | -0.021; -0.040 | -0.023; -0.023 |
| GSE84133 human | +0.239; +0.262 | +0.220; +0.192 |
| GSE86469 | +0.016; -0.027 | +0.024; -0.020 |
| GSE81547 | +0.109; +0.016 | +0.116; +0.025 |
| E-MTAB-5061 | +0.063; +0.035 | +0.066; +0.045 |
| GSE84133 mouse | -0.045; +0.002 | -0.011; +0.049 |

Adjustment changes the consensus-pair medians by at most 0.034 and E by at most 0.070; it changes the matched-null level more (GSE84133 human 0.052 to 0.091; mouse 0.078 to 0.029). The empirical P of the GSE84133 human median goes from 0.071 to 0.023 (E: 0.016 to 0.002).

## 5. Primary endpoint: meta-analysis (deliverable D, 4.4; `meta/`)
E (excess) = mean Fisher z of the evaluable consensus pairs minus the mean over 2000 matched random pair sets; SE = SD of those 2000 means (as pre-specified: it describes the variability of random pair sets, not the finite number of cells). Datasets qualify with at least 3 evaluable consensus pairs; all 5 human datasets qualified. Primary analysis: depth-adjusted, all frozen beta-cells (the donor with type 2 diabetes of GSE84133 stays in), tau2 by REML, Hartung-Knapp standard error with the t quantile on k - 1 df, margin +/-0.10 on the E scale (Fisher z; atanh(0.10) = 0.1003).

| dataset | beta-cells | evaluable pairs | E | SE of E |
|---|---|---|---|---|
| GSE81608 | 503 | 11 | -0.023 | 0.015 |
| GSE84133 human | 2525 | 3 | +0.192 | 0.087 |
| GSE86469 | 264 | 11 | -0.020 | 0.022 |
| GSE81547 | 348 | 9 | +0.025 | 0.022 |
| E-MTAB-5061 | 270 | 11 | +0.045 | 0.022 |
| **pooled (REML + Hartung-Knapp)** | | | **+0.013 [-0.052, +0.078]** | tau2 0.0011, I2 71% |

**Decision-rule outcome: the pooled 95% CI lies inside +/-0.10, so "no meaningful excess co-expression".** Sensitivity, never used for the rule (Table 3): the same verdict with the Knapp-Hartung factor floored at 1 (identical interval), DerSimonian-Laird with Hartung-Knapp [-0.053, +0.080], DerSimonian-Laird with the normal interval [-0.025, +0.053], fixed effect +0.003 [-0.016, +0.022]; not depth-adjusted -0.001 [-0.065, +0.063]; non-diabetic donors only (4 datasets, GSE84133 has no evaluable pair without its diabetic donor) +0.003 [-0.085, +0.090]; not adjusted and non-diabetic only -0.004 [-0.075, +0.067]. **Leave-one-out (Table 3b), primary:** without GSE81608 +0.026 [-0.062, +0.115] (inconclusive), without GSE84133 human +0.005 [-0.049, +0.059] (no excess), without GSE86469 +0.026 [-0.071, +0.123] (inconclusive), without GSE81547 +0.015 [-0.095, +0.125] (inconclusive), without E-MTAB-5061 +0.000 [-0.082, +0.083] (no excess). Four- and three-study Hartung-Knapp intervals are wide because of the t quantile (3.18 and 4.30), which is how the rule behaves with few datasets.

**Add-on, POST HOC (`gene_graph_null.py`, logged in protocol section 8 at 14:19 before its result was seen; `meta/gene_graph_*`).** The pre-specified null draws the two genes of each pair independently, but the evaluable consensus pairs share genes (up to 4 pairs per gene), so the SE of E might be too small. With a null that replaces every distinct gene of the evaluable pairs by one random gene of its own detection decile and keeps the pair graph (2000 draws, seed 2028), the SEs change by factors 1.05 (GSE81608), 1.06 (GSE84133 human), 1.07 (GSE86469), 1.21 (GSE81547), 0.98 (E-MTAB-5061) and 1.09 (mouse); the pooled estimate is **+0.012 [-0.050, +0.074]** (REML + Hartung-Knapp), the same verdict, and the five leave-one-out verdicts are the same as above (3 inconclusive). The rule stays defined on the pre-specified SE.

Single datasets (empirical P of E in the depth-adjusted analysis): GSE84133 human 0.002 and E-MTAB-5061 0.017 are below 0.05; GSE81547 0.126; GSE86469 0.819; GSE81608 0.943; mouse 0.176 (mouse is analysed alone: E +0.049, SE 0.052, 5 pairs).

## 6. Pair-level pooled rho (4.4; Table 4)
Eleven consensus pairs are evaluable in at least 3 human datasets (depth-adjusted rho, Fisher z, variance 1.06/(df - 1), REML + Hartung-Knapp). Pooled rho ranges from -0.032 (TCF7L2-ITPR1) to +0.163 (CACNA1C-CACNA1D, +0.163 [-0.107, +0.410], I2 98%: GSE84133 human 0.48 against 0.02 to 0.19 in the others). Only CACNA1D-PTBP1 has an interval that excludes zero (+0.095 [+0.001, +0.187]), without correction for 11 pairs. The pair behind the MLSS Synergy Bonus, CACNA1D-GLP1R, is evaluable in 4 datasets: +0.051 [-0.073, +0.173]. Ten pairs are never evaluable (Table 2b): WNT5A with CACNA1C, CACNA1D, PTBP1 and GLP1R; WNT5B with TCF7L2 and GLP1R; WNT9A-TCF7; TCF7L2-RYR3; RYR2-RYR3; RYR3-GLP1R.

## 7. Metacell sensitivity (deliverable E, 4.5; `metacell/`)
As written in the choices note (section C): within each donor with at least 40 beta-cells the cells are shuffled and cut into blocks of exactly 20 (sensitivity 10 and 50), the remainder is dropped, a donor needs at least 2 blocks; counts are summed (CP10K or CPM) or RPKM averaged, donor-centred, ranked, Spearman, no depth term; a gene is evaluable when nonzero in at least 50% of the metacells; the matched null is rebuilt on the metacell data; 20 random partitions per size (`default_rng(5000 + i)`). Table 5 of `WP4_summary_tables.md` has the full numbers (median and range over the partitions). Median over the 20 partitions; "pairs" = evaluable consensus pairs, cell level in brackets; rho = median consensus rho; E as in section 5:

| dataset | donors with at least 40 cells | size 20: metacells; pairs; rho; E | size 10: metacells; pairs; rho; E | size 50: metacells; pairs; rho; E |
|---|---|---|---|---|
| GSE81608 | 3 | 13; 15 (11); +0.022; +0.012 | 27; 11; -0.004; -0.005 | 2 metacells, 1 donor: not evaluable |
| GSE84133 human | 4 | 124; 7 (3); +0.032; +0.016 | 251; 7; +0.026; +0.018 | 48; 11; +0.022; +0.019 |
| GSE86469 | 3 | 7; 14 (11); +0.023; +0.053 | 17; 11; +0.020; -0.038 | none: not evaluable |
| GSE81547 | 4 | 11; 11 (9); +0.023; -0.037 | 22; 9; +0.028; -0.003 | none: not evaluable |
| E-MTAB-5061 | 2 | 5; 12 (11); +0.175; +0.105 | 10; 11; +0.097; +0.062 | none: not evaluable |
| GSE84133 mouse | 2 | 44; 7 (5); +0.031; +0.011 | 89; 5; +0.070; +0.055 | 17; 7; +0.029; -0.004 |

* **How many more pairs become evaluable.** At size 20, 1 to 4 more (GSE81608 +4, GSE84133 human +4, GSE86469 +3, GSE81547 +2, mouse +2, E-MTAB-5061 +1); at size 50 in GSE84133 human 11 pairs against 3 at the cell level. Genes that become evaluable (nonzero in at least half of the metacells in at least half of the 20 partitions; `metacell/metacell_gene_evaluability.csv`, recomputed with `metacell_genes.py` and identical in count to every row of `metacell_runs.csv`): GSE81608 WNT5A (18 of 20 partitions); GSE84133 human WNT9A, TCF7L2 and GLP1R at size 20 (and ITPR1, RYR2, SLC2A2 at size 50); GSE86469 RYR3; GSE81547 ITPR1; E-MTAB-5061 WNT9A and TCF7; mouse Hes1, Tcf7l2 and Itpr1. The four WNT5A pairs are therefore evaluable only in GSE81608 at size 20 (13 metacells from 3 donors) and WNT9A-TCF7 only in E-MTAB-5061 at size 20 (5 metacells from 2 donors); the two WNT5B pairs are evaluable in no dataset at any size (`metacell/metacell_pair_evaluability.csv`).
* **What the correlations show.** Where many metacells exist (GSE84133 human 124 to 251 at sizes 20 and 10, 48 at size 50; mouse 44 to 89 and 17), the median consensus rho is +0.02 to +0.03 (mouse size 10 +0.070) and E has a median of +0.02 (GSE84133 human; -0.06 to +0.10 across the partitions of all three sizes) and between -0.004 and +0.055 (mouse; -0.15 to +0.19 across partitions, 17 to 89 metacells), i.e. no excess over the matched null; the three human pairs that stood out at the cell level (adjusted median +0.220) shrink to +0.03 when random cells are aggregated. This is the pattern expected if part of the cell-level correlation came from per-cell detection effects, but the analysis cannot separate that from a loss of signal through averaging random cells (random metacells keep no cell-state structure). In GSE81608, GSE86469, GSE81547 and E-MTAB-5061 only 2 to 4 donors have at least 40 beta-cells, so 5 to 27 metacells enter and the median rho varies from -0.65 to +0.60 across partitions at size 20 (E-MTAB-5061: 5 metacells): these runs are not informative about co-expression, and they are shown because the protocol asked for them.
* Size 50 can be applied only to GSE84133 (human and mouse); in the other datasets no donor has 100 beta-cells, so no donor yields 2 blocks.

## 8. Perturbation concordance (deliverable F, 4.6; `perturbation/`)
Two fixed GEO DataSets queries (`perturbation_scout.py`; terms in `PRESPEC_WP4_choices.md` section D; `perturbation/queries.json`): Notch 19 hits, Wnt 27 hits, 43 unique series (3 in both). Every candidate is listed with its decision and reason in `perturbation/candidates_decisions.csv`, made from title, organism, assay type and summary only. Dropped 42: 17 stem-cell, progenitor or reprogramming designs; 10 decided from the title (other tissues, reporter, SAGE and microRNA assays, isogenic hESC); 10 read in full (other perturbations, the opposite direction or other models, for example the Notch gain-of-function series GSE193888, a mirror design taken up again in section 15, and GSE44045 and GSE44046, which are diet comparisons although the series titles refer to Sfrp5: an earlier version of this README called them Sfrp5 silencing, withdrawn in addendum H of `PRESPEC_WP4_choices.md`); 5 that are not expression profiles (genome binding and other assays). Kept 1: **GSE60803** (Lenz 2014; adult human beta-cell-derived cells, beta-catenin shRNA; Affymetrix HuGene 1.1 ST; 4 donors, 12 arrays). No Notch (gamma-secretase) inhibition design qualified, so the DAPT arm could not be compared.
Design and contrast were fixed before any value was read (choices note, addendum G): per donor, control = untreated array (S108, S109, S112) or nontarget shRNA array (S114); perturbed = mean of the donor's beta-catenin shRNA arrays; log2 fold change = mean over the 4 donors; probe set with the highest control expression per gene (Entrez id match; all 17 genes have a probe set). **The contrast is confounded**: the records describe the untreated and nontarget arrays as intact islet cells and the shRNA arrays as expanded, dedifferentiated cells, and only one donor has a nontarget control.

| contrast against the INS-1 DKK-1 arm (RT-qPCR, as reported) | genes | same sign | exact binomial P | Spearman rho | P of rho |
|---|---|---|---|---|---|
| primary: 4 donors, perturbed minus control | 17 | 7 | 0.63 | -0.037 | 0.89 |
| secondary, descriptive: donor S114 only (control n = 1) | 17 | 10 | 0.63 | -0.007 | 0.98 |

Same sign (primary): WNT5A, WNT5B, TCF7L2, KCNJ11, SLC2A2, INS, GLP1R. Opposite: HES1, WNT9A, LEF1, TCF7, ITPR1, RYR2, RYR3, CACNA1C, CACNA1D, PTBP1. The fold changes in GSE60803 are small (|log2 fold change| at most 0.28; `perturbation/gse60803_fold_changes.csv`). The result neither supports nor refutes the INS-1 DKK-1 changes: model systems differ (rat cell line against human cells), the contrast mixes the knockdown with expansion, and 17 genes give little power.

## 9. Figures (deliverable G; `make_fig7_v2.py`, `figures/`)
`figures/Figure_7_v2_islet_crosscheck.png/.pdf` (170 x 190 mm, 300 dpi): (a) detection of the 17 genes in the beta-cells of six datasets (numbers = percent; grey = below the 10% rule); (b) the 21 consensus pairs, depth-adjusted rho (dot) and matched random pairs (grey bar = null mean to 95th percentile), NE = not evaluable; panel titles carry the dataset and the number of evaluable pairs. `figures/Figure_8_draft_islet_meta.png/.pdf` (170 x 150 mm): (a) forest plot of E per dataset (square, +/-1.96 SE) and pooled estimates (diamond, Hartung-Knapp 95% CI) with the +/-0.10 margin shaded; (b) pair-level pooled rho (dots = datasets, diamond = pooled); (c) GSE60803 log2 fold changes against the INS-1 DKK-1 arm. No key strips; the symbols above go into the legends. Both are drafts; they plot saved results only.

## 10. Limits that must travel with the numbers
* Cell-to-cell co-variation in untreated beta-cells is not the response to DAPT or DKK-1; "no excess co-expression" does not contradict the INS-1 finding and does not test treatment-induced co-regulation.
* Only 11 of 21 consensus pairs are evaluable at the cell level in any dataset; the 7 Wnt pairs (the DKK-1 arm's pathway) are not (at the metacell level five of them, four WNT5A pairs and WNT9A-TCF7, become evaluable in one dataset each, from 5 to 13 metacells). Evaluable genes are the highly expressed ones; dropout attenuates all correlations (the matched null shares the attenuation through the detection deciles).
* The SE of E is the SD of random pair sets (pre-specified), not a sampling error of the cells; the consensus pairs share genes, which the post hoc gene-graph null shows to matter little (SE factors 0.98 to 1.21); with 3 to 5 datasets the Hartung-Knapp interval is wide and the verdict changes with the dataset left out (section 5).
* One dataset (GSE84133 human) contributes 3 pairs (CACNA1C, CACNA1D, PTBP1 among themselves) and a very different E (+0.19); the other four give -0.02 to +0.05.
* Donor identity: GSE81547 has no donor id in the record; donors are age-and-sex combinations (8), so two donors sharing both would be merged. E-MTAB-5061 and GSE86469 donors are record ids. Donors with type 2 diabetes stay in the primary analysis; the non-diabetic-only analysis is a sensitivity analysis.
* GSE60803 is confounded (section 8); no Notch series qualified.
* All of this is post hoc relative to the INS-1 data (one preparation) and descriptive for the mouse data.

## 11. Deviations, additions and code fixes (all logged in `../PROTOCOL_v2_prespecified.md` section 8)
Times are the machine clock (+0700); the clock jumped from 00:17 to 07:4x between two parts of the work.
1. 6 Nov 00:17: dataset table and beta-cell lists frozen from metadata; GSE85241 skipped (no per-cell labels).
2. 6 Nov 07:50-08:00: gate (section 2); GSE81608 through a lean copy of the loading block (statistics verbatim).
3. 6 Nov 13:40: `PRESPEC_WP4_choices.md` fixed (data layer per dataset, depth covariates, pooling estimators, margin scale, metacell rules, perturbation queries and rules) before any analysis of the new datasets; fingerprinted.
4. 6 Nov 13:43, addendum F: E-MTAB-5061 library size and detected genes over the 26179 RefSeq rows only (the IDF documents ERCC and eGFP rows); the table of the choices note had said "all rows". Fixed before any count was read.
5. 6 Nov about 13:45, addendum F2: the first parser run stopped on the identity check as first written (own type against ALL other types): the authors' "co-expression cell" and "unclassified endocrine cell" categories carry INS by definition. Criterion corrected to single-hormone types; no analysis statistic of E-MTAB-5061 existed. The first result file is kept (`logs/identity_check_E-MTAB-5061_first_version_all_types.json`).
6. 6 Nov 13:51, addendum G: GSE60803 contrast and its confounding fixed before any expression value was read; the S114-only contrast (control n = 1, below the 3-sample rule) is reported as secondary and descriptive.
7. Add-ons fixed in the choices note before the results (not in the protocol text): non-diabetic-only sensitivity for every dataset with donors who have type 2 diabetes (the protocol had it for GSE84133 only); alternative pooling estimators; leave-one-out; pair-level pooling variance with the Fieller factor 1.06.
8. 6 Nov about 14:05, metacell run: the first run (slow sampler) was stopped after GSE81608, GSE84133 human and part of GSE86469 and restarted with `fast_draws`, which draws `pool[rng.integers(0, n)]` instead of `rng.choice(pool)`. The two give identical values and leave the generator in the same state (tested for pool sizes 1 to 26000; `_selftest` in `metacell.py` checks it on every start); the restarted run reproduced the earlier GSE81608 and GSE84133 human output bit for bit (120 rows, max difference 0; `metacell/metacell_runs_slow_sampler_first_two_datasets.csv` keeps the earlier file). No effect on results.
8b. 6 Nov 14:19: post hoc add-on `gene_graph_null.py` (section 5), logged before its result was seen.
9. Code fixes without effect on results: a literal "%" in a report string (`meta_analysis.py`), the separator of Entrez ids in the GPL11532 annotation (`perturbation_gse60803.py`), the label margin of Fig. 8a, and `fetch_ensembl_ids.py` (an efetch XML request exceeded the 20 MB response limit of `wp4_http.get_text`; the NCBI Datasets gene report was used instead).
10. Python: the analysis stack is in Python 3.14 (`python`); the default `python` (3.13) has no numpy. Run with `-I`.

## 12. Numbers that can replace those in the manuscript's islet text (for WP6b)
* Methods subsection "Public islet single-cell data" and the aims: five human datasets and one mouse dataset instead of two; add the three datasets with their records and citations (section 3), the depth adjustment (section 4), the meta-analysis with the fixed margin, the metacell sensitivity and the perturbation search; state that the protocol, the dataset lists and the choices note were fingerprinted before the analyses (`../PROTOCOL_v2_SHA256.txt`).
* Results subsection (now titled "Public islet single-cell data neither confirm nor exclude the consensus pairs"): the GSE81608 and GSE84133 numbers stay as written (reproduced exactly); add Table 1 values; replace the closing paragraph ("The three primary analyses therefore do not agree ... neither confirm nor exclude") by the decision-rule sentence of section 1 with its limits (11 of 21 pairs evaluable, none of the Wnt pairs; I2 71%; leave-one-out in 3 of 5 "inconclusive"). Suggested title: "Public islet single-cell data show no excess co-expression of the evaluable consensus pairs".
* Abstract (Results sentence "could evaluate only 3 to 11 of the 21 consensus pairs, with inconsistent results"): "In five public human islet datasets, the 11 evaluable consensus pairs were not more co-expressed than expression-matched random pairs (pooled excess +0.01, 95% CI -0.05 to +0.08)"; keep the word count under 350.
* Fig. 7: replace by `figures/Figure_7_v2_islet_crosscheck.*` (legend text from section 9) and add Fig. 8 or move panel (a) of Fig. 8 into the figure; Tables S9 to S9c and the Additional file 4 sheets `Islet_*` need the new datasets (the per-pair tables are `results/*__cons.csv`, `results/*__pairs.csv`, `results/*__detection.csv`).
* Discussion, Limitations, Priorities, Conclusions: "11 of 21 in GSE81608, 3 in the human and 5 in the mouse GSE84133 data" becomes the evaluability per dataset (Table 1: 11, 3, 11, 9, 11 and 5 pairs); the Conclusions sentence "neither confirmed nor excluded" becomes "found no excess co-expression of the evaluable consensus pairs, which cannot test the treatment response".
* Declarations: the ethics statement and the data availability statement list only GSE81608 and GSE84133; add GSE86469, GSE81547, E-MTAB-5061 (ArrayExpress) and GSE60803. References: the six publications of section 3 (verify on PubMed before citing).

## 13. For the authors
* Whether the islet result is to be the headline in the abstract and how far to go with "no meaningful excess" given the leave-one-out results.
* The margin of +/-0.10 was set before the data (rho = 0.10 explains 1% of variance); a different margin changes the wording, not the numbers.
* Whether the confounded GSE60803 comparison should appear in the main text or only in an additional file.

## 14. Files and how to run
Scratch data root: `V2_DATA_ROOT` (default: the session scratch folder; `v2_data/wp4/` holds the downloads and `parsed/`). Python 3.14, run with `-I`; one heavy job at a time.

| step | script | writes | run time here (8 GB machine, OneDrive sync running) |
|---|---|---|---|
| gate | `gate_human_islet_check_lean.py`, `gate_report.py`, `gate_port.py` | `gate/`, `gate_output.txt` | minutes |
| freeze | `freeze_lists.py` (metadata only) | `frozen_*.csv` | seconds |
| download | `fetch_meta.py`, `fetch_expression.py GSE81547 \| E-MTAB-5061 \| GSE86469`, `fetch_ensembl_ids.py`, `fetch_gse60803.py` | scratch data root, `attempts.csv` | GSE81547 2 min, E-MTAB-5061 6 min |
| parse | `parse_gse81608.py`, `parse_new_datasets.py GSE86469 \| GSE81547 \| E-MTAB-5061` | `parsed/<tag>/` | 1 to 3 min each |
| analysis (4.2, 4.3) | `run_datasets.py [tag ...]` (uses `wp4_core.py`) | `results/` | 6 to 21 s per dataset |
| meta-analysis (4.4) | `meta_analysis.py` | `meta/` | seconds |
| post hoc add-on | `gene_graph_null.py` | `meta/gene_graph_*` | about 2 min |
| metacells (4.5) | `metacell.py [tag ...]` | `metacell/` | about 20 min for all six datasets |
| perturbation (4.6) | `perturbation_scout.py`, `perturbation_decisions.py`, `perturbation_gse60803.py` | `perturbation/` | under 1 min |
| citations | `fetch_dataset_citations.py` | `dataset_citations.csv` | seconds |
| tables and figures | `make_tables.py`, `make_fig7_v2.py` | `WP4_summary_tables.md`, `figures/` | under 1 min |

Seeds: matched null 2026 (main) and 2027 (pair level), 2000 draws; metacell partition i uses `default_rng(5000 + i)`; the meta-analysis is deterministic.

## 15. Extension requested after the primary result: positive control and wider search (6 Nov; full text `README_WP4b.md`, proposed wording `WORDING_islet_crosscheck_v3.md`)
Pre-specified and fingerprinted before the analyses (`PRESPEC_WP4b_positive_control_and_wider_search.md`), reported whichever way they pointed; the primary analysis and its decision rule are unchanged.
* **Positive control: the cross-check is sensitive.** Simulated coupling of nominal strength 0.3 among the evaluable consensus genes (conservative model: only where a gene is detected) would have changed the pooled verdict away from "no meaningful excess co-expression" in 97.5% of 200 simulations (70% at 0.2, 13% at 0.1, 0.5% at 0.05); ribosomal protein pairs give an excess of +0.101 to +0.444 in all five human datasets. Weaker coupling (0.1 or less) cannot be excluded; the smallest margin that gives "no meaningful excess" is 0.078.
* **Wider search: no further support, no contradiction.** 135 records (GEO, ArrayExpress, PubMed-linked series), 6 kept, 3 analysable, all Wnt-arm: sign agreement with the INS-1 DKK-1 arm 7 of 17 (GSE60803), 11 of 17 (GSE167073), 8 of 17 (GSE99503, mirror design); pooled Spearman rho +0.062 (95% CI -0.216 to +0.331), "inconclusive". No Notch series analysable; the closest model (TCF7L2 knockdown in rat INS-1 832/13 cells, E-MTAB-2034) has raw reads only.
* Correction: GSE44045 and GSE44046 are diet comparisons, not Sfrp5 silencing (sections 8 and 11, addendum H).
* Reading: consistent with the paper's conclusion that the pooled network reflects treatment structure; no evidence of co-regulation; not a test of the treatment response; no support for synergy.
