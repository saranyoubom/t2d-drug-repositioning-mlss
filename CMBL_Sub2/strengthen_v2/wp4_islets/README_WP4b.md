# WP4 extension (6 Nov 2026): positive control of the islet cross-check and a wider search for matched perturbation data

Requested by the user after the primary WP4 result was known, with the stated hope of a more supportive cross-check. Plan, queries, keep rules, reading rules and decision rules were fixed and fingerprinted BEFORE any of the analyses ran (`PRESPEC_WP4b_positive_control_and_wider_search.md`; `../PROTOCOL_v2_SHA256.txt`; protocol section 8 entries of 16:11 and 18:35). The primary analysis, the margin of +/-0.10, the datasets, the donors and the decision rule of protocol 4.4 are unchanged. Results are reported as they came out. Nothing here changes the manuscript or the packages; no commit or push.

## 1. Summary
1. **The cross-check is sensitive** (reading rule R1 met). If the evaluable consensus genes had been coupled at a nominal within-donor correlation of 0.3 (a conservative model: coupling only where a gene is detected), the pooled verdict would have been "no meaningful excess co-expression" in only 2.5% of 200 simulations (5 of 200; it would have read "inconclusive" in 97.5%); at 0.2 in 30%, at 0.1 in 87%, at 0.05 in 99.5% (199 of 200). Ribosomal protein pairs, a real co-expression module, give an excess of +0.101 to +0.444 in all five human datasets (reading rule met). So the null result is informative for coupling of about 0.2 or more and for ribosomal-level co-expression; weaker coupling (0.1 or less) cannot be excluded.
2. **The wider search found no further support and no contradiction.** 135 records from three sources, 6 kept, 3 analysable, all for the Wnt arm: sign agreement with the INS-1 DKK-1 arm 7 of 17 (GSE60803), 11 of 17 (GSE167073) and 8 of 17 genes (GSE99503, a mirror design); pooled Spearman rho +0.062 (95% CI -0.216 to +0.331), outcome "inconclusive" (3 series, 17 genes each). No Notch-arm series could be analysed. The closest model (TCF7L2 knockdown in rat INS-1 832/13 cells, E-MTAB-2034, 4 vs 4) is available only as raw reads and was not processed.
3. **Reading.** The cross-check supports the paper's own conclusion that the pooled network reflects treatment structure rather than baseline co-expression modules, and it gives no evidence of co-regulation of the consensus pairs. It does not test the treatment response, and it does not support synergy. Wording: `WORDING_islet_crosscheck_v3.md`.

## 2. A: positive control (`positive_control.py`; `positive_control/`)
* **A1, synthetic spike-in** into the real cells of each dataset (one latent factor per replicate added to the distinct genes of the evaluable consensus pairs on the donor-centred log scale; loading sqrt(rho / (1 - rho)) times the gene's SD; nominal rho 0.05, 0.1, 0.2, 0.3, 0.5; mode "detected" = only in cells where the gene is detected, zeros stay zeros (primary, conservative); mode "all" = every cell (optimistic upper bound); 200 replicates; the unchanged depth-adjusted pipeline and the matched null of the primary run). Gate: without injection the code reproduces E and SE of every dataset to 1e-9 (passed for all six).

| nominal rho | mode "detected": no meaningful excess / inconclusive / co-expressed more | pooled E (mean) | mode "all": no meaningful excess / inconclusive / co-expressed more | pooled E (mean) |
|---|---|---|---|---|
| 0.05 | 0.99 / 0.01 / 0.00 | +0.020 | 0.00 / 1.00 / 0.00 | +0.119 |
| 0.10 | 0.87 / 0.13 / 0.00 | +0.028 | 0.00 / 0.08 / 0.92 | +0.186 |
| 0.20 | 0.30 / 0.70 / 0.00 | +0.044 | 0.00 / 0.00 / 1.00 | +0.290 |
| 0.30 | 0.03 / 0.97 / 0.00 | +0.063 | 0.00 / 0.00 / 1.00 | +0.381 |
| 0.50 | 0.00 / 1.00 / 0.00 | +0.108 | 0.00 / 0.00 / 1.00 | +0.582 |

  The primary result (no injection) is pooled E +0.013 [-0.052, +0.078]. **R1 (false-reassurance rate at 0.3, mode "detected", at most 20%) is MET: 2.5% (5 of 200).** Per dataset in mode "detected" at nominal 0.3 (recovered E against the primary E): E-MTAB-5061 +0.121 (+0.045), GSE81547 +0.042 (+0.025), GSE81608 +0.030 (-0.023), GSE86469 +0.035 (-0.020), GSE84133 mouse +0.057 (+0.049), GSE84133 human +0.192 (+0.192); attenuation of the injected effect 0.24, 0.06, 0.17, 0.18, 0.03 and 0.00. Single-dataset power (empirical P of E below 0.05) at 0.3: GSE81608 0.655, GSE86469 0.465, GSE81547 0.750, mouse 0.010; GSE84133 human and E-MTAB-5061 are already below 0.05 without injection.
* **A2, ribosomal protein pairs** (100 random sets of 11 pairs of cytosolic ribosomal protein genes per dataset against matched random pairs): mean E_RP +0.162 (GSE81608), +0.215 (GSE84133 human), +0.101 (GSE86469), +0.444 (GSE81547), +0.153 (E-MTAB-5061), +0.177 (mouse); in units of the SE of the consensus E: 10.5, 2.5, 4.6, 20.5, 7.0, 3.4. Reading rule (above +0.10 in at least 4 of the 5 human datasets): MET (5 of 5; GSE86469 only just, 52% of its sets are above 0.10).
* **A3, margin sensitivity:** the primary interval [-0.052, +0.078] reads "no meaningful excess" for margins of 0.10 and above and "inconclusive" at 0.05 and 0.075; the smallest margin that gives "no meaningful excess" is 0.078. The pre-set margin stays +/-0.10.
* **Limits.** (i) R1 is a low bar by construction: the primary interval ends only 0.022 below the margin, so a modest rise of the pooled E turns the verdict into "inconclusive"; the stricter question, whether the cross-check could have said "co-expressed more", is answered "no" for the conservative model at every nominal rho up to 0.5 and "yes" from 0.1 for the optimistic model. (ii) The answer depends on how coupling is modelled; the two modes bracket it and the ribosomal control gives a real-data anchor. (iii) The mouse dataset is insensitive in mode "detected" (power 0.01 at 0.3). (iv) Nominal rho is not a prediction of the correlation of real regulatory modules. (v) The test code was run once on GSE86469 with 5 replicates before the full run; nothing was changed after it.

## 3. B: wider search for matched perturbation data (`wider_search*.py`; `wider_search/`)
* **Sources, each run once** (queries in `wider_search/queries.json`): GEO DataSets 88 records; PubMed 371 PMIDs linked to 65 GEO records (62 series); the ArrayExpress collection of BioStudies 23 hits (9 outside GEO); 135 unique records (`candidates_all.csv`).
* **Decisions** (`frozen_candidates.csv`, fingerprinted before any value was read; every record has a reason): 6 kept, 129 dropped: 39 stem-cell, progenitor or reprogramming designs, 22 not expression profiles, 12 organisms other than human, mouse or rat, 12 with fewer than 6 samples, 20 with no Notch or Wnt perturbation in title or summary, 24 for a specific reason after reading (title, summary or sample records). Mirror designs (Notch gain of function, Wnt activation) are kept with the expected sign reversed. The sample records corrected one earlier statement: GSE44045 and GSE44046 are diet comparisons, not Sfrp5 silencing (addendum H).
* **Kept series:**

| series | design | arm, expected sign | data | result |
|---|---|---|---|---|
| GSE60803 | beta-catenin shRNA in expanded human beta-cell-derived cells, 4 donors | DKK-1, +1 | series matrix | analysed; confounded (controls not expanded) |
| GSE167073 | Fzd9 knockout against wild type in islets of MycER-BclXL mice, 4 vs 4 | DKK-1, +1 | series matrix and platform table | analysed; receptor-level, tumour model |
| GSE99503 | GSK3 inhibitor (1-azakenpaullone) against vehicle in scramble-siRNA EndoC-betaH1 cells, 3 vs 3 | DKK-1, -1 (mirror) | NCBI-generated counts | analysed; inhibitor not pathway-specific |
| GSE193888 | Notch gain of function (NICD) in mouse islets, 4 vs 4 | DAPT, -1 (mirror) | none by plain download | not analysable (K4) |
| GSE249170 | Jagged1 gain of function in mouse islets, 3 vs 3 | DAPT, -1 (mirror) | none by plain download | not analysable (K4) |
| E-MTAB-2034 | TCF7L2 knockdown against control in rat INS-1 832/13 cells, 4 vs 4 | DKK-1, +1 | raw reads only | not analysable (K4) |

* **Concordance with the matching INS-1 arm** (17 genes; INS-1 log2 fold changes "as reported"; `concordance_results.csv`, `concordance_report.md`):

| series | same sign | exact binomial P | Spearman rho (P) | strong INS-1 genes (q < 0.05) same sign |
|---|---|---|---|---|
| GSE60803 | 7 of 17 | 0.629 | -0.037 (0.889) | 3 of 6 |
| GSE167073 | 11 of 17 | 0.332 | +0.184 (0.480) | 4 of 6 |
| GSE99503 (mirror) | 8 of 17 | 1.000 | +0.037 (0.889) | 1 of 6 |

  **Pooled (REML + Hartung-Knapp on Fisher z; 3 series):** rho +0.062, 95% CI -0.216 to +0.331, I2 0%: "inconclusive". The Notch series, the "mirror excluded" and the "GSE60803 excluded" sets have fewer than 3 analysable series, so the rule is not applied.
* **Limits.** Three analysable series of 17 genes each give almost no power; all three are weak matches (confounding, receptor-level loss in a tumour model, a non-specific inhibitor); the Notch arm is empty. The three series without data could be added by processing raw reads (E-MTAB-2034, ENA) or by obtaining count tables from the authors (GSE193888, GSE249170); aligning reads needs software that is not installed here, and nothing was installed or downloaded for it.

## 4. Deviations, additions and code fixes (all in protocol section 8)
Times are the machine clock (+0700). 16:10 the plan note fixed (with the framing paragraph); 16:11 protocol entry. Before the full run, a code test (GSE86469, 5 replicates, 5 ribosomal sets, separate folder, deleted) showed an attenuation of about 0.2 in mode "detected" and an E_RP of +0.09; no mode, level, gene list or reading rule was changed after it. The report writer crashed once on an unescaped percent sign after all six dataset runs had saved their files (fixed, no effect). The first ribosomal id request returned 20 of 77 symbols (default API page size); repeated with a larger page size. 16:40 the wider-search files and addendum H (Sfrp5 correction) fingerprinted before any value of a new series was read. The GPL17400 platform table was fetched through the GEO text export because no annotation file exists. 18:35 results and deviations logged.

## 5. Files and how to run
`positive_control.py [tag ...|report]` (about 1 min per dataset at full speed; outputs in `positive_control/`), `rp_list.py`, `fetch_rp_ids.py`, `wider_search.py` (three sources, under 1 min), `wider_search_meta.py ACC ...` (sample records), `wider_search_freeze.py`, `wider_search_analyze.py`, `check_wp4b_numbers.py` (every number of this README rebuilt from the result files). Seeds: `SeedSequence([2029, dataset, rho, mode, replicate])` for the spike-in, `default_rng(2030 + dataset)` and `SeedSequence([2030, dataset, decile, decile])` for the ribosomal control. Data root as in README.md section 14 (`v2_data/wp4/wider/`). Python 3.14, `-I`.
