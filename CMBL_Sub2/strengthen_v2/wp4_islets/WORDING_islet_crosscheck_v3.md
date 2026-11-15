# Proposed wording for the islet cross-check, version 3 (6 Nov 2026; for WP6b; numbers from the result files of `README.md`, `README_WP4b.md`)

Status: proposals, not manuscript text. They were written AFTER the results and follow the reading rules fixed before them. The numbers are quoted from files; the table at the end names the file for each.

## 0. Framing (decided with the user)
The network, the drug ranking and the cross-check are **computational, data-driven, hypothesis-generating** results from one preparation. Use: "identifies", "ranks", "prioritizes", "is consistent with", "suggests", "would have been detected". Do not use: "validates", "confirms", "supports co-regulation", "synergy", "mechanism". The wet-lab readouts are observations of a single experiment. Wet-lab confirmation is the next step and is described as PROPOSED, not done (section 8).

What the cross-check is: a test of whether the 21 consensus pairs are co-expressed from cell to cell in untreated human or mouse beta-cells. What it is not: a test of the treatment response, of co-regulation, or of synergy.

## 1. Abstract, Results sentence (replaces "Public islet single-cell data (GSE81608, GSE84133) could evaluate only 3 to 11 of the 21 consensus pairs, with inconsistent results.")
Long form (about 55 words): "In five public human islet datasets, the 11 evaluable consensus pairs were not more co-expressed than expression-matched random pairs (pooled excess +0.01, 95% CI −0.05 to +0.08); positive controls were detected, 10 pairs, including all Wnt pairs, could not be evaluated, and three perturbation datasets gave inconclusive concordance."
Short form (about 30 words): "In five public human islet datasets, evaluable consensus pairs showed no excess co-expression over matched random pairs (pooled +0.01, 95% CI −0.05 to +0.08); 10 of 21 pairs could not be evaluated."

## 2. Methods paragraph (new, in "Public islet single-cell data")
"To judge what the cross-check could detect, a latent coupling factor was added to the genes of the evaluable consensus pairs of each dataset (nominal within-donor correlation 0.05 to 0.5, added only to cells in which a gene was detected, or to all cells as an upper bound; 200 simulations per setting), and the pooled analysis was repeated; ribosomal protein pairs served as a real-gene positive control (100 random sets of 11 pairs per dataset). The reading rules were fixed before the runs. Perturbation datasets were searched once in GEO, ArrayExpress and PubMed-linked GEO series with a fixed query. Records were kept if they profiled islets or beta-cell lines after inhibition of Notch or Wnt signaling, or after the mirror manipulation with the sign reversed, with at least three samples per arm; decisions and contrasts were frozen before any expression value was read. The log2 fold changes of the 17 genes were compared with the matching INS-1 arm by sign agreement (exact binomial test) and Spearman correlation and pooled by random-effects meta-analysis."

## 3. Results paragraph (new, after the pooled result)
"The cross-check could detect co-expression when it was present. Ribosomal protein pairs showed an excess over matched random pairs of +0.10 to +0.44 in all five human datasets, and simulated coupling of nominal strength 0.3 changed the pooled verdict away from "no meaningful excess" in 97.5% of simulations (70% at 0.2, 13% at 0.1), so coupling of that size among the evaluable consensus genes would have been flagged; weaker coupling could not be excluded. Of 135 perturbation records, six met the rules and three had analyzable data, all for the Wnt arm (GSE60803, GSE167073, GSE99503); the fold changes agreed in sign with the INS-1 DKK-1 arm for 7, 11 and 8 of 17 genes (pooled Spearman ρ +0.06, 95% CI −0.22 to +0.33, inconclusive). No Notch inhibition series with processed data was found. The closest model, TCF7L2 knockdown in rat INS-1 832/13 cells, was available only as raw reads."
Note: "the pooled verdict away from 'no meaningful excess'" means "inconclusive" in the conservative model (97%, 70%, 13% are the shares not reading "no meaningful excess": 1 − 0.03, 1 − 0.30, 1 − 0.87); "co-expressed more" occurred only in the optimistic model.

## 4. Discussion paragraph (replaces the paragraph that says the data "neither confirm nor exclude")
"This cross-check is computational and tests baseline co-variation in untreated beta-cells, not the response to DAPT or DKK-1. That the evaluable consensus pairs show no excess co-expression is consistent with our conclusion that the pooled INS-1 network mainly reflects treatment structure; it neither confirms nor refutes co-regulation. Its sensitivity was shown by positive controls but depends on how coupling is modeled: under a conservative model that restricts coupling to detected values, it flags nominal coupling of 0.2 to 0.3 but not 0.1 or less, and it could not have returned "co-expressed more" even at 0.5. Ten of the 21 pairs, including all seven Wnt pairs, could not be evaluated because the genes are rarely detected in single beta-cells. The perturbation datasets were few, one was confounded, they cover one pathway, and the concordance was inconclusive. The candidate pairs and the potency-ranked drug combinations therefore remain hypotheses."

## 5. Limitations (one bullet)
"The cross-check could evaluate 11 of 21 consensus pairs and none of the seven Wnt pairs; it measures cell-to-cell rather than treatment-induced variation; its sensitivity depends on the coupling model (conservative model: nominal coupling of 0.2 or more); and only three analyzable perturbation datasets, all for the Wnt arm, were available."

## 6. Conclusions (one sentence)
"A cross-check in five public islet datasets found no excess co-expression of the evaluable consensus pairs, with sensitivity shown by positive controls, and the available perturbation datasets gave inconclusive concordance; it cannot test the treatment response, so the candidate pairs and the potency-ranked drug combinations remain hypotheses for independent wet-lab confirmation."

## 7. Figure legend insertions (for Fig. 8; WP6a)
* (c) "Log2 fold changes of the 17 genes in three perturbation datasets against the INS-1 DKK-1 arm: beta-catenin shRNA in expanded human beta-cell-derived cells (GSE60803), Fzd9 knockout in mouse islets (GSE167073) and GSK3 inhibition in EndoC-βH1 cells (GSE99503; Wnt activation, fold changes multiplied by −1)."
* (d) "Share of 200 simulations in which the pooled verdict was "no meaningful excess co-expression" after coupling of the stated nominal strength was added to the evaluable consensus genes; solid line, coupling added only to cells in which a gene was detected (conservative); dashed line, added to all cells (upper bound)."
* (e) "Excess co-expression of ribosomal protein pairs (mean of 100 random sets of 11 pairs) and of the consensus pairs in each dataset; dotted lines, the ±0.10 margin."

## 8. Wet-lab confirmation (proposed; for Discussion, Priorities)
"Wet-lab confirmation is the next step and has not been done. A design that can test the network at the level of single pairs needs independent biological replicates: in simulations with three groups and a six-gene module, 6 to 8 independent experiments per group detected module-level coupling at r = 0.5 with power 0.78 to 0.96, and 12 to 20 per group detected 52% to 95% of the true pairs (Additional file 5). Each independent experiment should contain all three groups, all 17 genes on the same plate or run with a randomized layout, at least two reference genes, and efficiency-corrected quantification reported according to MIQE. For the drug ranking, dose-response matrices for the top-ranked pair (liraglutide plus isradipine) and for low-ranked pairs, analyzed against Bliss and Loewe reference models, would test whether the potency ranking predicts combination effects."

## 9. What not to write
* "The public data support/confirm the consensus pairs or the network" (they show no excess co-expression).
* "The public data exclude co-expression" (weak coupling, 0.1 or less, cannot be excluded; Wnt pairs were not evaluable).
* "Perturbation data confirm the DKK-1 or DAPT effects" (3 weak-match series, inconclusive, no Notch series).
* Any statement that the potency ranking predicts synergy.

## 10. Numbers and sources
| number | file |
|---|---|
| pooled E +0.013, 95% CI −0.052 to +0.078; 11 of 21 pairs evaluable; 10 never evaluable | `meta/meta_report.md`; `WP4_summary_tables.md` Tables 2b and 3 |
| 97.5% / 70% / 13% / 0.5% not "no meaningful excess" at nominal 0.3 / 0.2 / 0.1 / 0.05 (shares 0.025, 0.300, 0.870, 0.995 reading "no meaningful excess"; 5, 60, 174, 199 of 200) | `positive_control/positive_control_pooled.csv` |
| ribosomal pairs +0.101 to +0.444 | `positive_control/rp_control.csv` |
| 135 records, 6 kept, 3 analysable; 7, 11, 8 of 17; pooled ρ +0.062 [−0.216, +0.331] | `wider_search/queries.json`, `frozen_candidates.csv`, `concordance_results.csv`, `concordance_meta.csv` |
| 6 to 8 and 12 to 20 independent experiments per group | `../wp2_benchmark2/design_table.md` |
