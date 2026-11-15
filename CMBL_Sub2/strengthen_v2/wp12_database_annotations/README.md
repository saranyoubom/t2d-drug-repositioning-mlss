# WP12: annotations from the updated databases (9 Nov 2026)

**Strength of this evidence.** A descriptive annotation layer with no hypothesis and no decision rule. The reported analysis (Fig. 5 and 6, the 153 combinations) stays primary and is not changed. Pre-specified in `PRESPEC_WP12_annotations.md` and fingerprinted (`FINGERPRINTS.txt`, `../PROTOCOL_v2_SHA256.txt`) on 9 Nov 2026, 00:33 +0700, before the openFDA retrieval (00:34 to 00:35). Independent check: `check_wp12_numbers.py`, result in `check_wp12_last_run.txt`.

**Why it exists.** The user answered the open author items on the druggability statement, the quartile classes and the market status of albiglutide, bepridil, nilvadipine and nitrendipine with "add it according to the updated database". The three items are therefore derived here from the DrugCentral files of the Database folder (25 Sep 2026) and from the FDA's own application database (Drugs@FDA through openFDA). The INS-1 source and authentication cannot come from a database and stay with the authors.

## 1. What the databases say (the answers in one place)
| item | result |
|---|---|
| Druggability of the hub genes | GLP1R and CACNA1D are Tclin targets (target development level of Pharos/IDG, as stored in the DrugCentral file) with 7 and 24 distinct drugs that have a human record in the file used for the reported run (8 and 24 in the file of 25 Sep 2026); TCF7L2 and PTBP1 have no drug-target record in either file. Hubs were chosen by degree alone, so the old sentence that druggability status chose them stays deleted. |
| Score-quartile classes | They exist only as a description of the ranking: the 77 cross-target combinations score 7.65 to 9.96, the 76 same-target combinations 2.70 to 4.05, and nothing lies between. The first quartile (MLSS 7.88 or higher) holds the 42 best cross-target combinations, including all 15 of Fig. 6. Classes Q3 and Q4 contain only same-target pairs. |
| Market status of the four drugs | Drugs@FDA (retrieved 9 Nov 2026): albiglutide (Tanzeum, BLA 125431, approved 15 April 2014) and bepridil (Vascor, NDA 019002, and Bepadin, NDA 019001, both approved 28 December 1990) have only discontinued products; nilvadipine and nitrendipine have no FDA application. This is the "discontinued products or not marketed" reason that the WP9 README gave as probable, now confirmed for the United States. |

## 2. Hub genes (part 1; `results/hub_annotation.csv`)
Exact gene-symbol match against the GENE column (split at "|"), as the MLSS port does. "Reported file" = `R_17genes/drug.target.interaction.csv` (19378 rows, the file of the reported run); "25 Sep 2026" = `Database/drug.target.interaction.tsv.gz` (22364 rows).
| gene | role | records, all species / human, reported file | distinct drugs with a human record | TDL | records, all species / human, 25 Sep 2026 | distinct drugs | drugs that entered the MLSS (reported / 25 Sep 2026) |
|---|---|---|---|---|---|---|---|
| GLP1R | hub (degree 6) | 8 / 7 | 7 | Tclin | 9 / 8 | 8 | 7 / 8 (tirzepatide added) |
| TCF7L2 | hub (degree 5) | 0 / 0 | 0 | no record | 0 / 0 | 0 | 0 / 0 |
| CACNA1D | hub (degree 5) | 25 / 25 | 24 | Tclin | 25 / 25 | 24 | 11 / 11 |
| PTBP1 | hub (degree 5) | 0 / 0 | 0 | no record | 0 / 0 | 0 | 0 / 0 |
| WNT5A | context (degree 4) | 0 / 0 | 0 | no record | 0 / 0 | 0 | not a hub |
| CACNA1C | context (degree 4) | 41 / 35 | 33 | Tclin | 41 / 35 | 34 | not a hub |
TDL is stored only for genes that have a record, so no TDL is given for TCF7L2 and PTBP1. The MLSS used fewer drugs than the records contain (for example 11 of 24 for CACNA1D): the MLSS rules (Methods, WP3 README section 8) need usable potency values; this table counts records, not usable potency values.

## 3. Score-quartile classes (part 2; `results/mlss_quartile_classes.csv`)
Definition (pre-specified): t75, t50, t25 are the 75th, 50th and 25th percentiles of the scores of the set; Q1 = score at or above t75, Q2 = at or above t50, Q3 = at or above t25, Q4 = below t25. Ties stay together.
| set | t75 / t50 / t25 | Q1 / Q2 / Q3 / Q4 | rank ranges |
|---|---|---|---|
| 153 reported combinations | 7.88 / 7.65 / 2.75 | 42, 35, 48 and 28 | 1-42, 43-77, 78-125, 126-153 |
| 171 combinations (25 Sep 2026) | 7.88 / 7.65 / 2.87 | 47, 41, 40 and 43 | |
* The sizes differ from 25% each because 47 distinct scores share the 153 places (14 combinations have exactly the score 7.65, the lowest cross-target score).
* The classes are not independent information: Q1 and Q2 are exactly the 77 cross-target combinations and Q3 and Q4 the 76 same-target combinations (scores 2.70 to 4.05, against 7.65 to 9.96 for the cross-target ones). The "complementarity" term of the MLSS separates the two groups completely.
* All 15 combinations of Fig. 6 are in Q1.
* With the 171 combinations the t25 threshold moves from 2.75 to 2.87; 14 of the 153 reported combinations change class (all from Q3 to Q4); the 18 new tirzepatide combinations fall into 5, 6, 6 and 1 of the classes Q1 to Q4. Scores of the shared combinations are identical (WP10).
* **Recommendation:** do not restore the four labels "Highest Confidence", "Very High", "High" and "Moderate Confidence": the MLSS is an uncalibrated potency-and-complementarity score, and those words would call same-target pairs "High" or "Moderate Confidence" although the score ranks them below every cross-target pair. State the cross-target/same-target split and Q1 instead (WORDING_annotations.md).

## 4. Regulatory and market status (part 3; `results/market_status_by_drug.csv`, `market_status_nda_bla_applications.csv`, raw responses in `results/raw/`)
openFDA endpoint `drug/drugsfda.json`, query `products.active_ingredients.name:"<name>"`, 19 saved responses (amlodipine needs two pages), retrieved 9 Nov 2026 00:34 to 00:35 +0700. DrugCentral columns are list membership (FDA, EMA, PMDA) from the files of the Database folder and are not a regulatory audit.
| drug | applications in Drugs@FDA | NDA/BLA application(s) | marketing status of the products | first US approval (ORIG, AP) | DrugCentral FDA / EMA / PMDA |
|---|---|---|---|---|---|
| albiglutide | 1 | BLA 125431 (Tanzeum, GlaxoSmithKline) | Discontinued | 15 April 2014 | yes / no / no |
| bepridil | 2 | NDA 019001 (Bepadin, MedPointe), NDA 019002 (Vascor, Johnson & Johnson) | Discontinued | 28 December 1990 | yes / no / no |
| nilvadipine | 0 | none | no FDA application | none | no / no / no |
| nitrendipine and (S)-nitrendipine | 0 | none | no FDA application | none | no / no / no |
The other 14 drug entries each have at least one product with the status Prescription (none is Over-the-counter). For isradipine both brand applications (Dynacirc, Dynacirc CR) are discontinued and only two generic applications (Elite Labs ANDA 077169, Watson/Teva ANDA 077317) are listed as Prescription.
* **Consequence for the screen:** 10 of the 15 top-ranked combinations, 35 of the 77 cross-target pairs and 75 of the 153 combinations contain one of albiglutide, bepridil, nilvadipine or nitrendipine. The first two numbers are the ones in the WP9 README (10 of 15, 35 of 77), whose probable reason is now confirmed. For the translational reading this means that 10 of the 15 best-scoring pairs include at least one drug that cannot currently be prescribed in the United States.
* **Limits:** United States only; "Discontinued" says that the application holder lists no marketed product, not why; "Prescription" means an active listing in Drugs@FDA, not a check of sales. EMA and PMDA are list membership only: the status of nitrendipine and nilvadipine in other countries was not assessed, the DrugCentral lists are not an audit, so "not approved anywhere" must NOT be written.

## 5. Deviations and limits
* The first version of part 1 counted only the file of 25 Sep 2026; the reported file was added so that the numbers of the Methods (which describe the reported run) are the primary ones. No rule changed.
* The quartile classes use 47 distinct score values; a different tie rule (rank-based quartiles) would split ties and change class sizes, not the cross/same-target structure.
* TDL is the annotation of the DrugCentral file (Pharos/IDG); the old three-level labels "Yes/Target/Potential" have no source and are not reconstructed.

## 6. Files and re-running
`PRESPEC_WP12_annotations.md`, `FINGERPRINTS.txt`, `wp12_annotations.py` (`retrieve` needs the network; `analyze` reads saved files only), `check_wp12_numbers.py`, `WORDING_annotations.md`, `results/` (`hub_annotation.csv`, `mlss_quartile_classes.csv`, `market_status_by_drug.csv`, `market_status_nda_bla_applications.csv`, `summary.json`, `analyze_stdout.txt`, `raw/`). Use `python -I`.
