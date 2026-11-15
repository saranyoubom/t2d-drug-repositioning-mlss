# WP12 pre-specification: annotations from the updated databases (written 9 Nov 2026, before the openFDA Drugs@FDA retrieval)

Request of the user (9 Nov 2026, in answer to the open author items "source of the druggability and quartile statements" and "market status of albiglutide and nitrendipine"): "add it according to the updated database". Reading used here: derive the three statements from the updated files and databases instead of leaving them for the authors, as a DESCRIPTIVE annotation layer. The reported analysis (Fig. 5 and 6, the 153 combinations) stays primary and is not changed. No decision rule, no hypothesis test. If the user meant something else, only the wording changes, not the tables.

## 1. Hub-gene target annotation (replaces the unsupported "druggability status (Yes/Target/Potential)")
* Input: `Database/drug.target.interaction.tsv.gz` (DrugCentral dump of 25 Sep 2026, sha256 19086849...). Hubs: the four genes with degree >= 5 in the 17-gene consensus network (GLP1R, TCF7L2, CACNA1D, PTBP1, as in `wp10_versions/results/drugcentral_summary.json`); Wnt5a and Cacna1c (degree 4) are listed for context.
* For each gene (exact gene symbol match against the GENE column, split at "|", as the MLSS port does; species counted separately, human records are the ones the MLSS uses): number of records, number of distinct drugs, number of human records, TDL values (Target Development Level of Pharos/IDG as stored in the file: Tclin, Tchem, Tbio, Tdark; empty = not stated). A gene without any record is reported as "no drug-target record in DrugCentral" (no TDL is invented).
* Statement allowed by this table: which hub genes have a DrugCentral drug-target record and what their TDL is. The old statement that hubs were chosen by druggability stays deleted (hubs are chosen by degree only).

## 2. Score-quartile classes of the MLSS combinations (replaces the unsupported "four classes by score quartile")
* Inputs: the 153 reported combinations (`wp10_versions/results/mlss_comparison.csv`, column MLSS_old; identical to the reported run) and the 171 combinations of the 25 Sep 2026 release (`mlss_new_all_combinations.csv`).
* Definition: thresholds t75, t50, t25 = 75th, 50th, 25th percentile of the MLSS scores of the set (numpy linear interpolation). Class Q1: score >= t75; Q2: t50 <= score < t75; Q3: t25 <= score < t50; Q4: score < t25. Ties at a boundary stay together, so the classes need not be equal in size. The labels Q1 to Q4 are rank labels of an uncalibrated score, not levels of confidence; the four words of the old text ("Highest Confidence", "Very High", "High", "Moderate Confidence") are NOT used unless the authors insist.
* Output: class of every combination in both sets, thresholds, class sizes, number of the 153 whose class changes when the 18 tirzepatide combinations are added, class of each of the reported top-15 combinations.

## 3. Regulatory and market status of the 19 MLSS drugs (the four drugs of the safety screen in particular)
* Sources: (a) DrugCentral approval lists of the Database folder (FDA, EMA, PMDA membership; `wp10_versions/results/approval_tiers_old_vs_new.csv`); (b) openFDA Drugs@FDA endpoint `https://api.fda.gov/drug/drugsfda.json`, public, no key, retrieved on the day of the run (date stored with every response).
* Query per drug name (lower case; "(S)-nitrendipine" is queried as nitrendipine): `search=products.active_ingredients.name:"<name>"`, limit 100. Per application: application number, sponsor, brand names, `marketing_status` of every product (Prescription, Discontinued, Over-the-counter, None (Tentative Approval)), and the earliest submission with type ORIG and status AP as the approval date.
* Per drug: number of applications, set of marketing statuses, whether any product is not discontinued ("currently listed as marketed in the US"), first US approval date. A drug with no application is reported as "no FDA application".
* Scope of the statement: United States only (the FAERS reports come mainly from the US); EMA and PMDA are given only as DrugCentral list membership and are not an audit of the European or Japanese market. A status of "Discontinued" in Drugs@FDA means that the application holder has no marketed product listed; it does not say why.
* Responses are stored in `results/raw/` (file name = sha1 of the URL) so that the analysis can be repeated offline.

## 4. Checks
`check_wp12_numbers.py`: independent parser of the raw responses and of the input files; every number of README.md and WORDING.md is rebuilt from them.
