# Wording proposals from WP12 (9 Nov 2026; for the WP8 fix pass, not applied to any manuscript yet)

Rules: American spelling, no thousands commas, every number is in `README.md` of this folder (checked by `check_wp12_numbers.py`), the reported analysis stays primary. Items that need a decision of the authors are marked **[AUTHOR]**.

## 1. Methods, MLSS: annotation of the hub genes (fills the place of the deleted "druggability status" sentence; yellow item 3 of `wp6b_author_items.md`)
"Hub genes (degree of 5 or more) were annotated with the target development level (TDL) stored in the DrugCentral interaction file. GLP1R and CACNA1D are Tclin targets, with 7 and 24 distinct drugs that have a human record in the file used here (8 and 24 in the release of 25 September 2026); TCF7L2 and PTBP1 have no drug-target record in either file. The annotation did not select the hub genes, which were chosen by degree alone."
(The sentence replaces the AUTHOR flag; the claim that druggability chose the hubs is not made.)

## 2. Methods or Results, MLSS: the ranking in classes (fills the place of the deleted "four classes by score quartile" sentence; yellow item 6)
"The 77 combinations of drugs with different targets scored 7.65 to 9.96 and the 76 combinations of drugs with the same target 2.70 to 4.05. The first score quartile (MLSS 7.88 or higher) holds the 42 best combinations, all of them with different targets, and includes the 15 combinations of Fig. 6; the quartile classes (42, 35, 48 and 28 combinations, unequal because many scores are tied) are rank labels of an uncalibrated score and not levels of confidence."
**[AUTHOR]** Option: restore the four old words ("Highest Confidence", "Very High", "High", "Moderate Confidence") as labels of Q1 to Q4. Not recommended: Q3 and Q4 are the same-target pairs that the score places below every cross-target pair, so "High" and "Moderate Confidence" would misdescribe them.
Sensitivity sentence (optional, with the other WP10 sentences): "With the release of 25 September 2026 (171 combinations) the thresholds are unchanged except the lowest (2.75 to 2.87), 14 of the 153 combinations move from the third to the fourth class, and the 18 tirzepatide combinations fall into 5, 6, 6 and 1 of the classes Q1 to Q4."

## 3. Discussion, limits of the safety screen: the market status (replaces the yellow item 8; the reason is now confirmed)
"Albiglutide and bepridil have only discontinued products in the Drugs@FDA database (Tanzeum, BLA 125431, approved 15 April 2014; Vascor and Bepadin, NDA 019002 and 019001, approved 28 December 1990), and nilvadipine and nitrendipine have no FDA application (openFDA, queried on 9 November 2026). Discontinued or never marketed products carry no harmonized drug name in the adverse event reports, which is why the primary search found no report for these four drugs; the secondary free-text search is described above."
Check before use: the last clause restates the WP9 observation (pre-specification section 9, README of `wp9_safety`); the regulatory status of nitrendipine and nilvadipine outside the United States was not assessed and must not be stated.

## 4. Discussion or Limitations, translational reading (optional; the authors decide)
"Ten of the 15 best-scoring combinations, and 35 of the 77 combinations of drugs with different targets, include at least one of these four drugs, which cannot currently be prescribed in the United States; the potency ranking therefore does not prioritize combinations that are available for a trial."
(The numbers 10 of 15 and 35 of 77 are already in the safety-screen text; 75 of the 153 combinations is new.)

## 5. Data availability and Methods sentence for the new source
"Approval and marketing status were taken from Drugs@FDA through the openFDA drugsfda endpoint (public, no key), queried on 9 November 2026; the responses are in Additional file 6 and the repository."
**[AUTHOR]** The raw responses are about 1 MB; confirm that they go into the Zenodo release together with the openFDA adverse event responses (yellow item 12).

## 6. Additional file 6 (proposal for a WP6a round 5, after WP6b has finished; not started)
New sheets: `Ann_hubs_TDL` (the table of README section 2), `Ann_MLSS_classes` (`results/mlss_quartile_classes.csv`: 171 rows, classes in the 153 and in the 171), `Ann_US_status` (`results/market_status_by_drug.csv` plus `market_status_nda_bla_applications.csv`). Each with a source line (file, date, retrieval time) and the caution that status is list membership or an application database, not a regulatory audit.
