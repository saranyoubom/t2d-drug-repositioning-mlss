# WORDING_safety_screen.md: proposed text for the manuscript (WP9, 7 Nov 2026)

**Status.** Proposed text: adapt it, do not paste it blindly. Every number is in `README.md` and traced in section 9; `check_wp9_numbers.py` rebuilds the numbers of the README from the raw responses. Framing: computational and hypothesis-generating; the screen is an annotation of spontaneous reports, not a prediction. Style of the project: American spelling, no thousands commas, gene symbols italic (done by the build script). Cite the screen as Additional file 5, part 5.5 (protocol, README) and Additional file 6 (sheets `PV_*`); the table of the top 15 pairs is Table 2.

## 1. Methods (new subsection, after the Methods of the MLSS)
**Spontaneous-report screen of the ranked drug pairs (post hoc).** After the ranking was fixed, we annotated the 18 drugs and the 77 cross-target pairs of the Multi-Layer Synergy Score with spontaneous adverse event reports. The screen describes public reports; it does not predict toxicity or adverse effects and does not test the ranking. We used the openFDA drug adverse event endpoint (FDA Adverse Event Reporting System; snapshot of 28 October 2026 with 21054159 reports, retrieved on 7 November 2026) and counted the reports that list a drug through the harmonized active-ingredient field. Eight event groups of 40 MedDRA preferred terms (hypoglycemia; hyperglycemia and loss of glycemic control; pancreatic events; gastrointestinal intolerance; biliary events; hypotension and syncope; edema; conduction and rhythm), the rules and the thresholds were fixed before any count was retrieved in a fingerprinted pre-specification (Additional file 5). For each drug with at least 100 reports, a signal was a reporting odds ratio against all other reports with a lower 95% limit above 1 and at least 3 reports. For each pair with at least 20 reports that list both drugs, a co-report excess flagged an event group whose proportion among those reports had an exact 95% lower limit above the proportion among the reports of each drug without the other; this flag is descriptive and uncorrected and is not a test of interaction. Drugs without harmonized records were also searched by name in the free-text drug field in a secondary analysis that is labeled post hoc. The saved responses, their fingerprint and the scripts are in Additional file 6 and the repository.

## 2. Results (new subsection after the MLSS subsection)
Suggested title: *Spontaneous reports hardly ever list the top-ranked pairs together* (9 words).
Of the 17 search terms that stand for the 18 drugs, 13 had at least 100 reports through the harmonized drug field; albiglutide, bepridil, nilvadipine and nitrendipine had none. Thirty-five of the 77 cross-target pairs, and 10 of the 15 top-ranked pairs, therefore contained a drug that the primary search could not see (Table 2). Among the 13 drugs, 70 of 104 drug-group tests were signals, and 10 drugs had signals in five or more of the eight event groups, also for groups with no obvious link to the drug class (hypoglycemia and pancreatic events for calcium channel blockers); glucagon, the treatment of hypoglycemia, had a reporting odds ratio of 27.3 (95% CI 24.8 to 30.0) for hypoglycemia. Signals therefore reflected the treated population and the class more than single drugs. No report listed both drugs for 48 of the 77 pairs, 1 to 19 reports did for 13 pairs, and 16 pairs had at least 20 (20 to 2921 reports, median 125). Of the 15 top-ranked pairs, 14 had fewer than 20 reports that list both drugs and 11 had none; only rank 10 (liraglutide with nifedipine, 167 reports) could be screened, and it had no co-report excess. The 8 co-report excesses, in 5 of the 16 evaluable pairs (128 pair-group tests), all lay outside the top 15 (Additional file 6). A secondary free-text search, labeled post hoc, added drug profiles for albiglutide (8899 reports) and nitrendipine (1523 reports) but no evaluable pair among the top 15.

## 3. Discussion (two short paragraphs; the second can sit under "Priorities")
**Limits of the safety screen.** The screen describes unvalidated spontaneous reports; it gives no incidence, risk, causality or interaction. Signals against the whole database were common (70 of 104) and mostly marked the class and the older patients who are treated for diabetes and hypertension, so the screen did not separate drugs within a class, and the glucagon result shows confounding by indication. Drugs that are discontinued or not marketed in the United States carry no harmonized records, which hid 5 of the 18 entries (albiglutide, bepridil, nilvadipine, nitrendipine and its enantiomer); the secondary free-text search is noisier. The co-report flags were uncorrected and conservative, and 48 of the 77 pairs had no co-report at all. The screen therefore neither supports nor weakens any pair.
**What the screen adds.** The ranking is blind to real-world use. Most of the top-ranked pairs combine drugs that spontaneous reports almost never list together and that include drugs without harmonized records, whereas the 16 pairs that were reported together often enough (ranks 10 to 71) combine drugs with many reports: amlodipine, nifedipine, felodipine or nisoldipine with liraglutide, semaglutide, exenatide, lixisenatide, glucagon or dulaglutide. A potency ranking from a network hub list is therefore not a clinical shortlist; it needs a market and safety filter before dose-response testing, and the formal assessment of drug-drug interactions remains necessary.

## 4. Limitations (one sentence to add to the list)
The spontaneous-report screen describes unvalidated reports; it gives no incidence, risk, causality or interaction, it could not see 5 of the 18 drugs (discontinued or not marketed in the United States), and its signals were not specific because the comparator was the whole database.

## 5. Abstract fragments (the abstract is limited to 350 words; one sentence is enough)
* Methods: A spontaneous-report screen of the ranked pairs (openFDA, 21054159 reports) was pre-specified before any count was retrieved.
* Results: Of the 15 top-ranked pairs, 14 were reported together fewer than 20 times, and signals against the whole database were not specific (70 of 104 drug-group tests).
* Conclusions: The ranking carries no safety information and favors pairs without real-world co-reporting.

## 6. Conclusions (one sentence)
A spontaneous-report screen showed that the potency ranking favors pairs that are hardly ever reported together and that it carries no safety information, so a market and safety filter should precede any testing of the pairs.

## 7. Declarations: data availability (add to the existing sentence)
The adverse event reports were retrieved from openFDA (https://open.fda.gov/) on 7 November 2026 from the snapshot of 28 October 2026; the 443 saved responses, their SHA-256 fingerprint and the scripts are in the repository ([[DOI to be inserted]]). openFDA states that its data are unvalidated and must not be used for decisions about medical care.

## 8. What not to write
Never: "toxicity prediction", "adverse effect prediction", "predicts", "safe", "unsafe", "risk" for what the screen found, "validated", "confirms", "interaction signal" (the co-report excess is a descriptive flag), any statement that a signal means harm, any statement that the screen supports or weakens the ranking, and any reason for the missing reports (market status) that the authors have not confirmed (README section 9, item 1). Do not describe the signals as specific to a drug. Do not call openFDA counts incidence, events or patients: they are reports.

## 9. Number-to-source table (all in `results/`; README section numbers in brackets)
| number in the text | meaning | source |
|---|---|---|
| 21054159; 28 October 2026 | reports in the database; snapshot | `summary.json`: N_reports, last_updated [2] |
| 443 | saved responses | `RAW_SHA256.txt`; file count of `raw/` [2] |
| 17; 13 | search terms; evaluable (100 reports or more) | `drug_counts.csv`; `summary.json`: terms, evaluable_drugs [3] |
| albiglutide, bepridil, nilvadipine, nitrendipine with 0 reports | not evaluable | `drug_counts.csv`; `summary.json`: not_evaluable_drugs [3] |
| 35 of 77; 10 of 15; 5 of 18 | pairs and entries with a drug that the primary search could not see | `pairs_all77.csv` (signals_drug_A or B = n/e); `drug_counts.csv` [3] |
| 70 of 104; 10 drugs | signals; drugs with signals in five or more groups | `summary.json`: single_signals_total, single_tests_total; `signals_single.csv` [4] |
| 27.3 (95% CI 24.8 to 30.0) | glucagon, hypoglycemia, reporting odds ratio | `signals_single.csv`: glucagon, HYPO [4] |
| 48; 13; 16 | pairs with no co-report; 1 to 19; 20 or more | `pairs_all77.csv`: reports_listing_both [5] |
| 20 to 2921; median 125 | co-reports of the 16 evaluable pairs | `pairs_all77.csv`: status evaluable [5] |
| 14 of 15; 11; rank 10; 167 | top 15 pairs with fewer than 20 co-reports; with none; the one evaluable pair and its co-reports | `pairs_top15.csv` [5] |
| 8; 5; 16; 128 | co-report excesses; pairs with an excess; evaluable pairs; pair-group tests | `summary.json`: co_report_flags, co_report_flagged_pair_count, pairs_evaluable, pair_group_tests [5] |
| 8899; 1523 | free-text report totals for albiglutide and nitrendipine | `fallback_summary.json`: drugs [6] |
| 40; 8 | preferred terms; event groups | `PRESPEC_WP9_safety_screen.md` section 5 |
| 100; 20; 3; 95% | thresholds: reports per drug; co-reports per pair; reports for a signal; interval | `PRESPEC_WP9_safety_screen.md` sections 4 and 6 |
