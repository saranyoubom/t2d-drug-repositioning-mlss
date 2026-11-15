# WP9: pharmacovigilance safety screen of the MLSS drugs and drug pairs (7 Nov 2026)

> Release copy: working document of the analysis. Sections that only listed points for the authors or proposed wording for the manuscript are omitted and marked; the reading key is in `../README.md`.

**Strength of this evidence.** An annotation of the 18 MLSS drugs and the 77 cross-target pairs with spontaneous reports (FDA Adverse Event Reporting System, through openFDA). It is **not a prediction** of toxicity or adverse effects, it is not a test of the MLSS ranking, and it gives no risk estimate: signals and co-report counts are hypotheses. Pre-specified in `PRESPEC_WP9_safety_screen.md` (fingerprints in `../PROTOCOL_v2_SHA256.txt`, protocol section 8 entry of 7 Nov); post hoc with respect to the MLSS (ranking fixed earlier, not changed); the free-text analysis of section 7 is secondary and post hoc (pre-specification section 9). Done in the main session; nothing committed, pushed or uploaded.

## 1. Files
| file | content |
|---|---|
| `PRESPEC_WP9_safety_screen.md` | rules fixed before any count of an MLSS drug was retrieved; section 9 = dated entries (clarifications, secondary free-text fallback) |
| `fda_fetch.py`, `wp9_screen.py`, `wp9_fallback.py`, `freeze_raw.py` | cached client; retrieval and analysis of the primary screen; secondary free-text analysis; fingerprint of the raw responses |
| `results/raw/` | 443 saved openFDA responses (JSON with URL and retrieval time); fingerprint in `results/RAW_SHA256.txt` |
| `results/drug_counts.csv`, `signals_single.csv`, `pairs_all77.csv`, `pairs_top15.csv`, `pairs_long.csv`, `summary.json` | primary results (`-` = none, `n/e` = drug not evaluable, `n/a` = pair not evaluable) |
| `results/fallback_signals_single.csv`, `fallback_pairs.csv`, `fallback_summary.json` | secondary, post hoc, free-text names |
| `results/probe_log.txt`, `retrieve_log.txt`, `fallback_retrieve_log.txt`, `check_log.txt` | logs; `check_wp9_numbers.py` and `check_wp9_last_run.txt` = independent rebuild of every number from the raw responses |
| `WORDING_safety_screen.md` | proposed Methods, Results, Discussion, Limitations and Declarations text with a number-to-source table |

## 2. Data and retrieval
* openFDA drug adverse event endpoint (`api.fda.gov/drug/event.json`), snapshot `last_updated` **2026-10-28** (the FDA updates it quarterly); **21,054,159** reports in all. Retrieved on 7 Nov 2026 from 15:19 to 16:04 (+0700): 443 responses. A DNS failure of the machine stopped the first run at pair 43 of 70; it was resumed and saved responses were not requested again. The raw responses were fingerprinted at 16:04:43, before the analysis (`d3635ab15c6b2ed0...`).
* openFDA states in every response: "Do not rely on openFDA to make decisions regarding medical care. While we make every effort to ensure that data is accurate, you should assume all results are unvalidated." Terms: `https://open.fda.gov/terms/`, license: `https://open.fda.gov/license/`.
* Counts are numbers of **reports**, not events, patients or incidence. Eight event groups of 40 MedDRA preferred terms (section 5 of the pre-specification).

## 3. Coverage: which drugs the screen can see
17 search terms stand for the 18 MLSS entries ((S)-nitrendipine is searched as nitrendipine). **13 of 17 terms are evaluable** (100 reports or more through `patient.drug.openfda.substance_name`). **Four have no report at all: albiglutide, bepridil, nilvadipine and nitrendipine** (so 5 of the 18 MLSS entries). The probable reason, to be confirmed by the authors: products that are discontinued or not marketed in the US carry no harmonized `openfda` fields. The effect on the ranking: **35 of the 77 pairs contain such a drug, and 10 of the top 15** (albiglutide in ranks 1, 2, 7 and 15; nitrendipine in ranks 1, 4, 7, 8, 11, 12, 13 and 14).

| drug | target | reports | drug | target | reports |
|---|---|---|---|---|---|
| amlodipine | CACNA1D | 434,732 | semaglutide | GLP1R | 78,023 |
| nifedipine | CACNA1D | 41,738 | dulaglutide | GLP1R | 89,369 |
| felodipine | CACNA1D | 11,420 | exenatide | GLP1R | 52,322 |
| nisoldipine | CACNA1D | 1,198 | liraglutide | GLP1R | 51,303 |
| nimodipine | CACNA1D | 1,122 | glucagon | GLP1R | 4,235 |
| isradipine | CACNA1D | 560 | lixisenatide | GLP1R | 3,660 |
| clevidipine | CACNA1D | 367 | albiglutide, bepridil, nilvadipine, nitrendipine | | 0 (not evaluable) |

## 4. Single-drug signals (primary)
Rule: reporting odds ratio against all other reports, lower 95% limit above 1 and 3 or more reports (Woolf interval). **70 of 104 drug-group tests are signals** (13 evaluable drugs x 8 groups). By drug: 8 of 8 groups for amlodipine, felodipine, glucagon and nifedipine; 5 of 8 for dulaglutide, exenatide, isradipine, liraglutide, nisoldipine and semaglutide; 3 of 8 for clevidipine and lixisenatide; 2 of 8 for nimodipine. By group, drugs with a signal: PANC 11, HYPO 10, HYPER 10, GI 10, BILI 8, HYPOT 8, CARD 7, OEDEMA 6.
Examples (`signals_single.csv`): liraglutide PANC 3,481 of 51,303 reports, ROR 18.4 (17.8 to 19.0); glucagon HYPO 484 of 4,235, ROR 27.3 (24.8 to 30.0), as expected because glucagon is the treatment for hypoglycemia (confounding by indication); amlodipine OEDEMA ROR 3.08 (3.04 to 3.13) and HYPOT 2.88 (2.84 to 2.92); isradipine (560 reports) HYPOT 2.9 (2.0 to 4.3) and OEDEMA 3.1 (2.0 to 4.6); liraglutide HYPOT 0.91 (0.85 to 0.98) and OEDEMA 0.79 (0.73 to 0.86) are not signals.
**Reading.** Signals are widespread and not specific: the comparator is the whole database, and these drugs are used by older patients with diabetes and hypertension, so even events outside the drug's labelled profile (hypoglycemia and pancreatic events for calcium channel blockers) reach the rule. Class-wide and population-wide patterns are what the screen mostly shows; it does not separate drugs within a class.

## 5. Pairs (primary)
* Of the 77 cross-target pairs, **48 have no report listing both drugs, 13 have 1 to 19, and 16 have 20 or more** (evaluable; ranks 10, 16, 19, 20, 23, 25, 26, 29, 30, 33, 45, 52, 55, 64, 67 and 71; 20 to 2,921 co-reports, median 125). **Of the top 15 pairs, 14 have fewer than 20 co-reports (11 have none; the others have 1, 2 or 3); only rank 10 (liraglutide + nifedipine, 167 co-reports) is evaluable.** The best-ranked pairs are therefore pairs that spontaneous reports almost never list together.
* Co-report excess (flag, descriptive): **8 flags in 5 of the 16 evaluable pairs, among 128 pair-group tests; none in the top 15** (rank 10 has none): rank 23 amlodipine + semaglutide (BILI), rank 25 amlodipine + exenatide (PANC, BILI), rank 30 amlodipine + glucagon (HYPOT, CARD), rank 52 exenatide + felodipine (PANC), rank 64 amlodipine + dulaglutide (GI, BILI). No multiplicity correction and no adjustment for age, indication or co-medication: some flags may be chance, and patients who receive two drugs are different from those who receive one.
* Pairs with at least one drug that has a signal in the group, top 15 / all 77: HYPO 13 / 70, HYPER 12 / 70, PANC 13 / 71, GI 13 / 65, BILI 11 / 61, HYPOT 7 / 53, OEDEMA 7 / 41, CARD 2 / 47 (descriptive; no test of the ranking was planned).

Top 15 pairs (`pairs_top15.csv`; signals = group codes with a signal; reports listing both = co-reports):
| rank | pair | MLSS | co-reports | signals drug A | signals drug B | co-report excess |
|---|---|---|---|---|---|---|
| 1 | (S)-nitrendipine + albiglutide | 9.9604 | 0 | n/e | n/e | n/a |
| 2 | albiglutide + isradipine | 9.9604 | 0 | n/e | HYPO PANC GI HYPOT OEDEMA | n/a |
| 3 | exenatide + isradipine | 9.9604 | 3 | HYPO HYPER PANC GI BILI | HYPO PANC GI HYPOT OEDEMA | n/a |
| 4 | (S)-nitrendipine + liraglutide | 9.9604 | 0 | n/e | HYPO HYPER PANC GI BILI | n/a |
| 5 | isradipine + liraglutide | 9.9604 | 1 | HYPO PANC GI HYPOT OEDEMA | HYPO HYPER PANC GI BILI | n/a |
| 6 | isradipine + lixisenatide | 9.9604 | 0 | HYPO PANC GI HYPOT OEDEMA | HYPO HYPER PANC | n/a |
| 7 | albiglutide + nitrendipine | 9.9604 | 0 | n/e | n/e | n/a |
| 8 | liraglutide + nitrendipine | 9.9604 | 0 | HYPO HYPER PANC GI BILI | n/e | n/a |
| 9 | isradipine + semaglutide | 9.9604 | 2 | HYPO PANC GI HYPOT OEDEMA | HYPO HYPER PANC GI BILI | n/a |
| 10 | liraglutide + nifedipine | 9.9469 | 167 | HYPO HYPER PANC GI BILI | all 8 groups | none |
| 11 | (S)-nitrendipine + semaglutide | 9.9361 | 0 | n/e | HYPO HYPER PANC GI BILI | n/a |
| 12 | nitrendipine + semaglutide | 9.9361 | 0 | n/e | HYPO HYPER PANC GI BILI | n/a |
| 13 | (S)-nitrendipine + exenatide | 9.9037 | 0 | n/e | HYPO HYPER PANC GI BILI | n/a |
| 14 | exenatide + nitrendipine | 9.9037 | 0 | HYPO HYPER PANC GI BILI | n/e | n/a |
| 15 | albiglutide + nifedipine | 9.8792 | 0 | n/e | all 8 groups | n/a |
(n/e = drug not evaluable through the harmonized field; n/a = fewer than 20 co-reports.)

## 6. Secondary, post hoc: free-text names for the drugs that the primary mapping cannot see
Search through `patient.drug.medicinalproduct` with names fixed in advance (pre-specification section 9; the numbers in brackets are reporting odds ratios): albiglutide (TANZEUM, ALBIGLUTIDE) **8,899 reports**, signals HYPO (2.88), HYPER (4.11), PANC (2.34) and GI (1.09); bepridil (VASCOR, BEPRIDIL) **356**, signals HYPOT (3.72) and CARD (20.4); nilvadipine **138**, signals OEDEMA (3.50) and CARD (6.80); nitrendipine **1,523**, signals in 7 of 8 groups (all but HYPO; CARD (10.7), HYPOT (6.15), BILI (6.03), OEDEMA (4.70)), which looks like reports with many co-medications rather than a drug effect. Of the 39 pairs that contain one of these four drugs, **one has 20 or more co-reports (albiglutide + amlodipine, rank 21, 86 co-reports; excess flag GI)**; none of the top 15 is evaluable (co-reports among them: nitrendipine + liraglutide 5, nitrendipine + semaglutide 4, albiglutide + nifedipine 3, all others 0). The free-text search adds drug profiles but not pair information for the top-ranked pairs. Free-text matching is noisier (spelling, brand names) and is not merged with the primary counts.

## 7. What this screen shows, and what it does not
* **Shows:** (1) spontaneous reports hardly ever list the top-ranked pairs together (14 of the 15 have fewer than 20 co-reports; 11 none), and the evaluable pairs are the common drugs of lower rank; the potency ranking therefore prioritizes combinations without real-world co-reporting. (2) The ranking carries no safety information: nearly every drug has signals in several groups, and the screen does not separate the pairs of the ranking. (3) 5 of the 18 MLSS entries (35 of 77 pairs, 10 of the top 15) are drugs that the harmonized FAERS fields cannot see (discontinued or not marketed in the US).
* **Does not show:** that any pair is safe or unsafe, any incidence or risk, any interaction (co-report excess is a descriptive flag), any causal link, or that the MLSS ranks safety. It cannot test the treatment response either.

## 8. Limits
Spontaneous reporting (under-reporting, stimulated reporting for newer drugs, duplicate reports, no denominators or exposure data); mapping through harmonized `openfda` fields (reports without them are not seen); confounding by indication, age and co-medication (glucagon for hypoglycemia is the clear example); disproportionality against the whole database is not risk and flags class-wide and population-wide patterns; pair counts are sparse (48 of 77 pairs have no co-report); the co-report flag is uncorrected and conservative; no information on dose or timing; one snapshot of a database that changes every quarter; the free-text analysis is post hoc and noisy.

## 9. Items for the authors (omitted)

Omitted from this copy: this section listed points for the authors that were settled while the article was prepared.

## 10. Re-running and checks
`python -I wp9_screen.py retrieve` (network; resumable; at most 850 new requests per run; the openFDA limit without a key is 1000 per day), `python -I wp9_fallback.py retrieve`, `python -I freeze_raw.py` (before the analysis), `python -I wp9_screen.py analyze`, `python -I wp9_fallback.py analyze`, `python -I check_wp9_numbers.py` (own parser, own formulas; compares with the output files; result in `check_wp9_last_run.txt`). Python 3.14 with numpy and scipy.
