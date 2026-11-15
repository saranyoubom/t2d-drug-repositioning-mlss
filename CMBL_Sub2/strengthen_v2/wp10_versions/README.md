# WP10: database-version sensitivity of the STRING classification and of the MLSS (8 Nov 2026)

> Release copy: working document of the analysis. Sections that only listed points for the authors or proposed wording for the manuscript are omitted and marked; the reading key is in `../README.md`.

**Strength of this evidence.** A sensitivity layer, descriptive, with no decision rule. The reported analysis stays primary and unchanged. Pre-specified in `PRESPEC_WP10_database_versions.md` and fingerprinted (`FINGERPRINTS.txt`: pre-specification, the new DrugCentral files, the reported STRING and DrugCentral inputs) before the STRING 12.5 network was retrieved and before the MLSS was run on the new file. Done in the main session; nothing committed, pushed or uploaded. One finding is a correction candidate for the reported analysis (section 3).

## 1. Files and versions
| item | content |
|---|---|
| STRING 12.5 | API at `https://version-12-5.string-db.org` (version endpoint answered 12.5), rat (10116), the 21 genes of the matrix (Glut2 queried as Slc2a2), minimum score 0.4, functional network: 59 distinct pairs (59 in the reported v12.0 export). Raw: `results/string_v12_5_network.tsv`, `results/string_version.json` |
| DrugCentral | `Database/drug.target.interaction.tsv.gz` (sha256 `19086849...`, downloaded by the user on 7 Nov 2026 from drugcentral.org, database dump of 25 Sep 2026, the "2027 Update"): 22364 rows against 19378 in the file of the reported run, 2897 distinct drugs against 2587. The 1.4 GB dump in the same folder was not used |
| code | `string_sensitivity.py` (gate, retrieve, analyze), `drugcentral_sensitivity.py` (uses `../wp3_mlss/mlss_port.py` unchanged), `check_wp10_numbers.py` and `check_wp10_last_run.txt` |
| results | `results/string_gate.json`, `string_summary.json`, `string_classification.csv` (120 pairs, both versions, both name treatments); `drugcentral_summary.json`, `mlss_comparison.csv`, `mlss_new_all_combinations.csv`, `pic50_changes.csv` |

## 2. Part A: STRING 12.0 versus 12.5
Gate: the Python port of the Fig. 5 classification reproduces the reported numbers from the 12.0 export: 120 pairs with |r| >= 0.7, Known Validated 5, Novel Discovery 16, STRING Only 29, Not Validated 70, **34 in STRING (28.3%)**: PASS.
| classification of the 120 pairs | Known Validated | Novel Discovery | STRING Only | Not Validated | in STRING |
|---|---|---|---|---|---|
| v12.0, names as in the R code (reported) | 5 | 16 | 29 (24.2%) | 70 (58.3%) | 34 (28.3%) |
| v12.5, names as in the R code | 5 | 16 | 30 (25.0%) | 69 (57.5%) | 35 (29.2%) |
| v12.0, with the Glut2 alias | 5 | 16 | 33 (27.5%) | 66 (55.0%) | 38 (31.7%) |
| v12.5, with the Glut2 alias | 5 | 16 | 34 (28.3%) | 65 (54.2%) | 39 (32.5%) |
Both releases have 59 distinct pairs among the 21 genes; among the 17 analyzed genes 38 in 12.0 and 39 in 12.5. Three pairs are gained in 12.5 (Kcnj11-Ryr2, Wnt5a-Wnt9a, Wnt5b-Wnt9a) and two are lost (Cacna1c-Itpr1, Cacna1d-Itpr1); each changes one pair between STRING Only and Not Validated; no consensus pair changes class. **Reading: the STRING release does not change the findings; the 5 known and 16 database-unregistered consensus pairs are the same in 12.0 and 12.5.**

## 3. A correction candidate in the reported analysis: Glut2 and Slc2a2
The correlation matrix names the gene Glut2; STRING (both releases) names it Slc2a2. The R code (`Fig_5abcd_Gene_Interaction_Network_Analysis.r`) matches names literally, so the 4 pairs of Glut2 that STRING holds (Glp1r-Glut2, Glut2-Ins1, Glut2-Kcnj11, Glut2-Tcf7l2; 15 of the 120 pairs involve Glut2) were counted as absent from STRING. With the alias, STRING Only is 33 instead of 29, Not Validated 66 instead of 70, and STRING coverage 38 of 120 (31.7%) instead of 34 of 120 (28.3%) in 12.0. Known Validated (5) and Novel Discovery (16) do not change. Affected in the manuscript: the numbers 29, 70, 24.2%, 58.3% and 28.3% (34 of 120) in the Results and the Discussion, and Fig. 5 panels a to c. This is a naming error, not a tuning choice; the author decision is in section 6.
**Decision (recorded 9 Nov 2026).** The user chose "correct it" on 8 Nov 2026 and approved the corrected Fig. 5 on 9 Nov 2026 at 00:26 +0700 (`fig5_corrected/APPROVED.txt`, which holds the checksums of the approved files). The corrected numbers are therefore the reported numbers of the manuscript from now on: STRING Only 33 (27.5%), Not analytically consistent 66 (55.0%), 38 of 120 pairs in STRING (31.7%) with STRING 12.0 as used, and 39 of 120 (32.5%) with STRING 12.5; Known Validated 5 and Novel Discovery 16 do not change. The run as originally reported (literal names; 29 / 70 / 34 of 120) stays in `fig5_corrected/as_reported/` and in `R_17genes/` as the record of the earlier version; the corrected figure and its data exports are in `fig5_corrected/composite_alias/` and `fig5_corrected/alias/`. The only change of the analysis is one line of the R script (the STRING node name Slc2a2 read as Glut2, marked `ALIAS FIX`); STRING's own identifier mapping confirms that Glut2 is the only one of the 17 names that differs (`results/string_id_mapping_rat_v12_5.tsv`).

## 4. Part B: MLSS with the DrugCentral release of 25 Sep 2026
Gate: the port on the old file reproduces the reported scores (153 of 153, rank order 153 of 153): PASS. Then the same port, the same rules and the same inputs, with the new interaction file:
* **Drugs:** 18 become 19; the only addition is **tirzepatide** (GLP1R drugs 7 become 8; CACNA1D drugs stay 11); no drug is removed; TCF7L2 and PTBP1 have no drug in either release.
* **Combinations:** 153 become 171; cross-target pairs 77 become 88 (tirzepatide with each of the 11 calcium channel blocker entries).
* **Shared combinations:** all 153 keep their score (0 changed; Spearman rho 1.0 for scores and for ranks). The only pIC50 that differs is the new value for tirzepatide at GLP1R (8.37); no existing value changed.
* **Top 15:** **14 of the 15 reported combinations stay in the top 15**; ranks 1 to 12 are unchanged; **isradipine + tirzepatide enters at rank 13 (MLSS 9.9158)** and the old ranks 13, 14 and 15 become 14, 15 and 16 (the old rank 15, albiglutide + nifedipine, drops to 16). No tie group straddles rank 15 in either ranking.
* **Tirzepatide combinations:** 18 (11 cross-target at ranks 13, 26, 27, 30, 33, 69, 70, 71, 72, 73 and 74; 7 same-target at rank 99 or lower).
**Reading: the ranking is robust to the DrugCentral release (the hub-target data of the existing drugs did not change); the reported run lacks tirzepatide, which adds one pair to the top 15.** The tier-1 approval list (`FDA+EMA+PMDA_Approved.csv`) is byte-identical to the one of the reported run; the FDA-only list (tier 3) was not downloaded and is not re-evaluated.

## 5. Limits
Only the interaction file and the triple-approval list were replaced; the INS-1 matrix, the hubs (degree >= 5), the correlation lookup and the weights are as reported, so the sensitivity says what the database release changes, nothing else. The reported file carries no version stamp, so which release it was is not known; the hub-target content differs from the 2026 file only by tirzepatide, which fits an older release. STRING: one species, one threshold (0.4), the functional network only; edges come from the API, not from the web export, and the reported export has no stamp either. The safety screen (WP9) covers the 18 reported drugs; tirzepatide (rank 13 with isradipine) is not screened (optional, post hoc, about 100 requests).

## 6. Items for the authors (omitted)

Omitted from this copy: this section listed points for the authors that were settled while the article was prepared.

## 7. Proposed wording (omitted)

Omitted from this copy: this section proposed wording for the article; the wording used is that of the article.

## 8. Deviations (dated)
* 8 Nov 2026, about 01:20: the first run of `drugcentral_sensitivity.py` stopped after the gate because of a pandas indexing error in the comparison code (frozenset labels); the comparison code was changed to dictionaries and re-run; no rule, input or output definition changed.
* 9 Nov 2026: the user's decision on the Glut2 / Slc2a2 correction and the approval of the corrected Fig. 5 were recorded (sections 3, 6 and 7 updated, this entry added); the user's own STRING 12.5 web export of 8 Nov 2026 was compared with the API network of part A (section 10); no number, rule, input or output of the analysis was changed.

## 9. Re-running and checks
`python -I string_sensitivity.py gate | retrieve | analyze`; `python -I drugcentral_sensitivity.py` (needs pandas and scipy; slow on a loaded machine); `python -I check_wp10_numbers.py` (own csv parser, no pandas: rebuilds the STRING classification from the raw files, the DrugCentral counts and the tirzepatide potency from the raw files, compares with the result files and with the numbers of this README; result in `check_wp10_last_run.txt`).

## 10. Cross-check with the user's own STRING 12.5 web export (8 Nov 2026, added by WP11)
The user's web export (`Database/STRING12.5/`, file times 8 Nov 2026 13:16 to 13:23; network `https://version-12-5.string-db.org/cgi/network?networkId=beSoaw2hyqcw`, saved in `STRING12.5.txt`) holds 118 rows = 59 distinct pairs among the 21 genes. They are the same 59 pairs as in the API network of part A, with identical combined scores (largest difference 0.000), and the Fig. 5 classification from the web export equals the one in section 2 (literal names 5 / 16 / 30 / 69 with 35 pairs in STRING, 29.2%; with the Glut2 alias 5 / 16 / 34 / 65 with 39, 32.5%). Against the old export of Fig. 5: 55 pairs are in both files and 54 of them have a different combined score, so the old file is not a 12.5 export; 4 pairs are gained in 12.5 (Kcnj11-Ryr2, Wnt2-Wnt9a, Wnt5a-Wnt9a, Wnt5b-Wnt9a) and 4 are lost (Cacna1c-Itpr1, Cacna1d-Itpr1, Hey1-Lef1, Itpr1-Ryr1); Wnt2, Hey1 and Ryr1 are not among the 17 analyzed genes, which is why section 2 counts 3 gained and 2 lost. That the old file is release 12.0 stays the authors' statement (no version stamp; download date not recorded, copy date 2 Nov 2026). Script `string_web_export_check.py`, result `results/string_web_export_check.json`; its numbers are also checked by `../wp11_plate_provenance/check_wp11_numbers.py`.
