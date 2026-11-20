# strengthen_v2: version-2 analyses for the revised manuscript (November 2026)

Scripts, pre-specification, fingerprints and results of the version-2 analyses that accompany the manuscript *Notch and Wnt perturbations account for pooled co-expression networks in INS-1 cells: an exploratory in vitro and computational study*.
They extend the 17-gene analysis of the folders next to this one (`../R_17genes/`, `../R_panels_authors/`, `../reanalysis_17genes/`, `../primer_check/`, `../additional_files/`). The INS-1 data are one cell preparation (12 observations): every INS-1 analysis here is post hoc and exploratory.
The public-data analyses carry the name of their dataset. Nothing here predicts toxicity, efficacy or synergy, and the numbers of the manuscript and its Additional files are authoritative.

## Pre-specification and fingerprints
* `PROTOCOL_v2_prespecified.md` is the written protocol; its section 8 is the log of every deviation. Its first version was fingerprinted on 4 November 2026 at 13:39:33 (+0700) before the analyses were run; every later version has its own line in `PROTOCOL_v2_SHA256.txt`.
* `PROTOCOL_v2_SHA256.txt` and the lists of the work packages (`FINGERPRINTS.txt`, `INPUT_FINGERPRINTS.txt`, `results/RAW_SHA256.txt`) record the SHA-256 of the pre-specification files, of frozen lists and inputs, and of the saved raw responses, with the clock time of the record.
* Files whose SHA-256 is recorded in these lists and still equal it are copied into this repository byte for byte. `python -I verify_fingerprints.py` recomputes them (MATCH / CHANGED / NOT SHIPPED). The protocol file itself changes after each deviation entry, so only its last recorded hash matches.

## Folders
| Path | Content |
|---|---|
| `PROTOCOL_v2_prespecified.md`, `PROTOCOL_v2_SHA256.txt` | protocol with its deviation log; fingerprints |
| `coexp_checks.py`, `test_coexp_checks.py` | calibrated within-group coupling tests with a permutation engine; their tests (about 3 minutes) |
| `recon_run_structure.py` | planning reconnaissance on the INS-1 data that the protocol discloses |
| `wp1_ins1_structure/` | WP1: run structure of the INS-1 data (gene-to-run map, run-structure statistic, calibrated tests), post hoc |
| `wp2_benchmark2/` | WP2: benchmark of the checks in public expression data (37 development series, 40 held-out arrays, 19 human RNA-seq series), with the frozen lists |
| `wp3_mlss/` | WP3: independent Python port of the MLSS (v4.0), scenarios, direction table |
| `wp4_islets/` | WP4: cross-check in public human and mouse islet single-cell data, positive control, wider search for perturbation data |
| `wp5_reporting/` | WP5a: effect sizes with intervals and across-gene Benjamini-Hochberg q; environment records (`requirements.txt`, `R_environment.md`, `README_v2.md`) |
| `wp9_safety/` | WP9: annotation of the ranked drugs and pairs with openFDA adverse event reports (not a prediction); saved responses in `results/raw_responses.zip` |
| `wp10_versions/` | WP10: database-version sensitivity (STRING 12.5, DrugCentral release of 25 September 2026) and the Glut2 / Slc2a2 correction of Fig. 5 (`fig5_corrected/`) |
| `wp11_plate_provenance/` | WP11: provenance of the June, August and November RT-qPCR values (plates, thresholds, sorted rows) |
| `wp12_database_annotations/` | WP12: hub-gene annotation, score-quartile classes, United States market status from Drugs@FDA; saved responses in `results/raw_responses.zip` |
| `wp13_gapdh_corrected/` | WP13: Gapdh reference sensitivity of the 26 August run |
| `wp14_table_s8/` | WP14: Table S8 recomputed from the Cq values of the analyzed matrix |
| `ars_review_191126/post_hoc/` | post hoc checks of the benchmark (19-20 November 2026; added in release v2.1.0): the share of the pooled covariance carried by group means against its expectation without within-group coupling, and T_r2 after removal of one shared axis, on the 37 development series, the INS-1 matrix and the 59 validation sets; the method notes were written and fingerprinted before the runs (`FINGERPRINTS.txt`). See its own `README.md`. |
| `fig10_201126/` | code and numbers of Fig. 10 (candidate checks for pooled co-expression networks), added in release v2.1.0; the script imports `wp6a_assets/figstyle.py` and, to check the INS-1 numbers, reads the manuscript file (environment variable `FIG10_MANUSCRIPT`; the manuscript is not part of this repository) |
| `wp6a_assets/figstyle.py` | the figure style module that `fig10_201126/make_fig10.py` imports; the rest of `wp6a_assets/` is not included |
| `../benchmark_route3/` | first part of the benchmark (route 3), which the version-2 benchmark builds on, with its own pre-specification and fingerprints |

## Reading key
The READMEs and pre-specifications are working documents written while the analyses ran, and some of them use working language: "the user" means the authors, "WP" a work package, "main session", "planner" and "agent" working sessions of the AI assistant (Claude Code, Anthropic) that the Methods of the article describe under *Use of AI tools*.
Items named A1, B3, C11 and the like refer to an internal list of manuscript fixes, and `briefs/` and `MANUSCRIPT_FIX_LIST.md` are internal files; none of them is part of this repository. Times are local time (UTC+7).
Sections of four READMEs that only listed points for the authors or proposed wording for the manuscript are omitted where they stood and marked "(omitted)".
The `WORDING_*.md` files are the proposals of wording that the `check_*_numbers.py` scripts verify against the results; the text of the article is authoritative.
In the working documents and in the code, `Analysis_outputs/` is the folder `CMBL_Sub2/` of this repository (the path components in the scripts were renamed accordingly).
`wp13_gapdh_corrected/README.md`, section 4, records the two options as they stood on 9 November 2026; the authors chose option A: the reported analysis stays primary, and the corrected version (C0) and the worst case (C1) are sensitivity analyses.
`ars_review_191126/` is named after a simulated review of the manuscript (written with the AI assistant) that raised the two questions of the post hoc checks; the review itself is not part of this repository.

## Running the code
* Python 3.14.4 with the packages of `wp5_reporting/requirements.txt` (numpy 2.4.4, scipy 1.17.1, pandas 3.0.3, statsmodels 0.14.6, matplotlib, openpyxl); run every script with `python -I`, from its own folder. R 4.4.0 was used for the R scripts.
* The scripts find `../../R_17genes/`, `../../reanalysis_17genes/` and `../../benchmark_route3/` relative to their own folder, so keep this layout. Some scripts write results into their own folder: run them in a copy.
* Inputs that are **not** included: the raw CFX Maestro exports and the calculation sheets of the RT-qPCR runs (environment variables `RTPCR_DIR` and `EXPORT_DIR`); the DrugCentral file of 25 September 2026 and the STRING 12.5 export (a folder `Database/` three levels above `wp10_versions/`); public data from GEO and ArrayExpress, which the `fetch_*.py` scripts download into `data_cache/` (environment variable `V2_DATA_ROOT`).
* The relative-expression matrix of the INS-1 experiment is `../R_17genes/Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx` (environment variable `INS1_MATRIX`).
* The validation mode of the post hoc checks (`ars_review_191126/post_hoc/post_hoc_validation_191126.py DATA_DIR`) needs 61 files from NCBI GEO (about 210 MB): `--list` shows them without network access, and nothing is downloaded unless `--allow-download` is given.
* Before running `check_wp9_numbers.py`, unpack `wp9_safety/results/raw_responses.zip` into `wp9_safety/results/raw/` (see `wp9_safety/results/README_raw_responses.md`).
* The `check_*numbers.py` scripts rebuild the numbers of the READMEs and WORDING files from the result files with their own parsers. `check_wp9_numbers.py` (after unpacking the raw responses; 120 checks), `check_wp13_numbers.py` (17 checks), `check_wp4b_numbers.py` and `test_coexp_checks.py` (56 checks) run on the files of this repository and passed in the repository layout (smoke test of 15 November 2026). The other checks need inputs that are not included: Additional file 1 of the manuscript (`check_wp14_numbers.py`), Additional file 4 (`check_effect_sizes.py`), the DrugCentral file of 25 September 2026 and the STRING 12.5 export in `Database/` (`check_wp10_numbers.py`, `check_wp12_numbers.py`) and the raw instrument exports (`check_wp11_numbers.py`); the first lines of each script say what it reads.

## Licenses and data
The code is MIT-licensed (see the LICENSE of the repository). Tables derived from DrugCentral and STRING are included for reproducibility under the licenses of those sources; please cite the original databases. The openFDA responses are public data, which openFDA describes as unvalidated and not to be used for decisions about medical care. GEO and ArrayExpress data are not redistributed here.
