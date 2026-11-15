# Reproducibility bundle for the v2 analyses (6 Nov 2026; for the repository folder, NOT pushed)

What this is: the scripts, cached inputs, outputs, seeds and fingerprints behind the version-2 work on the INS-1 co-expression pipeline (protocol `../PROTOCOL_v2_prespecified.md`, plan `../../../PLAN_v2_strengthen_JTM.md`): WP1 INS-1 structure and calibrated checks, WP2 benchmark part 2, WP3 MLSS port and scenarios, WP4 public islet cross-check version 2, WP5a effect sizes and MIQE checklist, WP5b this bundle and Table 1. Every package folder has its own README with the results; this file says how to run things and what is needed. Nothing has been committed, pushed, released or uploaded.

## 1. Environment
* Windows 11, 8 GB RAM (about 1 GB free during the runs), project inside a OneDrive folder (sync load slows runs by a factor of 2 to 4). One heavy job at a time.
* Python 3.14.4 with the packages in `requirements.txt` (exact versions read from the installed packages). Always run `python -I`; scripts add their own folder to `sys.path`. The default `python` on this machine is 3.13 without numpy; use `python` or `py -3.14`.
* R is not used by any v2 script (`R_environment.md` records the R side of the authors' scripts).
* `V2_DATA_ROOT`: folder with the downloaded data (default in the scripts: the session scratch folder). Sub-folders: `bench_data/` (development series, 270 MB), `v2_data/arrays_heldout/` (276 MB), `v2_data/rnaseq/` (56 MB), `v2_data/wp4/` (376 MB: islet downloads and `parsed/`, 34 MB), `gse81608/` (36 MB), `gse84133/` (69 MB, including the parsed per-donor files). Whole data root: about 1.1 GB. Downloaded files are untrusted data: only read, never executed.
* Inputs outside the project: the authors' calculation sheets and exports (`Undergrad_Glucose uptake signaling\RT-PCR\`, read only; environment variables `RTPCR_DIR`, `EXPORT_DIR`, `INS1_MATRIX`, `BENCH_DIR` override the default paths), and the INS-1 matrix in `SBA\Figures\Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx`.

## 2. Order and run times (this machine)
`python -I run_all_v2.py --dry-run` lists every step and checks every input and every script; `--only wp1,wp2,wp3,wp4,wp5a`, `--skip-heavy` and `--list` select; steps run in sequence, each in its own folder.

| order | package | step | writes | run time here | heavy |
|---|---|---|---|---|---|
| 1 | WP1 | `build_gene_run_map.py` (about 10 s), `check_sheet_vs_exports.py` (about 10 s), `wp1_ins1_analysis.py` | `gene_run_map.csv`, `results_ins1.json`, `tbl_2_*.csv` | about 4 min (README of WP1); test suite `test_coexp_checks.py` 56 checks | yes |
| 2 | WP3 | `mlss_port.py` gate, `direction_pic50.py`, `mlss_scenarios.py` | `gate_scores.csv`, `mlss_scenarios.csv`, `mlss_permutation_draws.csv` | scenarios about 5 min (the first run was stopped once for low memory) | yes |
| 3 | WP5a | `effect_sizes.py` | `effect_sizes.csv/.md`, `gate_reproduction.json` | about 1 min | no |
| 4 | WP2 | `gate_part1.py`, `calibrate_t2.py run 2` (25 min) then `summarize`, `apply_dev.py`, `select_heldout.py` (downloads), `run_heldout.py arrays dry` and `full`, `select_rnaseq.py`, `run_rnaseq.py`, `design_table.py run 2 300` (10 min), `summarize_validation.py`, `make_wp2_figure.py` | `calibration_*`, `validation_*`, `design_table.*`, figure draft | see `../wp2_benchmark2/README.md` section 9 | yes (calibration, held-out runs) |
| 5 | WP4 | downloads and parsers (README of WP4, section 14), `run_datasets.py`, `meta_analysis.py`, `gene_graph_null.py`, `metacell.py`, `perturbation_gse60803.py`, extension: `positive_control.py` (about 15 min) and `wider_search_analyze.py` (the live searches `wider_search.py`, `wider_search_meta.py`, `wider_search_freeze.py` are not re-run: their results are frozen and fingerprinted), `make_tables.py`, `make_fig7_v2.py` | `results/`, `meta/`, `metacell/`, `perturbation/`, `WP4_summary_tables.md`, `figures/` | cell-level analysis 6 to 21 s per dataset; meta-analysis seconds; gene-graph add-on about 2 min; metacell about 20 to 25 min; perturbation under 1 min | metacell only |
| 6 | WP5b | `make_inventory.py`, `run_all_v2.py` | `script_inventory.md` | seconds | no |

## 3. Large downloads (plain HTTPS from NCBI and EBI; at most 3 requests per second; nothing about the user is sent)
| need | source | size | used by |
|---|---|---|---|
| development series of the benchmark (37 GEO series matrices and platform files, cached) | NCBI GEO | 270 MB cached in `bench_data/` | WP2 (and `benchmark_route3/`) |
| held-out arrays (40 analysable of 300 screened) | NCBI GEO series matrices | 276 MB in `v2_data/arrays_heldout/` (659 files) | WP2 |
| held-out RNA-seq (19 analysable of 65 human series; mouse not available) | NCBI-generated raw counts (`geo/download/?type=rnaseq_counts`) | 56 MB in `v2_data/rnaseq/` (637 files) | WP2 |
| GSE81608 RPKM table and series matrix | GEO | 37.2 MB + 44 KB | WP4 gate and analysis |
| GSE84133 RAW.tar (parsed per donor by `../../islet_GSE84133/load_beta_counts.py`) | GEO | 29.1 MB | WP4 gate and analysis |
| GSE86469 RSEM expected counts | GEO | 12.7 MB | WP4 |
| GSE81547 RAW.tar (2544 per-cell files) | GEO | 176 MB | WP4 |
| E-MTAB-5061 processed.1.zip (the 723 MB text is streamed from the zip, never written out) | EBI | 159 MB | WP4 |
| GSE60803 series matrix and GPL11532 annotation (the 2 GB platform family file is not used) | GEO | 1.4 MB + 5.6 MB | WP4 |
| Ensembl ids of the 17 genes; publications attached to the dataset records; GEO DataSets queries | NCBI Datasets, E-utilities | under 1 MB | WP4 |
No dense single-cell matrix is ever loaded; the loaders stream and keep sparse matrices.

## 4. Random seeds (all fixed in the code; a different seed changes numbers in the last digits, never a verdict reported in the READMEs)
| package | seed | used for |
|---|---|---|
| `coexp_checks.py` (WP1, WP2) | 2026 (argument `seed`) | permutation nulls: 2000 permutations (consensus count, adjusted counts, coupling test), 20000 for the run-structure test, 10000 where WP1 says so |
| WP1 `wp1_ins1_analysis.py` | 2026 (`SEED`) | 20000 label permutations (2.2, 2.3), 2000 adjusted-count permutations (2.5), 10000 for 2.6 and 2.8, 5000 simulations of the 26 August layout |
| WP2 | 2026 + series number (`bench2_core.py`, `gate_part1.py`); `SeedSequence([2026, n, condition, replicate])` (`calibrate_t2.py`); `SeedSequence([2026, 37, n, ri, rep])` (`design_table.py`); `random.Random(2026)` shuffles the GEO accession lists (`select_heldout.py`, `select_rnaseq.py`); gene sets A and B drawn with `default_rng(2026 + n)` and `default_rng([2026 + n, 1])` | calibration, design simulation, held-out selection, random gene sets |
| WP3 `mlss_scenarios.py` | 2026 (`SEED`), 1000 draws (`NPERM`) | potency and target-mapping permutations |
| WP4 `wp4_core.py` | 2026 (main matched null), 2027 (pair-level null), 2000 draws | matched random pairs (`gate_human_islet_check_lean.py` uses the same two seeds) |
| WP4 `metacell.py`, `metacell_genes.py` | `default_rng(5000 + i)`, i = 0 to 19 | metacell partitions |
| WP4 `gene_graph_null.py` | 2028, 2000 draws | post hoc gene-graph null |
| WP4 `positive_control.py` | `SeedSequence([2029, dataset, rho, mode, replicate])`; ribosomal control `default_rng(2030 + dataset)` and `SeedSequence([2030, dataset, decile, decile])` | spike-in replicates; random ribosomal pair sets and null cells |
| WP4 `meta_analysis.py`, `perturbation_*.py` | none | deterministic |
Test files use their own small seeds (`test_coexp_checks.py`).

## 5. Pre-specification files and fingerprints (verified 6 Nov 2026 with sha256 of the files against `../PROTOCOL_v2_SHA256.txt`)
| file | fingerprints recorded | latest matches the current file |
|---|---|---|
| `../PROTOCOL_v2_prespecified.md` | 11 (4 Nov 13:39 frozen, then one per section 8 update by WP1, WP2 (three), WP4 (three) and the WP4 extension (three)) | YES |
| `../wp2_benchmark2/frozen_heldout_arrays.csv`, `frozen_heldout_rnaseq.csv` | 1 each (frozen before any expression value was read) | YES |
| `../wp4_islets/frozen_islet_datasets.csv`, `frozen_beta_cells.csv` | 1 each (frozen 6 Nov 00:17 from metadata only) | YES |
| `../wp4_islets/PRESPEC_WP4_choices.md` | 6 (13:40 fixed, addenda F, F2, a wording fix, G, H) | YES |
| `../wp4_islets/PRESPEC_WP4b_positive_control_and_wider_search.md` | 1 (6 Nov 16:10, before any analysis of the extension) | YES |
| `../wp4_islets/wider_search/frozen_candidates.csv`, `candidates_all.csv`, `queries.json` | 1 each (6 Nov 16:40, before any value of a new candidate was read) | YES |
| `../wp4_islets/ensembl_ids_17genes.csv` | 1 | YES |
Mismatches: none. The first protocol fingerprint (4 Nov 13:39:33) belongs to the text before any section 8 entry; that version was not kept as a file, so only the latest fingerprint can be re-checked; every earlier one is followed by the section 8 entry that explains the change. The benchmark part 1 protocol has its own fingerprint file (`../../benchmark_route3/PROTOCOL_SHA256.txt`).

## 6. Checks run on this bundle
* `run_all_v2.py --dry-run`: all steps listed, all inputs and scripts present, with and without `V2_DATA_ROOT`.
* Light steps for real (`--only wp5a,wp4 --skip-heavy`, log `run_all_light_test.log`, comparison `run_all_light_test_compare.txt`): exit 0 for all 7 steps. Times on the throttled CPU (about 22% of full speed at that moment): `effect_sizes.py` 80 s, `run_datasets.py` 556 s (72 s earlier the same day at full speed), `meta_analysis.py` 52 s, `gene_graph_null.py` 122 s, `perturbation_gse60803.py` 52 s, `make_tables.py` 17 s, `make_fig7_v2.py` 84 s. Of the 99 output files (WP4 results, meta, perturbation, figures, summary tables and the WP5a effect sizes) 97 are identical byte for byte to the files before the re-run; the 2 that differ are the figure PDFs (they embed their creation time). Not re-run: WP1 (about 4 min, `wp1_ins1_analysis.py`), WP3 scenarios (about 5 min), the WP2 steps and the WP4 metacell run (heavy; WP1 and WP3 have tests and gates of their own, see their READMEs).
* `../wp4_islets/check_readme_numbers.py`: 59 headline numbers of the WP4 README rebuilt from the result files, 0 failures.
* Script inventory: `script_inventory.md` (68 scripts, generated from their docstrings by `make_inventory.py`).
