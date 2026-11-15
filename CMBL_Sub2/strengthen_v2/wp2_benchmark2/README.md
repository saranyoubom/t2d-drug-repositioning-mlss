# WP2: benchmark part 2 (protocol v2, section 3)

Status 5 Nov 2026: complete except for the mouse RNA-seq arm (not done, see 6). Protocol and fingerprints: `../PROTOCOL_v2_prespecified.md`, `../PROTOCOL_v2_SHA256.txt`. Progress log with times: `logs/progress.txt`.

## 0. Who did what
The WP2 agent wrote the scripts `bench2_core.py`, `bench_data_io.py`, `gate_part1.py`, `calibrate_t2.py`, `select_heldout.py`, `select_rnaseq.py`, `endpoints_lib.py`, `apply_dev.py`, `design_table.py`, `summarize_validation.py`, `run_heldout.py`, `ncbi_http.py` and ran gate A, then stopped twice on the weekly limit. The main session (planner) fixed two import-order bugs (`calibrate_t2.py`, `design_table.py`), wrote `run_rnaseq.py` and `make_wp2_figure.py`, ran everything from the calibration on, and made the memory-lean loaders (section 8). Two background-job sets were stopped by the system for low memory (logged); nothing was lost.

## 1. Gate A (deliverable A): PASS
All 37 analysed development series reload (58 of 58 rows identical in status, DE counts, gene lists) and `bench_pipeline.analyse` reproduces part 1 exactly for gene sets A and B: 74 of 74 (observed count, null mean difference 0, null SD, P). `gate_part1_reproduction.txt` (section 3 compares the new module `coexp_checks`: identical observed values in 71-73 of 74 matrices; the differences are confined to the three matrices with degenerate genes, section 3).

## 2. Calibration (deliverable B, protocol 3.2): FAIL as pre-specified
37 series x 340 replicates (100 without injection; r in 0.3, 0.5, 0.7, 0.9 x m in 3, 6, 10, 20 each), 500 permutations; `calibration_verdict.txt`, `calibration_tables.csv`, `calibration_replicates.csv`.
* T_r2 type I error 276/3700 = 7.5% (Wilson 6.7-8.4%), limit 6%: FAIL. Power at r = 0.7, m = 6: 671/740 = 90.7% (88.4-92.6%): meets 80%. Overall FAIL; the secondary statistic named by the fixed rule is T_p (type I 5.8%, power 90.8%).
* Diagnosis (post hoc): one series, GSE20677, rejects 100/100 (13 of its 17 set-A genes are constant in all 12 samples; the part 1 ANOVA flags them DE because within-group SS is exactly 0 and between SS about 4e-29, F = inf; decoupling cannot change them). Without GSE20677 and GSE28477 (both with degenerate genes): T_r2 type I 175/3500 = 5.0% (4.3-5.8%), power 90.9%; T_p 3.3% / 91.1%; T_eig 4.9% / 90.7%; count test 4.3% / 29.3% (`calibration_sensitivity_no_degenerate.csv`).
* Power depends on how many genes are coupled: T_r2 with 3 genes at r = 0.9 only 35%; 6 genes: 13% (r 0.3), 43% (0.5), 91% (0.7), 100% (0.9); 10 genes: 30%, 95%, 100%, 100%.
* The consensus-count test: type I 4.1%; power 28.2% (r = 0.7, 6 genes), 67.6% (r = 0.9); 13.4% at r = 0.5. This agrees with part 1 (27%).

## 3. Development series (deliverable C, protocol 3.3, 3.4): `dev_application_summary.md`
Set A (17 random DE genes) / set B (17 random expressed genes), 37 series: T_r2 coupling in 25/37 = 68% (51-80%) / 29/37 = 78% (63-89%); consensus count not above null (E1) 32/37 = 86% / 24/37 = 65% (part 1: 81% / 68%); E2 on the raw scheme alone median 0%; group means carry 98% [96, 99] / 69% [42, 86] of covariance. Without degenerate series: 23/35 = 66% / 28/35 = 80%. Degenerate genes (no within-group variation) occur in 2 of 37 series: GSE20677 (A 13, B 3 genes) and GSE28477 (A 6, B 4); part 1 values for them differ from the new module (see fix list A17).

## 4. Held-out arrays (deliverable D)
`select_heldout.py`: same query and shuffle as part 1 (positions 1-150 reproduced exactly), window 151-450, series sharing a sample with the 37 development series excluded; 53 GPL570 and 82 GPL1261 eligible; list frozen 5 Nov 00:56:35, sha256 71ae9a1e... in the fingerprint file, before any expression value was read. Dry run on DE counts, then the full run with the same code path (0 mismatches): 20 analysable per platform (GPL570 to list position 377, GPL1261 to 245). Not analysable, with reasons in `heldout_arrays_results.jsonl`: fewer than 34 DE genes (GPL570 12, GPL1261 4), selected samples missing in the matrix (1), and two empty gene matrices: GSE294689 (series matrix without expression rows) and GSE39840 (upper-case probe ids that match none of the 54,675 platform ids exactly; the exact mapping of part 1 is kept). No held-out matrix has degenerate genes.

## 5. Validation endpoints (deliverable F, protocol 3.6): `validation_summary.md`, `WP2_summary_table.md/.csv`
All 59 validation sets (40 held-out arrays, 19 human RNA-seq), T_r2 coupling in 47/59 = 80% (68-88%) for set A and 48/59 = 81% (70-89%) for set B (binomial P vs 5% below 1e-49); held-out arrays versus development, Fisher P = 0.3 (A), 1 (B). Per arm, A: held-out GPL570 85%, GPL1261 75%, RNA-seq 79%. E1 (count not above null): arrays 85%, RNA-seq 58% (A). Group means carry 97-98% of covariance among DE genes. After treatment adjustment the raw scheme keeps a median 0% of consensus pairs (set A); ratio-scheme adjusted counts sit at their own null (median observed/null 0.96-1.42 for A across arms), so they cannot be read as evidence. Share of raw pooled-significant pairs among group-effect genes 0.73 [0.63, 0.88] against the geometric expectation 0.61.

## 6. RNA-seq arm (deliverable E, protocol 3.5)
Query as in the protocol (retrieved per n_samples because it returns more than 10,000 records), shuffled with `random.Random(2026)`, first 300 per organism screened on metadata: 65 human and 109 mouse eligible; list frozen 5 Nov 01:09:49 (sha256 39ec5dc8...), before any count was read. NCBI-generated counts: human, 10 of 20 probes outside the window returned a file; the list of 65 gave 19 analysable (28 without counts, 10 with selected samples absent from the count file, 8 with fewer than 34 DE genes) and the list was exhausted below the target of 20. **Mouse: HTTP 404 for 109 of 109 and for 20 of 20 probes: arm not done** (rule of 3.5). The GEO documentation page and the annotation file answer with a reCAPTCHA page; this was not worked around and is not needed (gene sets are random draws; Entrez ids are the gene names). `rnaseq_pattern_check.csv`.

## 7. Design table (deliverable F, protocol 3.7): `design_table.md/.csv`
G = 17, one module of 6 genes, 3 groups, 300 independent simulated experiments per cell. With 4 per group (this study) T_r2 power at r = 0.5 is 0.39 and at r = 0.7 0.87, but the share of true module pairs detected at BH q < 0.05 is 0.01 and 0.06: the best single pair needs |r| >= 0.90 (residual df 8). Pair detection of r = 0.5 needs 12 per group (0.52), 20 for 0.95. Wells from one preparation are not independent experiments.

## 8. Deviations, additions, code changes (all also in protocol section 8 or `logs/progress.txt`)
* ADD-ON (post hoc, logged before any held-out coupling statistic): `n_zero_var_genes` and a sensitivity analysis without degenerate series beside every result.
* Mouse RNA-seq arm not done; human arm 19 of the target 20.
* Code only: import-order fixes; `run_rnaseq.py`; lean loaders (`read_matrix`, `annotation` stream the table; verified bit-identical on all 74 development matrices, gene lists and DE counts; extraction 48 s instead of 431 s; earlier files in `compact/before_lean/`).
* The part 1 results for two series are affected by constant genes (fix list A17).

## 9. Re-running (one heavy job at a time on the 8 GB machine; data root = environment variable V2_DATA_ROOT)
`python -I gate_part1.py extract | compare 2 | consistency | report`; `calibrate_t2.py run 2` (25 min) then `summarize`; `apply_dev.py`; `select_heldout.py --download-only`, then no flag (freezes), `run_heldout.py arrays dry`, `arrays full`; `select_rnaseq.py lists | pattern-check | screen human,mouse`, `run_rnaseq.py dry | full`; `summarize_validation.py`; `design_table.py run 2 300` then `summarize`; `make_wp2_figure.py`. Seeds: draws 2026 + GSE number (set B [2026 + n, 1]); tests seeded 2026 + n; calibration SeedSequence([2026, n, condition, replicate]); design SeedSequence([2026, 37, n, r index, replicate]); selection shuffles random.Random(2026).

## 10. Numbers for the manuscript
Calibration table (section 2) with the verdict stated as failed-as-pre-specified plus the diagnosis; the count test's 28% (68% at r = 0.9); prevalence 80% (47/59) and 81% (48/59); covariance share 98%; ratio-scheme counts at their null; the design table; limits: log-scale arrays and CPM only, one design (3 x 4), random gene sets, mouse RNA-seq absent. Draft figure: `WP2_benchmark2_figure_draft.png/.pdf` (symbols to be explained in the legend, no key strip).

## 11. Author input needed
None for WP2.
