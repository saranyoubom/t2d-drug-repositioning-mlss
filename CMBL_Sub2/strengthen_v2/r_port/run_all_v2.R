## run_all_v2.R - R counterpart of wp5_reporting/run_all_v2.py: re-creates the v2 tables, checks and figures by calling each package's own R script (no logic is copied here). One heavy job at a time: steps run in sequence, never in parallel.
## Usage: Rscript run_all_v2.R [--source archive|R] [--dry-run] [--only wp1,wp3,...] [--skip-heavy] [--list]
##   --source archive (default): the verification run of the port. The summary steps of route 3 and WP2 re-summarise the archived PYTHON results and Additional file 6 is built from the
##                     archived results, so every R output can be compared with the archive by agreement_report.R (the gate of R_PORT_PLAN.md, calibrated on this mode).
##   --source R        the pipeline from the data, in R: downloads (GEO, EBI), the full route-3 run, the WP2 data runs (extract, part-1 gate, core analysis, held-out arrays, RNA-seq,
##                     calibration, design table), the post hoc checks and the islet analyses; the summary steps and Additional file 6 read the R results. About 28 h on the 8 GB laptop
##                     with one worker (calibration about 8 h and design table about 6 h of it; run those two by hand with 2 workers to halve them). Where R differs from the archive
##                     (the known differences of agreement_known_reasons.tsv, e.g. the BLAS-noise series), this mode carries the R value into the summaries.
## Paths are relative to this file. Each result goes to results_R/<package>/ (the archived Python results of strengthen_v2/ are never written); compare with agreement_report.R. --dry-run only checks
## the scripts and the repository inputs and prints what would run; inputs that an earlier step writes (R:results_R/...) are checked when their step starts.
## Inputs marked R: are relative to this folder; the others are files of the repository (relative to strengthen_v2/). Downloads (GEO, NCBI, EBI) are steps of their own and are cached in results_R/*_data/;
## the live queries of WP2 (selection), WP4 (scouts) and WP9/WP10/WP12 (retrieve) are NOT re-run: their results are frozen and fingerprinted, the analyses read the frozen lists and responses.
## Changed 27 Nov 2026 (fix F2 of audit_271126/AUDIT_REPORT_271126.md; the previous version is run_all_v2.R.bak_271126): a clean checkout no longer stops before step 1 (inputs written by earlier
## steps are checked when their step starts); WP10 reads the frozen STRING 12.5 retrieval of the repository (wp10_versions/results); the downloads are steps; the islet data step runs before the
## islet chain; --source R added.
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); V2 <- normalizePath(file.path(R_PORT, ".."))
RS <- Sys.getenv("RSCRIPT", file.path(R.home("bin"), "Rscript.exe")); if (!file.exists(RS)) RS <- file.path(R.home("bin"), "Rscript")
## each step: package, script (R or shell) and arguments, inputs that must exist, what it writes, rough minutes on the 8 GB laptop, heavy?, mode ("both", "archive" or "R"), environment variables, stdout log
S <- function(pkg, cmd, inputs, out, mins, heavy = FALSE, mode = "both", env = character(0), log = NULL, optional = FALSE) list(pkg = pkg, cmd = cmd, inputs = inputs, out = out, mins = mins, heavy = heavy, mode = mode, env = env, log = log, optional = optional)
## optional steps are skipped (with a message) when a repository input is missing: the Additional file 6 workbook of the article is not in the repository
WP10_FROZEN <- file.path(V2, "wp10_versions", "results")      # the STRING 12.5 retrieval of 8 Nov 2026 (string_v12_5_network.tsv, string_version.json); the live retrieve of the R port gave the identical file
## saved openFDA responses: a folder results/raw in the project, a packed raw_responses.zip in the repository (unpack_openfda.R unpacks it into results_R/openfda_raw)
RAW9_REPO <- file.path(V2, "wp9_safety", "results", "raw"); RAW12_REPO <- file.path(V2, "wp12_database_annotations", "results", "raw")
RAW9 <- if (dir.exists(RAW9_REPO)) RAW9_REPO else file.path(R_PORT, "results_R", "openfda_raw", "wp9_safety"); RAW12 <- if (dir.exists(RAW12_REPO)) RAW12_REPO else file.path(R_PORT, "results_R", "openfda_raw", "wp12_database_annotations")
IN9 <- if (dir.exists(RAW9_REPO)) "wp9_safety/results/raw" else "R:results_R/openfda_raw/wp9_safety"; IN_OPENFDA <- if (dir.exists(RAW9_REPO)) character(0) else c("wp9_safety/results/raw_responses.zip", "wp12_database_annotations/results/raw_responses.zip")
AF6_TEMPLATE <- "wp6a_assets/Additional file 6_JTM_working_201126b.xlsx"
DRUGCENTRAL <- "../../Database/drug.target.interaction.tsv.gz"      # DrugCentral dump of 25 Sep 2026 (SHA-256 1908684983a7..., fingerprinted in Additional file 5): not in the repository; WP10 DrugCentral and WP12 need it
R_CORE <- "R:results_R/wp2/compact/dev_core_results.jsonl"
ENV_WP2_R <- c(WP2_DEV_CORE = file.path(R_PORT, "results_R", "wp2", "compact", "dev_core_results.jsonl"), WP2_DEV_CSV = file.path(R_PORT, "results_R", "wp2", "dev_application_per_series.csv"),
               WP2_HELDOUT = file.path(R_PORT, "results_R", "wp2", "heldout_arrays_results.jsonl"), WP2_RNASEQ = file.path(R_PORT, "results_R", "wp2", "heldout_rnaseq_results.jsonl"))
STEPS <- list(
  S("g0", c("rng/test_numpy_rng.R"), "R:rng/numpy_rng.R", "gate G0: numpy streams bit for bit", 2), S("g0", c("rng/test_pyformat.R"), "R:rng/reference/pyformat_reference.tsv", "gate: Python number parsing and formatting (331 values, 5 formats)", 1), S("g0", c("rng/test_pyrandom.R"), "R:rng/pyrandom.R", "gate G0b: Python random.Random", 1),
  S("g1", c("test_coexp_checks.R"), c("coexp_checks.py", "R:coexp_checks.R"), "gate 1.1: 52 checks of the coexpression tool", 30),
  S("wp1", c("wp1_ins1_analysis.R"), c("wp1_ins1_structure/ins1_sheet_cq.csv", "wp1_ins1_structure/gene_run_map.csv"), "results_R/wp1: tbl_2_*.csv, results_ins1.json", 8, TRUE),
  S("wp3", c("wp3_mlss.R"), c("wp3_mlss/direction_table.csv", "wp3_mlss/pic50_sources.csv"), "results_R/wp3: mlss_scenarios.csv, mlss_permutation_draws.csv", 15, TRUE),
  S("wp5a", c("effect_sizes.R"), character(0), "results_R/wp5: effect_sizes.csv/.md", 1), S("wp5a", c("check_effect_sizes.R"), "R:results_R/wp5/effect_sizes.csv", "independent check of effect_sizes.csv", 1),
  S("openfda", c("unpack_openfda.R"), IN_OPENFDA, "results_R/openfda_raw: the saved openFDA responses of WP9 and WP12, unpacked when the repository holds them packed (WP9 folder digest checked)", 1),
  S("wp9", c("wp9_safety.R", "analyze"), c(IN9, "wp3_mlss/direction_table.csv"), "results_R/wp9: signals, pairs, summary.json (offline, frozen openFDA responses)", 2, env = c(WP9_RAW_ARCHIVE = RAW9)), S("wp9", c("wp9_safety.R", "fallback-analyze"), c(IN9), "results_R/wp9: fallback_signals_single.csv, fallback_pairs.csv, fallback_summary.json (secondary free-text fallback, post hoc)", 1, env = c(WP9_RAW_ARCHIVE = RAW9)),
  S("wp9", c("wp9_check_numbers.R"), c("R:results_R/wp9/summary.json", "R:results_R/wp9/fallback_summary.json"), "120 checks of the WP9 numbers", 2, env = c(WP9_RAW_ARCHIVE = RAW9)),
  S("wp10", c("wp10_string_sensitivity.R", "analyze"), c("wp10_versions/results/string_v12_5_network.tsv", "wp10_versions/results/string_version.json"), "results_R/wp10: string_summary.json, string_classification.csv (frozen STRING 12.5 retrieval)", 1,
    env = c(WP10_RETRIEVED_DIR = WP10_FROZEN)),
  S("wp10", c("wp10_drugcentral_sensitivity.R"), DRUGCENTRAL, "results_R/wp10: drugcentral_summary.json, mlss_comparison.csv (needs the DrugCentral dump of 25 Sep 2026)", 2, optional = TRUE),
  S("wp10", c("wp10_string_web_export_check.R"), c("R:results_R/wp10/string_summary.json", "../../Database/STRING12.5/string_interactions_21_nodes.tsv"), "results_R/wp10/string_web_export_check.json (the STRING 12.5 web export against the frozen retrieval)", 1, env = c(WP10_RETRIEVED_DIR = WP10_FROZEN)), S("wp10", c("wp10_check_numbers.R"), c("R:results_R/wp10/string_summary.json", DRUGCENTRAL), "54 checks of the WP10 numbers (needs the DrugCentral dump)", 1, optional = TRUE, env = c(WP10_RETRIEVED_DIR = WP10_FROZEN)),
  S("wp11", c("wp11_plate_provenance/reconstruct_threshold.R"), c("wp1_ins1_structure/ins1_sheet_cq.csv"), "results_R/wp11/threshold_fit.txt", 5), S("wp11", c("wp11_plate_provenance/check_threshold150.R"), character(0), "results_R/wp11/threshold150_*", 1),
  S("wp11", c("wp11_plate_provenance/check_refA_vs_plate1.R"), character(0), "results_R/wp11/refA_vs_plate1.txt", 1), S("wp11", c("wp11_plate_provenance/check_aug26_nov09_order.R"), character(0), "results_R/wp11/aug26_nov09_order_check.*", 1),
  S("wp11", c("wp11_plate_provenance/export_tables.R"), c("R:results_R/wp11/threshold_fit.txt", "R:results_R/wp11/refA_vs_plate1.txt"), "results_R/wp11/*.csv", 1), S("wp11", c("wp11_plate_provenance/fingerprint_inputs.R"), character(0), "results_R/wp11/INPUT_FINGERPRINTS.txt", 1),
  S("wp11", c("wp11_plate_provenance/check_wp11_numbers.R"), "R:results_R/wp11/export_tables_done", "31 checks of the WP11 numbers", 1),
  S("wp12", c("wp12_annotations.R", "analyze"), DRUGCENTRAL, "results_R/wp12 (offline, frozen openFDA responses; needs the DrugCentral dump)", 1, env = c(WP12_RAW_ARCHIVE = RAW12), optional = TRUE), S("wp13", c("wp13_gapdh_versions.R"), character(0), "results_R/wp13: gapdh_versions.json", 1), S("wp14", c("wp14_effsens_analyzed_matrix.R"), "wp1_ins1_structure/ins1_sheet_cq.csv", "results_R/wp14: Table S8 recomputed", 1),
  S("reanalysis", c("reanalysis_17genes/pipeline.R"), character(0), "results_R/reanalysis_17genes: pipeline_20genes.json, pipeline_17genes.json (the other scripts of that folder are run the same way, see run_log_*.txt)", 5),
  S("wp2data", c("wp2_benchmark2/wp2_prefetch.R"), "wp2_benchmark2/frozen_heldout_arrays.csv", "download of the GEO files of the development, held-out and RNA-seq series into results_R/wp2_data (cached files are skipped)", 60, TRUE,
    log = "R:results_R/wp2_data/prefetch_log_retry.txt"),
  S("route3", c("benchmark_route3/validate_pipeline.R"), character(0), "validation of the route-3 pipeline on the INS-1 data (16 checks)", 4),
  S("route3", c("benchmark_route3/summarize_benchmark.R"), "../benchmark_route3/results_per_dataset.jsonl", "results_R/benchmark_route3: summary files from the archived per-dataset records", 2, mode = "archive"),
  S("route3", c("run_route3_full.sh"), "R:results_R/wp2_data/prefetch_log_retry.txt", "results_R/benchmark_route3_full: the 58 frozen series analysed in R, summarised and compared with the archive", 95, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/gate_part1.R", "extract"), "R:results_R/wp2_data/prefetch_log_retry.txt", "results_R/wp2/compact: dev_sets.json, dev_prepared.json (the 74 development matrices read in R)", 20, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/gate_part1.R", "compare"), "R:results_R/wp2/compact/dev_sets.json", "results_R/wp2/compact/gate_bench_pipeline.jsonl (part-1 gate, 2000 permutations per matrix)", 45, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/gate_part1.R", "consistency"), "R:results_R/wp2/compact/dev_sets.json", "results_R/wp2/compact/dev_core_results.jsonl (coexpression checks on the 74 matrices)", 30, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/gate_part1.R", "report"), R_CORE, "results_R/wp2/gate_part1_reproduction.txt", 1, mode = "R"),
  S("wp2", c("wp2_benchmark2/run_heldout.R", "arrays", "dry"), "R:results_R/wp2_data/prefetch_log_retry.txt", "results_R/wp2/heldout_arrays_dry.jsonl", 20, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/run_heldout.R", "arrays", "full"), "R:results_R/wp2/heldout_arrays_dry.jsonl", "results_R/wp2/heldout_arrays_results.jsonl", 55, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/run_rnaseq.R", "dry"), "R:results_R/wp2_data/prefetch_log_retry.txt", "results_R/wp2/heldout_rnaseq_dry.jsonl", 60, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/run_rnaseq.R", "full"), "R:results_R/wp2/heldout_rnaseq_dry.jsonl", "results_R/wp2/heldout_rnaseq_results.jsonl", 60, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/calibrate_t2.R", "run", "1", "1"), "R:results_R/wp2/compact/dev_sets.json", "results_R/wp2/compact/calib_series (one worker; resumes)", 480, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/design_table.R", "run", "1", "1"), "R:results_R/wp2/compact/dev_sets.json", "results_R/wp2/compact/design_cells (one worker; resumes)", 360, TRUE, mode = "R"),
  S("wp2", c("wp2_benchmark2/apply_dev.R"), "wp2_benchmark2/compact/dev_core_results.jsonl", "results_R/wp2: dev_application_summary.md (from the archived core results)", 1, mode = "archive"),
  S("wp2", c("wp2_benchmark2/apply_dev.R"), R_CORE, "results_R/wp2: dev_application_summary.md, dev_application_per_series.csv (from the R core results)", 1, mode = "R", env = ENV_WP2_R),
  S("wp2", c("wp2_benchmark2/summarize_validation.R"), "wp2_benchmark2/heldout_arrays_results.jsonl", "results_R/wp2: validation_summary.md, WP2_summary_table.md (from the archived results)", 1, mode = "archive"),
  S("wp2", c("wp2_benchmark2/summarize_validation.R"), c("R:results_R/wp2/dev_application_per_series.csv", "R:results_R/wp2/heldout_arrays_results.jsonl", "R:results_R/wp2/heldout_rnaseq_results.jsonl"),
    "results_R/wp2: validation_summary.md, WP2_summary_table.md (from the R results)", 1, mode = "R", env = ENV_WP2_R),
  S("wp2", c("wp2_benchmark2/calibrate_t2.R", "summarize"), c("R:results_R/wp2/compact/calib_series", "wp2_benchmark2/compact/dev_core_results.jsonl"), "results_R/wp2: calibration_tables.csv, calibration_verdict.txt (calibration runs of the R port, made by calibrate_t2.R run; archived core results)", 1, mode = "archive"),
  S("wp2", c("wp2_benchmark2/calibrate_t2.R", "summarize"), c("R:results_R/wp2/compact/calib_series", R_CORE), "results_R/wp2: calibration_tables.csv, calibration_verdict.txt (from the R runs)", 1, mode = "R", env = ENV_WP2_R),
  S("wp2", c("wp2_benchmark2/design_table.R", "summarize"), "R:results_R/wp2/compact/design_cells", "results_R/wp2: design_table.md (from the R design runs, made by design_table.R run)", 1, mode = "both"),
  S("wp2run", c("run_wp2_validation.sh"), c("R:wp2_benchmark2/gate_part1.R", "R:results_R/wp2_data/prefetch_log_retry.txt"), "data runs of WP2 against the archive (extract, held-out arrays and RNA-seq, about 3 h)", 200, TRUE, mode = "archive"),
  S("posthoc", c("run_post_hoc.sh"), c("R:post_hoc_191126.R", "R:results_R/wp2/compact/dev_sets.json"), "results_R/post_hoc: dev (14 min), INS-1 (2 min), validation (60 min)", 80, TRUE),
  S("wp4data", c("wp4_islets/wp4_http.R"), "wp4_islets/frozen_islet_datasets.csv", "download of the raw data of the six islet datasets into results_R/wp4_data (about 450 MB; size checked first)", 30, TRUE),
  S("wp4data", c("wp4_islets/run_datasets.R"), c("wp4_islets/frozen_beta_cells.csv", "wp4_islets/ensembl_ids_17genes.csv"), "results_R/wp4/results: cell-level analysis of the six islet datasets (about 90 min)", 90, TRUE),
  S("wp4", c("run_wp4_chain.sh"), c("R:results_R/wp4/results/index.csv"), "results_R/wp4: meta, gene-graph null, metacells (about 2 h), positive control, wider search, tables", 180, TRUE),
  S("wp4", c("wp4_islets/perturbation_gse60803.R"), "wp4_islets/perturbation/candidates_decisions.csv", "results_R/wp4/perturbation", 5), S("wp4", c("wp4_islets/wider_search_freeze.R"), "wp4_islets/wider_search/candidates_all.csv", "results_R/wp4/wider_search_freeze", 1),
  S("wp4", c("wp4_islets/perturbation_decisions.R"), "wp4_islets/perturbation/candidates.csv", "results_R/wp4/perturbation_decisions", 1), S("wp4", c("wp4_islets/freeze_lists.R"), character(0), "results_R/wp4/frozen", 1),
  S("check4", c("wp4_islets/check_readme_numbers.R"), "R:results_R/wp4/results/index.csv", "numbers of the islet README (CHECK_SOURCE=R)", 1), S("check4", c("wp4_islets/check_wp4b_numbers.R"), "R:results_R/wp4/positive_control/positive_control_pooled.csv", "numbers of README_WP4b and the wording (CHECK_SOURCE=R)", 1),
  S("figures", c("figures/run_all_figures.sh"), c("R:figures/figstyle.R"), "results_R/figures: Figs 7, 8, 9, 10, A5.1, S1 and the route-3 benchmark figure (PNG 300 dpi and PDF)", 10),
  S("af6", c("af6/build_af6.R"), AF6_TEMPLATE, "results_R/af6/Additional file 6.xlsx: the results workbook from the archived results", 3, mode = "archive", optional = TRUE),
  S("af6", c("af6/build_af6.R"), AF6_TEMPLATE, "results_R/af6/Additional file 6.xlsx: the results workbook from the R results", 3, mode = "R", env = c(AF6_SOURCE = "R"), optional = TRUE),
  S("af6", c("compare_af6.R"), c(AF6_TEMPLATE, "R:results_R/af6/Additional file 6.xlsx"), "cell by cell comparison with the submitted workbook", 2, optional = TRUE),
  S("report", c("agreement_report.R"), character(0), "results_R/agreement_report.md: every R result against the saved Python result (calibrated on --source archive)", 15))
NOTES <- c("two modes: --source archive (default; verification of the port against the archived Python results, the mode the agreement gate is calibrated on) and --source R (the pipeline from the data in R)",
           "live queries (not re-run): wp4_islets/perturbation_scout.R, wider_search.R, wider_search_meta.R (GEO and ArrayExpress queries), wp9_safety.R retrieve, wp12_annotations.R retrieve, wp10_string_sensitivity.R retrieve, benchmark_route3/select_datasets.R, wp2_benchmark2/select_*.R",
           "heavy steps (one at a time): the downloads (wp2_prefetch.R, wp4_http.R), run_route3_full.sh, gate_part1.R extract/compare/consistency, run_heldout.R, run_rnaseq.R, calibrate_t2.R run, design_table.R run, run_post_hoc.sh, wp4_islets/run_datasets.R, run_wp4_chain.sh")
input_path <- function(rel) if (startsWith(rel, "R:")) file.path(R_PORT, substring(rel, 3)) else file.path(V2, rel)
made_by_a_step <- function(rel) startsWith(rel, "R:results_R/")                                  # written by an earlier step: checked when its step starts
check <- function(steps) {                                                                         # before step 1: every script, and every input that is a file of the repository or of this folder
  missing <- character(0); later <- character(0); skip <- logical(length(steps))
  for (k in seq_along(steps)) { s <- steps[[k]]; sp <- file.path(R_PORT, s$cmd[1]); if (!file.exists(sp)) { missing <- c(missing, sp); cat("MISSING SCRIPT [", s$pkg, "]: ", sp, "\n", sep = "") }
    absent <- Filter(function(rel) !made_by_a_step(rel) && !file.exists(input_path(rel)), s$inputs)
    if (s$optional && length(absent)) { skip[k] <- TRUE; cat("OPTIONAL STEP SKIPPED [", s$pkg, "] ", paste(s$cmd, collapse = " "), ": not in this checkout: ", paste(absent, collapse = ", "), "\n", sep = ""); next }
    for (rel in s$inputs) { p <- input_path(rel); if (file.exists(p)) next
      if (made_by_a_step(rel)) { later <- c(later, rel); next }
      missing <- c(missing, p); cat("MISSING INPUT [", s$pkg, "]: ", p, "\n", sep = "") } }
  list(missing = missing, later = unique(later), skip = skip)
}
with_env <- function(env, expr) {                                                                  # set environment variables for one step and restore them afterwards
  if (!length(env)) return(expr); old <- Sys.getenv(names(env), unset = NA, names = TRUE); do.call(Sys.setenv, as.list(env))
  on.exit(for (k in names(env)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, setNames(list(old[[k]]), k)), add = TRUE); expr
}
main <- function() {
  a <- commandArgs(trailingOnly = TRUE); dry <- "--dry-run" %in% a; skip_heavy <- "--skip-heavy" %in% a; only <- if ("--only" %in% a) strsplit(a[match("--only", a) + 1L], ",", fixed = TRUE)[[1]] else NULL
  src <- if ("--source" %in% a) a[match("--source", a) + 1L] else "archive"; if (!(src %in% c("archive", "R"))) stop("--source must be archive or R")
  if ("--list" %in% a) { cat(paste0(NOTES, collapse = "\n"), "\n"); return(invisible()) }
  steps <- Filter(function(s) s$mode %in% c("both", src) && (is.null(only) || s$pkg %in% only), STEPS); ck <- check(steps); steps <- steps[!ck$skip]
  for (s in steps) cat(if (dry) "would run" else "running", sprintf(" [%s] %s  -> %s (about %d min%s)\n", s$pkg, paste(s$cmd, collapse = " "), s$out, s$mins, if (s$heavy) ", heavy" else ""), sep = "")
  cat(sprintf("source: %s; %d steps, about %.1f h; %d inputs are written by earlier steps and checked when their step starts\n", src, length(steps), sum(vapply(steps, function(s) s$mins, 0)) / 60, length(ck$later)))
  if (dry) { cat("dry run finished;", if (!length(ck$missing)) "all scripts and repository inputs present" else sprintf("%d scripts or inputs missing", length(ck$missing)), "\n"); return(invisible()) }
  if (length(ck$missing)) stop(sprintf("stopping: %d scripts or inputs are missing (listed above)", length(ck$missing)))
  for (s in steps) {
    if (skip_heavy && s$heavy) { cat(sprintf("[%s] %s skipped (heavy step; run it alone, see --list)\n", s$pkg, s$cmd[1])); next }
    gone <- Filter(function(rel) made_by_a_step(rel) && !file.exists(input_path(rel)) && !grepl("export_tables_done", rel, fixed = TRUE), s$inputs)
    if (length(gone)) stop(sprintf("step [%s] %s: input not written yet (an earlier step did not run or failed): %s", s$pkg, s$cmd[1], paste(gone, collapse = ", ")))
    t0 <- Sys.time(); sc <- file.path(R_PORT, s$cmd[1]); lg <- if (is.null(s$log)) "" else input_path(s$log); if (nzchar(lg)) dir.create(dirname(lg), showWarnings = FALSE, recursive = TRUE)
    rc <- with_env(s$env, if (endsWith(sc, ".sh")) system2("bash", c(sc, s$cmd[-1]), stdout = lg) else system2(RS, c("--vanilla", sc, s$cmd[-1]), stdout = lg))
    if (is.null(rc) || !is.numeric(rc)) rc <- 0L
    cat(sprintf("[%s] %s exit %d in %.0f s\n", s$pkg, paste(s$cmd, collapse = " "), rc, as.numeric(difftime(Sys.time(), t0, units = "secs")))); if (rc != 0) stop("step failed: ", paste(s$cmd, collapse = " "))
  }
}
if (sys.nframe() == 0L) main()
