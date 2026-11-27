## build_af6.R - R port of wp6a_assets/build_af6.py together with the three later changes of the submitted workbook, so that one run gives the final Additional file 6:
##   apply_move_tables_181126.py (one phrase of the Index: the PV_pairs_top15 line now points to Supplementary Table S11), add_af6_sheets_201126.py (sheets PostHoc_per_matrix, PostHoc_summary and Fig_10_numbers, their Index rows and a sentence of the last note)
##   and apply_audit_fixes_201126.py (nine source cells of Fig_10_numbers and a second sentence of the last note).
## Additional file 6: results tables of the version-2 analyses per dataset and series, as one workbook. A NEW workbook written from the result files named in the sheet 'Index'; nothing is analyzed here. One sheet per table; sheet 'Index' lists every sheet
## with its source file and a one-line description. Header row bold and frozen, column widths set, no merged cells, numbers stored as numbers. Large source files are not copied (calibration_replicates.csv, 1.9 MB: only its summary tables are included).
## Usage: Rscript build_af6.R [--no-extension]
##   AF6_SOURCE = archive (default: the saved result files of strengthen_v2/, the numbers cited in the manuscript) | R (the results of the R run in r_port/results_R/, where a result exists)
##   AF6_OUT    = output file (default results_R/af6/Additional file 6.xlsx; the Python files are never touched)
## Compare with the submitted workbook: Rscript compare_af6.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "xlsx_io.R")); suppressMessages(library(jsonlite))
V2 <- normalizePath(file.path(R_PORT, "..")); SOURCE <- Sys.getenv("AF6_SOURCE", "archive"); OUT <- Sys.getenv("AF6_OUT", file.path(R_PORT, "results_R", "af6", "Additional file 6.xlsx"))
EXT <- !("--no-extension" %in% commandArgs(TRUE)); V2_CUTOFF <- "2026-11-06 16:00:00"      # attempts logged later belong to the WP4 extension (positive control, wider search)
## archived result path (relative to strengthen_v2/) -> the same file of the R run (relative to r_port/results_R/); files that are not listed are inputs (frozen lists, records of live queries, pre-specifications)
RMAP <- c("wp1_ins1_structure/tbl_" = "wp1/tbl_", "wp5_reporting/effect_sizes.csv" = "wp5/effect_sizes.csv", "wp3_mlss/gate_scores.csv" = "wp3/gate_scores.csv", "wp3_mlss/mlss_scenarios.csv" = "wp3/mlss_scenarios.csv", "wp3_mlss/mlss_scenario_scores.csv" = "wp3/mlss_scenario_scores.csv",
          "wp2_benchmark2/calibration_" = "wp2/calibration_", "wp2_benchmark2/dev_application_per_series.csv" = "wp2/dev_application_per_series.csv", "wp2_benchmark2/validation_" = "wp2/validation_", "wp2_benchmark2/WP2_summary_table.csv" = "wp2/WP2_summary_table.csv", "wp2_benchmark2/design_table.csv" = "wp2/design_table.csv",
          "wp4_islets/results/" = "wp4/results/", "wp4_islets/meta/" = "wp4/meta/", "wp4_islets/metacell/" = "wp4/metacell/", "wp4_islets/perturbation/gse60803_" = "wp4/perturbation/gse60803_", "wp4_islets/perturbation/candidates_decisions.csv" = "wp4/perturbation_decisions/candidates_decisions.csv",
          "wp4_islets/positive_control/" = "wp4/positive_control/", "wp4_islets/wider_search/frozen_candidates.csv" = "wp4/wider_search_freeze/frozen_candidates.csv", "wp4_islets/wider_search/concordance_" = "wp4/wider_search/concordance_", "wp4_islets/wider_search/fold_changes_" = "wp4/wider_search/fold_changes_",
          "wp9_safety/results/" = "wp9/", "wp10_versions/results/" = "wp10/", "wp11_plate_provenance/results/" = "wp11/", "wp12_database_annotations/results/" = "wp12/", "wp13_gapdh_corrected/results/" = "wp13/", "ars_review_191126/post_hoc/results_" = "post_hoc/results_", "fig10_201126/fig10_numbers.csv" = "figures_fromR/fig10_numbers.csv")
no_r_file <- character(0)
src <- function(rel) {
  if (SOURCE == "R") for (i in seq_along(RMAP)) if (startsWith(rel, names(RMAP)[i])) { p <- file.path(R_PORT, "results_R", paste0(RMAP[[i]], substring(rel, nchar(names(RMAP)[i]) + 1L))); if (file.exists(p)) return(p); no_r_file <<- union(no_r_file, rel) }
  file.path(V2, rel)
}
csv <- function(rel) pd_read_csv(src(rel))
csv_text <- function(rel, sep = ",") pd_read_csv(src(rel), sep, na_values = "")      # text such as n/a and n/e stays as written; truly empty cells stay empty
txt_lines <- function(rel) readLines(src(rel), warn = FALSE, encoding = "UTF-8")
jload <- function(rel) fromJSON(src(rel), simplifyVector = FALSE)

islet_tables <- function() {
  idx <- csv("wp4_islets/results/index.csv"); tags <- chr_of(idx, "tag"); vars <- chr_of(idx, "variant"); modes <- chr_of(idx, "mode"); pairs <- list(); det <- list(); nullcols <- c("null_mean", "null_95pct", "empirical_p_ge_observed")
  for (i in seq_len(nrow_t(idx))) {
    stem <- paste0(tags[i], "__", vars[i], "__", modes[i]); p <- csv(paste0("wp4_islets/results/", stem, "__pairs.csv")); cn <- csv(paste0("wp4_islets/results/", stem, "__cons.csv"))
    cons <- filter_rows(p, vapply(col_of(p, "consensus_invitro"), isTRUE, TRUE))
    if (nrow_t(cn) && all(nullcols %in% cn$names)) { k <- match(c("gene1", "gene2", nullcols), cn$names); cons <- merge_left(cons, tbl(cn$names[k], cn$cols[k], cn$dtype[k]), c("gene1", "gene2")) }
    else { stopifnot(sum(vapply(col_of(cons, "evaluable"), isTRUE, TRUE)) == 0); for (x in nullcols) cons <- tbl(c(cons$names, x), c(cons$cols, list(rep(list(NA), nrow_t(cons)))), c(cons$dtype, "float")) }      # no evaluable pair: the cons file is empty
    cons <- insert_col(cons, "mode", modes[i], "obj"); cons <- insert_col(cons, "variant", vars[i], "obj"); cons <- insert_col(cons, "dataset", tags[i], "obj"); pairs[[length(pairs) + 1L]] <- drop_col(cons, "consensus_invitro")
    if (modes[i] == "adjusted") { d <- csv(paste0("wp4_islets/results/", stem, "__detection.csv")); d <- insert_col(d, "variant", vars[i], "obj"); d <- insert_col(d, "dataset", tags[i], "obj"); det[[length(det) + 1L]] <- d }
  }
  list(idx = idx, pairs = concat(pairs), det = concat(det))
}
meta_tables <- function() {
  ms <- jload("wp4_islets/meta/meta_summary.json"); est <- list(); loo <- list()
  for (v in ms$variants) {
    for (e in v$estimators) est[[length(est) + 1L]] <- c(list(variant = v$variant, mode = v$mode, datasets = paste(unlist(v$qualifying), collapse = ";")), e)
    for (e in v$leave_one_out) loo[[length(loo) + 1L]] <- c(list(variant = v$variant, mode = v$mode), e)
  }
  list(est = df_from_records(est), loo = df_from_records(loo))
}
concordance_table <- function() {
  cc <- jload("wp4_islets/perturbation/gse60803_concordance.json"); rows <- list()
  for (kl in list(c("all_donors", "primary: 4 donors, perturbed mean minus control"), c("S114_only", "secondary, descriptive: donor S114 only (nontarget control, n = 1)"))) {
    x <- cc[[kl[1]]]; rows[[length(rows) + 1L]] <- list(contrast = kl[2], series = cc$series, platform = cc$platform, n_genes = x$n_genes, n_same_sign = x$n_same_sign, binomial_p = x$binomial_p, spearman_rho = x$spearman_rho, spearman_p = x$spearman_p,
                                                       same_sign_genes = paste(unlist(x$same_sign_genes), collapse = ";"), opposite_sign_genes = paste(unlist(x$opposite_sign_genes), collapse = ";"), ins1_arm = cc$ins1_arm, confounded = cc$confounded)
  }
  df_from_records(rows)
}
## str() and repr() of the values of a JSON file, as Python writes them
py_repr <- function(o) if (is.null(o)) "None" else if (is.list(o)) { if (!is.null(names(o))) paste0("{", paste(sprintf("%s: %s", vapply(names(o), py_repr, ""), vapply(o, py_repr, "")), collapse = ", "), "}") else paste0("[", paste(vapply(o, py_repr, ""), collapse = ", "), "]") } else if (is.character(o)) paste0("'", gsub("'", "\\'", o, fixed = TRUE), "'") else pstr(o)
py_str <- function(o) if (is.null(o)) "None" else if (is.list(o)) py_repr(o) else if (is.character(o)) o else pstr(o)
flat_row <- function(path, o) { num <- !is.null(o) && is.numeric(o); list(item = paste(path, collapse = "."), value_number = if (num) as.numeric(o) else NA, value_text = if (num) NA else py_str(o)) }
queries_table <- function() {      # queries.json of the wider search flattened to one row per entry (dictionaries are walked; lists are written as Python writes them)
  q <- jload("wp4_islets/wider_search/queries.json"); rows <- list()
  walk <- function(o, path) { if (is.list(o) && !is.null(names(o))) for (k in names(o)) walk(o[[k]], c(path, k)) else rows[[length(rows) + 1L]] <<- flat_row(path, o) }
  walk(q, character(0)); df_from_records(rows)
}
flat_json <- function(files) {      # json files flattened to one row per entry (item = path, the prefix of the file first); a list of plain values is joined with '; ' (an empty list is written as []); a list of lists or dicts gets one row per element
  rows <- list(); add <- function(path, o) rows[[length(rows) + 1L]] <<- flat_row(path, o)
  walk <- function(o, path) {
    if (is.list(o) && !is.null(names(o))) { for (k in names(o)) walk(o[[k]], c(path, k)) }
    else if (is.list(o)) { if (!length(o)) add(path, "[]") else if (!any(vapply(o, is.list, TRUE))) add(path, paste(vapply(o, py_str, ""), collapse = "; ")) else for (i in seq_along(o)) walk(o[[i]], c(path, as.character(i - 1L))) }
    else add(path, o)
  }
  for (f in files) walk(jload(f[1]), if (nzchar(f[2])) f[2] else character(0))
  df_from_records(rows)
}
pv_event_groups <- function() {      # the table of section 5 of the WP9 pre-specification (code | theme | preferred terms)
  md <- paste(txt_lines("wp9_safety/PRESPEC_WP9_safety_screen.md"), collapse = "\n"); sec <- substr(md, regexpr("## 5. Event groups", md, fixed = TRUE), regexpr("## 6. Statistics", md, fixed = TRUE) - 1L)
  lines <- trimws(strsplit(sec, "\n", fixed = TRUE)[[1]]); lines <- lines[startsWith(lines, "|")]
  rows <- lapply(lines, function(l) trimws(strsplit(paste0(gsub("^[|]+|[|]+$", "", l), "|"), "|", fixed = TRUE)[[1]])); hdr <- rows[[1]]; body <- rows[-(1:2)]
  stopifnot(identical(hdr, c("code", "theme", "preferred terms")), length(body) == 8, all(vapply(body, length, 1L) == 3L))
  tbl(hdr, lapply(1:3, function(j) lapply(body, function(r) r[j])), rep("obj", 3))
}
## ---------------------------------------------------------------- sheets with several blocks (WP11, WP12, WP13): source lines, then per block a title row, a header row and the data rows; cells are typed one by one
TEXT_COLS <- c("first_US_NDA_BLA_approval", "first_ORIG_AP_date")
typed <- function(t, text_col = FALSE) {
  if (t == "") return(NA)
  if (text_col) return(t)
  if (grepl("^[+-]?[0-9]+$", t)) { s <- sub("^[+-]+", "", t); return(if (nchar(s) > 1 && startsWith(s, "0")) t else pint(as.numeric(t))) }
  if (grepl("^[+-]?(?:[0-9]+[.][0-9]*|[.][0-9]+|[0-9]+)(?:[eE][+-]?[0-9]+)?$", t, perl = TRUE)) return(py_float(t))
  t
}
csv_block <- function(title, rel, sep = ",") {
  x <- read_rows(src(rel), sep); hdr <- unlist(x[1, ], use.names = FALSE); body <- x[-1, , drop = FALSE]; nb <- nrow(body); nh <- length(hdr)
  boolc <- vapply(seq_len(nh), function(j) any(body[[j]] != "") && all(body[[j]] %in% c("", "True", "False")), TRUE)
  list(title = title, header = hdr, rows = lapply(seq_len(nb), function(i) lapply(seq_len(nh), function(j) { v <- body[[j]][i]; if (boolc[j] && v != "") v == "True" else typed(v, hdr[j] %in% TEXT_COLS) })))
}
blocks_sheet <- function(lines, blocks) structure(list(lines = lines, blocks = blocks), class = "af6_blocks")
dims <- function(t) if (inherits(t, "af6_blocks")) c(sum(vapply(t$blocks, function(b) length(b$rows), 1L)), max(vapply(t$blocks, function(b) length(b$header), 1L))) else c(nrow_t(t), length(t$cols))
grab <- function(pat, s) regmatches(s, regexec(pat, s, perl = TRUE))[[1]]
plate1_block <- function() {      # the plate-1 lines of threshold_fit.txt as a table of thresholds
  txt <- txt_lines("wp11_plate_provenance/results/threshold_fit.txt"); seg <- txt[which(startsWith(txt, "=== 24 June plate 1"))[1]:length(txt)]
  a <- grab("threshold of the export = ([0-9.]+) \\(mean \\|diff\\| ([0-9.]+), max ([0-9.]+)\\)", seg[1])
  b <- grab("T = ([0-9.]+) RFU, mean \\|recomputed - sheet\\| = ([0-9.]+) cycles over ([0-9]+) values \\(blocks skipped because a well never reaches T: ([0-9]+)\\)", seg[startsWith(seg, "threshold that best reproduces the sheet values")][1])
  lc <- seg[startsWith(seg, "mean |recomputed - sheet| at other thresholds:")][1]; cm <- lapply(regmatches(lc, gregexpr("T [0-9]+: [0-9.]+", lc))[[1]], function(z) grab("T ([0-9]+): ([0-9.]+)", z)[-1])
  d <- grab("Gapdh blocks T = ([0-9.]+) \\(mean ([0-9.]+)\\); target blocks T = ([0-9.]+) \\(mean ([0-9.]+)\\); the five analyzed genes only \\(([^)]*)\\) T = ([0-9.]+) \\(mean ([0-9.]+)\\)", seg[startsWith(seg, "thresholds fitted separately:")][1])
  stopifnot(length(a) == 4, length(b) == 5, length(cm) == 8, length(d) == 8); n <- function(x) py_float(x)
  rows <- list(list("threshold of the 2025 re-export (the export itself)", n(a[2]), n(a[3]), n(a[4]), NA, NA), list("best single threshold for the workbook, all blocks", n(b[2]), n(b[3]), NA, pint(n(b[4])), pint(n(b[5]))),
               list("Gapdh blocks alone", n(d[2]), n(d[3]), NA, NA, NA), list("target blocks alone", n(d[4]), n(d[5]), NA, NA, NA), list(sprintf("the five analyzed genes alone (%s)", d[6]), n(d[7]), n(d[8]), NA, NA, NA))
  rows <- c(rows, lapply(cm, function(z) list("mean difference at a fixed threshold", pint(n(z[1])), n(z[2]), NA, NA, NA)))
  list(title = "2. 24 June 2020, plate 1 (file without \" 2\"): thresholds, from the plate-1 lines of threshold_fit.txt (RFU = relative fluorescence units; differences in cycles)",
       header = c("fit", "threshold_RFU", "mean_abs_difference_cycles", "max_abs_difference_cycles", "n_values", "blocks_skipped"), rows = rows)
}
SH <- function(name, group, srcfiles, desc, data, label = NULL) list(name = name, group = group, src = srcfiles, desc = desc, data = data, label = label)

round5_sheets <- function() {
  G8 <- "WP11 plate provenance (descriptive, read only, post hoc)"; G9 <- "WP12 database annotations (descriptive; the reported analysis is unchanged)"; p11 <- "wp11_plate_provenance/results/"; p12 <- "wp12_database_annotations/results/"; out <- list()
  ## ---- Prov_aug26_rows
  b <- csv_block("26 August run: the rows of the calculation workbook against the instrument export (aug26_workbook_vs_instrument.csv)", paste0(p11, "aug26_workbook_vs_instrument.csv")); h <- b$header; R <- b$rows
  ix <- setNames(match(c("target", "group_in_workbook", "workbook_Cq_equals_instrument_Cq_of_named_well", "block_value_set_equals_instrument"), h), c("target", "group", "eqc", "same"))
  val <- function(r, k) r[[ix[[k]]]]; eq <- function(rs) sum(vapply(rs, function(r) isTRUE(val(r, "eqc") == 1), TRUE)); tg <- Filter(function(r) val(r, "target") != "Gapdh", R); gd <- Filter(function(r) val(r, "target") == "Gapdh", R); off <- Filter(function(r) isTRUE(val(r, "same") == 0), R)
  lines <- c("Source: wp11_plate_provenance/results/aug26_workbook_vs_instrument.csv (WP11, 8 Nov 2026). Descriptive provenance check, read only, post hoc; it changes no number of the reported analysis.",
             "One row per Cq row of the calculation workbook of the 26 August run (Cacna1c, Cacna1d, Kcnj11, Ptbp1, Tcf7l2 and Gapdh; 4 wells per gene and group), compared with the instrument export: the Cq of the well named in the row, the difference, and whether the block holds the same four values as the instrument and is ascending down the workbook (1 = yes, 0 = no).",
             sprintf("Counts of these rows: %d of %d carry the Cq of the well named in the row (target rows %d of %d, Gapdh rows %d of %d); %d rows lie in a block whose four values differ from the instrument (%s).",
                     eq(R), length(R), eq(tg), length(tg), eq(gd), length(gd), length(off), paste(py_sort(unique(vapply(off, function(r) sprintf("%s %s", val(r, "target"), val(r, "group")), ""))), collapse = ", ")))
  out[[length(out) + 1L]] <- SH("Prov_aug26_rows", G8, paste0(p11, "aug26_workbook_vs_instrument.csv"), "WP11: the Cq rows of the calculation workbook of the 26 August run compared with the instrument export: the Cq of the well named in the row, the difference, and whether the block holds the same four values and is ascending (the counts are in the first rows of the sheet); descriptive, read only, post hoc.", blocks_sheet(lines, list(b)))
  ## ---- Prov_threshold
  lines <- c("Source: wp11_plate_provenance/results/threshold_validation.csv, threshold_fit.txt, refA_dcq_comparison.csv and refA_block_residuals.csv (WP11, 8 Nov 2026). Reconstruction of the Cq threshold from the baseline-subtracted amplification curves; descriptive, read only, post hoc.",
             "Blocks below: 1, validation of the reconstruction against the thresholds of the exports (the known threshold must come back); 2, plate 1 of 24 June 2020 (thresholds that reproduce the 2025 re-export and the workbook); 3, dCq (target minus Gapdh) of the workbook against the 2025 re-export and against the curves at the fitted threshold; 4, largest absolute difference of the four sorted values of each block at the fitted threshold (cycles).")
  blocks <- list(csv_block("1. Validation: the fitted threshold must come back as the known one (threshold_validation.csv)", paste0(p11, "threshold_validation.csv")), plate1_block(),
                 csv_block("3. dCq (target minus Gapdh refA) of the workbook (sheet), the 2025 re-export at 45.70 RFU and the curves at 209.91 RFU (refA_dcq_comparison.csv)", paste0(p11, "refA_dcq_comparison.csv")),
                 csv_block("4. Largest absolute difference of the four sorted values of each block at the fitted threshold, in cycles; skipped = a well never reaches the threshold (refA_block_residuals.csv)", paste0(p11, "refA_block_residuals.csv")))
  out[[length(out) + 1L]] <- SH("Prov_threshold", G8, paste(paste0(p11, c("threshold_validation.csv", "threshold_fit.txt", "refA_dcq_comparison.csv", "refA_block_residuals.csv")), collapse = "; "), "WP11: reconstruction of the Cq threshold from the baseline-subtracted curves: validation against the known thresholds of the exports, the thresholds of plate 1 of 24 June 2020, the dCq comparison of the workbook with the 2025 re-export and with the curves, and the block differences (four blocks, each with its own header row).", blocks_sheet(lines, blocks))
  ## ---- Ann_hubs_TDL
  lines <- c("Source: wp12_database_annotations/results/hub_annotation.csv (WP12, 9 Nov 2026). The hub genes (degree 5 or more) and two context genes of the network, annotated from the DrugCentral drug-target files of the Database folder: the file of the reported run (R_17genes/drug.target.interaction.csv) and the file of 25 Sep 2026.",
             "TDL = target development level (Pharos/IDG) as stored in the DrugCentral file; \"no record in the file\" means that the gene has no drug-target record, so no TDL is stored. Counts are records in the file, not usable potency values: the MLSS used fewer drugs than the records contain. Descriptive; the hubs were chosen by degree alone.")
  out[[length(out) + 1L]] <- SH("Ann_hubs_TDL", G9, paste0(p12, "hub_annotation.csv"), "WP12: the four hub genes and two context genes with their DrugCentral records (reported file and file of 25 Sep 2026), distinct drugs, target development level (TDL) and MLSS drugs; counts are records, not usable potency values.", blocks_sheet(lines, list(csv_block("Hub and context genes (hub_annotation.csv)", paste0(p12, "hub_annotation.csv")))))
  ## ---- Ann_MLSS_classes
  sm <- jload("wp12_database_annotations/results/summary.json"); q <- sm$quartiles
  lines <- c("Source: wp12_database_annotations/results/mlss_quartile_classes.csv and summary.json (WP12, 9 Nov 2026). Score-quartile classes Q1 to Q4 of the MLSS ranking: a description of the ranking, not calibrated labels; descriptive, the reported analysis is unchanged.",
             "Definition (pre-specified): t75, t50 and t25 are the 75th, 50th and 25th percentiles of the scores of the set; Q1 = score at or above t75, Q2 = at or above t50, Q3 = at or above t25, Q4 = below t25; ties stay together. Sets: the 153 combinations of the reported run and the 171 combinations with the DrugCentral file of 25 Sep 2026 (tirzepatide added).")
  t1 <- list(title = "Thresholds and class sizes (summary.json)", header = c("set", "n_combinations", "t75", "t50", "t25", "Q1", "Q2", "Q3", "Q4"),
             rows = list(c(list("153 combinations of the reported run", q$n_reported), unname(q$thresholds_reported_153[c("t75", "t50", "t25")]), unname(q$class_sizes_reported_153[c("Q1", "Q2", "Q3", "Q4")])),
                         c(list("171 combinations with the DrugCentral file of 25 Sep 2026", q$n_171), unname(q$thresholds_171[c("t75", "t50", "t25")]), unname(q$class_sizes_171[c("Q1", "Q2", "Q3", "Q4")]))))
  t2 <- list(title = "Other values of summary.json", header = c("item", "value"),
             rows = list(list("distinct score values among the 153 reported combinations", q$distinct_scores_reported_153), list("largest absolute score difference of the 153 shared combinations (file of 25 Sep 2026 against the reported file)", q$max_abs_score_difference_shared),
                         list("reported combinations that change class when the 171 are used (all from Q3 to Q4)", length(q$reported_combinations_changing_class_when_171_used))))
  out[[length(out) + 1L]] <- SH("Ann_MLSS_classes", G9, paste0(p12, "mlss_quartile_classes.csv; ", p12, "summary.json"), "WP12: score-quartile classes Q1 to Q4 of the 153 reported and the 171 combinations (thresholds and class sizes of summary.json in the first blocks, then every combination with rank, score, class and cross-/same-target status); a description of the ranking, not calibrated labels.",
                                blocks_sheet(lines, list(t1, t2, csv_block("Classes of the combinations (mlss_quartile_classes.csv; in_reported_153 = 1 for the 153 reported combinations)", paste0(p12, "mlss_quartile_classes.csv")))))
  ## ---- Ann_US_status and Ann_US_applications
  mr <- sm$market_retrieval; caution <- "Caution: list membership or application database, not a regulatory audit. United States only; \"Discontinued\" says that the application holder lists no marketed product, not why; \"Prescription\" means an active listing in Drugs@FDA, not a check of sales; nothing is stated about other countries."
  ap <- csv_block("NDA and BLA applications (market_status_nda_bla_applications.csv)", paste0(p12, "market_status_nda_bla_applications.csv")); dcol <- match("drug", ap$header)
  lines <- c(sprintf("Source: wp12_database_annotations/results/market_status_by_drug.csv (WP12; annotation of 9 Nov 2026). openFDA drug/drugsfda.json (Drugs@FDA), query products.active_ingredients.name:\"<name>\", %d saved responses retrieved from %s to %s; the DrugCentral columns are list membership (FDA, EMA, PMDA) in the files of the Database folder (25 Sep 2026).", mr$responses_saved, mr$first_retrieved, mr$last_retrieved), caution)
  out[[length(out) + 1L]] <- SH("Ann_US_status", G9, paste0(p12, "market_status_by_drug.csv; ", p12, "summary.json"), "WP12: market status of the drug entries in Drugs@FDA (openFDA, 9 Nov 2026) with the DrugCentral FDA, EMA and PMDA list membership; list membership or application database, not a regulatory audit; United States only.",
                                blocks_sheet(lines, list(csv_block("Market status by drug (market_status_by_drug.csv)", paste0(p12, "market_status_by_drug.csv")))))
  lines <- c(sprintf("Source: wp12_database_annotations/results/market_status_nda_bla_applications.csv (WP12; annotation of 9 Nov 2026). The %d NDA and BLA application rows of the %d drugs that have a US application in Drugs@FDA (openFDA, %d saved responses retrieved from %s to %s).", length(ap$rows), length(unique(lapply(ap$rows, function(r) r[[dcol]]))), mr$responses_saved, mr$first_retrieved, mr$last_retrieved), caution)
  out[[length(out) + 1L]] <- SH("Ann_US_applications", G9, paste0(p12, "market_status_nda_bla_applications.csv"), "WP12: the NDA and BLA applications behind Ann_US_status (sponsor, brands, marketing status of the products, date of the first original approval); list membership or application database, not a regulatory audit; United States only.", blocks_sheet(lines, list(ap)))
  out
}
round6_sheets <- function() {      # WP13, the Gapdh reference sensitivity of the 26 August run: the sheet Gapdh_versions (source lines first; five blocks)
  G10 <- "WP13 Gapdh reference sensitivity of the 26 August run (descriptive, post hoc; the reported analysis stays primary)"; w13 <- "wp13_gapdh_corrected/results/"
  J <- jload("wp13_gapdh_corrected/results/gapdh_versions.json"); VER <- c("reported", "C0", "C1"); G5 <- c("Tcf7l2", "Kcnj11", "Cacna1d", "Cacna1c", "Ptbp1")
  b1 <- csv_block("1. The three versions side by side (gapdh_versions_summary.csv; the fold changes and r are rounded as in the file)", paste0(w13, "gapdh_versions_summary.csv"))
  b2 <- list(title = "2. DAPT fold change against the control mean of the five genes of the 26 August run, per version (gapdh_versions.json, key fc_dapt; per-replicate 2^-dCt relative to the control mean)",
             header = c("gene", "reported", "C0", "C1"), rows = lapply(G5, function(g) c(list(g), lapply(VER, function(v) J[[v]]$fc_dapt[[g]]))))
  rows3 <- list(); for (v in c("C0", "C1")) for (sk in list(c("lost", "lost"), c("new", "new"))) for (p in J[[v]][[sk[2]]]) rows3[[length(rows3) + 1L]] <- list(v, sk[1], p)
  b3 <- list(title = sprintf("3. Consensus pairs lost and new against the 21 consensus pairs of the reported version (gapdh_versions.json, keys lost and new; the reported version has none; C0 keeps %d of the 21 and C1 keeps %d)", J$C0$kept, J$C1$kept), header = c("version", "status", "pair"), rows = rows3)
  W <- J$qc_gapdh_wells
  b4 <- list(title = "4. The 12 Gapdh wells of the 26 August run: sample group, Cq and plateau fluorescence in RFU (gapdh_versions.json, key qc_gapdh_wells)", header = c("well", "sample", "cq", "plateau_RFU"), rows = lapply(W, function(w) list(w$well, w$sample, w$cq, w$plateau_RFU)))
  med <- as.numeric(np_median(vapply(Filter(function(w) !(w$well %in% c("A10", "B10")), W), function(w) as.numeric(w$plateau_RFU), 1)))
  b5 <- list(title = sprintf("5. Plateau fluorescence of wells A10 and B10 as a share of the median of the other ten Gapdh wells (derived from block 4; the median of the ten wells other than A10 and B10 is %s RFU)", py_fmt_g(med, 6L)),
             header = c("well", "plateau_RFU", "median_plateau_of_the_other_ten_wells_RFU", "share_of_that_median_percent"),
             rows = lapply(Filter(function(w) w$well %in% c("A10", "B10"), W), function(w) list(w$well, w$plateau_RFU, med, py_round(100 * w$plateau_RFU / med, 1L))))
  lines <- c("Source: wp13_gapdh_corrected/results/gapdh_versions_summary.csv and gapdh_versions.json (WP13, 9 Nov 2026; independent check of 10 Nov 2026: 17 checks, 0 failures). Sensitivity analysis of the Gapdh reference of the DAPT group of the 26 August run; descriptive, post hoc, one preparation; the reported analysis stays primary (decision of 9 Nov 2026).",
             "Versions: \"reported\" = DAPT reference rows 22.5685, 22.7983, 22.5685, 22.7983 (wells C10 and D10, each entered twice); C0 = wells A10 and B10 excluded without duplication (every DAPT replicate divided by the mean of C10 and D10, 22.6834); C1 = all four wells (mean 23.6215), the worst case. Only the DAPT group of the 26 August run changes; the other runs and groups are as reported.",
             "Blocks below: 1, the three versions side by side; 2, DAPT fold changes of the five genes of the run per version; 3, consensus pairs lost and new against the 21 of the reported version; 4, the 12 Gapdh wells of the run with Cq and plateau fluorescence (RFU); 5, the plateau of wells A10 and B10 as a share of the median of the other ten wells. The report is part 5.9 of Additional file 5.")
  list(SH("Gapdh_versions", G10, paste0(w13, "gapdh_versions_summary.csv; ", w13, "gapdh_versions.json"), "WP13: the DAPT reference of the 26 August run in three versions (reported; C0, wells A10 and B10 excluded without duplication; C1, all four wells): the versions side by side, the DAPT fold changes of the five genes, the consensus pairs lost and new, and the 12 Gapdh wells with Cq and plateau fluorescence; descriptive, post hoc, the reported analysis stays primary.",
          blocks_sheet(lines, list(b1, b2, b3, b4, b5))))
}

## ---------------------------------------------------------------- the 82 sheets of build_af6.py
sheets <- function() {
  S <- list(); add <- function(name, group, srcfiles, desc, data, label = NULL) S[[length(S) + 1L]] <<- SH(name, group, srcfiles, desc, data, label)
  G1 <- "INS-1 (post hoc, one preparation)"; w1 <- "wp1_ins1_structure/"
  for (z in list(
    c("INS1_gene_run_map", "gene_run_map.csv", "Verified gene-to-run map: the calculation-sheet run that supplied each matrix column (exact value match), with the evidence for each gene."),
    c("INS1_run_structure", "tbl_2_2_run_structure.csv", "Run-structure test (protocol 2.2): mean within-group r of same-run against between-run pairs, per grouping and run, with label-permutation P."),
    c("INS1_degree_by_run", "tbl_2_3_degree_by_run.csv", "Degree of the eligible genes by run and hub counts against expectation (protocol 2.3)."),
    c("INS1_eligible_pairs", "tbl_2_3_eligible_pairs.csv", "The 91 eligible gene pairs with run labels, same-run flag and consensus status (protocol 2.3)."),
    c("INS1_hubs_vs_run", "tbl_2_3_hubs_vs_run.csv", "Each gene with its run, degree, hub status and same-run edges against expectation (protocol 2.3)."),
    c("INS1_pairs_involving_run", "tbl_2_3_pairs_involving_run.csv", "Eligible and consensus pairs that involve each run, with the hypergeometric P (protocol 2.3)."),
    c("INS1_pairs_reference", "tbl_2_4_pairs.csv", "Per pair: within-group r on target Cq and on dCq, and the r expected from the shared reference alone (protocol 2.4)."),
    c("INS1_shared_reference", "tbl_2_4_shared_reference_summary.csv", "Mean within-group r of same-run and between-run pairs on target Cq, on dCq and expected from the shared reference (protocol 2.4)."),
    c("INS1_sensitivity_versions", "tbl_2_4_sensitivity_versions.csv", "Consensus network and checks when the Gapdh reference is replaced by run-and-group means or run-wide means (protocol 2.4)."),
    c("INS1_adjusted_calibrated", "tbl_2_5_adjusted_calibrated.csv", "Adjusted q < 0.05 counts per scheme against their permutation null (protocol 2.5)."),
    c("INS1_q_per_pair", "tbl_2_5_q_per_pair.csv", "Treatment-adjusted correlation, p and q per pair and scheme; ratio schemes are artifact-prone (protocol 2.5)."),
    c("INS1_coupling", "tbl_2_6_coupling.csv", "Within-group coupling statistics and consensus count against the independent-gene and the run-block null (protocol 2.6)."),
    c("INS1_geometry", "tbl_2_7_geometry.csv", "Share of pooled-significant raw pairs against the geometric expectation (protocol 2.7)."),
    c("INS1_scale_sensitivity", "tbl_2_8_scale_sensitivity.csv", "Consensus count, null, covariance share, adjusted counts and hubs on the linear, log and Spearman scales (protocol 2.8)."))) add(z[1], G1, paste0(w1, z[2]), z[3], csv(paste0(w1, z[2])))
  add("Effect_sizes", "Effect sizes (RT-qPCR panel, one preparation)", "wp5_reporting/effect_sizes.csv", "Fold changes with 95% intervals for RT-qPCR, DNA content, 2-NBDG and C-peptide, Gapdh variants, BH across the 17 genes; within-experiment values.", csv("wp5_reporting/effect_sizes.csv"))
  G3 <- "MLSS (recomputation, WP3)"
  for (z in list(
    c("MLSS_scenarios", "mlss_scenarios.csv", "The 40 scenarios of protocol section 5: agreement with the reported ranking, top-15 overlap, permutation summaries."),
    c("MLSS_scenario_scores", "mlss_scenario_scores.csv", "MLSS of each of the 153 combinations under every deterministic scenario."),
    c("MLSS_gate_scores", "gate_scores.csv", "Gate: the 153 scores of the R script beside the scores of the Python port (largest difference 5e-15)."),
    c("MLSS_direction", "direction_table.csv", "Direction table: drug action against the mRNA change of its target under DAPT and DKK-1."),
    c("MLSS_pic50_sources", "pic50_sources.csv", "Source records of the pIC50 values: type, database, assay, species and reason when not used."))) add(z[1], G3, paste0("wp3_mlss/", z[2]), z[3], csv(paste0("wp3_mlss/", z[2])))
  G2 <- "Benchmark part 2 (public data, WP2)"; w2 <- "wp2_benchmark2/"
  vl <- txt_lines(paste0(w2, "calibration_verdict.txt")); vl <- vl[nzchar(trimws(vl))]
  add("Bench_calib_verdict", G2, paste0(w2, "calibration_verdict.txt"), "Verdict on the pre-specified calibration requirement (protocol 3.2): type I error 7.5% (fail), power 90.7% (pass); post hoc sensitivity without the 2 series with constant genes.", tbl("line", list(as.list(vl)), "obj"))
  for (z in list(
    c("Bench_calibration", "calibration_tables.csv", "Calibration on decoupled development data: type I error and power per statistic, injected r and number of coupled genes (37 series), Wilson intervals."),
    c("Bench_calib_by_platform", "calibration_by_platform.csv", "Type I error and power at r = 0.7 with 6 coupled genes by array platform."),
    c("Bench_calib_no_degenerate", "calibration_sensitivity_no_degenerate.csv", "Post hoc: calibration without the 2 series that have genes with no within-group variation."),
    c("Bench_calib_per_series", "calibration_per_condition.csv", "Rejections per development series, statistic and condition (summary of calibration_replicates.csv, which is not copied)."),
    c("Bench_dev_per_series", "dev_application_per_series.csv", "Development series: endpoints of protocol 3.3 and 3.4 per series and gene set (A: 17 random differentially expressed genes; B: 17 random expressed genes)."),
    c("Bench_validation_per_series", "validation_per_series.csv", "All 96 analyzed series (37 development, 40 held-out arrays, 19 human RNA-seq): endpoints per series and gene set."),
    c("Bench_validation_endpoints", "validation_endpoints_long.csv", "Prevalence and median endpoints per arm and gene set with Wilson intervals (long format)."),
    c("Bench_summary", "WP2_summary_table.csv", "Summary table of the validation endpoints by arm (the table printed in Additional file 5, part 4)."),
    c("Bench_design", "design_table.csv", "Design table (protocol 3.7): power of T_r2 and pair detection by experiments per group and r (300 simulated experiments per cell)."),
    c("Bench_frozen_arrays", "frozen_heldout_arrays.csv", "Frozen list of held-out microarray series (fingerprinted before any expression value was read)."),
    c("Bench_frozen_rnaseq", "frozen_heldout_rnaseq.csv", "Frozen list of held-out RNA-seq series (fingerprinted before any count was read)."))) add(z[1], G2, paste0(w2, z[2]), z[3], csv(paste0(w2, z[2])))
  G4 <- "Public islet cross-check, version 2 (WP4)"; it <- islet_tables(); mt <- meta_tables(); w4 <- "wp4_islets/"
  add("Islet_dataset_index", G4, paste0(w4, "results/index.csv"), "One row per dataset, variant (primary, non-diabetic donors only) and mode (adjusted, unadjusted): β-cells, donors, evaluable genes and pairs, median rho, matched null, E and SE of E.", it$idx)
  add("Islet_consensus_pairs", G4, paste0(w4, "results/<dataset>__<variant>__<mode>__pairs.csv and __cons.csv"), "The 21 consensus pairs in every dataset, variant and mode: evaluable flag, partial Spearman rho, p, q and the matched null (mean, 95th percentile, empirical P); not evaluable pairs have empty values.", it$pairs)
  add("Islet_detection", G4, paste0(w4, "results/<dataset>__<variant>__adjusted__detection.csv"), "Detection of the 17 genes in the β-cells of each dataset and variant (percentage of β-cells with a value above 0; evaluable = at least 10%).", it$det)
  add("Islet_meta_study_effects", G4, paste0(w4, "meta/study_effects__<variant>__<mode>.csv"), "Excess co-expression E and its SE per dataset for the four variants of the meta-analysis.",
      concat(unlist(lapply(c("primary", "nondiabetic"), function(v) lapply(c("adjusted", "unadjusted"), function(m) csv(sprintf("%smeta/study_effects__%s__%s.csv", w4, v, m)))), recursive = FALSE)))
  add("Islet_meta_estimates", G4, paste0(w4, "meta/meta_summary.json"), "Random-effects pooled E per variant, mode and estimator (primary: REML + Hartung-Knapp, margin +/-0.10) with the verdict of the fixed decision rule.", mt$est)
  add("Islet_meta_leave_one_out", G4, paste0(w4, "meta/meta_summary.json"), "Leave-one-dataset-out pooled E (REML + Hartung-Knapp) per variant and mode with the verdict.", mt$loo)
  add("Islet_gene_graph_datasets", G4, paste0(w4, "meta/gene_graph_by_dataset.csv"), "POST HOC add-on: E and SE per dataset under a null that keeps the gene graph of the evaluable pairs, beside the pre-specified SE.", csv(paste0(w4, "meta/gene_graph_by_dataset.csv")))
  add("Islet_gene_graph_meta", G4, paste0(w4, "meta/gene_graph_meta.csv"), "POST HOC add-on: pooled E under the gene-graph null for the five estimators.", csv(paste0(w4, "meta/gene_graph_meta.csv")))
  add("Islet_pair_pooled_rho", G4, paste0(w4, "meta/pair_level_pooled_rho.csv"), "Pair-level pooled rho (REML + Hartung-Knapp on Fisher z) for the 11 pairs evaluable in at least 3 human datasets, with the rho of each dataset.", csv(paste0(w4, "meta/pair_level_pooled_rho.csv")))
  add("Islet_metacell_summary", G4, paste0(w4, "metacell/metacell_summary.csv"), "Metacell sensitivity (20 random partitions per size): metacells, evaluable genes and pairs, median rho and E as median and range over partitions.", csv(paste0(w4, "metacell/metacell_summary.csv")))
  add("Islet_metacell_pairs", G4, paste0(w4, "metacell/metacell_pair_evaluability.csv"), "Number of partitions in which each consensus pair is evaluable at the metacell level, per dataset and size.", csv(paste0(w4, "metacell/metacell_pair_evaluability.csv")))
  add("Islet_metacell_genes", G4, paste0(w4, "metacell/metacell_gene_evaluability.csv"), "Number of partitions in which each gene is evaluable at the metacell level, per dataset and size, beside the cell-level status.", csv(paste0(w4, "metacell/metacell_gene_evaluability.csv")))
  add("Pert_candidates", G4, paste0(w4, "perturbation/candidates_decisions.csv"), if (EXT) "First perturbation search (protocol 4.6): the 43 candidate GEO series from the two fixed queries with the decision and its reason (the wider search of the WP4 extension is in the Ext_ sheets)." else "Perturbation search (protocol 4.6): the 43 candidate GEO series from the two fixed queries with the decision and its reason.", csv(paste0(w4, "perturbation/candidates_decisions.csv")))
  add("Pert_GSE60803_fold_changes", G4, paste0(w4, "perturbation/gse60803_fold_changes.csv"), "GSE60803 (β-catenin shRNA): log2 fold change of the 17 genes per donor and overall beside the INS-1 DKK-1 arm and the sign agreement.", csv(paste0(w4, "perturbation/gse60803_fold_changes.csv")))
  add("Pert_GSE60803_concordance", G4, paste0(w4, "perturbation/gse60803_concordance.json"), "Sign agreement (exact binomial) and Spearman correlation of GSE60803 with the INS-1 DKK-1 arm; the contrast is confounded (controls are not expanded cells).", concordance_table())
  add("Islet_dataset_citations", G4, paste0(w4, "dataset_citations.csv"), "Publications attached to the dataset records (PubMed via the records; to be verified before citing).", csv(paste0(w4, "dataset_citations.csv")))
  att <- csv(paste0(w4, "attempts.csv"))
  if (EXT) add("Islet_attempts", G4, paste0(w4, "attempts.csv"), "Every dataset and source attempted, with its outcome (rows up to 14:14 on 6 Nov belong to the version-2 analyses, rows from 16:12 on to the WP4 extension).", att)
  else add("Islet_attempts", G4, paste0(w4, "attempts.csv"), "Every dataset and source attempted in the version-2 analyses with its outcome (rows logged up to 14:14 on 6 Nov; later rows of the file belong to the WP4 extension and are not copied).", filter_rows(att, chr_of(att, "when") < V2_CUTOFF))
  if (EXT) {
    G5 <- "WP4 extension: positive control and wider search (post hoc relative to the primary result)"; pc <- paste0(w4, "positive_control/"); ws <- paste0(w4, "wider_search/")
    add("Ext_posctrl_by_dataset", G5, paste0(pc, "positive_control_by_dataset.csv"), "Positive control A1: coupling of nominal rho spiked into the real β-cells (200 replicates per dataset, rho and mode): recovered excess E, attenuation and share of replicates with empirical P < 0.05.", csv(paste0(pc, "positive_control_by_dataset.csv")))
    add("Ext_posctrl_pooled", G5, paste0(pc, "positive_control_pooled.csv"), "Positive control A1: verdict of the pooled five human datasets per nominal rho and mode (shares of replicates: no meaningful excess, inconclusive, co-expressed more).", csv(paste0(pc, "positive_control_pooled.csv")))
    add("Ext_posctrl_margin", G5, paste0(pc, "margin_sensitivity.csv"), "Positive control A3: verdict of the primary pooled interval for margins of 0.05 to 0.20 (descriptive; the pre-set margin stays +/-0.10).", csv(paste0(pc, "margin_sensitivity.csv")))
    add("Ext_posctrl_ribosomal", G5, paste0(pc, "rp_control.csv"), "Positive control A2: excess E of ribosomal protein gene pairs (a real co-expression module) per dataset, with the share of random sets at or above +0.10.", csv(paste0(pc, "rp_control.csv")))
    add("Ext_search_queries", G5, paste0(ws, "queries.json"), "Wider search: the three term groups (cells, Notch, Wnt) and the query of each source with the number of records it returned (S1 GEO DataSets, S2 ArrayExpress collection of BioStudies, S3 PubMed linked to GEO series), the 135 unique records and their split by source combination; one row per entry of queries.json (numbers in value_number, text in value_text).", queries_table())
    add("Ext_search_candidates", G5, paste0(ws, "frozen_candidates.csv"), "Wider search for perturbation data (GEO, BioStudies, PubMed-linked): all 135 unique records with decision, reason, expected sign, mirror flag, contrast and data route, frozen before any expression value was read.", csv(paste0(ws, "frozen_candidates.csv")))
    add("Ext_search_results", G5, paste0(ws, "concordance_results.csv"), "Wider search: per kept series (6) the status (analyzable or not), the weak-match note and the sign agreement and Spearman correlation with the matching INS-1 arm.", csv(paste0(ws, "concordance_results.csv")))
    add("Ext_search_meta", G5, paste0(ws, "concordance_meta.csv"), "Wider search: pooled Spearman correlation (random effects, Hartung-Knapp) with the outcome of the pre-specified decision rule; sets with fewer than 3 analyzable series are not pooled.", csv(paste0(ws, "concordance_meta.csv")))
    for (z in list(c("GSE60803", "β-catenin shRNA in expanded human β-cell-derived cells; controls are not expanded cells"), c("GSE167073", "Fzd9 knockout in islets of MycER-BclXL mice"), c("GSE99503", "GSK3 inhibitor in EndoC-βH1 cells; mirror design, expected sign −1")))
      add(paste0("Ext_FC_", z[1]), G5, paste0(ws, "fold_changes_", z[1], ".csv"), sprintf("Wider search, %s (%s): log2 fold change of the 17 genes (perturbed minus control), expected sign, sign-adjusted fold change (mirror designs multiplied by −1) beside the INS-1 DKK-1 arm.", z[1], z[2]), csv(paste0(ws, "fold_changes_", z[1], ".csv")))
  }
  ## ---- WP9: pharmacovigilance safety screen
  G6 <- "WP9 safety screen: spontaneous-report annotation (post hoc to the MLSS; not a prediction)"; w9 <- "wp9_safety/"; r9 <- paste0(w9, "results/")
  add("PV_event_groups", G6, paste0(w9, "PRESPEC_WP9_safety_screen.md"), "WP9: the 8 event groups of the screen (code, theme, MedDRA preferred terms matched exactly; 40 terms), fixed in section 5 of the pre-specification before any count was retrieved.", pv_event_groups())
  add("PV_drug_counts", G6, paste0(r9, "drug_counts.csv"), "WP9: the 17 search terms for the 18 MLSS entries ((S)-nitrendipine is searched as nitrendipine) with target gene, reports through the harmonized drug field and whether the drug is evaluable (100 reports or more).", csv_text(paste0(r9, "drug_counts.csv")))
  add("PV_signals_single", G6, paste0(r9, "signals_single.csv"), "WP9: reporting odds ratio (Woolf 95% interval) of each of the 13 evaluable drugs for each of the 8 event groups against all other reports; signal = at least 3 reports and a lower limit above 1 (104 drug-group tests).", csv_text(paste0(r9, "signals_single.csv")))
  add("PV_pairs_all77", G6, paste0(r9, "pairs_all77.csv"), "WP9: the 77 cross-target pairs in MLSS rank order with the reports that list both drugs, the status (evaluable at 20 or more), the signals of each drug (n/e = drug not evaluable) and the co-report excess flags (n/a = fewer than 20 co-reports; - = none).", csv_text(paste0(r9, "pairs_all77.csv")))
  add("PV_pairs_top15", G6, paste0(r9, "pairs_top15.csv"), "WP9: the 15 top-ranked pairs (the primary values behind Table 2 of the manuscript; n/e = drug not evaluable).", csv_text(paste0(r9, "pairs_top15.csv")))
  add("PV_pairs_long", G6, paste0(r9, "pairs_long.csv"), "WP9: event-group detail of the 16 evaluable pairs (128 pair-group tests): reports in the group, proportion among the reports that list both drugs with exact 95% interval, proportions among the reports of each drug without the other, and the co-report excess flag.", csv_text(paste0(r9, "pairs_long.csv")))
  add("PV_summary", G6, paste0(r9, "summary.json"), "WP9: summary.json flattened to key and value rows (snapshot date, number of reports, evaluable drugs, signals, pair counts, co-report flags, per-group lists); a list is joined with \"; \" and an empty list is shown as [].", flat_json(list(c(paste0(r9, "summary.json"), ""))))
  lab <- "SECONDARY, POST HOC, FREE-TEXT NAMES (patient.drug.medicinalproduct): "
  add("PV_fallback_signals", G6, paste0(r9, "fallback_signals_single.csv"), paste0(lab, "signals of the four drugs with fewer than 100 reports through the harmonized field (albiglutide, bepridil, nilvadipine, nitrendipine), searched by the free-text names fixed in section 9 of the pre-specification; not merged with the primary counts."),
      csv_text(paste0(r9, "fallback_signals_single.csv")), label = paste0(lab, "not merged with the primary counts; the four drugs without harmonized records (albiglutide, bepridil, nilvadipine, nitrendipine); same signal rule as the primary screen"))
  add("PV_fallback_pairs", G6, paste0(r9, "fallback_pairs.csv"), paste0(lab, "co-report counts of the pairs that contain one of these four drugs (the partner drug keeps its primary counts; a pair with two such drugs has one row per drug); not merged with the primary counts."),
      csv_text(paste0(r9, "fallback_pairs.csv")), label = paste0(lab, "not merged with the primary counts; the partner drug keeps its primary counts (n/e = partner not evaluable); same pair rule as the primary screen"))
  ## ---- WP10: database-version sensitivity
  G7 <- "WP10 database-version sensitivity (descriptive; the reported analysis is unchanged)"; w10 <- "wp10_versions/results/"
  add("Ver_string_classification", G7, paste0(w10, "string_classification.csv"), "WP10: the 120 pairs with |r| >= 0.7 of the Fig. 5 classification with the four-method consensus flag, the STRING scores (releases 12.0 and 12.5, names as in the R code) and the category in four settings (12.0 as reported, 12.5, and both with the Glut2 alias for Slc2a2); descriptive.", csv_text(paste0(w10, "string_classification.csv")))
  add("Ver_string_summary", G7, paste0(w10, "string_gate.json; ", w10, "string_summary.json"), "WP10: string_gate.json (the port reproduces the reported classification from the 12.0 export; the Glut2 alias counts) and string_summary.json (counts per class and release with and without the alias, pairs gained and lost, category changes) flattened to key and value rows (item with the prefix gate or summary; value_number; value_text; a list is joined with \"; \" and an empty list is shown as []).",
      flat_json(list(c(paste0(w10, "string_gate.json"), "gate"), c(paste0(w10, "string_summary.json"), "summary"))))
  add("Ver_string_v125_network", G7, paste0(w10, "string_v12_5_network.tsv"), "WP10: the 59 distinct pairs of the STRING 12.5 functional network (rat, the 21 genes of the matrix, Glut2 queried as Slc2a2, minimum score 0.4) with the evidence scores as returned by the API.", csv_text(paste0(w10, "string_v12_5_network.tsv"), sep = "\t"))
  add("Ver_mlss_summary", G7, paste0(w10, "drugcentral_summary.json"), "WP10: drugcentral_summary.json flattened to key and value rows: rows and drugs of the DrugCentral file of the reported run and of the file of 25 Sep 2026, MLSS drugs per hub, combinations, shared combinations (scores unchanged), the reported top 15 in the new ranking and the combinations with tirzepatide.", flat_json(list(c(paste0(w10, "drugcentral_summary.json"), ""))))
  add("Ver_mlss_comparison", G7, paste0(w10, "mlss_comparison.csv"), "WP10: every combination of the union of the two MLSS runs (171) with old and new rank and score and the status (both, new only); the 153 reported combinations keep their scores.", csv_text(paste0(w10, "mlss_comparison.csv")))
  add("Ver_pic50_changes", G7, paste0(w10, "pic50_changes.csv"), "WP10: the pIC50 values that differ between the two DrugCentral files: only the new value of tirzepatide at GLP1R (8.37); no value of an existing drug changed.", csv_text(paste0(w10, "pic50_changes.csv")))
  c(S, round5_sheets(), round6_sheets())
}

NOT_INCLUDED <- list(
  c("wp2_benchmark2/calibration_replicates.csv", "1.9 MB, one row per simulated replicate; the summary tables are in Bench_calibration, Bench_calib_by_platform, Bench_calib_no_degenerate and Bench_calib_per_series."),
  c("wp3_mlss/mlss_permutation_draws.csv", "2 x 1000 permutation draws of the MLSS scenarios; their summaries are in MLSS_scenarios."),
  c("wp2_benchmark2/*.jsonl, eligibility_*.csv, dev_endpoints_long.csv", "per-series raw run records and eligibility lists; the endpoints are in Bench_validation_per_series, Bench_validation_endpoints and Bench_dev_per_series."),
  c("wp4_islets/results/*__pairs.csv (all 136 pairs)", "only the 21 consensus pairs are copied (Islet_consensus_pairs); the 136-pair files stay in the project folder."),
  c("wp9_safety/results/raw/*.json, probe_log.txt, retrieve_log.txt, fallback_retrieve_log.txt, check_log.txt", "the 443 saved openFDA responses (about 1 MB; their combined fingerprint is in part 5.5 of Additional file 5), the request logs and the check log of the safety screen; the tables behind the analysis are in the PV_ sheets."),
  c("wp11_plate_provenance/results/threshold150_matches.csv, threshold150_summary.txt, refA_vs_plate1.txt, aug26_nov09_order_check.txt and .json", "the exploration outputs of the provenance check as text or json (block matches at 150 RFU, the group comparison of plate 1 and the order check of the 26 August and 9 November runs); the tables behind the answers are in the Prov_ sheets."),
  c("wp12_database_annotations/results/raw/*.json, analyze_stdout.txt", "the 19 saved openFDA Drugs@FDA responses and the print of the analysis; the tables behind the annotations are in the Ann_ sheets."),
  c("wp13_gapdh_corrected/README.md, wp13_gapdh_versions.py, check_wp13_numbers.py", "the report (copied word for word into part 5.9 of Additional file 5), the script and the independent check of the Gapdh reference sensitivity of the 26 August run; the tables behind the report are in the sheet Gapdh_versions."),
  c("wp10_versions/results/mlss_new_all_combinations.csv, string_version.json, drugcentral_log.txt", "all 171 combinations of the MLSS with the new DrugCentral file with the four score components, the answer of the STRING server about its version and the run log; the comparison of the two rankings is in Ver_mlss_comparison."),
  if (EXT) c("wp4_islets/positive_control/inj_*.npz, rp_*.json; wider_search/candidates_all.csv, sample_records_*.tsv", "raw simulation arrays, per-dataset ribosomal sets, the unfiltered list of every record returned by the wider search and the sample records read for the keep decisions; the tables behind the decisions are in the Ext_ sheets.")
  else c("wp4_islets/positive_control/, wp4_islets/wider_search/ and the attempts logged from 16:12 on 6 Nov", "the WP4 extension (positive control, wider search for perturbation data); not included in this build (--no-extension)."))
NOTES <- c("The INS-1 data come from one cell preparation (12 observations) and every INS-1 analysis here is post hoc. Public-data tables carry the name of their dataset.",
           "Every table is copied from the result file named in the column \"Source file\" without change of values; where a source is a JSON file it is flattened to rows. Numbers are stored as numbers; empty cells are missing values.",
           "The sheets Prov_* (WP11), Ann_* (WP12) and Gapdh_versions (WP13) start with rows that state the source, the date and, for the status sheets, the caution, and have several blocks (a title row, a header row and the data rows of each; blank rows between them). In them text cells (n/a, skipped, \"no record in the file\") stay text, dates written YYYYMMDD stay text, True and False are booleans and other numbers are numbers.",
           "Companion text: the parts 5.6 (database versions, sheets Ver_*), 5.7 (plate provenance, sheets Prov_*), 5.8 (database annotations, sheets Ann_*) and 5.9 (Gapdh reference sensitivity of the 26 August run, sheet Gapdh_versions) of Additional file 5 are the reports behind these sheets. In Additional file 5 the items that the authors must supply are listed in part 7; those that the authors have already decided are marked resolved there.",
           "Workbook built 6 to 10 Nov 2026 by build_af6.py (wp6a_assets).")
SRCCOL <- "Source file (relative to Analysis_outputs/strengthen_v2/)"

## ---------------------------------------------------------------- the post hoc sheets and the later edits (20 Nov 2026)
posthoc_sheets <- function() {
  PH <- "ars_review_191126/post_hoc/"; J <- function(...) jload(paste0(PH, paste(c(...), collapse = "/")))
  num <- function(v) {      # '' -> None, True / False -> bool, digits -> int, else float, else text
    if (is.na(v) || v == "") return(NA)
    if (v %in% c("True", "False")) return(v == "True")
    if (grepl("^-*[0-9]+$", v)) return(pint(as.numeric(v)))
    if (grepl("^\\s*[+-]?(([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?|[Ii][Nn][Ff]([Ii][Nn][Ii][Tt][Yy])?|[Nn][Aa][Nn])\\s*$", v)) return(py_float(v))
    v
  }
  dict_rows <- function(rel) { x <- read_rows(src(rel)); h <- unlist(x[1, ], use.names = FALSE); x <- x[-1, , drop = FALSE]; lapply(seq_len(nrow(x)), function(i) setNames(as.list(unlist(x[i, seq_along(h)], use.names = FALSE)), h)) }
  FIELDS <- c("gse", "platform", "set", "obs_median", "obs_mean", "obs_n_pairs", "null_n_pairs_mean", "null_no_sig_frac", "null_mean", "null_q025", "null_q975", "excess", "above_q975", "T_r2", "T_r2_p", "T_r2_pc1", "T_r2_pc1_p", "pc1_share", "null_pc1_mean")
  rows <- c(lapply(dict_rows(paste0(PH, "results_dev/per_matrix.csv")), function(r) c(list("development", "development arrays"), lapply(FIELDS, function(k) num(r[[k]])))),
            lapply(dict_rows(paste0(PH, "results_validation/per_matrix.csv")), function(r) c(list("validation", r$arm), lapply(FIELDS, function(k) num(r[[k]])))))
  I <- J("results_ins1", "summary.json"); comb <- c(I$share, I$coupling); comb <- comb[!duplicated(names(comb), fromLast = TRUE)]
  rows <- c(rows, list(c(list("INS-1", "RT-qPCR, 3 groups x 4 wells", "INS-1", "RT-qPCR", "17 genes"), lapply(FIELDS[-(1:3)], function(k) if (is.null(comb[[k]])) NA else comb[[k]]))))
  HEAD_M <- c("series_set", "arm", FIELDS)
  SUM_H <- c("series_set", "arm", "gene_set", "n_matrices", "share_observed_q1", "share_observed_median", "share_observed_q3", "share_expected_q1", "share_expected_median", "share_expected_q3", "excess_q1", "excess_median", "excess_q3", "n_above_null_q975", "n_share_defined",
             "T_r2_flagged_original", "original_ci_lo", "original_ci_hi", "T_r2_flagged_after_PC1_removal", "pc1_ci_lo", "pc1_ci_hi", "paired_both", "paired_only_before", "paired_only_after", "paired_neither", "pc1_variance_share_q1", "pc1_variance_share_median", "pc1_variance_share_q3")
  srows <- list()
  add_summary <- function(label, arm, S) for (gs in c("A", "B", "A+B")) { v <- S$sets[[gs]]
    srows[[length(srows) + 1L]] <<- c(list(label, arm, gs, v$n_matrices), v$obs_share_quartiles, v$null_share_quartiles, v$excess_quartiles, list(v$n_above_null_q975, v$n_share_defined, v$flagged_original), v$flagged_original_ci, list(v$flagged_pc1), v$flagged_pc1_ci,
                                      list(v$paired$both, v$paired$only_original, v$paired$only_after_removal, v$paired$neither), v$pc1_share_quartiles) }
  add_summary("development", "development arrays (37 series)", J("results_dev", "summary.json")); add_summary("validation", "all validation sets (59)", J("results_validation", "summary.json"))
  for (af in list(c("arrays GPL1261", "summary_arrays_GPL1261.json"), c("arrays GPL570", "summary_arrays_GPL570.json"), c("rnaseq human", "summary_rnaseq_human.json"))) add_summary("validation", af[1], J("results_validation", af[2]))
  f10 <- read_rows(src("fig10_201126/fig10_numbers.csv")); f10h <- unlist(f10[1, ], use.names = FALSE); f10b <- lapply(seq_len(nrow(f10) - 1L) + 1L, function(i) as.list(unlist(f10[i, seq_along(f10h)], use.names = FALSE)))
  ## audit fix of the same day: nine 'source' cells of Fig_10_numbers
  PS <- function(a, b) sprintf("sheet PostHoc_summary: series_set validation, arm all validation sets (59), gene_set %s, column %s", a, b)
  NEW <- c("share carried by group means, gene set A (median of series)" = PS("A", "share_observed_median"), "expected share without coupling, gene set A" = PS("A", "share_expected_median"), "share carried by group means, gene set B" = PS("B", "share_observed_median"),
           "expected share without coupling, gene set B" = PS("B", "share_expected_median"), "matrices above the 97.5th null percentile" = paste0(PS("A+B", "n_above_null_q975"), " (2 of 118 matrices)"), "T_r2 flagged after removal of PC1, gene set A" = paste0(PS("A", "T_r2_flagged_after_PC1_removal"), " (32 of 59)"),
           "T_r2 flagged, original, same sets (gene set A)" = paste0(PS("A", "T_r2_flagged_original"), " (47 of 59)"), "INS-1 share carried by group means" = "sheet PostHoc_per_matrix: series_set INS-1, column obs_median", "INS-1 expected share" = "sheet PostHoc_per_matrix: series_set INS-1, column null_mean")
  done <- character(0); f10b <- lapply(f10b, function(r) { lab <- r[[1]]; if (lab %in% names(NEW)) { r[[3]] <- NEW[[lab]]; done <<- c(done, lab) }; r }); stopifnot(setequal(done, names(NEW)))
  srcrel <- "ars_review_191126/post_hoc/"; G <- "Post hoc checks of the benchmark (20 Nov 2026)"
  entries <- list(
    list("PostHoc_per_matrix", G, paste0(srcrel, "results_dev/per_matrix.csv; ", srcrel, "results_validation/per_matrix.csv; ", srcrel, "results_ins1/summary.json"),
         "One row per matrix: the share of the pooled covariance carried by group means (observed, expected without within-group coupling, 95% null interval, excess) and T_r2 before and after removal of the first principal component (P, share of the residual variance in PC1); 74 development matrices, 118 validation matrices, the INS-1 matrix. Method: Additional file 5, part 4.4.", length(rows), length(HEAD_M)),
    list("PostHoc_summary", G, paste0(srcrel, "results_dev/summary.json; ", srcrel, "results_validation/summary*.json"),
         "Summary of the same checks: quartiles of the observed and expected share and of their difference, matrices above the 97.5th null percentile, T_r2 flagged before and after removal of PC1 with Wilson 95% intervals and the paired table, by gene set (A, B, both), for the development series, all validation sets and each validation arm.", length(srows), length(SUM_H)),
    list("Fig_10_numbers", "Fig. 10 (20 Nov 2026)", "fig10_201126/fig10_numbers.csv", "Every number of Fig. 10 with the result file or the manuscript text it comes from.", length(f10b), length(f10h)))
  list(entries = entries, sheets = list(list(name = "PostHoc_per_matrix", header = HEAD_M, rows = rows, widths = c(14, 30, 12, 12, 8, rep(13, 16))),
                                        list(name = "PostHoc_summary", header = SUM_H, rows = srows, widths = c(12, 30, 9, rep(13, 25))),
                                        list(name = "Fig_10_numbers", header = f10h, rows = f10b, widths = c(60, 16, 90))))
}

## ---------------------------------------------------------------- writing
width_of <- function(header, col) { vals <- c(header, vapply(Filter(function(x) !missing_v(x), utils::head(col, 200)), pstr, "")); max(9, min(60, floor(max(nchar(vals, type = "chars")) * 1.05) + 2)) }
put <- function(s, t, label = NULL) {
  r0 <- 1L; if (!is.null(label)) { sheet_cell(s, 1L, 1L, label, TRUE); r0 <- 2L }
  sheet_cells(s, rep(r0, length(t$names)), seq_along(t$names), as.list(t$names), TRUE)
  n <- nrow_t(t); if (n) { v <- unlist(t$cols, recursive = FALSE); keep <- !vapply(v, missing_v, TRUE); sheet_cells(s, (r0 + rep(seq_len(n), times = length(t$cols)))[keep], rep(seq_along(t$cols), each = n)[keep], v[keep], FALSE) }
  for (j in seq_along(t$names)) sheet_width(s, j, width_of(t$names[j], t$cols[[j]]))
  s$freeze <- r0
}
put_blocks <- function(s, b) {      # lines in rows 1..k (row 1 bold), then per block a blank row, its title row (bold), its header row (bold) and its data rows; row 1 frozen; no merged cells
  r <- 1L; ncol <- max(vapply(b$blocks, function(x) length(x$header), 1L)); widths <- rep(9, ncol)
  for (k in seq_along(b$lines)) { sheet_cell(s, r, 1L, b$lines[k], k == 1L); r <- r + 1L }
  for (bl in b$blocks) {
    r <- r + 1L; sheet_cell(s, r, 1L, bl$title, TRUE); r <- r + 1L
    for (j in seq_along(bl$header)) { h <- bl$header[j]; sheet_cell(s, r, j, h, TRUE); widths[j] <- max(widths[j], min(60, floor(nchar(h) * 1.05) + 2)) }
    r <- r + 1L
    for (row in bl$rows) { for (j in seq_along(row)) { v <- row[[j]]; if (missing_v(v)) next; sheet_cell(s, r, j, v, FALSE); widths[j] <- max(widths[j], min(60, floor(nchar(pstr(v)) * 1.05) + 2)) }; r <- r + 1L }
  }
  for (j in seq_along(widths)) sheet_width(s, j, widths[j]); s$freeze <- 1L
}
main <- function() {
  t0 <- Sys.time(); S <- sheets(); names_ <- vapply(S, function(s) s$name, ""); stopifnot(!anyDuplicated(names_), all(nchar(names_) <= 31))
  rows <- lapply(S, function(s) { d <- dims(s$data); setNames(list(s$name, s$group, s$src, s$desc, pint(d[1]), pint(d[2])), c("Sheet", "Group", SRCCOL, "Description", "Data rows", "Columns")) })
  note <- function(sheet, group, srcf, desc) setNames(list(sheet, group, srcf, desc, NA, NA), c("Sheet", "Group", SRCCOL, "Description", "Data rows", "Columns"))
  notes <- c(lapply(NOT_INCLUDED, function(z) note("(not included)", "", z[1], z[2])), lapply(NOTES, function(t) note("(note)", "", "", t)))
  full <- concat(list(df_from_records(rows), df_from_records(notes)))
  base_widths <- vapply(seq_along(full$names), function(j) width_of(full$names[j], full$cols[[j]]), 1)      # the widths of the first build (the later edits keep them)
  ## apply_move_tables_181126.py: one phrase of the Index
  dj <- match("Description", full$names); full$cols[[dj]] <- lapply(full$cols[[dj]], function(v) if (is.character(v)) sub("the primary values behind Table 2 of the manuscript; n/e = drug not evaluable", "the primary values behind Supplementary Table S11 (Additional file 1); n/e = drug not evaluable", v, fixed = TRUE) else v)
  ## add_af6_sheets_201126.py: three Index rows before the notes and a sentence of the last note; apply_audit_fixes_201126.py: another sentence of the last note
  ph <- posthoc_sheets(); first_note <- which(vapply(full$cols[[1]], function(v) identical(v, "(note)"), TRUE))[1]
  for (j in seq_along(full$names)) full$cols[[j]] <- append(full$cols[[j]], lapply(ph$entries, function(e) if (is.numeric(e[[j]])) pint(e[[j]]) else e[[j]]), after = first_note - 1L)
  last <- length(full$cols[[dj]]); stopifnot(startsWith(full$cols[[dj]][[last]], "Workbook built 6 to 10 Nov 2026"))
  full$cols[[dj]][[last]] <- paste0(full$cols[[dj]][[last]], " The sheets PostHoc_per_matrix, PostHoc_summary and Fig_10_numbers were added on 20 Nov 2026 (post hoc checks and Fig. 10; Additional file 5, part 4.4).",
                                    " In the column Group and in the notes, WP followed by a number is a work package of the analysis (explained in the reading key of Additional file 5).")
  book <- list(); idx <- sheet_new("Index"); put(idx, full)
  for (j in seq_along(base_widths)) sheet_width(idx, j, base_widths[j]); sheet_width(idx, 4L, 90); sheet_width(idx, 3L, 60); book[[1]] <- idx
  for (s in S) { sh <- sheet_new(s$name); if (inherits(s$data, "af6_blocks")) put_blocks(sh, s$data) else put(sh, s$data, s$label); book[[length(book) + 1L]] <- sh }
  for (p in ph$sheets) {
    sh <- sheet_new(p$name); sheet_cells(sh, rep(1L, length(p$header)), seq_along(p$header), as.list(p$header), TRUE)
    for (i in seq_along(p$rows)) { r <- p$rows[[i]]; keep <- !vapply(r, function(v) is.null(v) || (length(v) == 1 && is.logical(v) && is.na(v)), TRUE); if (any(keep)) sheet_cells(sh, rep(i + 1L, sum(keep)), seq_along(r)[keep], r[keep], FALSE) }
    for (j in seq_along(p$widths)) sheet_width(sh, j, p$widths[j]); sh$freeze <- 1L; book[[length(book) + 1L]] <- sh
  }
  write_workbook(book, OUT)
  cat(sprintf("saved %s (%d sheets; %d KB) from the %s results in %.0f s\n", OUT, length(book), file.size(OUT) %/% 1024, SOURCE, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (SOURCE == "R" && length(no_r_file)) cat("no R result (yet) for these files; the archived files were used:\n ", paste(no_r_file, collapse = "\n  "), "\n")
}
if (sys.nframe() == 0L) main()
