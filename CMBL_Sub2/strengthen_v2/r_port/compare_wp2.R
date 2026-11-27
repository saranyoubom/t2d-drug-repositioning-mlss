## compare_wp2.R - agreement of the R port of wp2_benchmark2 with the archived Python results.
## Usage: Rscript compare_wp2.R <section ...>   sections: post (apply_dev + summarize_validation outputs), calib, design, select, heldout, rnaseq, core (dev_core_results.jsonl), gate
## WP2_OUT = folder of the R run (default results_R/wp2), WP2_ARCHIVE = archived Python results.
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(R_PORT, "wp2_benchmark2", "wp2_common.R"))
secs <- commandArgs(trailingOnly = TRUE); if (!length(secs)) secs <- "post"
rd_lines <- function(f) { x <- readBin(f, "raw", file.size(f)); strsplit(rawToChar(x), "\r\n", fixed = TRUE)[[1]] }
cmp_md <- function(f, quiet = FALSE) {                       # text files line by line (CRLF as written by Python); reports the number of differing lines
  a <- rd_lines(file.path(ARC, f)); b <- rd_lines(file.path(OUT, f)); n <- max(length(a), length(b)); a <- c(a, rep(NA, n - length(a))); b <- c(b, rep(NA, n - length(b))); bad <- which(is.na(a) | is.na(b) | a != b)
  cat(sprintf("%-6s %-42s %d of %d lines identical%s\n", if (!length(bad)) "AGREE" else "DIFFER", f, n - length(bad), n, if (length(bad)) sprintf("; first difference at line %d:\n         Python: %s\n         R     : %s", bad[1], substr(a[bad[1]], 1, 220), substr(b[bad[1]], 1, 220)) else ""))
  invisible(length(bad) == 0)
}
cmp_csv <- function(f, tol = 1e-12, key = NULL) agree_csv(file.path(ARC, f), file.path(OUT, f), tol = tol, key = key, label = f)
cmp_json <- function(f, tol = 1e-12) agree_json(file.path(ARC, f), file.path(OUT, f), tol = tol, label = f)
cmp_jsonl <- function(py_path, r_path, label, tol = 1e-9, keyf = function(x) paste(x$gse, x$tag)) {      # records matched by key, every number compared (relative tolerance tol)
  a <- read_jsonl(py_path); b <- read_jsonl(r_path); bm <- setNames(b, vapply(b, keyf, "")); n <- 0; bad <- 0; worst <- 0; msg <- character(0)
  for (r in a) {
    k <- keyf(r); if (is.null(bm[[k]])) { bad <- bad + 1; msg <- c(msg, paste("missing in R:", k)); next }
    res <- agree_json(r, bm[[k]], tol = tol, label = k, quiet = TRUE); n <- n + res$n; worst <- max(worst, res$worst); if (!res$ok) { bad <- bad + 1; msg <- c(msg, paste0(k, ": ", paste(head(res$msgs, 3), collapse = " | "))) }
  }
  cat(sprintf("%-6s %-42s %d of %d records compared (R has %d), %d values, largest relative difference %.2g, %d records differ%s\n", if (!bad && length(a) == length(b)) "AGREE" else "DIFFER", label, length(a) - bad, length(a), length(b), n, worst, bad,
              if (length(msg)) paste0("\n         ", paste(head(msg, 6), collapse = "\n         ")) else ""))
  invisible(list(ok = !bad, n = n, worst = worst))
}
for (s in secs) {
  if (s == "gate") {
    cat("== gate A step 1: reload of the development series from GEO (gate_part1.R extract) ==\n")
    a <- fromJSON(file.path(ARC, "compact", "extract_check.json"), simplifyVector = FALSE); b <- fromJSON(file.path(OUT, "compact", "extract_check.json"), simplifyVector = FALSE)
    cat(sprintf("rows processed: Python %d, R %d; identical rows (status, DE counts, expressed counts, scale, q99, gene lists A and B): Python %d, R %d; analysed per platform: Python %s, R %s\n", a$n_rows_processed, b$n_rows_processed, a$n_all_ok, b$n_all_ok,
                paste(sprintf("%s %d", names(a$analysed_per_platform), unlist(a$analysed_per_platform)), collapse = ", "), paste(sprintf("%s %d", names(b$analysed_per_platform), unlist(b$analysed_per_platform)), collapse = ", ")))
    if (length(b$series_not_readable_from_GEO_archived_record_used)) cat("series that GEO no longer serves (archived record used):", paste(unlist(b$series_not_readable_from_GEO_archived_record_used), collapse = ", "), "\n")
    bad <- Filter(function(c) !c$all_ok, b$checks); for (c in bad) cat(sprintf("  mismatch %s: %s\n", c$gse, paste(names(c$ok)[!unlist(c$ok)], collapse = ", ")))
    pa <- fromJSON(file.path(ARC, "compact", "dev_prepared.json"), simplifyVector = FALSE); pb <- fromJSON(file.path(OUT, "compact", "dev_prepared.json"), simplifyVector = FALSE)
    cat(sprintf("%-6s dev_prepared.json: %d series in Python, %d in R; same series, same order: %s\n", if (identical(names(pa), names(pb))) "AGREE" else "DIFFER", length(pa), length(pb), identical(names(pa), names(pb))))
    agree_json(pa, pb, tol = 1e-12, label = "dev_prepared.json (genes, DE counts)")
    Za <- fromJSON(file.path(R_PORT, "oracle", "dev_sets_archive.json"), simplifyVector = FALSE); Zb <- fromJSON(file.path(OUT, "compact", "dev_sets.json"), simplifyVector = FALSE)
    agree_json(Za[intersect(names(Za), names(Zb))], Zb[intersect(names(Za), names(Zb))], tol = 1e-12, label = sprintf("12 x 17 matrices (%d of %d)", length(intersect(names(Za), names(Zb))), length(Za)))
  }
  if (s == "core") { cat("== analyse_matrix on the 74 development matrices (gate_part1.R consistency) ==\n"); cmp_jsonl(file.path(ARC, "compact", "dev_core_results.jsonl"), file.path(OUT, "compact", "dev_core_results.jsonl"), "dev_core_results.jsonl") }
  if (s == "heldout") { cat("== held-out arrays (run_heldout.R) ==\n"); for (f in c("heldout_arrays_dry.jsonl", "heldout_arrays_results.jsonl")) if (file.exists(file.path(OUT, f))) cmp_jsonl(file.path(ARC, f), file.path(OUT, f), f, keyf = function(x) x$gse) }
  if (s == "rnaseq") { cat("== held-out RNA-seq (run_rnaseq.R) ==\n"); for (f in c("heldout_rnaseq_dry.jsonl", "heldout_rnaseq_results.jsonl")) if (file.exists(file.path(OUT, f))) cmp_jsonl(file.path(ARC, f), file.path(OUT, f), f, keyf = function(x) x$gse) }
  if (s == "select") {
    cat("== selection (select_heldout.R, select_rnaseq.R) ==\n")
    for (f in c("eligibility_heldout.csv", "frozen_heldout_arrays.csv", "eligibility_heldout_rnaseq.csv", "frozen_heldout_rnaseq.csv")) if (file.exists(file.path(OUT, f))) {
      a <- readBin(file.path(ARC, f), "raw", file.size(file.path(ARC, f))); b <- readBin(file.path(OUT, f), "raw", file.size(file.path(OUT, f)))
      cat(sprintf("%-6s %-42s %s (sha256 archive %s, R %s)\n", if (identical(a, b)) "AGREE" else "DIFFER", f, if (identical(a, b)) "byte-identical" else "bytes differ", substr(sha256_file(file.path(ARC, f)), 1, 12), substr(sha256_file(file.path(OUT, f)), 1, 12)))
    }
  }
  if (s == "calib") {
    cat("== calibration (calibrate_t2.R summarize) ==\n")
    for (f in c("calibration_replicates.csv", "calibration_per_condition.csv", "calibration_tables.csv", "calibration_by_platform.csv", "calibration_sensitivity_no_degenerate.csv")) if (file.exists(file.path(OUT, f))) cmp_csv(f)
    if (file.exists(file.path(OUT, "calibration_verdict.txt"))) cmp_md("calibration_verdict.txt")
  }
  if (s == "design") {
    cat("== design table (design_table.R summarize) ==\n")
    if (file.exists(file.path(OUT, "design_table.csv"))) cmp_csv("design_table.csv"); if (file.exists(file.path(OUT, "design_table.md"))) cmp_md("design_table.md")
  }
  if (s == "post") {
    cat("== post-processing (apply_dev.R, summarize_validation.R) ==\n")
    for (f in c("dev_application_per_series.csv", "dev_endpoints_long.csv", "validation_per_series.csv", "validation_endpoints_long.csv", "WP2_summary_table.csv")) if (file.exists(file.path(OUT, f))) cmp_csv(f)
    for (f in c("dev_application_summary.json", "validation_summary.json")) if (file.exists(file.path(OUT, f))) cmp_json(f)
    for (f in c("dev_application_summary.md", "validation_summary.md", "WP2_summary_table.md")) if (file.exists(file.path(OUT, f))) cmp_md(f)
  }
}
