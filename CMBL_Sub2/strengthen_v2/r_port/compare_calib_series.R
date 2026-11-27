## compare_calib_series.R - the per-series results of calibrate_t2.R run (results_R/wp2/compact/calib_series/*.rds) against the archived Python files (wp2_benchmark2/compact/calib_series/*.json), value by value.
## Usage: Rscript compare_calib_series.R [series ...]   (default: every series that has an RDS file)
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); suppressMessages(library(jsonlite))
ARC <- normalizePath(file.path(R_PORT, "..", "wp2_benchmark2", "compact", "calib_series")); OUT <- Sys.getenv("WP2_OUT", file.path(R_PORT, "results_R", "wp2")); SER <- file.path(OUT, "compact", "calib_series")
ids <- commandArgs(trailingOnly = TRUE); if (!length(ids)) ids <- sub("[.]rds$", "", list.files(SER, pattern = "[.]rds$"))
tot_n <- 0; tot_bad <- 0
for (gse in ids) {
  rows <- readRDS(file.path(SER, paste0(gse, ".rds"))); arch <- fromJSON(file.path(ARC, paste0(gse, ".json")), simplifyVector = FALSE); bad <- 0; nc <- 0; worst <- 0; first <- NULL
  if (length(rows) != length(arch)) { cat(sprintf("DIFFER %-10s rows R %d, Python %d\n", gse, length(rows), length(arch))); next }
  for (i in seq_along(rows)) for (k in names(arch[[i]])) {
    a <- arch[[i]][[k]]; b <- rows[[i]][[k]]; nc <- nc + 1
    if (is.null(a)) { if (!is.na(b)) { bad <- bad + 1; if (is.null(first)) first <- sprintf("row %d %s: Python null, R %s", i, k, format(b)) } }
    else if (is.character(a)) { if (!identical(a, b)) { bad <- bad + 1; if (is.null(first)) first <- sprintf("row %d %s: Python '%s', R '%s'", i, k, a, b) } }
    else { d <- abs(a - b) / max(1, abs(b)); worst <- max(worst, d); if (is.na(d) || d > 1e-12) { bad <- bad + 1; if (is.null(first)) first <- sprintf("row %d %s: Python %s, R %s", i, k, format(a, digits = 15), format(b, digits = 15)) } }
  }
  tot_n <- tot_n + nc; tot_bad <- tot_bad + bad
  cat(sprintf("%-6s %-10s %d rows, %d values compared, %d differ, largest relative difference %.2g%s\n", if (bad == 0) "AGREE" else "DIFFER", gse, length(rows), nc, bad, worst, if (is.null(first)) "" else paste0("; first: ", first)))
}
cat(sprintf("\n%d series, %d values compared, %d differ\n", length(ids), tot_n, tot_bad))
