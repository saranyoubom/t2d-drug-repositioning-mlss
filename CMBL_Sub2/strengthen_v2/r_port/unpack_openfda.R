## unpack_openfda.R - the saved openFDA responses of WP9 and WP12 for run_all_v2.R (added 27 Nov 2026).
## In the repository the responses are packed (<wp>/results/raw_responses.zip, see README_raw_responses.md there); in the project folder they are a folder <wp>/results/raw.
## For each of wp9_safety and wp12_database_annotations: if the folder raw/ exists nothing is done; otherwise the archive is unpacked into results_R/openfda_raw/<wp>/
## (the archived folders of strengthen_v2/ are never written). For WP9 the folder digest is checked against the one written before the analysis (results/RAW_SHA256.txt:
## SHA-256 over the lines "file name:SHA-256 of the file", files sorted by name); a mismatch stops the run.
## Usage: Rscript unpack_openfda.R
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); V2 <- normalizePath(file.path(R_PORT, ".."))
sha256_file <- function(p) { if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_); digest::digest(file = p, algo = "sha256") }
sha256_text <- function(s) { if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_); digest::digest(charToRaw(s), algo = "sha256", serialize = FALSE) }
for (wp in c("wp9_safety", "wp12_database_annotations")) {
  res <- file.path(V2, wp, "results"); raw <- file.path(res, "raw"); zip <- file.path(res, "raw_responses.zip"); out <- file.path(R_PORT, "results_R", "openfda_raw", wp)
  if (dir.exists(raw)) { cat(sprintf("%s: folder %s exists, nothing to unpack\n", wp, raw)); next }
  if (!file.exists(zip)) stop(sprintf("%s: neither %s nor %s exists", wp, raw, zip))
  if (dir.exists(out) && length(list.files(out))) { cat(sprintf("%s: already unpacked in %s (%d files)\n", wp, out, length(list.files(out)))) } else {
    dir.create(out, showWarnings = FALSE, recursive = TRUE); utils::unzip(zip, exdir = out, junkpaths = TRUE); cat(sprintf("%s: unpacked %d files into %s\n", wp, length(list.files(out)), out)) }
  if (wp == "wp9_safety" && file.exists(file.path(res, "RAW_SHA256.txt"))) {
    want <- strsplit(readLines(file.path(res, "RAW_SHA256.txt"), n = 1, warn = FALSE), " ", fixed = TRUE)[[1]][1]
    fs <- sort(list.files(out, pattern = "[.]json$"), method = "radix"); h <- vapply(file.path(out, fs), sha256_file, "")
    if (anyNA(h)) { cat("wp9_safety: package digest not installed, folder digest not checked\n"); next }
    got <- sha256_text(paste0(paste0(fs, ":", h, "\n"), collapse = ""))
    if (!identical(got, want)) stop(sprintf("wp9_safety: folder digest %s differs from RAW_SHA256.txt %s", got, want))
    cat(sprintf("wp9_safety: folder digest of the %d responses equals RAW_SHA256.txt (%s)\n", length(fs), substr(want, 1, 12)))
  }
}
