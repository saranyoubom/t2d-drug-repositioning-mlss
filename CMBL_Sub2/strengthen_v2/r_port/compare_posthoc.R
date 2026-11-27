## compare_posthoc.R - agreement of the R port of the post hoc checks (post_hoc_191126.R) with the archived Python results (ars_review_191126/post_hoc/results_<mode>).
## Usage: Rscript compare_posthoc.R <mode ...>   modes: dev, ins1, validation (default: all three). POSTHOC_OUT = folder of the R run (default results_R/post_hoc), POSTHOC_ARCHIVE = archived Python results.
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(R_PORT, "wp2_benchmark2", "wp2_common.R"))
PA <- Sys.getenv("POSTHOC_ARCHIVE", normalizePath(file.path(R_PORT, "..", "ars_review_191126", "post_hoc"))); PR <- Sys.getenv("POSTHOC_OUT", file.path(R_PORT, "results_R", "post_hoc"))
modes <- commandArgs(trailingOnly = TRUE); if (!length(modes)) modes <- c("dev", "ins1", "validation")
rd_lines <- function(f) { x <- readBin(f, "raw", file.size(f)); strsplit(rawToChar(x), "\r\n", fixed = TRUE)[[1]] }
cmp_md <- function(a_path, b_path, label) {
  if (!file.exists(a_path) || !file.exists(b_path)) { cat(sprintf("MISSING %-52s %s\n", label, if (!file.exists(a_path)) "(archive)" else "(R run)")); return(invisible(FALSE)) }
  a <- rd_lines(a_path); b <- rd_lines(b_path); n <- max(length(a), length(b)); a <- c(a, rep(NA, n - length(a))); b <- c(b, rep(NA, n - length(b))); bad <- which(is.na(a) | is.na(b) | a != b)
  cat(sprintf("%-6s %-52s %d of %d lines identical%s\n", if (!length(bad)) "AGREE" else "DIFFER", label, n - length(bad), n,
              if (length(bad)) sprintf("; first difference at line %d:\n         Python: %s\n         R     : %s", bad[1], substr(a[bad[1]], 1, 230), substr(b[bad[1]], 1, 230)) else ""))
  invisible(!length(bad))
}
cmp <- function(mode, f, kind, tol = 1e-9) {
  a <- file.path(PA, paste0("results_", mode), f); b <- file.path(PR, paste0("results_", mode), f); lab <- paste0(mode, "/", f)
  if (!file.exists(a) || !file.exists(b)) { cat(sprintf("MISSING %-52s %s\n", lab, if (!file.exists(a)) "(archive)" else "(R run)")); return(invisible()) }
  if (kind == "csv") agree_csv(a, b, tol = tol, label = lab) else if (kind == "json") agree_json(a, b, tol = tol, label = lab) else cmp_md(a, b, lab)
}
for (m in modes) {
  cat(sprintf("== post hoc %s ==\n", m))
  if (m == "dev") for (f in c("per_matrix.csv", "summary.json", "summary.md")) cmp(m, f, sub(".*[.]", "", f))
  if (m == "ins1") for (f in c("summary.json", "summary.md")) cmp(m, f, sub(".*[.]", "", f))
  if (m == "validation") for (f in c("per_matrix.csv", "summary.json", "summary.md", "summary_arrays_GPL1261.json", "summary_arrays_GPL1261.md", "summary_arrays_GPL570.json", "summary_arrays_GPL570.md", "summary_rnaseq_human.json", "summary_rnaseq_human.md")) cmp(m, f, sub(".*[.]", "", f))
}
