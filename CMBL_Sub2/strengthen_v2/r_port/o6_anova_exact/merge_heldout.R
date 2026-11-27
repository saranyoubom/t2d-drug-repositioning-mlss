## merge_heldout.R - 27 Nov 2026, impact check O6, step 6: GSE68424 is the only held-out array series whose DE list changes under ANOVA_EXACT (10,332 -> 10,326 DE genes; o6_summary.md).
## Gate: the GSE68424 line of the numpy-order rerun (heldout68424_numpy) must equal the line of the R run (results_R/wp2/heldout_arrays_results.jsonl) in every field (numbers within 1e-9 relative, the rest identical).
## Then writes results_R/o6_anova_exact/heldout_arrays_results_exact.jsonl = the R run with the GSE68424 line replaced by the exact-rule line (all other series are unchanged under the exact rule: same
## DE count, so the same draws; o6_heldout_series.csv) and prints sets A and B of GSE68424 under both rules. Read-only on results_R/wp2.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
O <- file.path(R_PORT, "results_R", "o6_anova_exact"); RUN <- file.path(R_PORT, "results_R", "wp2", "heldout_arrays_results.jsonl"); G <- "GSE68424"
nan_null <- function(l) gsub("(?<=[:\\[,] )(NaN|-?Infinity)(?=[,}\\]])", "null", l, perl = TRUE)
pj <- function(l) fromJSON(nan_null(l), simplifyVector = FALSE)
raw <- readBin(RUN, "raw", file.info(RUN)$size); eol <- if (any(raw == as.raw(13))) "\r\n" else "\n"
L <- readLines(RUN, warn = FALSE, encoding = "UTF-8"); L <- L[nzchar(L)]; gse <- vapply(L, function(l) pj(l)$gse, "", USE.NAMES = FALSE)
if (sum(gse == G) != 1) stop("GSE68424 must appear once in the R run")
one <- function(f) { x <- readLines(f, warn = FALSE, encoding = "UTF-8"); x <- x[nzchar(x)]; if (length(x) != 1) stop("expected one line in ", f); x }
ln_np <- one(file.path(O, "heldout68424_numpy", "heldout_arrays_results.jsonl")); ln_x <- one(file.path(O, "heldout68424_exact", "heldout_arrays_results.jsonl"))
diffs <- character(0)
cmp <- function(a, b, path) {
  if (is.list(a) || is.list(b)) {
    if (!is.list(a) || !is.list(b) || length(a) != length(b) || !identical(names(a), names(b))) { diffs <<- c(diffs, sprintf("%s: structure differs", path)); return(invisible()) }
    for (i in seq_along(a)) cmp(a[[i]], b[[i]], paste0(path, "/", if (is.null(names(a))) i else names(a)[i])); return(invisible())
  }
  if (is.null(a) && is.null(b)) return(invisible())
  if (is.numeric(a) && is.numeric(b)) { if (abs(a - b) > 1e-9 * max(1, abs(a), abs(b))) diffs <<- c(diffs, sprintf("%s: %s vs %s", path, format(a, digits = 17), format(b, digits = 17))); return(invisible()) }
  if (!identical(a, b)) diffs <<- c(diffs, sprintf("%s: '%s' vs '%s'", path, paste(a, collapse = ","), paste(b, collapse = ",")))
}
cmp(pj(L[gse == G]), pj(ln_np), G)
cat(sprintf("GATE (numpy-order rerun of %s = line of the R run, every field): %s%s\n", G, if (!length(diffs)) "PASS" else "FAIL", if (length(diffs)) paste0("\n  ", paste(head(diffs, 30), collapse = "\n  ")) else ""))
if (length(diffs)) quit(status = 1)
out <- L; out[gse == G] <- ln_x
con <- file(file.path(O, "heldout_arrays_results_exact.jsonl"), "wb"); writeBin(charToRaw(paste0(paste(enc2utf8(out), collapse = eol), eol)), con); close(con)
a <- pj(ln_np); b <- pj(ln_x)
cat(sprintf("%s: status %s -> %s; DE %s -> %s\n", G, a$status, b$status, a$n_de, b$n_de))
for (tag in c("A", "B")) {
  fa <- a[[tag]]$flat; fb <- b[[tag]]$flat; ga <- unlist(a[[tag]]$genes); gb <- unlist(b[[tag]]$genes)
  cat(sprintf("set %s: genes shared %d of 17; n_de_in_set %s -> %s; n_zero_var_genes %s -> %s\n", tag, length(intersect(ga, gb)), a[[tag]]$n_de_in_set, b[[tag]]$n_de_in_set, fa$n_zero_var_genes, fb$n_zero_var_genes))
  for (k in c("observed", "z", "perm_p", "T_r2_p", "T_p_p", "T_eig_p", "share_between", "adj_any", "adj_all4", "retained_any", "retained_all4"))
    cat(sprintf("  %-14s %s -> %s\n", k, format(if (is.null(fa[[k]])) NA else fa[[k]], digits = 6), format(if (is.null(fb[[k]])) NA else fb[[k]], digits = 6)))
}
cat(sprintf("written: heldout_arrays_results_exact.jsonl (%d lines, %s line ends)\n", length(out), if (eol == "\r\n") "CRLF" else "LF"))
