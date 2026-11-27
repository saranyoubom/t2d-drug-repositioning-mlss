## compare_o6.R - 26 Nov 2026, impact check O6 (audit_261126/AUDIT_REPORT.md): compares the dry runs of run_o6_route3.sh / run_o6_heldout.sh.
## Gate: the R dry run with numpy's order (archived behaviour) must give the archived status and DE count for every series of the archived dry run; otherwise the exact-rule numbers are not trusted.
## Then: per series DE count numpy vs exact, the series whose analysability changes, series screened only under the exact rule, and the resulting analysed sets.
## Output: results_R/o6_anova_exact/o6_summary.md (and o6_route3_series.csv, o6_heldout_series.csv). Read-only on the archive.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
O <- file.path(R_PORT, "results_R", "o6_anova_exact"); ARC <- normalizePath(file.path(R_PORT, ".."))
## The gate compares with the archived FULL run (status "analysed" read as "analysable (dry run)"): the archived route 3 dry run predates the duplicate rule (protocol deviation D1, 4 Nov 2026),
## so it calls GSE53185 analysable, while the full run and the R dry runs call it a duplicate of GSE53183.
rj <- function(p) { if (!file.exists(p)) return(NULL); L <- readLines(p, warn = FALSE, encoding = "UTF-8"); L <- L[nzchar(trimws(L))]
  L <- gsub("(?<=[:\\[,] )(NaN|-?Infinity)(?=[,}\\]])", "null", L, perl = TRUE)   # Python's json writes NaN and Infinity as bare tokens (only the identifying fields and n_de are read here)
  do.call(rbind, lapply(L, function(l) { j <- fromJSON(l, simplifyVector = FALSE); data.frame(gse = j$gse, platform = j$platform, order = as.integer(j$order), status = if (identical(j$status, "analysed")) "analysable (dry run)" else j$status,
    n_de = if (is.null(j$n_de)) NA_integer_ else as.integer(j$n_de), stringsAsFactors = FALSE) })) }
ok_status <- function(s) s %in% c("analysable (dry run)", "analysed")
KNOWN <- c("GSE2634: archived 'not analysable: ValueError samples missing: ['GSM50690', 'GSM50691', 'GSM50692']' DE NA; R 'not analysable: 0 DE genes (< 34)' DE 0" =
  "bug of the archived Python reader (pandas skiprows; audit_261126 O1); not analysable in both, no number moves")
md <- c("# Impact check O6: route 3 ANOVA artefact (constant genes), exact arithmetic against numpy's order", "", sprintf("Written %s by o6_anova_exact/compare_o6.R. Dry runs only (DE counts and the frozen selection loop).", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
section <- function(title, arc_path, numpy_path, exact_path, csv_name) {
  A <- rj(arc_path); N <- rj(numpy_path); X <- rj(exact_path)
  out <- c(paste("##", title), "")
  if (is.null(A) || is.null(N) || is.null(X)) return(c(out, "missing input: run the dry runs first", ""))
  m <- merge(A[, c("gse", "status", "n_de")], N[, c("gse", "status", "n_de")], by = "gse", all = TRUE, suffixes = c("_arc", "_np"))
  bad <- m[is.na(m$status_arc) | is.na(m$status_np) | m$status_arc != m$status_np | xor(is.na(m$n_de_arc), is.na(m$n_de_np)) | (!is.na(m$n_de_arc) & !is.na(m$n_de_np) & m$n_de_arc != m$n_de_np), ]
  key <- sprintf("%s: archived '%s' DE %s; R '%s' DE %s", bad$gse, bad$status_arc, bad$n_de_arc, bad$status_np, bad$n_de_np); kn <- key %in% names(KNOWN)
  gate <- !any(!kn) && nrow(A) == nrow(N)
  out <- c(out, sprintf("Gate (R numpy order = archived full run, status and DE count of every series): %s (%d archived series, %d R series, %d differing, %d of them known by exact key)", if (gate) "PASS" else "FAIL", nrow(A), nrow(N), nrow(bad), sum(kn)))
  if (nrow(bad)) out <- c(out, "", "Differing series:", sprintf("- %s%s", key, ifelse(kn, paste0(" | KNOWN: ", KNOWN[key]), " | UNEXPLAINED")))
  t <- merge(N[, c("gse", "platform", "order", "status", "n_de")], X[, c("gse", "platform", "order", "status", "n_de")], by = c("gse", "platform", "order"), all = TRUE, suffixes = c("_numpy", "_exact"))
  t <- t[order(t$platform, t$order, t$gse), ]; write.csv(t, file.path(O, csv_name), row.names = FALSE)
  t$a_np <- ok_status(t$status_numpy); t$a_x <- ok_status(t$status_exact)
  out <- c(out, "", "| platform | analysed (numpy) | analysed (exact) | screened (numpy) | screened (exact) | DE count changed | lost analysability | gained (incl. newly screened) |", "|---|---|---|---|---|---|---|---|")
  for (g in sort(unique(t$platform))) { s <- t[t$platform == g, ]
    out <- c(out, sprintf("| %s | %d | %d | %d | %d | %d | %d | %d |", g, sum(s$a_np), sum(s$a_x), sum(!is.na(s$status_numpy)), sum(!is.na(s$status_exact)),
      sum(!is.na(s$n_de_numpy) & !is.na(s$n_de_exact) & s$n_de_numpy != s$n_de_exact), sum(s$a_np & !s$a_x), sum(!s$a_np & s$a_x))) }
  ch <- t[!is.na(t$n_de_numpy) & !is.na(t$n_de_exact) & t$n_de_numpy != t$n_de_exact, ]
  out <- c(out, "", sprintf("Series whose DE count changes (%d):", nrow(ch)), if (nrow(ch)) sprintf("- %s %s (#%d): DE %d -> %d; %s -> %s", ch$platform, ch$gse, ch$order, ch$n_de_numpy, ch$n_de_exact, ch$status_numpy, ch$status_exact) else "- none")
  sw <- t[xor(t$a_np, t$a_x) | xor(is.na(t$status_numpy), is.na(t$status_exact)), ]
  out <- c(out, "", sprintf("Series whose role changes (analysed in one run only, or screened in one run only) (%d):", nrow(sw)),
    if (nrow(sw)) sprintf("- %s %s (#%d): numpy '%s' (DE %s); exact '%s' (DE %s)", sw$platform, sw$gse, sw$order, ifelse(is.na(sw$status_numpy), "not screened", sw$status_numpy), sw$n_de_numpy,
      ifelse(is.na(sw$status_exact), "not screened", sw$status_exact), sw$n_de_exact) else "- none")
  out <- c(out, "", sprintf("Analysed set identical under both rules: %s", identical(sort(t$gse[t$a_np]), sort(t$gse[t$a_x]))), "")
  out
}
md <- c(md, section("Route 3 / development set (frozen_datasets.csv, 20 analysable per platform)", file.path(ARC, "..", "benchmark_route3", "results_per_dataset.jsonl"),
  file.path(O, "route3_dry_numpy", "dry_run_per_dataset.jsonl"), file.path(O, "route3_dry_exact", "dry_run_per_dataset.jsonl"), "o6_route3_series.csv"))
md <- c(md, section("WP2 held-out arrays (frozen_heldout_arrays.csv as archived, 20 analysable per platform)", file.path(ARC, "wp2_benchmark2", "heldout_arrays_results.jsonl"),
  file.path(O, "heldout_numpy", "heldout_arrays_dry.jsonl"), file.path(O, "heldout_exact", "heldout_arrays_dry.jsonl"), "o6_heldout_series.csv"))
rc <- file.path(O, "rnaseq_constant_genes.csv")
if (file.exists(rc)) { d <- read.csv(rc); md <- c(md, "## WP2 held-out RNA-seq", "", sprintf("Cached count files read: %d; series with at least one gene constant in all 12 samples: %d (largest number %s). The exact rule changes only such genes, so the RNA-seq DE counts are unchanged where this is 0.",
  sum(!is.na(d$n_constant)), sum(d$n_constant > 0, na.rm = TRUE), max(d$n_constant, na.rm = TRUE)), "") }
writeLines(md, file.path(O, "o6_summary.md"), useBytes = TRUE); cat(md, sep = "\n")
