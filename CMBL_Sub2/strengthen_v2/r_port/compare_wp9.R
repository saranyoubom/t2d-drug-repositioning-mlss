## compare_wp9.R - agreement of wp9_safety.R with the archived Python results (wp9_safety/results). Usage: Rscript compare_wp9.R   (WP9_OUT = folder of the R run, default results_R/wp9)
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(R_PORT, "agree.R"))
ARC <- normalizePath(file.path(R_PORT, "..", "wp9_safety", "results")); OUT <- Sys.getenv("WP9_OUT", file.path(R_PORT, "results_R", "wp9")); have <- function(f) file.exists(file.path(OUT, f))
for (f in c("drug_counts.csv", "signals_single.csv", "pairs_all77.csv", "pairs_top15.csv", "pairs_long.csv", "fallback_signals_single.csv", "fallback_pairs.csv")) if (have(f)) agree_csv(file.path(ARC, f), file.path(OUT, f), tol = 1e-12, label = f)
for (f in c("summary.json", "fallback_summary.json")) if (have(f)) agree_json(file.path(ARC, f), file.path(OUT, f), tol = 1e-12, label = f)
for (f in c("drug_counts.csv", "signals_single.csv", "pairs_all77.csv", "pairs_top15.csv", "pairs_long.csv", "summary.json", "fallback_signals_single.csv", "fallback_pairs.csv", "fallback_summary.json", "RAW_SHA256.txt")) if (have(f) && file.exists(file.path(ARC, f))) {
  a <- readBin(file.path(ARC, f), "raw", file.size(file.path(ARC, f))); b <- readBin(file.path(OUT, f), "raw", file.size(file.path(OUT, f)))
  cat(sprintf("%-6s %-42s %s\n", if (identical(a, b)) "SAME" else "BYTES", f, if (identical(a, b)) "byte-identical" else sprintf("bytes differ (%d vs %d)", length(a), length(b))))
}
