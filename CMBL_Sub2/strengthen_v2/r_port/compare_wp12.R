## compare_wp12.R - agreement of wp12_annotations.R with the archived Python results (wp12_database_annotations/results). Usage: Rscript compare_wp12.R   (WP12_OUT = folder of the R run, default results_R/wp12)
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(R_PORT, "agree.R"))
ARC <- normalizePath(file.path(R_PORT, "..", "wp12_database_annotations", "results")); OUT <- Sys.getenv("WP12_OUT", file.path(R_PORT, "results_R", "wp12"))
for (f in c("hub_annotation.csv", "mlss_quartile_classes.csv", "market_status_by_drug.csv", "market_status_nda_bla_applications.csv")) agree_csv(file.path(ARC, f), file.path(OUT, f), tol = 1e-12, label = f)
agree_json(file.path(ARC, "summary.json"), file.path(OUT, "summary.json"), tol = 1e-12, label = "summary.json")
for (f in c("hub_annotation.csv", "mlss_quartile_classes.csv", "market_status_by_drug.csv", "market_status_nda_bla_applications.csv", "summary.json")) {   # byte comparison (CRLF, quoting, number formats)
  a <- readBin(file.path(ARC, f), "raw", file.size(file.path(ARC, f))); b <- readBin(file.path(OUT, f), "raw", file.size(file.path(OUT, f)))
  cat(sprintf("%-6s %-42s %s\n", if (identical(a, b)) "SAME" else "BYTES", f, if (identical(a, b)) "byte-identical" else sprintf("bytes differ (%d vs %d)", length(a), length(b))))
}
