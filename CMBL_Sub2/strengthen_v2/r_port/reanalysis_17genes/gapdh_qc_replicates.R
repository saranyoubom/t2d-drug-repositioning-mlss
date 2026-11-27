## gapdh_qc_replicates.R - R port of reanalysis_17genes/gapdh_qc_replicates.py: plateau RFU of the 26 August target wells by plate row (qc26_wells.csv of gapdh_qc_amplification.R), the target wells
## with a plateau below 700 RFU, and the SD of the Gapdh Cq of the replicate wells per sample group in every run. Prints only (the Python script writes no file).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "cfx.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes"))
txt <- function(x) if (is.null(x)) NA_character_ else as.character(x)
d <- read.csv(file.path(OUT, "qc26_wells.csv"), stringsAsFactors = FALSE); d$row <- substr(d$well, 1, 1); t <- d[d$gene != "Gapdh", ]
cat("plateau RFU by plate row (target genes, 26 Aug):\n")
by_row <- do.call(rbind, lapply(py_sort(unique(t$row)), function(r) { v <- t$plateau[t$row == r]; data.frame(row = r, count = length(v), median = round(median(v)), min = round(min(v))) })); print(by_row, row.names = FALSE)
cat("\ntarget wells with plateau < 700 RFU:\n"); low <- t[t$plateau < 700, c("gene", "well", "sample", "cq", "plateau")]; low$cq <- round(low$cq, 2); low$plateau <- round(low$plateau, 2); print(low, row.names = FALSE)
out_rows <- list()
for (run in cfx_runs()) {
  sh <- tryCatch(cfx_sheet_rows(run, "Quantification Cq Results"), error = function(e) NULL); rows <- sh[["0"]]; if (is.null(rows)) next
  h <- vapply(rows[[1]], txt, ""); iT <- match("Target", h); iS <- match("Sample", h); iC <- match("Cq", h)
  recs <- Filter(function(r) { x <- r[[iT]]; !is.null(x) && !(is.character(x) && !nzchar(x)) && !(is.numeric(x) && x == 0) }, rows[-1])
  df <- data.frame(t = vapply(recs, function(r) txt(r[[iT]]), ""), s = vapply(recs, function(r) txt(r[[iS]]), ""), cq = vapply(recs, function(r) if (is.numeric(r[[iC]])) r[[iC]] else suppressWarnings(as.numeric(txt(r[[iC]]))), 1), stringsAsFactors = FALSE)
  df <- df[grepl("gapdh", tolower(df$t), fixed = TRUE), ]
  for (s in py_sort(unique(df$s[!is.na(df$s)]))) {
    v <- df$cq[!is.na(df$s) & df$s == s]; vv <- v[!is.na(v)]
    out_rows[[length(out_rows) + 1]] <- data.frame(run = substr(run, 1, 28), sample = s, n = length(vv), mean_Cq = round(mean(vv), 2), SD = if (length(vv) > 1) round(sd(vv), 2) else NA_real_, range = round(max(vv) - min(vv), 2), stringsAsFactors = FALSE)
  }
}
o <- do.call(rbind, out_rows); cat("\nGapdh replicate wells per group, SD (cycles):\n"); print(o, row.names = FALSE)
cat("\nmedian SD of all groups:", median(o$SD, na.rm = TRUE), "| groups with SD > 0.5:", sum(o$SD > 0.5, na.rm = TRUE), "of", nrow(o), "\n")
