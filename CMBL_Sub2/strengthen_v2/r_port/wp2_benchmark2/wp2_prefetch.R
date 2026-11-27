## wp2_prefetch.R - download (paced, <= 3 requests per second, cached) the GEO / NCBI files the R run of WP2 reads, so that the analysis scripts run from the local cache.
## What is fetched: the series matrices of the 58 part-1 series (development arms), the 59 held-out array series, the platform annotations GPL570 and GPL1261, and the NCBI-generated human RNA-seq count
## files of the series that had counts on 5 Nov 2026 (the series without counts answered HTTP 404 then; the R run of run_rnaseq.R asks again). The series list is read from the archived Python logs.
## Files land in <WP2_DATA_ROOT>/bench_data/{matrices,platforms}, v2_data/arrays_heldout/matrices and v2_data/rnaseq/counts (untrusted data: read as text/numbers only).
## Usage: Rscript wp2_prefetch.R [dev] [heldout] [rnaseq]   (default: all three)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "bench_pipeline.R")); source(file.path(HERE, "bench_data_io.R"))
what <- commandArgs(trailingOnly = TRUE); if (!length(what)) what <- c("dev", "heldout", "rnaseq")
plat <- file.path(DATA, "bench_data", "platforms"); failed <- character(0)
try_get <- function(label, expr) { r <- tryCatch({ force(expr); "ok" }, error = function(e) conditionMessage(e)); cat(sprintf("%s %s\n", format(Sys.time(), "%H:%M:%S"), if (identical(r, "ok")) paste("ok  ", label) else paste("FAIL", label, "-", substr(r, 1, 100)))); if (!identical(r, "ok")) failed <<- c(failed, label); invisible(r) }
if ("dev" %in% what) {
  io_configure(file.path(DATA, "bench_data", "matrices"), plat); recs <- read_jsonl(file.path(P1, "results_per_dataset.jsonl"))
  for (r in recs) if (!grepl("duplicate", r$status)) try_get(paste(r$gse, r$platform), series_path(r$gse, r$platform))
  for (g in c("GPL570", "GPL1261")) try_get(g, annotation(g))
}
if ("heldout" %in% what) {
  io_configure(file.path(DATA, "v2_data", "arrays_heldout", "matrices"), plat); recs <- read_jsonl(file.path(ARC, "heldout_arrays_results.jsonl"))
  for (r in recs) if (!grepl("duplicate", r$status)) try_get(paste(r$gse, r$platform), series_path(r$gse, r$platform))
}
if ("rnaseq" %in% what) {
  recs <- read_jsonl(file.path(ARC, "heldout_rnaseq_results.jsonl"))
  for (r in recs) if (r$organism == "human" && !grepl("HTTP 404|duplicate", r$status)) try_get(paste(r$gse, "counts"), { f <- fetch_counts(DATA, r$gse, "human"); if (is.null(f$path)) stop(f$why) })
}
cat(sprintf("prefetch finished; %d failures%s\n", length(failed), if (length(failed)) paste0(": ", paste(failed, collapse = ", ")) else ""))
