## compare_route3.R - benchmark_route3: the R results (results_R/benchmark_route3/) against the archived Python results (Analysis_outputs/benchmark_route3/).
## Usage: Rscript compare_route3.R   (RB_OUT = folder of the R run, default results_R/benchmark_route3; the agreement report uses the full run results_R/benchmark_route3_full; RB_ARCHIVE = archive folder)
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "agree.R"))
ARC <- Sys.getenv("RB_ARCHIVE", normalizePath(file.path(HERE, "..", "..", "benchmark_route3"))); OUT <- Sys.getenv("RB_OUT", file.path(HERE, "results_R", "benchmark_route3"))
res <- list(); have <- function(f) file.exists(file.path(OUT, f))
for (f in c("processing_log.csv", "per_dataset_results.csv")) if (have(f)) res[[f]] <- agree_csv(file.path(ARC, f), file.path(OUT, f))$ok
if (have("summary_endpoints.json")) res[["summary_endpoints.json"]] <- agree_json(file.path(ARC, "summary_endpoints.json"), file.path(OUT, "summary_endpoints.json"), tol = 2e-9)$ok
if (have("results_per_dataset.jsonl")) {                       # the R benchmark run: compared dataset by dataset (only the series that were run)
  a <- lapply(readLines(file.path(ARC, "results_per_dataset.jsonl"), warn = FALSE), fromJSON, simplifyVector = FALSE); names(a) <- vapply(a, function(r) r$gse, "")
  b <- lapply(readLines(file.path(OUT, "results_per_dataset.jsonl"), warn = FALSE), fromJSON, simplifyVector = FALSE)
  for (r in b) res[[paste0("series ", r$gse)]] <- agree_json(a[[r$gse]], r, tol = 2e-9, label = paste0("results_per_dataset.jsonl ", r$gse))$ok
}
if (have("archived_records_used.txt")) for (l in readLines(file.path(OUT, "archived_records_used.txt"), warn = FALSE)) if (nzchar(trimws(l))) {   # RB_ARCHIVED fallback of run_benchmark.R (26 Nov 2026): such a record was copied from the archive, not computed by R, so it can never count as agreement
  cat(sprintf("MISSING results_per_dataset.jsonl %s not computed by R: the archived record was copied (%s)\n", sub(":.*$", "", l), l)); res[[paste0("archived ", l)]] <- FALSE }
cat(sprintf("\nbenchmark_route3 so far (%d comparisons): %s\n", length(res), if (length(res) && all(unlist(res))) "ALL AGREE" else "DIFFERENCES FOUND"))
