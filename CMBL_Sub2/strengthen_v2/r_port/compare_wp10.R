## compare_wp10.R - WP10 part A (string_sensitivity): the R results (results_R/wp10/) against the archived Python results (wp10_versions/results/).
## Usage: Rscript compare_wp10.R   (R_OUT_DIR = folder of the R run, WP10_ARCHIVE = archive folder). Also compares the retrieved STRING file when `retrieve` was run in R.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "agree.R"))
ARC <- Sys.getenv("WP10_ARCHIVE", normalizePath(file.path(HERE, "..", "wp10_versions", "results"))); OUT <- Sys.getenv("R_OUT_DIR", file.path(HERE, "results_R", "wp10"))
res <- list(); have <- function(f) file.exists(file.path(OUT, f))
for (f in c("string_gate.json", "string_summary.json")) if (have(f)) res[[f]] <- agree_json(file.path(ARC, f), file.path(OUT, f), tol = 2e-9)$ok
if (have("string_classification.csv")) res[["string_classification.csv"]] <- agree_csv(file.path(ARC, "string_classification.csv"), file.path(OUT, "string_classification.csv"))$ok
if (have("string_v12_5_network.tsv")) { a <- trimws(readLines(file.path(ARC, "string_v12_5_network.tsv"), warn = FALSE)); b <- trimws(readLines(file.path(OUT, "string_v12_5_network.tsv"), warn = FALSE)); same <- identical(a, b)
  cat(sprintf("%-6s %-42s %d lines (archive) / %d lines (R)\n", if (same) "AGREE" else "DIFFER", "string_v12_5_network.tsv (retrieved again)", length(a), length(b))); res[["string_v12_5_network.tsv"]] <- same }
## part B (wp10_drugcentral_sensitivity.R): DrugCentral file of 25 Sep 2026 against the reported file
if (have("drugcentral_summary.json")) res[["drugcentral_summary.json"]] <- agree_json(file.path(ARC, "drugcentral_summary.json"), file.path(OUT, "drugcentral_summary.json"), tol = 2e-9)$ok
for (f in c("mlss_comparison.csv", "mlss_new_all_combinations.csv", "pic50_changes.csv")) if (have(f)) res[[f]] <- agree_csv(file.path(ARC, f), file.path(OUT, f), tol = 2e-9)$ok
cat(sprintf("\nWP10 so far (%d comparisons): %s\n", length(res), if (length(res) && all(unlist(res))) "ALL AGREE" else "DIFFERENCES FOUND"))
