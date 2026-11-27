## compare_reanalysis.R - reanalysis_17genes: the R results (results_R/reanalysis_17genes/) against the archived Python results (Analysis_outputs/reanalysis_17genes/).
## Usage: Rscript compare_reanalysis.R   (RA_OUT = folder of the R run, RA_ARCHIVE = archive folder). Lists that Python wrote in set-iteration order (hash order, differs between runs) are compared as sets.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(HERE, "agree.R"))
ARC <- Sys.getenv("RA_ARCHIVE", normalizePath(file.path(HERE, "..", "..", "reanalysis_17genes"))); OUT <- Sys.getenv("RA_OUT", file.path(HERE, "results_R", "reanalysis_17genes"))
res <- list(); have <- function(f) file.exists(file.path(OUT, f))
set_sort <- function(l) l[order(vapply(l, function(p) paste(unlist(p), collapse = "|"), ""))]
for (n in c(20, 17)) {                                                    # pipeline.R
  f <- sprintf("pipeline_%dgenes.json", n)
  if (have(f)) {
    a <- fromJSON(file.path(ARC, f), simplifyVector = FALSE); b <- fromJSON(file.path(OUT, f), simplifyVector = FALSE)
    for (k in c("adj_all_four", "sign_inconsistent")) { a[[k]] <- set_sort(a[[k]]); b[[k]] <- set_sort(b[[k]]) }
    res[[f]] <- agree_json(a, b, tol = 2e-9, label = f)$ok
  }
}
for (f in c("perm_detail.json", "table_counts.json", "denom_sens.json", "gapdh_exact.json")) if (have(f)) res[[f]] <- agree_json(file.path(ARC, f), file.path(OUT, f), tol = 2e-9)$ok
if (have("gapdh_exact_pairs.csv")) res[["gapdh_exact_pairs.csv"]] <- agree_csv(file.path(ARC, "gapdh_exact_pairs.csv"), file.path(OUT, "gapdh_exact_pairs.csv"))$ok
if (have("perm_arrays_17.csv")) res[["perm_arrays_17"]] <- agree_csv(file.path(HERE, "oracle", "perm_arrays_17_archive.csv"), file.path(OUT, "perm_arrays_17.csv"), tol = 0, label = "perm_arrays_17.npz (as CSV)")$ok
if (have("network_metrics.rds")) res[["network_metrics_17genes"]] <- agree_csv(file.path(ARC, "network_metrics_17genes.csv"), readRDS(file.path(OUT, "network_metrics.rds"))[["17|consensus"]], label = "network_metrics_17genes.csv")$ok
if (have("cfx_summary.json")) res[["cfx_summary.json"]] <- agree_json(file.path(ARC, "cfx_summary.json"), file.path(OUT, "cfx_summary.json"), tol = 2e-9)$ok
if (have("cfx_long.json")) res[["cfx_long.json"]] <- agree_json(file.path(HERE, "oracle", "cfx_long_python.json"), file.path(OUT, "cfx_long.json"), tol = 2e-9, label = "cfx_long.json (Python melt.py, rerun 24 Nov)")$ok   # cfx_long.json is not archived: oracle = a fresh run of melt.py
if (have("qc26_wells.csv")) res[["qc26_wells.csv"]] <- agree_csv(file.path(ARC, "qc26_wells.csv"), file.path(OUT, "qc26_wells.csv"), tol = 1e-6)$ok      # logistic fits: MINPACK in both, tolerance 1e-6
if (have("efficiency_sensitivity_group_level.csv")) res[["efficiency_sensitivity_group_level.csv"]] <- agree_csv(file.path(ARC, "efficiency_sensitivity_group_level.csv"), file.path(OUT, "efficiency_sensitivity_group_level.csv"))$ok
PC <- Sys.getenv("PRIMER_ARCHIVE", normalizePath(file.path(HERE, "..", "..", "primer_check")))                  # primer_evidence.csv is archived in Analysis_outputs/primer_check
if (have("primer_evidence.csv") && file.exists(file.path(PC, "primer_evidence.csv"))) res[["primer_evidence.csv"]] <- agree_csv(file.path(PC, "primer_evidence.csv"), file.path(OUT, "primer_evidence.csv"))$ok
cat(sprintf("\nreanalysis_17genes so far (%d comparisons): %s\n", length(res), if (length(res) && all(unlist(res))) "ALL AGREE" else "DIFFERENCES FOUND"))
