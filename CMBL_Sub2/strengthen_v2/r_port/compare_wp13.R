## compare_wp13.R - WP13: the R results (results_R/wp13/) against the archived Python results (wp13_gapdh_corrected/results/). Usage: Rscript compare_wp13.R   (R_OUT_DIR, WP13_ARCHIVE)
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "agree.R"))
ARC <- Sys.getenv("WP13_ARCHIVE", normalizePath(file.path(HERE, "..", "wp13_gapdh_corrected", "results"))); OUT <- Sys.getenv("R_OUT_DIR", file.path(HERE, "results_R", "wp13")); res <- list()
res[["gapdh_versions.json"]] <- agree_json(file.path(ARC, "gapdh_versions.json"), file.path(OUT, "gapdh_versions.json"), tol = 2e-6)$ok       # the QC plateau values come from the logistic fits (agree to 2.6e-8)
res[["gapdh_versions_summary.csv"]] <- agree_csv(file.path(ARC, "gapdh_versions_summary.csv"), file.path(OUT, "gapdh_versions_summary.csv"))$ok
cat(sprintf("\nWP13 so far (%d comparisons): %s\n", length(res), if (all(unlist(res))) "ALL AGREE" else "DIFFERENCES FOUND"))
