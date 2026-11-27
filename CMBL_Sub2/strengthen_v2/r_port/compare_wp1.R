## compare_wp1.R - agreement of the R port of WP1 with the saved Python results (wp1_ins1_structure/results_ins1.json and tbl_2_*.csv)
args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "agree.R"))
PY <- Sys.getenv("WP1_DIR", file.path(here, "..", "wp1_ins1_structure")); RR <- Sys.getenv("R_OUT_DIR", file.path(here, "results_R", "wp1"))
files <- c("tbl_2_2_run_structure.csv", "tbl_2_3_hubs_vs_run.csv", "tbl_2_3_degree_by_run.csv", "tbl_2_3_pairs_involving_run.csv", "tbl_2_3_eligible_pairs.csv", "tbl_2_4_pairs.csv", "tbl_2_4_shared_reference_summary.csv",
           "tbl_2_4_sensitivity_versions.csv", "tbl_2_5_adjusted_calibrated.csv", "tbl_2_5_q_per_pair.csv", "tbl_2_6_coupling.csv", "tbl_2_7_geometry.csv", "tbl_2_8_scale_sensitivity.csv")
res <- list()
for (f in files) res[[f]] <- agree_csv(file.path(PY, f), file.path(RR, f), key = if (f == "tbl_2_5_q_per_pair.csv") c("gene_a", "gene_b") else NULL)
## known benign difference: raw q values tie at the 5th place of 'smallest_raw_q_pairs'; numpy's unstable argsort picked another pair of the tie than R's stable order
res[["results_ins1.json"]] <- agree_json(file.path(PY, "results_ins1.json"), file.path(RR, "results_ins1.json"), skip = c("/runtime_seconds", "/2_5_adjusted_calibrated/smallest_raw_q_pairs[5]"))
ok <- vapply(res, function(r) r$ok, TRUE)
cat(sprintf("\nWP1 agreement: %d of %d files agree\n", sum(ok), length(ok)))
