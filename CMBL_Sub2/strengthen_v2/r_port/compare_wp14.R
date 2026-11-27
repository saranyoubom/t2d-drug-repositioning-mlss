args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "agree.R"))
r <- agree_csv(file.path(here, "..", "wp14_table_s8", "results", "efficiency_sensitivity_analyzed_matrix.csv"), file.path(here, "results_R", "wp14", "efficiency_sensitivity_analyzed_matrix.csv"), tol = 2e-6)
