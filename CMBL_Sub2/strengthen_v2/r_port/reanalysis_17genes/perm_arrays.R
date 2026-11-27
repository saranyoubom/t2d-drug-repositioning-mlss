## perm_arrays.R - R port of reanalysis_17genes/perm_arrays.py: for the 17 genes, per permutation (seed 2026, 2000 residual permutations) the number of significant pairs (p < 0.05) in each scheme and
## the number of consensus pairs (any sign, sign-consistent). The Python script saved perm_arrays_17.npz; this one saves perm_arrays_17.csv (one column per array, same names and order).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "pipeline.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
M <- P_load(P_GENES17); pn <- P_permutation_null(M, 2000, 2026)
rec <- list(Raw = pn$sig[["Raw"]], `Hes1-normalized` = pn$sig[["Hes1-normalized"]], `Kcnj11-normalized` = pn$sig[["Kcnj11-normalized"]], `Ins1-normalized` = pn$sig[["Ins1-normalized"]], cons = pn$ca, cons_sign = pn$cs)
py_write_csv(as.data.frame(rec, check.names = FALSE), file.path(OUT, "perm_arrays_17.csv"))
cat("saved", paste0("{", paste(sprintf("'%s': %s", names(rec), vapply(rec, function(v) format(mean(v), digits = 15), "")), collapse = ", "), "}"), "\n")
