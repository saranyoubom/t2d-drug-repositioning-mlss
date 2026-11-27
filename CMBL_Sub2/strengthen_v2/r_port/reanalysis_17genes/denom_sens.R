## denom_sens.R - R port of reanalysis_17genes/denom_sens.py: consensus pairs when Hes1 as normalization gene is replaced by another gene (Tcf7, Cacna1c, and Hey1 for the 20 genes):
## pairs evaluable (not involving the new denominator), retained, new consensus size and pairs gained. Writes denom_sens.json.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "pipeline.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run <- function(genes, alts) {
  M <- P_load(genes); base <- P_sig_sets(P_schemes(M)); cons <- P_consensus(base, use_sign = FALSE); consk <- .pkey(cons); out <- list(n_consensus = nrow(cons))
  for (D in alts) {
    sch <- list(Raw = M)
    for (dd in c("Ins1", "Kcnj11")) sch[[paste0(dd, "-normalized")]] <- M[, colnames(M) != dd, drop = FALSE] / M[, dd]
    sch[[paste0(D, "-normalized")]] <- M[, colnames(M) != D, drop = FALSE] / M[, D]
    newc <- .pkey(P_consensus(P_sig_sets(sch), use_sign = FALSE))
    evaluable <- consk[!(cons[, 1] == D | cons[, 2] == D)]; kept <- evaluable[evaluable %in% newc]
    out[[D]] <- list(evaluable = length(evaluable), retained = length(kept), new_consensus = length(newc), gained = length(setdiff(newc, consk)))
  }
  out
}
r17 <- run(P_GENES17, c("Tcf7", "Cacna1c")); r20 <- run(P_GENES20, c("Hey1", "Tcf7", "Cacna1c"))
cat("20 genes (check vs manuscript: Hey1/Tcf7 retained all evaluable, Cacna1c 12 of 22):", py_json_typed(r20, indent = NULL), "\n")
cat("17 genes:", py_json_typed(r17, indent = NULL), "\n")
write_json_typed(list(`20` = r20, `17` = r17), file.path(OUT, "denom_sens.json"))
