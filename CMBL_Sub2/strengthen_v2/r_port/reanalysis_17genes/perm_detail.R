## perm_detail.R - R port of reanalysis_17genes/perm_detail.py: mean and 2.5 / 97.5 percentiles of the permutation null (seed 2026, 2000 permutations) of the number of significant pairs per scheme and of
## the consensus counts, with the observed counts and P(null >= observed), for the 20 and the 17 genes. Writes perm_detail.json (indent 1).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "pipeline.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run <- function(genes, n = 2000, seed = 2026) {
  M <- P_load(genes); pn <- P_permutation_null(M, n, seed)
  rec <- list(Raw = pn$sig[["Raw"]], `Ins1-normalized` = pn$sig[["Ins1-normalized"]], `Kcnj11-normalized` = pn$sig[["Kcnj11-normalized"]], `Hes1-normalized` = pn$sig[["Hes1-normalized"]], cons = pn$ca, cons_sign = pn$cs)
  obs <- P_sig_sets(P_schemes(M)); out <- lapply(rec, function(a) list(mean = mean(a), lo = np_percentile(a, 2.5), hi = np_percentile(a, 97.5)))
  out$observed <- lapply(obs, function(s) as.integer(sum(s$p[upper.tri(s$p)] < 0.05)))
  out$observed$cons <- as.integer(nrow(P_consensus(obs, use_sign = FALSE))); out$observed$cons_sign <- as.integer(nrow(P_consensus(obs, use_sign = TRUE)))
  out$P_ge_obs <- list(cons = mean(rec$cons >= out$observed$cons), cons_sign = mean(rec$cons_sign >= out$observed$cons_sign))
  out
}
res <- list(`20` = run(P_GENES20), `17` = run(P_GENES17))
for (k in names(res)) {
  cat("==", k, "genes\n")
  for (kk in names(res[[k]])) cat("  ", kk, paste0("{", paste(sprintf("'%s': %s", names(res[[k]][[kk]]), vapply(res[[k]][[kk]], function(x) format(x, digits = 15), "")), collapse = ", "), "}"), "\n")
}
write_json_typed(res, file.path(OUT, "perm_detail.json"))
