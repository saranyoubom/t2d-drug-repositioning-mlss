## bench2_core.R - R port of wp2_benchmark2/bench2_core.py: per-matrix analysis of protocol v2 section 3 for one 12 x 17 matrix (3 groups x 4 replicates, linear scale, the first three genes are the
## denominators), shared by the development series, the held-out arrays and the RNA-seq series (same code path). Every statistic comes from coexp_checks.R (WP1); this module only organizes the calls.
## Permutation counts: within-group coupling tests 2000 (3.1), consensus count 2000 (as in part 1), calibrated adjusted counts 500 (3.4). The random stream is seeded with 2026 + GSE number.
## Needs rng/numpy_rng.R and coexp_checks.R sourced first. Gene indices are 1-based in R (the denominators are genes 1 to 3).
BC_GROUPS <- rep(0:2, each = 4); BC_DES <- cc_design(BC_GROUPS); BC_DEN <- 1:3; BC_SCHEMES <- c("raw", "den1", "den2", "den3")

analyse_matrix <- function(X, genes, n, nperm_cp = 2000L, nperm_cons = 2000L, nperm_adj = 500L, eig = TRUE) {   # returns list(flat = scalars, detail = ...)
  X <- as.matrix(X); storage.mode(X) <- "double"; seed <- 2026 + as.integer(n); names <- as.character(genes); G <- ncol(X); up <- upper.tri(matrix(0, G, G))
  obs <- cc_observed_summary(X, BC_DES, BC_DEN, names); S0 <- obs$.S[[1]] & up
  cons <- obs$.cons_any; Q <- obs$.Q; n_cons <- as.integer(sum(cons))
  q_raw <- Q[[1]]; cons_raw <- as.integer(sum(cons & !is.na(q_raw) & q_raw < 0.05))                                   # E2 on the raw scheme alone
  Qs <- simplify2array(Q); present <- !is.na(Qs); qall <- apply(ifelse(present, Qs, -Inf), c(1, 2), max); qall[!is.finite(qall)] <- 1
  pairs_all4 <- as.integer(sum((qall < 0.05) & up))                                                                      # part 1 definition: significant in every scheme that contains the pair
  sh <- cc_group_mean_share(X, BC_DES, S0); geo <- cc_geometry_check(X, BC_DES)
  cp <- cc_coupling_test(X, BC_DES, nperm = nperm_cp, seed = seed, eig = eig)
  cn <- cc_consensus_null(X, BC_DES, BC_DEN, nperm = nperm_cons, seed = seed, genes = names)
  an <- cc_adjusted_null(X, BC_DES, BC_DEN, nperm = nperm_adj, seed = seed, genes = names)
  pooled <- obs$pooled_sig_per_scheme
  flat <- list(observed = n_cons, observed_same_sign = obs$consensus_same_sign, null_mean = cn$null_mean, null_sd = cn$null_sd, null_q95 = cn$null_q95, z = cn$z, perm_p = cn$p_ge, perm_p_plus1 = cn$p_ge_plus1,
               perm_p_same_sign = cn$same_sign$p_ge, pooled_raw = pooled[1], pooled_den1 = pooled[2], pooled_den2 = pooled[3], pooled_den3 = pooled[4],
               adj_any = obs$cons_adj_any, adj_all4 = obs$cons_adj_all_schemes, adj_pairs_all4 = pairs_all4, adj_cons_raw = cons_raw,
               retained_any = if (n_cons) obs$cons_adj_any / n_cons else NaN, retained_all4 = if (n_cons) obs$cons_adj_all_schemes / n_cons else NaN, retained_raw = if (n_cons) cons_raw / n_cons else NaN,
               share_between = sh$median, share_sig_all = geo$share_sig_all, share_sig_group_effect = geo$share_sig_group_effect, n_group_effect_genes = geo$n_group_effect_genes, geo_expected = geo$expected)
  sdw <- sqrt(colSums((BC_DES$Q %*% log(X))^2) / (BC_DES$N - BC_DES$K))                                              # pooled within-group SD (natural log) per gene
  flat$n_zero_var_genes <- as.integer(sum(sdw < 1e-6))                                                                  # degenerate genes (ADD-ON, 5 Nov 01:35): no within-group variation
  for (k in c("T_r2", "T_p", "T_eig")) if (!is.null(cp[[k]])) { flat[[paste0(k, "_obs")]] <- cp[[k]]$observed; flat[[paste0(k, "_null_mean")]] <- cp[[k]]$null_mean; flat[[paste0(k, "_p")]] <- cp[[k]]$p_ge_plus1 }
  for (i in seq_along(BC_SCHEMES)) {
    s <- BC_SCHEMES[i]; v <- an$schemes[[an$scheme_names[i]]]
    flat[[paste0("adjq_", s, "_obs")]] <- v$observed; flat[[paste0("adjq_", s, "_null_mean")]] <- v$null_mean; flat[[paste0("adjq_", s, "_null_sd")]] <- v$null_sd; flat[[paste0("adjq_", s, "_null_q95")]] <- v$null_q95
    flat[[paste0("adjq_", s, "_p")]] <- v$p_ge_plus1; flat[[paste0("adjq_", s, "_p_plain")]] <- v$p_ge
  }
  deg <- as.list(obs$degree); if (!length(deg)) deg <- setNames(list(), character(0))
  detail <- list(genes = as.list(names), consensus_pairs = lapply(obs$consensus_pairs, as.list), degree = deg, adj_q05_per_scheme = as.list(obs$adj_q_lt_alpha_per_scheme))
  list(flat = flat, detail = detail)
}
