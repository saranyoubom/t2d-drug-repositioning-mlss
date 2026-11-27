## test_coexp_checks.R - R port of test_coexp_checks.py (gate 1.1). Same checks, same order, same names and details as the Python output
## (wp1_ins1_structure/test_output.txt), so the two outputs can be compared line by line. The two command-line checks of the Python test have no R
## counterpart (the R port has no command line). Run: Rscript test_coexp_checks.R  (env CC_SECTIONS=1,2,3 runs only some sections; INS1_MATRIX overrides the matrix path)
options(warn = 1)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "rng", "numpy_rng.R")); source(file.path(here, "coexp_checks.R")); source(file.path(here, "bench_pipeline.R"))
if (Sys.getenv("CC_FAST", "0") == "1") { source(file.path(here, "cc_fast.R")); cat("batched permutation statistics (cc_fast.R) in use\n") }      # CC_FAST=1: run the 52 checks with the batched engine
suppressMessages({ library(jsonlite); library(readxl) })
SECT <- as.integer(strsplit(Sys.getenv("CC_SECTIONS", "1,2,3,4,5,6"), ",")[[1]])
root <- file.path(here, "..", "..")
SRC <- Sys.getenv("INS1_MATRIX", file.path(root, "R_17genes", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"))
REF <- fromJSON(file.path(root, "reanalysis_17genes", "pipeline_17genes.json"), simplifyVector = FALSE)
G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
FAILED <- character(0); NCHECK <- 0
section <- function(t) cat("\n== ", t, "\n", sep = "")
check <- function(name, cond, detail = "") {
  NCHECK <<- NCHECK + 1; cond <- isTRUE(cond)
  cat(if (cond) "PASS  " else "FAIL  ", name, if (nzchar(detail)) paste0(" | ", detail) else "", "\n", sep = "")
  if (!cond) FAILED <<- c(FAILED, name)
}
fi <- function(v) paste0("[", paste(as.integer(v), collapse = ", "), "]")
wilson <- function(k, n, z = 1.96) { p <- k / n; d <- 1 + z * z / n; cc <- p + z * z / (2 * n); h <- z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)); c((cc - h) / d, (cc + h) / d) }
cat("coexp_checks", CC_VERSION, "| R", R.version.string, "\n")
pairset <- function(pl) sort(vapply(pl, function(p) paste(sort(unlist(p)), collapse = "|"), ""))
arr3 <- function(z, d) aperm(array(z, dim = rev(d)), 3:1)         # numpy C-order (d1, d2, d3) from a flat vector of draws

X <- as.matrix(as.data.frame(read_excel(SRC))[, G17]); storage.mode(X) <- "double"
grp <- rep(c("CTRL", "DAPT", "DKK-1"), each = 4); den <- c("Ins1", "Kcnj11", "Hes1"); denx <- match(den, G17); des <- cc_design(grp)

## ------------------------------------------------------------------------------------------------------------------ 1. INS-1 regression
if (1 %in% SECT) {
  section("1. INS-1 regression: 17 genes, 3 groups x 4 wells, denominators Ins1, Kcnj11, Hes1")
  o <- cc_observed_summary(X, grp, den, G17)
  check("consensus pairs, any sign = 21", o$consensus_any == 21 && 21 == length(REF$consensus), sprintf("got %d", o$consensus_any))
  check("consensus pairs, same sign in all four schemes = 21", o$consensus_same_sign == 21 && 21 == length(REF$consensus_signcons), sprintf("got %d", o$consensus_same_sign))
  check("same 21 consensus pairs as pipeline_17genes.json", identical(pairset(o$consensus_pairs), pairset(REF$consensus)))
  refsig <- as.integer(vapply(REF$counts, function(c) c$sig, 1))
  check("pooled significant pairs raw/Ins1/Kcnj11/Hes1 = 66/80/84/120", identical(as.integer(o$pooled_sig_per_scheme), c(66L, 80L, 84L, 120L)) && identical(as.integer(o$pooled_sig_per_scheme), refsig), fi(o$pooled_sig_per_scheme))
  refq <- as.integer(vapply(REF$adj_sig_q05, function(v) as.integer(v), 1L))
  check("adjusted q < 0.05 raw/Ins1/Kcnj11/Hes1 = 0/51/11/15", identical(as.integer(o$adj_q_lt_alpha_per_scheme), c(0L, 51L, 11L, 15L)) && identical(as.integer(o$adj_q_lt_alpha_per_scheme), refq), fi(o$adj_q_lt_alpha_per_scheme))
  check("13 consensus pairs with q < 0.05 in >= 1 scheme", o$cons_adj_any == 13 && 13 == length(REF$cons_adj_any), sprintf("got %d", o$cons_adj_any))
  check("0 consensus pairs with q < 0.05 in all four schemes; 0 pairs significant after adjustment in all four", o$cons_adj_all_schemes == 0 && o$adj_pairs_all_schemes == 0 && 0 == length(REF$adj_all_four))
  check("91 eligible consensus pairs (14 non-denominator genes)", o$n_eligible_pairs == 91)
  check("hub degrees Glp1r 6; Tcf7l2, Cacna1d, Ptbp1 5; Wnt5a, Cacna1c 4", identical(as.integer(unname(o$degree[c("Glp1r", "Tcf7l2", "Cacna1d", "Ptbp1", "Wnt5a", "Cacna1c")])), c(6L, 5L, 5L, 5L, 4L, 4L)))
  bp <- B_pooled(X, denx); ba <- B_adjusted(X, denx)
  check("pooled significance matrices identical to bench_pipeline.pooled", all(mapply(function(a, b) identical(a, b), bp$S, o$.S)))
  nz <- function(m) ifelse(is.na(m), -1, m)
  check("adjusted q matrices equal bench_pipeline.adjusted (statsmodels BH) to 1e-12", all(mapply(function(a, b) all(abs(nz(a) - nz(b)) <= 1e-12), ba$Q, o$.Q)))
  sh <- cc_group_mean_share(X, grp, o$.S[[1]] & upper.tri(o$.S[[1]])); sh_old <- B_between_share(X, bp$S[[1]] & upper.tri(bp$S[[1]]))
  check("group-mean covariance share, median over 66 raw significant pairs = 0.97", round(sh$median, 2) == 0.97 && abs(sh$median - sh_old) < 1e-12, sprintf("median %.4f mean %.4f pooled %.4f", sh$median, sh$mean, sh$pooled))
  cn <- cc_consensus_null(X, grp, den, nperm = 10000, seed = 2026, genes = G17)
  check("consensus null mean 22.1 (10000 permutations; manuscript 22.1)", abs(cn$null_mean - REF$perm$mean) <= 0.4, sprintf("mean %.2f sd %.2f (manuscript %.2f, %.2f)", cn$null_mean, cn$null_sd, REF$perm$mean, REF$perm$sd))
  check("consensus null SD within 0.3 of the manuscript", abs(cn$null_sd - REF$perm$sd) <= 0.3)
  check("P(null >= 21) about 0.69", abs(cn$p_ge - REF$perm$p_ge_obs) <= 0.04, sprintf("got %.3f (manuscript %.3f)", cn$p_ge, REF$perm$p_ge_obs))
  check("same-sign null mean within 0.4 of the manuscript", abs(cn$same_sign$null_mean - REF$perm_signcons$mean) <= 0.4, sprintf("got %.2f (manuscript %.2f)", cn$same_sign$null_mean, REF$perm_signcons$mean))
  t0 <- Sys.time(); old <- B_permutation_counts(X, denx, 300, np_rng(11)); t_old <- as.numeric(Sys.time() - t0, units = "secs"); cn2 <- cc_consensus_null(X, grp, den, nperm = 3000, seed = 7, genes = G17)
  se <- sqrt(var(old$ca) / 300 + var(cn2$.null) / 3000)
  check("vectorized null agrees with the loop null (300 loop permutations): means within 3 SE, SDs within 25%", abs(mean(old$ca) - cn2$null_mean) <= 3 * se && abs(.pop_sd(old$ca) / cn2$null_sd - 1) < 0.25,
        sprintf("loop %.2f (sd %.2f), vectorized %.2f (sd %.2f)", mean(old$ca), .pop_sd(old$ca), cn2$null_mean, cn2$null_sd))
  an <- cc_adjusted_null(X, grp, den, nperm = 2000, seed = 2026, genes = G17); a <- an$schemes; nm <- an$scheme_names
  check("calibrated adjusted counts: observed 0/51/11/15 from the threshold counter equals the exact BH counts", identical(as.integer(vapply(nm, function(n) a[[n]]$observed, 1)), c(0L, 51L, 11L, 15L)))
  fl <- function(v, d) paste0("[", paste(sprintf(paste0("%.", d, "f"), v), collapse = ", "), "]")
  cat(sprintf("      adjusted q<0.05 [raw, Ins1, Kcnj11, Hes1]: observed %s | null mean %s | P(null >= obs) %s | excess %s\n", fi(vapply(nm, function(n) a[[n]]$observed, 1)),
              fl(vapply(nm, function(n) a[[n]]$null_mean, 1), 1), fl(vapply(nm, function(n) a[[n]]$p_ge, 1), 3), fl(vapply(nm, function(n) a[[n]]$excess, 1), 1)))
  check("adjusted nulls close to the planning reconnaissance (1000 permutations: 0.1, 29.0, 11.6, 6.4; P 1.0, 0.087, 0.465, 0.074)", abs(a$raw$null_mean - 0.1) < 0.3 && abs(a[["Ins1-normalized"]]$null_mean - 29.0) < 1.5 &&
          abs(a[["Kcnj11-normalized"]]$null_mean - 11.6) < 1.0 && abs(a[["Hes1-normalized"]]$null_mean - 6.4) < 1.0 && abs(a[["Hes1-normalized"]]$p_ge - 0.074) < 0.04 && a$raw$p_ge == 1.0)
  ct <- cc_coupling_test(X, grp, nperm = 10000, seed = 2026)
  cat("      INS-1 coupling (10000 permutations): ", paste(vapply(names(ct), function(k) sprintf("%s %.4g (null %.4g, P %.4f)", k, ct[[k]]$observed, ct[[k]]$null_mean, ct[[k]]$p_ge_plus1), ""), collapse = " | "), "\n", sep = "")
  rc <- cc_run_checks(X, grp, G17, den, nperm = 2000); w <- unlist(rc$warnings)
  check("run_checks warns: consensus count within null, ratio schemes not calibrated, group means carry > 90% of covariance", any(grepl("consensus count", w)) && any(grepl("Ins1-normalized", w)) && any(grepl("group means carry", w)), sprintf("%d warnings", length(w)))
}

## ------------------------------------------------------------------------------------------------------------------ 2. primitives
if (2 %in% SECT) {
  section("2. vectorized primitives against reference implementations")
  rng <- np_rng(1); good <- TRUE
  for (i in 1:200) { p <- round(np_random(rng, 60)^2, 3); good <- good && isTRUE(all.equal(.bh(p), p.adjust(p, "BH"), tolerance = 1e-12)) }
  check("BH adjusted p equals statsmodels fdr_bh (ties and zeros included)", good)
  tot <- 0; good <- TRUE
  for (i in 1:150) {
    sm <- cc_simulate(c(4, 4, 4), 17, rng); Xs <- cc_inject_factor(sm$X, sm$codes, 1:10, 0.8, rng); r <- cc_corr(t(Xs) %*% des$Q); iu <- upper.tri(r)
    exact <- sum(.bh(.pval(r[iu], des$df_resid)) < 0.05); fast <- .bh_count(abs(r[iu]), des$df_resid); good <- good && exact == fast; tot <- tot + exact
  }
  check("threshold BH counter equals exact BH count on 150 data sets with many significant pairs", good, sprintf("total significant pairs %d", tot))
  Rb <- arr3(np_standard_normal(rng, 6 * 17 * 12), c(6, 17, 12)); for (b in 1:6) Rb[b, , ] <- Rb[b, , ] - rowMeans(Rb[b, , ]); st <- .coupling_stats(Rb, des); ok <- TRUE
  for (b in 1:6) { r <- cor(t(Rb[b, , ])); iu <- upper.tri(r)
    ok <- ok && abs(st$T_r2[b] - mean(r[iu]^2)) < 1e-12 && st$T_p[b] == sum(abs(r[iu]) >= des$rc_resid) && abs(st$T_eig[b] - eigen(r, symmetric = TRUE, only.values = TRUE)$values[1]) < 1e-9 }
  Rb2 <- arr3(np_standard_normal(rng, 4 * 6 * 12), c(4, 6, 12)); st2 <- .coupling_stats(Rb2, des)
  ok <- ok && all(vapply(1:4, function(b) abs(st2$T_eig[b] - eigen(cor(t(Rb2[b, , ])), symmetric = TRUE, only.values = TRUE)$values[1]) < 1e-9, TRUE))
  check("T_r2, T_p, T_eig equal explicit np.corrcoef / eigvalsh computations (G > N via the Gram matrix, G < N directly)", ok)
  rho <- cor(X, method = "spearman"); osp <- cc_observed_summary(X, grp, den, G17, kind = "spearman"); n_sp <- sum(abs(rho[upper.tri(rho)]) >= des$rc_pooled)
  for (d in denx) { keep <- setdiff(1:17, d); rr <- cor(X[, keep] / X[, d], method = "spearman"); n_sp <- c(n_sp, sum(abs(rr[upper.tri(rr)]) >= des$rc_pooled)) }
  check("Spearman pooled counts per scheme equal scipy.stats.spearmanr", identical(as.integer(osp$pooled_sig_per_scheme), as.integer(n_sp)), fi(n_sp))
  olg <- cc_observed_summary(X, grp, den, G17, kind = "log"); rl <- cor(log(X)); n_lg <- sum(abs(rl[upper.tri(rl)]) >= des$rc_pooled)
  for (d in denx) { keep <- setdiff(1:17, d); rr <- cor(log(X[, keep] / X[, d])); n_lg <- c(n_lg, sum(abs(rr[upper.tri(rr)]) >= des$rc_pooled)) }
  check("log-scale Pearson pooled counts per scheme equal np.corrcoef of log values", identical(as.integer(olg$pooled_sig_per_scheme), as.integer(n_lg)), fi(n_lg))
  check("critical values: r_crit(10) = 0.576, r_crit(8) = 0.632 for 3 x 4", abs(des$rc_pooled - 0.5760) < 5e-4 && abs(des$rc_resid - 0.6319) < 5e-4 && des$df_resid == 8)
  d2 <- cc_design(c(rep("a", 3), rep("b", 5))); d4 <- cc_design(c(rep(0, 4), rep(1, 5), rep(2, 6), rep(3, 3)))
  check("Design handles K = 2 unequal sizes (df_resid 5) and K = 4 (df_resid 13)", d2$df_resid == 5 && d4$df_resid == 13 && identical(dim(d2$Q), c(8L, 8L)) && all(abs(colSums(d4$Q)) < 1e-12))
  sa <- np_rng(3)
  a2 <- matrix(np_standard_normal(sa, 2e6), 1e6, 2, byrow = TRUE); b2 <- matrix(np_standard_normal(sa, 2e6), 1e6, 2, byrow = TRUE)
  c2 <- abs(rowSums(a2 * b2) / sqrt(rowSums(a2^2)) / sqrt(rowSums(b2^2)))
  a3 <- matrix(np_standard_normal(sa, 3e6), 1e6, 3, byrow = TRUE); b3 <- matrix(np_standard_normal(sa, 3e6), 1e6, 3, byrow = TRUE)
  c3 <- abs(rowSums(a3 * b3) / sqrt(rowSums(a3^2)) / sqrt(rowSums(b3^2))); e3 <- cc_geometric_expectation(0.576, 3); e4 <- cc_geometric_expectation(0.576, 4)
  check("geometric expectation K = 3: (2/pi) arccos(0.576) = 0.609 and Monte Carlo", abs(e3 - 2 / pi * acos(0.576)) < 1e-12 && abs(e3 - mean(c2 >= 0.576)) < 0.003, sprintf("formula %.4f MC %.4f", e3, mean(c2 >= 0.576)))
  check("geometric expectation K = 4: 1 - r_crit = 0.424 and Monte Carlo; K = 2: 1", abs(e4 - 0.424) < 1e-12 && abs(e4 - mean(c3 >= 0.576)) < 0.003 && cc_geometric_expectation(0.576, 2) == 1.0, sprintf("MC %.4f", mean(c3 >= 0.576)))
}

## ------------------------------------------------------------------------------------------------------------------ 3. permutation engine
if (3 %in% SECT) {
  section("3. permutation engine")
  M <- matrix(rep(0:11, 6), 12, 6); gr <- c(rep("a", 4), rep("b", 4), rep("c", 4)); bl <- c("x", "x", "y", "y", "y", NA)
  bat <- cc_permuted_all(M, gr, 500, seed = 5, blocks = bl)               # entries are the source row index of each permuted sample
  within <- TRUE
  for (g in 1:6) for (k in 0:2) { sl <- bat[, g, 4 * k + 1:4]; within <- within && all(sl %in% (4 * k + 0:3)) && all(t(apply(sl, 1, sort)) == matrix(4 * k + 0:3, nrow(sl), 4, byrow = TRUE)) }
  check("every permutation stays inside its group (each gene a full within-group permutation)", within)
  eqg <- function(i, j) all(bat[, i, ] == bat[, j, ])
  check("genes of one block share the permutation, different blocks differ", eqg(1, 2) && eqg(3, 4) && eqg(4, 5) && !eqg(1, 3) && !eqg(5, 6))
  a1 <- cc_permuted_all(M, gr, 700, seed = 9, blocks = bl); a2_ <- cc_permuted_all(M, gr, 700, seed = 9, blocks = bl, max_elems = 150)
  check("result does not depend on the memory chunking (max_elems)", identical(a1, a2_))
  one <- cc_permuted_all(matrix(0:7, 8, 1), c(rep(0, 4), rep(1, 4)), 24000, seed = 1); codes <- apply(one[, 1, 1:4], 1, paste, collapse = ",")
  cnts <- as.vector(table(codes)); chi <- if (length(cnts) == 24) chisq.test(cnts)$p.value else NA
  check("group of 4: all 24 permutations occur equally often (chi-square P > 0.001)", !is.na(chi) && chi > 0.001, sprintf("P %.3f", chi))
  gl <- c(rep(0, 10), rep(1, 10)); M10 <- matrix(rep(0:19, 2), 20, 2); b10 <- cc_permuted_all(M10, gl, 20000, seed = 4)
  pos <- vapply(0:9, function(k) vapply(0:9, function(s) mean(b10[, 1, k + 1] == s), 1), numeric(10))
  check("group of 10 (random-key path): each source row lands in each position with probability 0.1", max(abs(pos - 0.1)) < 0.01 && all(t(apply(b10[, 1, 1:10], 1, sort)) == matrix(0:9, 20000, 10, byrow = TRUE)), sprintf("max deviation %.4f", max(abs(pos - 0.1))))
  Xr <- exp(matrix(np_standard_normal(np_rng(8), 72), 12, 6, byrow = TRUE)); one_block <- rep("r", 6)     # lognormal(0, 1, (12, 6)) = exp(normal)
  cn1 <- cc_consensus_null(Xr, gr, NULL, nperm = 40, seed = 3, blocks = one_block); ct1 <- cc_coupling_test(Xr, gr, nperm = 40, seed = 3, blocks = one_block)
  check("one block holding all genes leaves every correlation unchanged (null equals observed, P = 1)", cn1$null_sd == 0 && cn1$null_mean == cn1$observed && ct1$T_r2$null_sd < 1e-12 && ct1$T_r2$p_ge == 1.0)
  pairs_ok <- TRUE; Mbs <- cc_permuted_all(Xr, gr, 50, seed = 3, blocks = c("p", "p", "p", "q", "q", "q"))
  for (b in 1:50) pairs_ok <- pairs_ok && abs(cor(Mbs[b, 1, ], Mbs[b, 2, ]) - cor(Xr[, 1], Xr[, 2])) < 1e-12 && abs(cor(Mbs[b, 4, ], Mbs[b, 6, ]) - cor(Xr[, 4], Xr[, 6])) < 1e-12
  check("run-block permutation keeps same-block gene pairs exactly and breaks between-block pairs", pairs_ok && abs(cor(Mbs[1, 1, ], Mbs[1, 4, ]) - cor(Xr[, 1], Xr[, 4])) > 1e-6)
}

## ------------------------------------------------------------------------------------------------------------------ 4. synthetic known-answer tests
if (4 %in% SECT) {
  section("4. synthetic known-answer tests (3 groups x 4, 17 genes, group effects SD 1.0, noise SD 0.3 on the log scale)")
  rs <- np_rng(2026); ND <- 250; P <- list(T_r2 = numeric(0), T_p = numeric(0), T_eig = numeric(0), count = numeric(0)); t0 <- Sys.time()
  for (i in 0:(ND - 1)) {
    sm <- cc_simulate(c(4, 4, 4), 17, rs); Xs <- sm$X; g <- sm$codes; cc_ <- cc_coupling_test(Xs, g, nperm = 199, seed = 1000 + i); k <- cc_consensus_null(Xs, g, 15:17, nperm = 199, seed = 1000 + i)
    for (s in c("T_r2", "T_p", "T_eig")) P[[s]] <- c(P[[s]], cc_[[s]]$p_ge_plus1)
    P$count <- c(P$count, k$p_ge_plus1)
  }
  cat(sprintf("      %d independent data sets, 199 permutations each (%.0fs)\n", ND, as.numeric(Sys.time() - t0, units = "secs")))
  for (s in c("T_r2", "T_p", "T_eig", "count")) {
    lim <- 0.09; p <- P[[s]]; ks <- suppressWarnings(ks.test(p, "punif", exact = TRUE)); err <- mean(p <= 0.05); w <- wilson(sum(p <= 0.05), ND)
    check(sprintf("independent data, %s: P approximately uniform (mean 0.5 +- 0.07, type I error at 0.05 within 0.015-%.2f)", s, lim),
          abs(mean(p) - 0.5) < (if (s %in% c("T_r2", "T_eig")) 0.07 else 0.12) && err >= 0.015 && err <= lim,
          sprintf("mean P %.3f, KS D %.3f, P(P <= 0.05) = %.3f (Wilson %.3f-%.3f)", mean(p), unname(ks$statistic), err, w[1], w[2]))
  }
  ksr <- suppressWarnings(ks.test(P$T_r2, "punif", exact = TRUE)); check("independent data, T_r2: Kolmogorov-Smirnov against uniform P > 0.01", ksr$p.value > 0.01, sprintf("KS P %.3f", ksr$p.value))
  NR <- 100; power <- list()
  for (cfg in list(c(0.7, 6), c(0.9, 10), c(0.5, 3))) {
    r <- cfg[1]; m <- cfg[2]; hit <- c(T_r2 = 0, T_p = 0, T_eig = 0, count = 0)
    for (i in 0:(NR - 1)) {
      sm <- cc_simulate(c(4, 4, 4), 17, rs); Xs <- cc_inject_factor(cc_decouple(sm$X, sm$codes, rs), sm$codes, seq_len(m), r, rs)
      cc_ <- cc_coupling_test(Xs, sm$codes, nperm = 199, seed = 5000 + i); k <- cc_consensus_null(Xs, sm$codes, 15:17, nperm = 199, seed = 5000 + i)
      for (s in c("T_r2", "T_p", "T_eig")) hit[s] <- hit[s] + (cc_[[s]]$p_ge_plus1 < 0.05)
      hit["count"] <- hit["count"] + (k$p_ge_plus1 < 0.05)
    }
    power[[paste(r, m)]] <- hit
    cat(sprintf("      injected factor r = %.1f in m = %d genes (decoupled data, %d replicates): detection T_r2 %.2f, T_p %.2f, T_eig %.2f, consensus count %.2f\n", r, m, NR, hit["T_r2"] / NR, hit["T_p"] / NR, hit["T_eig"] / NR, hit["count"] / NR))
  }
  check("injected coupling r = 0.7 in 6 genes is detected by T_r2 in >= 80% of data sets (protocol 3.2 requirement, here on synthetic data)", power[["0.7 6"]]["T_r2"] / NR >= 0.80, sprintf("%.2f", power[["0.7 6"]]["T_r2"] / NR))
  check("injected coupling r = 0.9 in 10 genes is detected by T_r2 in >= 95% of data sets", power[["0.9 10"]]["T_r2"] / NR >= 0.95)
  check("the consensus count is less sensitive than T_r2 to the same injected coupling (r = 0.7, m = 6), as in benchmark part 1", power[["0.7 6"]]["count"] < power[["0.7 6"]]["T_r2"], sprintf("count %.2f vs T_r2 %.2f", power[["0.7 6"]]["count"] / NR, power[["0.7 6"]]["T_r2"] / NR))
  sm <- cc_simulate(c(30, 30, 30), 10, rs, noise_sd = 0.3); Xl <- sm$X; gl_ <- sm$codes
  Xi <- cc_inject_factor(cc_decouple(Xl, gl_, rs), gl_, 1:6, 0.7, rs); rl <- cc_residual_corr(Xi, gl_); iu6 <- upper.tri(matrix(0, 6, 6)); rn <- cc_residual_corr(cc_decouple(Xl, gl_, rs), gl_)
  check("inject_factor: mean residual correlation among the 6 injected genes = 0.7 +- 0.08 (3 x 30 samples); decoupled genes about 0", abs(mean(rl[1:6, 1:6][iu6]) - 0.7) < 0.08 && abs(mean(rn[upper.tri(rn)])) < 0.03,
        sprintf("injected %.3f, decoupled %.3f", mean(rl[1:6, 1:6][iu6]), mean(rn[upper.tri(rn)])))
  dd <- cc_decouple(Xl, gl_, rs); check("decouple keeps every group mean of the log values exactly", all(vapply(cc_design(gl_)$rows, function(r) isTRUE(all.equal(colMeans(log(Xl)[r, , drop = FALSE]), colMeans(log(dd)[r, , drop = FALSE]))), TRUE)))
  NS <- 60; raw_c <- numeric(0); ratio_c <- numeric(0); nullm <- numeric(0); pr_ <- NULL
  for (i in 0:(NS - 1)) {
    sm <- cc_simulate(c(4, 4, 4), 17, rs, noise_sd = 0.3, ref_noise_sd = 0.6, ref_gene = 15:17); a <- cc_adjusted_null(sm$X, sm$codes, 15:17, nperm = 199, seed = 7000 + i)
    s <- a$schemes; rn_ <- a$scheme_names[-1]; raw_c <- c(raw_c, s$raw$observed); ratio_c <- c(ratio_c, mean(vapply(rn_, function(n) s[[n]]$observed, 1)))
    nullm <- c(nullm, mean(vapply(rn_, function(n) s[[n]]$null_mean, 1))); pr_ <- rbind(pr_, vapply(rn_, function(n) s[[n]]$p_ge_plus1, 1))
  }
  cat(sprintf("      noisy shared denominators (extra log SD 0.6): adjusted q<0.05 count mean raw %.2f, ratio schemes %.1f (of 120 pairs); null mean of the ratio schemes %.1f; share of ratio-scheme P <= 0.05: %.3f\n", mean(raw_c), mean(ratio_c), mean(nullm), mean(pr_ <= 0.05)))
  check("independent genes with noisy shared denominators: ratio-scheme adjusted counts are inflated (>= 10x the raw count, > 20 pairs)", mean(ratio_c) > 20 && mean(ratio_c) >= 10 * max(mean(raw_c), 0.1))
  check("the calibrated null reproduces that inflation: P(null >= obs) of the ratio schemes is not small (share <= 0.05 is at most 0.10)", mean(pr_ <= 0.05) <= 0.10 && abs(mean(ratio_c) - mean(nullm)) < 0.15 * mean(ratio_c))
  Dp <- numeric(0); Ds <- numeric(0)
  for (i in 0:59) {
    sm <- cc_simulate(c(4, 4, 4), 12, rs); Xs <- sm$X; g <- sm$codes; runs <- c(rep("A", 5), rep("B", 4), rep("C", 3)); Dp <- c(Dp, cc_run_structure_test(Xs, g, runs, nperm = 999, seed = i)$r$p_two_sided)
    Xc <- cc_decouple(Xs, g, rs)
    for (idx in list(1:5, 6:9)) Xc <- cc_inject_factor(Xc, g, idx, 0.6, rs)
    Ds <- c(Ds, cc_run_structure_test(Xc, g, runs, nperm = 999, seed = i)$r$p_two_sided)
  }
  check("run-structure test: P approximately uniform on independent genes (mean 0.5 +- 0.1, P <= 0.05 in at most 12%)", abs(mean(Dp) - 0.5) < 0.1 && mean(Dp <= 0.05) <= 0.12, sprintf("mean %.3f, share %.3f", mean(Dp), mean(Dp <= 0.05)))
  check("run-structure test: detects within-run coupling (P <= 0.05 in >= 70% of data sets with r = 0.6 inside runs A and B)", mean(Ds <= 0.05) >= 0.70, sprintf("share %.2f", mean(Ds <= 0.05)))
}

## ------------------------------------------------------------------------------------------------------------------ 5. other designs
if (5 %in% SECT) {
  section("5. other designs (the command-line checks of the Python test have no R counterpart)")
  rs <- if (4 %in% SECT) rs else np_rng(99)
  so <- cc_simulate(c(3, 5), 8, rs); ro <- cc_run_checks(so$X, c(rep("a", 3), rep("b", 5)), nperm = 99, run_nperm = 99)
  sk <- cc_simulate(c(4, 5, 6, 3), 9, rs); rk <- cc_run_checks(sk$X, sk$codes, den = c(1, 2), runs = strsplit("aabbbccdd", "")[[1]], nperm = 99, run_nperm = 99)
  check("run_checks runs for K = 2 (unequal sizes) and K = 4 with run labels; geometry expectation 1.0 and 0.42 at the pooled r_crit", ro$meta$df_resid == 5 && rk$meta$df_resid == 13 && ro$geometry$expected == 1.0 &&
          abs(rk$geometry$expected - (1 - rk$geometry$r_crit)) < 1e-9 && !is.null(rk$run_structure))
  msg <- tryCatch({ cc_run_checks(matrix(c(1, 1, 3, 5, 0, 2, 4, 6), 4, 2), c(0, 0, 1, 1)); "" }, error = function(e) conditionMessage(e))
  check("non-positive input is refused with a clear message", grepl("positive", msg), substr(msg, 1, 60))
}

## ------------------------------------------------------------------------------------------------------------------ 6. speed
if (6 %in% SECT) {
  section("6. speed (this machine, G = 17, N = 12, 3 x 4, single process; informational in R)")
  sm <- cc_simulate(c(4, 4, 4), 17, np_rng(5)); Xs <- sm$X; g <- sm$codes; cc_coupling_test(Xs, g, nperm = 300)
  rate <- function(f, n) { t <- system.time(f())[["elapsed"]]; n / max(t, 1e-9) }
  sp <- c(T_r2_T_p = rate(function() cc_coupling_test(Xs, g, nperm = 4000, eig = FALSE), 4000), T_r2_T_p_T_eig = rate(function() cc_coupling_test(Xs, g, nperm = 4000), 4000),
          consensus = rate(function() cc_consensus_null(Xs, g, 15:17, nperm = 4000), 4000), adjusted_4_schemes = rate(function() cc_adjusted_null(Xs, g, 15:17, nperm = 4000), 4000),
          permutation_generation = rate(function() cc_perm_apply(Xs, g, 20000, 2026, fun = function(s) 0), 20000),
          spearman_consensus = rate(function() cc_consensus_null(Xs, g, 15:17, nperm = 1000, kind = "spearman"), 1000),
          run_label_permutations = rate(function() cc_run_structure_test(Xs, g, strsplit("aaaabbbbbcccccdde", "")[[1]], nperm = 20000), 20000))
  for (k in names(sp)) cat(sprintf("      %-26s %8.0f permutations/s\n", k, sp[[k]]))
}

cat(sprintf("\n%d checks, %d failed\n", NCHECK, length(FAILED)))
if (length(FAILED)) { cat("FAILED: ", paste(FAILED, collapse = "; "), "\nTESTS FAILED\n", sep = ""); quit(status = 1) }
cat("ALL TESTS PASSED\n")
