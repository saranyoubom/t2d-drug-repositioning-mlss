## bench_pipeline.R - R port of benchmark_route3/bench_pipeline.py (generalized pooled co-expression pipeline of the manuscript for any 12 x G matrix,
## 3 groups x 4 replicates in file order, any three denominator genes). Needs numpy_rng.R sourced first. Gene indices are 1-based.
B_N <- 12; B_K <- 3
B_GROUP <- rep(1:3, each = 4)
B_DESIGN <- cbind(1, as.numeric(B_GROUP == 2), as.numeric(B_GROUP == 3))
B_T <- qt(0.975, B_N - 2); B_RC10 <- B_T / sqrt(B_N - 2 + B_T^2)       # |r| that gives p = 0.05 at df 10 (pooled analysis)
B_DF_ADJ <- B_N - (B_K - 1) - 2                                           # 8: treatment-adjusted partial correlation

B_scheme_list <- function(G, den) c(list(list(keep = seq_len(G), d = NULL)), lapply(den, function(d) list(keep = setdiff(seq_len(G), d), d = d)))
B_scheme_matrix <- function(X, keep, d) if (is.null(d)) X[, keep, drop = FALSE] else X[, keep, drop = FALSE] / matrix(X[, d], nrow(X), length(keep))

B_pooled <- function(X, den, rc = B_RC10) {                  # per scheme: significance and sign matrices in the full G x G gene space
  G <- ncol(X); S <- list(); SG <- list(); cnt <- integer(0)
  for (s in B_scheme_list(G, den)) {
    keep <- s$keep; M <- B_scheme_matrix(X, keep, s$d); r <- suppressWarnings(cor(M)); r[is.na(r)] <- 0; diag(r) <- 0
    sig <- matrix(FALSE, G, G); sgn <- matrix(0, G, G); sig[keep, keep] <- abs(r) >= rc; sgn[keep, keep] <- sign(r)
    S[[length(S) + 1]] <- sig; SG[[length(SG) + 1]] <- sgn; cnt <- c(cnt, sum(sig & upper.tri(sig)))
  }
  list(S = S, SG = SG, cnt = cnt)
}
B_consensus_from <- function(S, SG) {
  allsig <- Reduce(`&`, S); same <- allsig & (abs(Reduce(`+`, SG)) == length(SG)); up <- upper.tri(allsig)
  list(any = allsig & up, same = same & up)
}
B_consensus_counts <- function(X, den) { p <- B_pooled(X, den); cs <- B_consensus_from(p$S, p$SG); c(sum(cs$any), sum(cs$same)) }

## null distribution of the consensus count: keep each gene's group means, permute its within-group log-scale residuals independently across genes
## (rng.permuted(R[ix], axis=0) shuffles every column with the same Fisher-Yates as Generator.shuffle, columns in order)
B_permutation_counts <- function(X, den, nperm, rng) {
  L <- log(X); G <- ncol(L); gm <- t(vapply(1:B_K, function(g) colMeans(L[B_GROUP == g, , drop = FALSE]), numeric(G))); resid <- L - gm[B_GROUP, , drop = FALSE]
  ca <- integer(nperm); cs <- integer(nperm); idx <- lapply(1:B_K, function(g) which(B_GROUP == g))
  for (b in seq_len(nperm)) {
    R <- resid
    for (ix in idx) { blk <- R[ix, , drop = FALSE]; for (j in seq_len(G)) blk[, j] <- np_shuffle(rng, blk[, j]); R[ix, ] <- blk }
    v <- B_consensus_counts(exp(gm[B_GROUP, , drop = FALSE] + R), den); ca[b] <- v[1]; cs[b] <- v[2]
  }
  list(ca = ca, cs = cs)
}

## treatment-adjusted partial correlations per scheme: q-value matrices in the full gene space (NA where absent) and q < 0.05 counts
B_adjusted <- function(X, den) {
  G <- ncol(X); Q <- list(); cnt <- integer(0)
  for (s in B_scheme_list(G, den)) {
    keep <- s$keep; M <- B_scheme_matrix(X, keep, s$d); R <- qr.resid(qr(B_DESIGN), M)
    r <- suppressWarnings(cor(R)); t <- r * sqrt(B_DF_ADJ / (1 - r^2)); p <- 2 * pt(-abs(t), B_DF_ADJ); m <- length(keep); iu <- upper.tri(p)
    pv <- p[iu]; pv[is.na(pv)] <- 1; q <- matrix(NA_real_, m, m); q[iu] <- p.adjust(pv, "BH")
    full <- matrix(NA_real_, G, G); full[keep, keep] <- q; Q[[length(Q) + 1]] <- full; cnt <- c(cnt, sum(q[iu] < 0.05))
  }
  list(Q = Q, cnt = cnt)
}

## median share of the pooled covariance carried by the group means, over the pairs in sigmask (raw linear scale)
B_between_share <- function(X, sigmask) {
  xm <- t(vapply(1:B_K, function(g) colMeans(X[B_GROUP == g, , drop = FALSE]), numeric(ncol(X))))[B_GROUP, , drop = FALSE]
  mu <- colMeans(X); xc <- X - matrix(mu, nrow(X), ncol(X), byrow = TRUE); bc <- xm - matrix(mu, nrow(X), ncol(X), byrow = TRUE)
  cov_tot <- crossprod(xc) / B_N; cov_btw <- crossprod(bc) / B_N; w <- which(sigmask, arr.ind = TRUE)
  if (nrow(w) == 0) return(NaN)
  median(cov_btw[w] / cov_tot[w])
}

## all endpoints of the protocol for one 12 x G matrix
B_analyse <- function(X, den, nperm = 2000, seed = 2026) {
  X <- as.matrix(X); storage.mode(X) <- "double"; p <- B_pooled(X, den); cs <- B_consensus_from(p$S, p$SG)
  rng <- np_rng(seed); pc <- B_permutation_counts(X, den, nperm, rng); ca <- pc$ca; csn <- pc$cs
  aj <- B_adjusted(X, den); obs <- sum(cs$any)
  Qs <- simplify2array(aj$Q); present <- !is.na(Qs)
  qmin <- apply(ifelse(present, Qs, Inf), c(1, 2), min); qall <- apply(ifelse(present, Qs, -Inf), c(1, 2), max)
  qall[!is.finite(qall)] <- 1          # numpy: nan_to_num(nanmax(...), nan = 1): a pair is judged on the schemes that contain it; only a pair in no scheme (the diagonal) gets 1 (fixed 24 Nov 2026: this line also set 1 for pairs missing from some scheme)
  adj_any <- cs$any & (qmin < 0.05); adj_all4 <- cs$any & (qall < 0.05); allq4 <- (qall < 0.05) & upper.tri(qall)
  sdv <- sqrt(mean((ca - mean(ca))^2)); sdc <- sum(cs$same)
  list(observed = obs, observed_same_sign = sdc, pooled_sig_per_scheme = p$cnt, null_mean = mean(ca), null_sd = sdv, null_q05 = unname(quantile(ca, 0.05, type = 7)),
       null_q95 = unname(quantile(ca, 0.95, type = 7)), perm_p = (1 + sum(ca >= obs)) / (nperm + 1), perm_p_plain = mean(ca >= obs), z = if (sdv > 0) (obs - mean(ca)) / sdv else NaN,
       null_mean_same_sign = mean(csn), perm_p_same_sign = (1 + sum(csn >= sdc)) / (nperm + 1), adj_q05_per_scheme = aj$cnt, cons_adj_any = sum(adj_any),
       cons_adj_all4 = sum(adj_all4), adj_pairs_all4 = sum(allq4), share_between_raw = B_between_share(X, cs$any & FALSE | (p$S[[1]] & upper.tri(p$S[[1]]))))
}
