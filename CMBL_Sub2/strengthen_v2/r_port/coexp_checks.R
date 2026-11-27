## coexp_checks.R - R port of strengthen_v2/coexp_checks.py (calibrated checks for pooled co-expression analyses of small designed experiments)
##
## Same definitions, same random streams: every random function draws from numpy_rng.R (numpy-compatible PCG64), so permutation counts, p-values
## and simulated data equal the Python results. Conventions (differences from Python are only the indexing):
##   X       samples x genes matrix on the LINEAR scale, strictly positive; rows in the order of `groups`.
##   groups  one label per sample, or a design made by cc_design().
##   den     denominator genes: names (with `genes`) or 1-BASED column indices. Gene indices (idx, ref_gene, keep) are 1-based everywhere.
##   runs    one label per gene (NA or "" = unresolved).
##   Stacks of permuted matrices are arrays (c, G, N) = permutations x genes x samples, as in numpy.
## Needs numpy_rng.R to be sourced first (rng/numpy_rng.R).

CC_VERSION <- "1.0 (WP1, 4 Nov 2026), R port"
## CC_CONSTANT_FIX=1 (added 27 Nov 2026, fix F3 of audit_271126/AUDIT_REPORT_271126.md; default off, so the port keeps reproducing the archived module and its
## numbers): the rules of coexp-checks 1.1.1 (tool_release_271126): a gene with the same value in every sample gets p = 1 in cc_anova_group_effect (F NaN),
## and .bh counts a NaN p as 1 instead of turning every q into NaN. Read once when this file is sourced; a test may set the variable directly.
CC_CONSTANT_FIX <- identical(Sys.getenv("CC_CONSTANT_FIX", ""), "1")
CC_MAX_ELEMS <- 1e6
CC_ALPHA <- 0.05

## ---------------------------------------------------------------------------------------------------------------- design
cc_rcrit <- function(df, alpha = CC_ALPHA) { t <- qt(1 - alpha / 2, df); t / sqrt(df + t * t) }

cc_design <- function(groups) {
  g <- as.vector(groups); labels <- unique(g); codes <- match(g, labels)
  N <- length(g); K <- length(labels)
  if (K < 2) stop("need at least two groups")
  sizes <- tabulate(codes, K)
  if (any(sizes < 2)) stop("every group needs at least two samples")
  rows <- lapply(seq_len(K), function(k) which(codes == k))
  df_pooled <- N - 2; df_resid <- N - K - 1
  if (df_resid < 3) stop(sprintf("residual df = N - K - 1 = %d is below 3", df_resid))
  A <- matrix(0, N, N); for (r in rows) A[r, r] <- 1 / length(r)
  structure(list(labels = labels, codes = codes, N = N, K = K, sizes = sizes, rows = rows, df_pooled = df_pooled, df_resid = df_resid,
                 rc_pooled = cc_rcrit(df_pooled), rc_resid = cc_rcrit(df_resid), Q = diag(N) - A), class = "cc_design")
}
as_design <- function(groups) if (inherits(groups, "cc_design")) groups else cc_design(groups)

.check_X <- function(X) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  if (!all(is.finite(X)) || any(X <= 0)) stop("X must be finite and strictly positive (linear scale); add a pseudocount for zeros")
  X
}
.names <- function(G, genes) if (!is.null(genes)) as.character(genes) else paste0("G", seq_len(G))
.den_idx <- function(den, genes) {
  if (is.null(den) || length(den) == 0) return(integer(0))
  out <- integer(0)
  for (d in den) {
    if (is.numeric(d)) out <- c(out, as.integer(d))
    else { if (is.null(genes) || !(d %in% genes)) stop(sprintf("denominator %s not among the genes", d)); out <- c(out, match(d, genes)) }
  }
  out
}
.unresolved <- function(b) is.null(b) || length(b) == 0 || is.na(b) || identical(as.character(b), "")

## ---------------------------------------------------------------------------------------------------------------- building blocks
cc_corr <- function(M, return_z = FALSE) {        # Pearson correlation of the rows of M (m x N); constant rows give r = 0
  N <- ncol(M); mu <- rowMeans(M); Z <- M - mu; ss <- rowSums(Z * Z); ok <- ss > 1e-24 * (ss + N * mu * mu)
  Z <- Z * ifelse(ok, 1 / sqrt(ifelse(ok, ss, 1)), 0); r <- tcrossprod(Z)
  if (return_z) list(r = r, Z = Z) else r
}
## r is clipped to [-1, 1] before the t transformation: noise rows of constant genes (residuals of 1e-17) are normalized to unit length and can correlate at 1 + 2e-16 with R's matrix product, which would give
## sqrt(negative) = NaN and (Python's _bh, below) turn every q of the scheme into NaN; |r| = 1 gives p = 0 as in numpy (added 24 Nov 2026, found on the degenerate series GSE20677 of WP2)
.pval <- function(r, df) { r <- pmin(pmax(r, -1), 1); t <- r * sqrt(df / (1 - r * r)); 2 * pt(-abs(t), df) }
.bh <- function(p) if (CC_CONSTANT_FIX) p.adjust(ifelse(is.na(p), 1, p), "BH") else if (anyNA(p)) rep(NaN, length(p)) else p.adjust(p, "BH")      # Python _bh: np.minimum.accumulate propagates a NaN p to every q; CC_CONSTANT_FIX: a NaN p counts as 1 (tool 1.1.1)
.bh_thresholds <- function(df, m, alpha) { t <- qt(alpha * seq_len(m) / m / 2, df, lower.tail = FALSE); t / sqrt(df + t * t) }
.bh_count <- function(absr, df, alpha = CC_ALPHA) {
  m <- length(absr); ok <- sort(absr, decreasing = TRUE) > .bh_thresholds(df, m, alpha)
  if (any(ok)) max(which(ok)) else 0L
}
.rank_rows <- function(M) { R <- M; for (i in seq_len(nrow(M))) R[i, ] <- rank(M[i, ]); R }
.transform <- function(M, kind) {
  if (kind == "pearson") M else if (kind == "log") log(M) else if (kind == "spearman") .rank_rows(M)
  else stop("kind must be 'pearson', 'log' or 'spearman'")
}
.schemes <- function(G, den) c(list(list(keep = seq_len(G), d = NULL)), lapply(den, function(d) list(keep = setdiff(seq_len(G), d), d = d)))
.scheme_matrix <- function(Xg, keep, d) {         # Xg: genes x samples
  M <- Xg[keep, , drop = FALSE]
  if (is.null(d)) M else M / matrix(Xg[d, ], nrow(M), ncol(M), byrow = TRUE)
}
.scheme_names <- function(genes, den) c("raw", if (length(den)) paste0(genes[den], "-normalized"))
## numpy sums over axis 0 of an (N, G) array add the rows one after the other (R's colSums/colMeans accumulate in long double): the difference matters for constant genes, where numpy's mean of
## identical values can differ from the value by an ulp (3v is rounded) and so decides between F = 0/0 = NaN and F = Inf (added 24 Nov 2026 with the degenerate-gene series of WP2)
.seq_colsums <- function(M) { s <- M[1, ]; for (i in seq_len(nrow(M))[-1]) s <- s + M[i, ]; s }
.seq_colmeans <- function(M) .seq_colsums(M) / nrow(M)
.group_means <- function(M, des) { out <- M; for (rows in des$rows) out[rows, ] <- matrix(.seq_colmeans(M[rows, , drop = FALSE]), length(rows), ncol(M), byrow = TRUE); out }
.blocks <- function(G, blocks) {
  if (is.null(blocks)) return(list(idx = seq_len(G), n = G))
  out <- integer(G); key <- list(); nxt <- 0L
  for (i in seq_len(G)) {
    b <- blocks[i]
    if (.unresolved(b)) { nxt <- nxt + 1L; out[i] <- nxt }
    else { k <- as.character(b); if (is.null(key[[k]])) { nxt <- nxt + 1L; key[[k]] <- nxt }; out[i] <- key[[k]] }
  }
  list(idx = out, n = nxt)
}
.perm_lex <- function(n) {                        # all permutations of 1..n in lexicographic order (itertools.permutations)
  if (n == 1) return(matrix(1L, 1, 1))
  sub <- .perm_lex(n - 1); out <- matrix(0L, factorial(n), n); r <- 0
  for (first in seq_len(n)) {
    rest <- setdiff(seq_len(n), first); blk <- matrix(rest[sub], nrow(sub), n - 1)
    out[r + seq_len(nrow(sub)), ] <- cbind(first, blk); r <- r + nrow(sub)
  }
  out
}
.perm_cache <- new.env()
cc_perm_table <- function(n) { k <- as.character(n); if (is.null(.perm_cache[[k]])) .perm_cache[[k]] <- .perm_lex(n); .perm_cache[[k]] }
.nz_rowmajor <- function(mask) {                  # (i, j) of the TRUE entries, row-major order like numpy.nonzero, 1-based
  w <- which(t(mask), arr.ind = TRUE)
  if (nrow(w) == 0) matrix(integer(0), 0, 2) else cbind(w[, 2], w[, 1])
}
.up <- function(G) upper.tri(matrix(0, G, G))

## The null engine. Python yields (c, G, N) stacks; here fun(stack) is called for every chunk and the results are returned as a list.
## Genes of one block share the permutation (blocks = NULL: every gene its own block). Groups of up to 8 samples draw from the table of all
## permutations, larger groups sort random keys. The random stream does not depend on the chunk size.
cc_perm_apply <- function(M, groups, nperm, seed = 2026, blocks = NULL, max_elems = NULL, fun) {
  des <- as_design(groups); M <- as.matrix(M); N <- nrow(M); G <- ncol(M)
  if (N != des$N) stop("groups and M disagree on the number of samples")
  bl <- .blocks(G, blocks); bog <- bl$idx; nb <- bl$n
  W <- sum(ifelse(des$sizes <= 8, nb, des$sizes * nb))
  rng <- if (inherits(seed, "environment")) seed else np_rng(seed)
  per <- max(1, (if (is.null(max_elems)) CC_MAX_ELEMS else max_elems) %/% (N * max(G, nb))); done <- 0; res <- list()
  while (done < nperm) {
    cn <- min(per, nperm - done)
    u <- matrix(np_random(rng, cn * W), nrow = cn, ncol = W, byrow = TRUE); off <- 0
    out <- array(0, dim = c(cn, G, N))
    for (k in seq_len(des$K)) {
      rows <- des$rows[[k]]; n <- des$sizes[k]
      if (n <= 8) {
        tab <- cc_perm_table(n); pidx <- floor(u[, off + seq_len(nb), drop = FALSE] * nrow(tab)) + 1
        perm <- array(tab[pidx, ], dim = c(cn, nb, n)); off <- off + nb
      } else {
        keys <- array(u[, off + seq_len(n * nb), drop = FALSE], dim = c(cn, n, nb))
        Kmat <- matrix(aperm(keys, c(2, 1, 3)), nrow = n)
        ord <- order(rep(seq_len(ncol(Kmat)), each = n), as.vector(Kmat))
        pm <- matrix(ord - rep((seq_len(ncol(Kmat)) - 1) * n, each = n), nrow = n)
        perm <- aperm(array(pm, dim = c(n, cn, nb)), c(2, 3, 1)); off <- off + n * nb
      }
      permg <- perm[, bog, , drop = FALSE]
      Mr <- M[rows, , drop = FALSE]
      idx <- permg + array(rep((seq_len(G) - 1) * n, each = cn), dim = c(cn, G, n))
      out[, , rows] <- Mr[idx]
    }
    res[[length(res) + 1]] <- fun(out); done <- done + cn
  }
  res
}
cc_permuted_all <- function(M, groups, nperm, seed = 2026, blocks = NULL, max_elems = NULL) {   # all stacks glued along the permutation axis
  ch <- cc_perm_apply(M, groups, nperm, seed, blocks, max_elems, function(s) s)
  d <- dim(ch[[1]]); tot <- sum(vapply(ch, function(s) dim(s)[1], 1)); out <- array(0, dim = c(tot, d[2], d[3])); p <- 0
  for (s in ch) { n <- dim(s)[1]; out[p + seq_len(n), , ] <- s; p <- p + n }
  out
}
.pvals <- function(null, obs) {
  n <- length(null); k <- sum(null >= obs - 1e-12 * max(1, abs(obs)))
  list(p_ge = k / n, p_ge_plus1 = (1 + k) / (n + 1))
}
.pop_sd <- function(x) sqrt(mean((x - mean(x))^2))
.desc <- function(null) list(null_mean = mean(null), null_sd = .pop_sd(null), null_q05 = unname(quantile(null, 0.05, type = 7)), null_q95 = unname(quantile(null, 0.95, type = 7)))
cc_residual_corr <- function(X, groups, log = TRUE) { X <- .check_X(X); des <- as_design(groups); cc_corr(t(if (log) base::log(X) else X) %*% des$Q) }

## ---------------------------------------------------------------------------------------------------------------- observed pipeline
cc_observed_summary <- function(X, groups, den = NULL, genes = NULL, kind = "pearson", alpha = CC_ALPHA) {
  X <- .check_X(X); des <- as_design(groups); N <- nrow(X); G <- ncol(X)
  if (N != des$N) stop("groups and X disagree on the number of samples")
  names <- .names(G, genes); den <- .den_idx(den, genes); sch <- .schemes(G, den); Xg <- t(X)
  S <- list(); SG <- list(); cnt <- integer(0); Q <- list(); qcnt <- integer(0)
  for (s in sch) {
    keep <- s$keep; d <- s$d; Tm <- .transform(.scheme_matrix(Xg, keep, d), kind); r <- cc_corr(Tm); diag(r) <- 0; m <- length(keep)
    sig <- matrix(FALSE, G, G); sgn <- matrix(0, G, G)
    sig[keep, keep] <- abs(r) >= des$rc_pooled; sgn[keep, keep] <- sign(r)
    S[[length(S) + 1]] <- sig; SG[[length(SG) + 1]] <- sgn; cnt <- c(cnt, sum(sig & .up(G)))
    iu <- .up(m); q <- matrix(NA_real_, m, m); rr <- cc_corr(Tm %*% des$Q); q[iu] <- .bh(.pval(rr[iu], des$df_resid))
    full <- matrix(NA_real_, G, G); full[keep, keep] <- q; Q[[length(Q) + 1]] <- full; qcnt <- c(qcnt, sum(q[iu] < alpha, na.rm = TRUE))
  }
  allsig <- Reduce(`&`, S); same <- allsig & (abs(Reduce(`+`, SG)) == length(SG)); up <- .up(G); cons <- allsig & up; cons_same <- same & up
  Qs <- simplify2array(Q); present <- !is.na(Qs)
  qmin <- apply(ifelse(present, Qs, Inf), c(1, 2), min); qmax <- apply(ifelse(present, Qs, -Inf), c(1, 2), max); in_all <- apply(present, c(1, 2), all)
  adj_any <- cons & (qmin < alpha); adj_all <- cons & in_all & (qmax < alpha); adj_pairs_all <- in_all & (qmax < alpha)
  ij <- .nz_rowmajor(cons); deg <- tabulate(c(ij[, 1], ij[, 2]), G)
  pr <- function(mask) { w <- .nz_rowmajor(mask); lapply(seq_len(nrow(w)), function(k) c(names[w[k, 1]], names[w[k, 2]])) }
  o <- order(-deg, seq_len(G)); o <- o[deg[o] > 0]
  list(scheme_names = .scheme_names(names, den), kind = kind, pooled_sig_per_scheme = cnt, consensus_any = sum(cons), consensus_same_sign = sum(cons_same),
       consensus_pairs = pr(cons), consensus_pairs_same_sign = pr(cons_same), degree = setNames(deg[o], names[o]),
       adj_q_lt_alpha_per_scheme = qcnt, cons_adj_any = sum(adj_any), cons_adj_all_schemes = sum(adj_all), adj_pairs_all_schemes = sum(adj_pairs_all),
       cons_adj_any_pairs = pr(adj_any), n_eligible_pairs = sum(in_all), .S = S, .SG = SG, .Q = Q, .cons_any = cons, .cons_same = cons_same, .den = den, .genes = names)
}

cc_group_mean_share <- function(X, groups, mask = NULL, kind = "pearson") {
  X <- .check_X(X); des <- as_design(groups); Y <- t(.transform(t(X), kind)); N <- nrow(Y)
  if (is.null(mask)) mask <- abs(cc_corr(t(Y))) >= des$rc_pooled
  mask <- as.matrix(mask) != 0 & .up(ncol(Y)); ij <- .nz_rowmajor(mask)
  if (nrow(ij) == 0) return(list(median = NaN, mean = NaN, pooled = NaN, n_pairs = 0L))
  mu <- colMeans(Y); xc <- Y - matrix(mu, N, ncol(Y), byrow = TRUE); bc <- .group_means(Y, des) - matrix(mu, N, ncol(Y), byrow = TRUE)
  tot <- crossprod(xc) / N; btw <- crossprod(bc) / N; sh <- btw[ij] / tot[ij]
  list(median = median(sh), mean = mean(sh), pooled = sum(btw[ij]) / sum(tot[ij]), n_pairs = nrow(ij))
}

## the share of cc_group_mean_share against what independence would give (tool version 1.1, 19 Nov 2026). Null: every gene's values are permuted independently
## within its group; the pooled-significant pairs are recomputed in every permuted matrix and the median share is taken again.
cc_group_mean_share_expected <- function(X, groups, nperm = 2000, seed = 2026, kind = "pearson") {
  X <- .check_X(X); des <- as_design(groups); Y <- t(.transform(t(X), kind)); N <- nrow(Y); G <- ncol(Y)
  obs <- cc_group_mean_share(X, des, kind = kind)
  bc <- .group_means(Y, des) - matrix(colMeans(Y), N, G, byrow = TRUE); btw <- crossprod(bc) / N; up <- .up(G); btwu <- btw[up]
  ch <- cc_perm_apply(Y, des, nperm, seed, NULL, fun = function(stack) {
    cn <- dim(stack)[1]; meds <- numeric(cn); nsig <- integer(cn)
    for (b in seq_len(cn)) {
      Yb <- matrix(stack[b, , ], G, N); r <- cc_corr(Yb); m <- abs(r[up]) >= des$rc_pooled
      xc <- Yb - rowMeans(Yb); tot <- tcrossprod(xc) / N; sh <- suppressWarnings(btwu / tot[up]); v <- sh[m]
      nsig[b] <- length(v); meds[b] <- if (length(v)) median(v) else NaN
    }
    list(meds = meds, nsig = nsig)
  })
  meds <- unlist(lapply(ch, function(x) x$meds)); nsig <- unlist(lapply(ch, function(x) x$nsig)); ok <- !is.na(meds)
  out <- list(observed_median = obs$median, observed_n_pairs = obs$n_pairs, null_n_pairs_mean = mean(nsig), null_no_sig_frac = 1 - mean(ok), nperm = nperm, kind = kind)
  if (any(ok)) {
    nm <- mean(meds[ok]); out$null_mean <- nm; out$null_q025 <- unname(quantile(meds[ok], 0.025, type = 7)); out$null_q975 <- unname(quantile(meds[ok], 0.975, type = 7))
    out$excess <- if (obs$n_pairs) obs$median - nm else NaN; out$above_q975 <- if (obs$n_pairs) obs$median > out$null_q975 else NA
  } else { out$null_mean <- NaN; out$null_q025 <- NaN; out$null_q975 <- NaN; out$excess <- NaN; out$above_q975 <- NA }
  out
}

## ---------------------------------------------------------------------------------------------------------------- within-group coupling
.coupling_stats <- function(Rb, des, masks = NULL, eig = TRUE) {      # Rb: (c, G, N) stack of log-scale residual matrices
  cn <- dim(Rb)[1]; G <- dim(Rb)[2]; N <- dim(Rb)[3]; up <- .up(G)
  T_r2 <- numeric(cn); T_p <- numeric(cn); T_eig <- if (eig) numeric(cn) else NULL
  mk <- lapply(masks, function(m) { m <- as.matrix(m) != 0; m & up }); Tm <- lapply(mk, function(m) numeric(cn)); nmk <- vapply(mk, function(m) max(sum(m), 1), 1)
  for (b in seq_len(cn)) {
    cz <- cc_corr(matrix(Rb[b, , ], G, N), TRUE); r <- cz$r; r2 <- r * r
    T_r2[b] <- (sum(r2) - sum(diag(r2))) / (G * (G - 1)); T_p[b] <- sum((abs(r) >= des$rc_resid) & up)
    if (eig) T_eig[b] <- (if (N < G) eigen(crossprod(cz$Z), symmetric = TRUE, only.values = TRUE)$values else eigen(r, symmetric = TRUE, only.values = TRUE)$values)[1]
    for (nm in names(mk)) Tm[[nm]][b] <- sum(r2 * mk[[nm]]) / nmk[[nm]]
  }
  out <- list(T_r2 = T_r2, T_p = T_p); if (eig) out$T_eig <- T_eig
  for (nm in names(mk)) out[[paste0("T_r2_", nm)]] <- Tm[[nm]]
  out
}
cc_coupling_test <- function(X, groups, nperm = 2000, seed = 2026, blocks = NULL, pair_masks = NULL, eig = TRUE) {
  X <- .check_X(X); des <- as_design(groups); R <- des$Q %*% log(X)
  obs <- lapply(.coupling_stats(array(t(R), dim = c(1, ncol(R), nrow(R))), des, pair_masks, eig), function(v) v[1])
  ch <- cc_perm_apply(R, des, nperm, seed, blocks, fun = function(Rb) .coupling_stats(Rb, des, pair_masks, eig))
  setNames(lapply(names(obs), function(k) { nul <- unlist(lapply(ch, function(x) x[[k]])); c(list(observed = obs[[k]]), .desc(nul), .pvals(nul, obs[[k]])) }), names(obs))
}

## ---------------------------------------------------------------------------------------------------------------- consensus and adjusted counts
.consensus_one <- function(Xb, elig, den, rc, kind) {                  # Xb: G x N; counts (any sign, same sign) of pooled-significant pairs in all schemes
  allsig <- NULL; agree <- NULL; pos0 <- NULL
  for (d in c(list(NULL), as.list(den))) {
    M <- if (is.null(d)) Xb[elig, , drop = FALSE] else Xb[elig, , drop = FALSE] / matrix(Xb[d, ], length(elig), ncol(Xb), byrow = TRUE)
    r <- cc_corr(.transform(M, kind)); s <- abs(r) >= rc; pos <- r > 0
    allsig <- if (is.null(allsig)) s else allsig & s
    if (is.null(pos0)) { pos0 <- pos; agree <- matrix(TRUE, nrow(s), ncol(s)) } else agree <- agree & (pos == pos0)
  }
  up <- .up(nrow(allsig)); c(sum(allsig & up), sum(allsig & agree & up))
}
.consensus_stack <- function(stack, elig, den, rc, kind) {
  cn <- dim(stack)[1]; G <- dim(stack)[2]; N <- dim(stack)[3]; a <- integer(cn); s <- integer(cn)
  for (b in seq_len(cn)) { v <- .consensus_one(matrix(stack[b, , ], G, N), elig, den, rc, kind); a[b] <- v[1]; s[b] <- v[2] }
  list(a = a, s = s)
}
cc_consensus_null <- function(X, groups, den = NULL, nperm = 2000, seed = 2026, blocks = NULL, kind = "pearson", genes = NULL) {
  X <- .check_X(X); des <- as_design(groups); G <- ncol(X); den <- .den_idx(den, genes); elig <- setdiff(seq_len(G), den)
  ov <- .consensus_one(t(X), elig, den, des$rc_pooled, kind); oa <- ov[1]; os <- ov[2]
  ch <- cc_perm_apply(X, des, nperm, seed, blocks, fun = function(Xb) .consensus_stack(Xb, elig, den, des$rc_pooled, kind))
  ca <- unlist(lapply(ch, function(x) x$a)); cs <- unlist(lapply(ch, function(x) x$s)); sdv <- .pop_sd(ca)
  c(list(observed = oa), .desc(ca), list(z = if (sdv > 0) (oa - mean(ca)) / sdv else NaN), .pvals(ca, oa),
    list(n_eligible_pairs = (length(elig) * (length(elig) - 1)) %/% 2, same_sign = c(list(observed = os), .desc(cs), .pvals(cs, os)), nperm = nperm, kind = kind, .null = ca, .null_same = cs))
}
.adjusted_stack <- function(stack, des, sch, alpha = CC_ALPHA) {
  cn <- dim(stack)[1]; G <- dim(stack)[2]; N <- dim(stack)[3]; out <- matrix(0L, cn, length(sch))
  for (b in seq_len(cn)) {
    Xb <- matrix(stack[b, , ], G, N)
    for (s in seq_along(sch)) { r <- cc_corr(.scheme_matrix(Xb, sch[[s]]$keep, sch[[s]]$d) %*% des$Q); out[b, s] <- .bh_count(abs(r[.up(nrow(r))]), des$df_resid, alpha) }
  }
  out
}
cc_adjusted_null <- function(X, groups, den = NULL, nperm = 2000, seed = 2026, blocks = NULL, genes = NULL, alpha = CC_ALPHA) {
  X <- .check_X(X); des <- as_design(groups); G <- ncol(X); den <- .den_idx(den, genes); sch <- .schemes(G, den); names <- .scheme_names(.names(G, genes), den)
  obs <- .adjusted_stack(array(t(X), dim = c(1, G, nrow(X))), des, sch, alpha)[1, ]
  nul <- do.call(rbind, cc_perm_apply(X, des, nperm, seed, blocks, fun = function(Xb) .adjusted_stack(Xb, des, sch, alpha)))
  schemes <- setNames(lapply(seq_along(names), function(s) c(list(observed = obs[s]), .desc(nul[, s]), list(excess = obs[s] - mean(nul[, s])), .pvals(nul[, s], obs[s]))), names)
  list(scheme_names = names, nperm = nperm, alpha = alpha, schemes = schemes, .null = nul)
}

## ---------------------------------------------------------------------------------------------------------------- run structure
cc_same_run_mask <- function(runs) {
  G <- length(runs); ok <- !vapply(seq_len(G), function(i) .unresolved(runs[i]), TRUE); lab <- as.character(runs)
  eq <- outer(lab, lab, "=="); both <- outer(ok, ok, "&"); up <- .up(G); list(same = eq & both & up, between = !eq & both & up)
}
cc_run_structure_test <- function(X, groups, runs, nperm = 20000, seed = 2026, genes = NULL, R = NULL) {
  X <- .check_X(X); des <- as_design(groups); names <- .names(ncol(X), genes); r <- if (is.null(R)) cc_residual_corr(X, des) else as.matrix(R)
  keep <- which(!vapply(seq_along(runs), function(i) .unresolved(runs[i]), TRUE)); lab <- runs[keep]; uniq <- unique(lab); code <- match(lab, uniq); g <- length(keep); K <- length(uniq)
  rr <- r[keep, keep, drop = FALSE]; iu <- which(.up(g), arr.ind = TRUE); P <- nrow(iu); nk <- tabulate(code, K); ni <- (nk * (nk - 1)) %/% 2; nc <- nk * (g - nk); nS <- sum(ni); nB <- P - nS
  vals <- list(r = rr[iu], abs_r = abs(rr[iu])); rng <- np_rng(seed); per <- max(1, 2000000 %/% max(P, 1))
  stat <- function(cp) {                                              # cp: c x g matrix of permuted run codes
    si <- cp[, iu[, 1], drop = FALSE]; sj <- cp[, iu[, 2], drop = FALSE]; same <- (si == sj) + 0; out <- list()
    for (tag in names(vals)) {
      v <- vals[[tag]]; Dk <- matrix(NA_real_, nrow(cp), K)
      for (k in seq_len(K)) if (ni[k] > 0 && nc[k] > 0) Dk[, k] <- drop(((si == k) & (sj == k)) %*% v) / ni[k] - drop(xor(si == k, sj == k) %*% v) / nc[k]
      out[[tag]] <- list(D = drop(same %*% v) / max(nS, 1) - drop((1 - same) %*% v) / max(nB, 1), Dk = Dk)
    }
    out
  }
  obs <- stat(matrix(code, 1)); nul <- lapply(vals, function(.) list(D = list(), Dk = list())); done <- 0
  while (done < nperm) {
    cn <- min(per, nperm - done); u <- matrix(np_random(rng, cn * g), cn, g, byrow = TRUE)
    cp <- matrix(0L, cn, g); for (b in seq_len(cn)) cp[b, ] <- code[order(u[b, ])]
    st <- stat(cp)
    for (tag in names(vals)) { nul[[tag]]$D[[length(nul[[tag]]$D) + 1]] <- st[[tag]]$D; nul[[tag]]$Dk[[length(nul[[tag]]$Dk) + 1]] <- st[[tag]]$Dk }
    done <- done + cn
  }
  res <- list(n_genes = g, runs = setNames(as.integer(nk), uniq), n_same_pairs = nS, n_between_pairs = nB, nperm = nperm, excluded = names[setdiff(seq_along(names), keep)])
  same0 <- code[iu[, 1]] == code[iu[, 2]]
  for (tag in names(vals)) {
    v <- vals[[tag]]; d0 <- obs[[tag]]$D[1]; dn <- unlist(nul[[tag]]$D); dkn <- do.call(rbind, nul[[tag]]$Dk)
    ent <- list(mean_same = if (nS) mean(v[same0]) else NaN, mean_between = if (nB) mean(v[!same0]) else NaN, D = d0, null_sd = .pop_sd(dn),
                p_two_sided = (1 + sum(abs(dn) >= abs(d0) - 1e-12)) / (nperm + 1), per_run = list())
    for (k in seq_len(K)) {
      dk0 <- obs[[tag]]$Dk[1, k]; inside <- code[iu[, 1]] == k & code[iu[, 2]] == k; cross <- xor(code[iu[, 1]] == k, code[iu[, 2]] == k)
      ent$per_run[[uniq[k]]] <- list(n_genes = nk[k], n_pairs_inside = ni[k], mean_inside = if (ni[k]) mean(v[inside]) else NaN, mean_cross = if (nc[k]) mean(v[cross]) else NaN,
                                     D = dk0, p_two_sided = if (!is.na(dk0)) (1 + sum(abs(dkn[, k]) >= abs(dk0) - 1e-12)) / (nperm + 1) else NaN)
    }
    res[[tag]] <- ent
  }
  res
}

## ---------------------------------------------------------------------------------------------------------------- geometry
cc_geometric_expectation <- function(r_crit, K) if (K < 3) 1.0 else pbeta(1 - r_crit^2, (K - 2) / 2, 0.5)
cc_anova_group_effect <- function(X, groups, log = TRUE) {
  X <- .check_X(X); des <- as_design(groups); Y <- if (log) base::log(X) else X; gm <- .group_means(Y, des); N <- des$N; K <- des$K
  ssw <- .seq_colsums((Y - gm)^2); ssb <- .seq_colsums((gm - matrix(.seq_colmeans(Y), N, ncol(Y), byrow = TRUE))^2)
  Fv <- (ssb / (K - 1)) / (ssw / (N - K)); p <- suppressWarnings(pf(Fv, K - 1, N - K, lower.tail = FALSE))
  if (CC_CONSTANT_FIX) { const <- colSums(Y != matrix(Y[1, ], N, ncol(Y), byrow = TRUE)) == 0; Fv[const] <- NaN; p[const] <- 1; return(list(F = Fv, p = p, q = .bh(p))) }   # tool 1.1.1
  list(F = Fv, p = p, q = if (anyNA(p)) rep(NaN, length(p)) else .bh(p))      # Python _bh: one NaN p (a constant gene with exactly zero between-group SS: 0/0) turns every q into NaN (np.minimum.accumulate propagates NaN)
}
cc_geometry_check <- function(X, groups, fdr = CC_ALPHA, log = TRUE) {
  X <- .check_X(X); des <- as_design(groups); G <- ncol(X); up <- .up(G); sig <- (abs(cc_corr(t(X))) >= des$rc_pooled) & up
  de <- cc_anova_group_effect(X, des, log)$q < fdr; de[is.na(de)] <- FALSE; both <- outer(de, de, "&") & up       # NaN < fdr is False in numpy
  list(r_crit = des$rc_pooled, expected = cc_geometric_expectation(des$rc_pooled, des$K), n_pairs_all = sum(up), n_sig_all = sum(sig), share_sig_all = sum(sig) / sum(up),
       n_group_effect_genes = sum(de), n_pairs_group_effect = sum(both), n_sig_group_effect = sum(sig & both),
       share_sig_group_effect = if (sum(both)) sum(sig & both) / sum(both) else NaN)
}

## ---------------------------------------------------------------------------------------------------------------- simulation helpers (protocol 3.2 and the tests)
cc_decouple <- function(X, groups, rng) {       # permute each gene's log values independently within every group; returns linear scale
  X <- .check_X(X); des <- as_design(groups); L <- log(X)
  for (rows in des$rows) {
    n <- length(rows); u <- matrix(np_random(rng, n * ncol(L)), n, ncol(L), byrow = TRUE); Lr <- L[rows, , drop = FALSE]
    for (j in seq_len(ncol(L))) L[rows, j] <- Lr[order(u[, j]), j]
  }
  exp(L)
}
cc_inject_factor <- function(X, groups, idx, r, rng) {   # idx: 1-based gene indices
  X <- .check_X(X); des <- as_design(groups); L <- log(X); f <- np_standard_normal(rng, des$N)
  for (rows in des$rows) f[rows] <- f[rows] - mean(f[rows])
  f <- f / sqrt(sum(f^2) / (des$N - des$K)); sd <- sqrt(colSums((des$Q %*% L)^2) / (des$N - des$K))
  for (j in idx) L[, j] <- L[, j] + sqrt(r / (1 - r)) * sd[j] * f
  exp(L)
}
cc_simulate <- function(n_per_group, G, rng, effect_sd = 1.0, noise_sd = 0.3, ref_noise_sd = 0.0, ref_gene = NULL) {   # ref_gene: 1-based indices
  K <- length(n_per_group); codes <- rep(seq_len(K), n_per_group); Nt <- length(codes)
  mu <- matrix(np_normal(rng, 0, effect_sd, K * G), K, G, byrow = TRUE)
  L <- mu[codes, , drop = FALSE] + matrix(np_normal(rng, 0, noise_sd, Nt * G), Nt, G, byrow = TRUE)
  if (!is.null(ref_gene) && ref_noise_sd > 0) for (j in ref_gene) L[, j] <- L[, j] + np_normal(rng, 0, ref_noise_sd, Nt)
  list(X = exp(L), codes = codes)
}

## ---------------------------------------------------------------------------------------------------------------- all-in-one
CC_ALL <- c("observed", "share", "coupling", "consensus_null", "adjusted_null", "geometry", "runs")
cc_run_checks <- function(X, groups, genes = NULL, den = NULL, runs = NULL, nperm = 2000, run_nperm = 20000, seed = 2026, which = NULL) {
  X <- .check_X(X); des <- as_design(groups); N <- nrow(X); G <- ncol(X); names <- .names(G, genes); denx <- .den_idx(den, genes)
  which <- if (is.null(which)) CC_ALL else which; res <- list(version = CC_VERSION)
  res$meta <- list(N = N, K = des$K, group_sizes = des$sizes, G = G, genes = names, denominators = names[denx], df_pooled = des$df_pooled, df_resid = des$df_resid,
                   r_crit_pooled = des$rc_pooled, r_crit_resid = des$rc_resid, nperm = nperm, run_nperm = run_nperm, seed = seed)
  W <- character(0)
  if (any(des$sizes < 3)) W <- c(W, "a group has fewer than 3 samples: within-group correlations rest on very few degrees of freedom")
  if (des$df_resid < 6) W <- c(W, sprintf("residual df = %d: the within-group tests have little power (r_crit = %.2f)", des$df_resid, des$rc_resid))
  if (G < 5) W <- c(W, sprintf("only %d genes: T_r2 and the network statistics are unstable", G))
  if (!is.null(runs) && length(runs) != G) stop("runs needs one label per gene")
  use_runs <- !is.null(runs) && "runs" %in% which
  if ("observed" %in% which || "share" %in% which) {
    obs <- cc_observed_summary(X, des, denx, names)
    if ("observed" %in% which) res$observed <- obs[!startsWith(names(obs), ".")]
    if ("share" %in% which) {
      gm <- cc_group_mean_share(X, des, obs$.S[[1]] & .up(G)); res$group_mean_share <- gm
      if (gm$median > 0.9) W <- c(W, sprintf("group means carry %.0f%% of the pooled covariance of significant raw pairs: pooled networks mostly restate the design", 100 * gm$median))
    }
  }
  if ("geometry" %in% which) res$geometry <- cc_geometry_check(X, des)
  if ("coupling" %in% which) {
    pm <- NULL; if (use_runs) { sm <- cc_same_run_mask(runs); pm <- list(same_run = sm$same, between_run = sm$between) }
    res$coupling <- ct <- cc_coupling_test(X, des, nperm, seed, pair_masks = pm)
    if (ct$T_r2$p_ge_plus1 < CC_ALPHA) W <- c(W, sprintf("T_r2 = %.3f is above the within-group permutation null (P = %.4f): within-group coupling exists; it does not by itself identify biology", ct$T_r2$observed, ct$T_r2$p_ge_plus1))
  }
  if ("consensus_null" %in% which) {
    cn <- cc_consensus_null(X, des, denx, nperm, seed, genes = names); res$consensus_null <- cn[!startsWith(names(cn), ".")]
    if (cn$p_ge_plus1 >= CC_ALPHA) W <- c(W, sprintf("consensus count %d is within the permutation null (mean %.1f, P = %.2f): weak evidence of no coupling, because the count test is insensitive to injected coupling in public benchmark data; read T_r2", as.integer(cn$observed), cn$null_mean, cn$p_ge_plus1))
    if (use_runs) { cb <- cc_consensus_null(X, des, denx, nperm, seed, blocks = runs, genes = names); res$consensus_null_run_block <- cb[!startsWith(names(cb), ".")] }
  }
  if ("adjusted_null" %in% which) {
    an <- cc_adjusted_null(X, des, denx, nperm, seed, genes = names); res$adjusted_calibrated <- an[!startsWith(names(an), ".")]
    for (n in names(an$schemes)) { v <- an$schemes[[n]]
      if (n != "raw" && v$observed > 0 && v$p_ge_plus1 >= CC_ALPHA) W <- c(W, sprintf("%s scheme: %d pairs with BH q < 0.05 is not above its null (mean %.1f, P = %.2f): BH is not calibrated for a ratio scheme (shared denominator)", n, as.integer(v$observed), v$null_mean, v$p_ge_plus1)) }
  }
  if (use_runs) {
    rs <- cc_run_structure_test(X, des, runs, run_nperm, seed, names); res$run_structure <- rs
    if (rs$r$p_two_sided < CC_ALPHA) W <- c(W, sprintf("within-group r differs between same-run and between-run pairs (D = %.2f, P = %.4f): compatible with shared reference, cDNA or plate effects, not co-regulation", rs$r$D, rs$r$p_two_sided))
  }
  res$warnings <- as.list(W); res
}
