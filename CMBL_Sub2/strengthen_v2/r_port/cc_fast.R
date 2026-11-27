## cc_fast.R - batched (vectorized over permutations) versions of the permutation statistics of coexp_checks.R. Source it AFTER coexp_checks.R: it replaces .coupling_stats, .consensus_stack and
## .adjusted_stack (same arguments, same return values, same random stream: cc_perm_apply is unchanged) by array operations that need a few hundred R calls per chunk of permutations instead of tens of
## thousands. The statistics are the same quantities; sums are formed in a different order (rowSums in long double), so values can differ from the loop versions in the last bits (1e-16 relative);
## ties at the 1e-12 tolerance of .pvals are handled as before. The loop versions stay in coexp_checks.R as the reference (CC_SLOW = TRUE restores them).
## Stacks are arrays (c, G, N) = permutations x genes x samples. Only the Pearson kind and the cases without pair masks are batched; everything else falls back to the loop versions.
.slow <- list(coupling = .coupling_stats, consensus = .consensus_stack, adjusted = .adjusted_stack)
.pair_ij <- local({ cache <- list(); function(G) { k <- as.character(G); if (is.null(cache[[k]])) { up <- upper.tri(matrix(0, G, G)); cache[[k]] <<- list(i = row(up)[up], j = col(up)[up]) }; cache[[k]] } })
.batch_Z <- function(S) {              # S (c, G, N): rows centred and scaled to unit length, 0 for constant rows (the rule of cc_corr)
  d <- dim(S); N <- d[3]; mu <- rowMeans(S, dims = 2); Zc <- S - array(mu, d); ss <- rowSums(Zc * Zc, dims = 2)
  ok <- ss > 1e-24 * (ss + N * mu * mu); f <- ifelse(ok, 1 / sqrt(ifelse(ok, ss, 1)), 0); Zc * array(f, d)
}
.batch_pairs <- function(Z) {          # Z (c, G, N) from .batch_Z -> matrix (G(G-1)/2, c) of the correlations of all gene pairs (upper triangle, column-major pair order)
  d <- dim(Z); c <- d[1]; G <- d[2]; N <- d[3]; ij <- .pair_ij(G); Zl <- lapply(seq_len(G), function(g) matrix(Z[, g, ], c, N))
  out <- matrix(0, length(ij$i), c); for (k in seq_along(ij$i)) out[k, ] <- rowSums(Zl[[ij$i[k]]] * Zl[[ij$j[k]]]); out
}
.coupling_stats_fast <- function(Rb, des, masks = NULL, eig = TRUE) {
  if (!is.null(masks) && length(masks)) return(.slow$coupling(Rb, des, masks, eig))
  d <- dim(Rb); cn <- d[1]; G <- d[2]; N <- d[3]; Z <- .batch_Z(Rb); P <- .batch_pairs(Z)
  out <- list(T_r2 = 2 * colSums(P * P) / (G * (G - 1)), T_p = as.numeric(colSums(abs(P) >= des$rc_resid)))
  if (eig) out$T_eig <- vapply(seq_len(cn), function(b) { Zb <- matrix(Z[b, , ], G, N); if (N < G) eigen(crossprod(Zb), symmetric = TRUE, only.values = TRUE)$values[1] else eigen(tcrossprod(Zb), symmetric = TRUE, only.values = TRUE)$values[1] }, 1)
  out
}
.consensus_stack_fast <- function(stack, elig, den, rc, kind) {
  if (kind != "pearson") return(.slow$consensus(stack, elig, den, rc, kind))
  d <- dim(stack); cn <- d[1]; N <- d[3]; m <- length(elig); base <- stack[, elig, , drop = FALSE]; allsig <- NULL; agree <- NULL; pos0 <- NULL
  for (dd in c(list(NULL), as.list(den))) {
    M <- if (is.null(dd)) base else base / aperm(array(matrix(stack[, dd, ], cn, N), c(cn, N, m)), c(1, 3, 2))       # ratio to the denominator gene, broadcast over the genes
    P <- .batch_pairs(.batch_Z(M)); s <- abs(P) >= rc; pos <- P > 0
    allsig <- if (is.null(allsig)) s else allsig & s
    if (is.null(pos0)) { pos0 <- pos; agree <- matrix(TRUE, nrow(P), cn) } else agree <- agree & (pos == pos0)
  }
  list(a = as.integer(colSums(allsig)), s = as.integer(colSums(allsig & agree)))
}
.adjusted_stack_fast <- function(stack, des, sch, alpha = CC_ALPHA) {
  d <- dim(stack); cn <- d[1]; G <- d[2]; N <- d[3]; out <- matrix(0L, cn, length(sch))
  for (s in seq_along(sch)) {
    keep <- sch[[s]]$keep; dd <- sch[[s]]$d; m <- length(keep); M <- stack[, keep, , drop = FALSE]
    if (!is.null(dd)) M <- M / aperm(array(matrix(stack[, dd, ], cn, N), c(cn, N, m)), c(1, 3, 2))
    Mr <- array(matrix(M, cn * m, N) %*% des$Q, c(cn, m, N)); A <- abs(.batch_pairs(.batch_Z(Mr))); thr <- .bh_thresholds(des$df_resid, nrow(A), alpha)
    for (b in seq_len(cn)) { ok <- sort(A[, b], decreasing = TRUE) > thr; out[b, s] <- if (any(ok)) max(which(ok)) else 0L }
  }
  out
}
cc_use_fast <- function(on = TRUE) {
  e <- globalenv(); if (on) { assign(".coupling_stats", .coupling_stats_fast, e); assign(".consensus_stack", .consensus_stack_fast, e); assign(".adjusted_stack", .adjusted_stack_fast, e) }
  else { assign(".coupling_stats", .slow$coupling, e); assign(".consensus_stack", .slow$consensus, e); assign(".adjusted_stack", .slow$adjusted, e) }
  invisible(on)
}
cc_use_fast(!isTRUE(as.logical(Sys.getenv("CC_SLOW", "FALSE"))))
