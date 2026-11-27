## pipeline.R - R port of reanalysis_17genes/pipeline.py: pooled co-expression pipeline of the manuscript for the INS-1 matrix (12 observations x genes, 2^-dCt;
## rows 1-4 CTRL, 5-8 DAPT, 9-12 DKK-1), four schemes (Raw and the Ins1-, Kcnj11- and Hes1-normalized matrices), consensus pairs, treatment-adjusted partial
## correlations, permutation null (numpy-compatible generator, seed 2026). Library plus main: `Rscript pipeline.R` writes pipeline_20genes.json and pipeline_17genes.json.
## The caller sets R_PORT (folder of numpy_rng.R and pycompat.R) and sources rng/numpy_rng.R and pycompat.R first. Pairs are keyed by the C-locale sorted gene names (Python tuples).
suppressMessages(library(readxl))
P_GROUP <- rep(1:3, each = 4); P_DEN <- c("Ins1", "Kcnj11", "Hes1")
P_GENES20 <- c("Wnt2", "Wnt2b", "Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Hey1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
P_EXCLUDED <- c("Hey1", "Wnt2", "Wnt2b", "Ryr1"); P_GENES17 <- P_GENES20[!(P_GENES20 %in% P_EXCLUDED)]
## the authors' original matrix (Submission/SBA/Figures, as the Python script), else the copy shipped in R_17genes (same values rounded to 16 digits); INS1_MATRIX overrides
P_matrix_path <- function() {
  cand <- c(Sys.getenv("INS1_MATRIX", ""), file.path(R_PORT, "..", "..", "..", "..", "SBA", "Figures", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"),
            file.path(R_PORT, "..", "..", "R_17genes", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"))
  cand <- cand[nzchar(cand) & file.exists(cand)]; if (!length(cand)) stop("INS-1 matrix not found; set INS1_MATRIX"); normalizePath(cand[1])
}
P_load <- function(genes) { df <- as.data.frame(read_excel(P_matrix_path())); M <- as.matrix(df[, genes, drop = FALSE]); storage.mode(M) <- "double"; rownames(M) <- NULL; M }
P_schemes <- function(M) {
  out <- list(Raw = M)
  for (d in P_DEN) if (d %in% colnames(M)) out[[paste0(d, "-normalized")]] <- M[, colnames(M) != d, drop = FALSE] / M[, d]
  out
}
P_corr_p <- function(M, df_resid = NULL) {
  r <- cor(M); dfree <- if (is.null(df_resid)) nrow(M) - 2 else df_resid
  t <- r * sqrt(dfree / (1 - r^2)); p <- 2 * pt(abs(t), dfree, lower.tail = FALSE); diag(p) <- NA
  list(r = r, p = p)
}
P_sig_sets <- function(sch) lapply(sch, function(m) { cp <- P_corr_p(m); list(genes = colnames(m), r = cp$r, p = cp$p) })
## the sub-matrices of every scheme over the genes common to all schemes, in C-locale name order (so i < j gives the sorted pair)
.P_common <- function(objs) py_sort(Reduce(intersect, lapply(objs, function(s) s$genes)))
.P_sub <- function(objs, field, cg) lapply(objs, function(s) { ix <- match(cg, s$genes); s[[field]][ix, ix, drop = FALSE] })
.P_pairs_rowmajor <- function(keep, cg) {            # pairs (name1, name2) of the upper triangle of a logical matrix, sorted by first then second name
  w <- which(keep & upper.tri(keep), arr.ind = TRUE); if (!nrow(w)) return(matrix(character(0), 0, 2)); w <- w[order(w[, 1], w[, 2]), , drop = FALSE]
  cbind(cg[w[, 1]], cg[w[, 2]])
}
P_consensus <- function(res, alpha = 0.05, use_sign = TRUE) {
  cg <- .P_common(res); P <- .P_sub(res, "p", cg); Rm <- .P_sub(res, "r", cg)
  ok <- Reduce(`&`, lapply(P, function(p) p <= alpha)); ok[is.na(ok)] <- FALSE
  if (use_sign) { s1 <- sign(Rm[[1]]); ok <- ok & Reduce(`&`, lapply(Rm[-1], function(r) sign(r) == s1)) }
  ok[is.na(ok)] <- FALSE
  .P_pairs_rowmajor(ok, cg)
}
## counts only (permutation loop): number of consensus pairs, any sign and sign-consistent
P_consensus_n <- function(res, alpha = 0.05) {
  cg <- .P_common(res); P <- .P_sub(res, "p", cg); Rm <- .P_sub(res, "r", cg)
  ok <- Reduce(`&`, lapply(P, function(p) p <= alpha)); ok[is.na(ok)] <- FALSE; s1 <- sign(Rm[[1]]); same <- Reduce(`&`, lapply(Rm[-1], function(r) sign(r) == s1)); same[is.na(same)] <- FALSE
  up <- upper.tri(ok); c(sum(ok & up), sum(ok & same & up))
}
P_counts <- function(res) lapply(res, function(s) {
  up <- upper.tri(s$p); p <- s$p[up]; r <- s$r[up]; sig <- p < 0.05
  list(pairs = length(p), sig = sum(sig), strong = sum(abs(r[sig]) >= 0.7), very = sum(abs(r[sig]) >= 0.8))
})
P_adjusted <- function(sch) {
  X <- cbind(1, as.numeric(P_GROUP == 2), as.numeric(P_GROUP == 3)); out <- list()
  for (k in names(sch)) {
    m <- sch[[k]]; Rm <- qr.resid(qr(X), m); cp <- P_corr_p(Rm, df_resid = 12 - 2 - 2); up <- upper.tri(cp$p); q <- matrix(NA_real_, ncol(m), ncol(m)); q[up] <- p.adjust(cp$p[up], "BH")
    q[lower.tri(q)] <- t(q)[lower.tri(q)]                       # symmetric (NA diagonal), so that re-indexing to another gene order keeps every pair
    out[[k]] <- list(genes = colnames(m), r = cp$r, p = cp$p, q = q)
  }
  out
}
## residual permutation null of the consensus count: group means kept, log-scale residuals of every gene permuted within group, columns in order (numpy Generator.permutation)
## Returns ca / cs (consensus counts, any sign / sign-consistent) and sig (per scheme, the number of pairs with p < 0.05) of every permutation. Within one R session the result is cached:
## pipeline.py, perm_arrays.py and perm_detail.py each redraw the same seeded streams, so the cached draws are the ones each script would make again (a standalone run computes them afresh).
.P_perm_cache <- new.env()
P_permutation_null <- function(M, n = 2000, seed = 2026) {
  key <- paste(paste(colnames(M), collapse = ","), n, seed, sprintf("%.17g", sum(M)), sep = "|"); if (!is.null(.P_perm_cache[[key]])) return(.P_perm_cache[[key]])
  rng <- np_rng(seed); L <- log(M); G <- ncol(L); gm <- t(vapply(1:3, function(g) colMeans(L[P_GROUP == g, , drop = FALSE]), numeric(G))); resid <- L - gm[P_GROUP, , drop = FALSE]
  idx <- lapply(1:3, function(g) which(P_GROUP == g)); ca <- integer(n); cs <- integer(n); gmx <- gm[P_GROUP, , drop = FALSE]; sig <- NULL
  for (b in seq_len(n)) {
    R <- resid
    for (ix in idx) for (cc in seq_len(G)) R[ix, cc] <- np_shuffle(rng, R[ix, cc])
    Pm <- exp(gmx + R); colnames(Pm) <- colnames(M); rs <- P_sig_sets(P_schemes(Pm))
    if (is.null(sig)) sig <- lapply(rs, function(s) integer(n))
    for (k in names(rs)) sig[[k]][b] <- sum(rs[[k]]$p[upper.tri(rs[[k]]$p)] < 0.05)
    v <- P_consensus_n(rs); ca[b] <- v[1]; cs[b] <- v[2]
  }
  .P_perm_cache[[key]] <- list(ca = ca, cs = cs, sig = sig)
}
P_rcrit <- function(df_, alpha) { t <- qt(1 - alpha / 2, df_); t / sqrt(df_ + t^2) }
.pairs_list <- function(m) lapply(seq_len(nrow(m)), function(i) list(m[i, 1], m[i, 2]))
.pkey <- function(m) paste(m[, 1], m[, 2], sep = "|")
## everything the pipeline reports for one gene set; numpy integers that the Python json.dump(default=str) wrote as strings are character here
P_analyze <- function(genes, nperm = 2000, seed = 2026, label = "") {
  M <- P_load(genes); sch <- P_schemes(M); res <- P_sig_sets(sch); cons <- P_consensus(res, use_sign = FALSE); cons_s <- P_consensus(res, use_sign = TRUE); cn <- P_counts(res); adj <- P_adjusted(sch)
  out <- list(label = label, genes = as.list(genes), n_genes = length(genes),
              counts = lapply(cn, function(c) list(pairs = c$pairs, sig = c$sig, strong = as.character(c$strong), very = as.character(c$very))), consensus = .pairs_list(cons), consensus_signcons = .pairs_list(cons_s))
  out$adj_sig_q05 <- lapply(adj, function(a) as.character(sum(a$q[upper.tri(a$q)] < 0.05))); out$adj_sig_p05 <- lapply(adj, function(a) as.character(sum(a$p[upper.tri(a$p)] < 0.05)))
  cg <- .P_common(adj); Qs <- .P_sub(adj, "q", cg); allq <- Reduce(`&`, lapply(Qs, function(q) q < 0.05)); allq[is.na(allq)] <- FALSE
  out$adj_all_four <- .pairs_list(.P_pairs_rowmajor(allq, cg))                       # Python: set order; compared as sets
  qlook <- function(a, pr) { i <- match(pr[1], a$genes); j <- match(pr[2], a$genes); a$q[min(i, j), max(i, j)] }
  nq <- function(pr) sum(vapply(adj, function(a) !is.na(qlook(a, pr)) && qlook(a, pr) < 0.05, TRUE))
  cnt_adj <- if (nrow(cons)) vapply(seq_len(nrow(cons)), function(i) nq(cons[i, ]), 1L) else integer(0)
  out$cons_adj_any <- .pairs_list(cons[cnt_adj > 0, , drop = FALSE]); out$cons_adj_counts <- setNames(lapply(cnt_adj, as.character), .pkey(cons))
  m_norm <- max(vapply(cn[names(cn) != "Raw"], function(c) c$pairs, 1)); m_raw <- cn$Raw$pairs
  out$rcrit_df8_q05_best <- list(raw = P_rcrit(8, 0.05 / m_raw), norm = P_rcrit(8, 0.05 / m_norm)); out$rcrit_df8_p05 <- P_rcrit(8, 0.05)
  pn <- P_permutation_null(M, nperm, seed); ca <- pn$ca; cs <- pn$cs
  summ <- function(a, obs) list(n = as.integer(nperm), mean = mean(a), sd = np_std(a), p_ge_obs = mean(a >= obs), observed = as.integer(obs), q05 = np_percentile(a, 5), q95 = np_percentile(a, 95))
  out$perm <- summ(ca, nrow(cons)); out$perm_signcons <- summ(cs, nrow(cons_s))
  P <- .P_sub(res, "p", cg <- .P_common(res)); Rm <- .P_sub(res, "r", cg); s1 <- sign(Rm[[1]]); allp <- Reduce(`&`, lapply(P, function(p) p <= 0.05)); allp[is.na(allp)] <- FALSE
  diff_sign <- !Reduce(`&`, lapply(Rm[-1], function(r) sign(r) == s1)); diff_sign[is.na(diff_sign)] <- FALSE
  out$sign_inconsistent <- .pairs_list(.P_pairs_rowmajor(allp & diff_sign, cg))      # Python: set order; compared as sets
  pget <- function(s, field, pr) { i <- match(pr[1], s$genes); j <- match(pr[2], s$genes); s[[field]][i, j] }
  out$cons_table <- setNames(lapply(seq_len(nrow(cons)), function(i) lapply(res, function(s) list(r = pget(s, "r", cons[i, ]), p = pget(s, "p", cons[i, ])))), .pkey(cons))
  out$adj_table <- setNames(lapply(seq_len(nrow(cons)), function(i) lapply(adj, function(a) list(r = pget(a, "r", cons[i, ]), q = qlook(a, cons[i, ])))), .pkey(cons))
  list(out = out, res = res, adj = adj, cons = cons, cons_s = cons_s)
}
## Python-style repr of a small named count list, for the printed summary
.py_counts <- function(l) paste0("{", paste(sprintf("'%s': %s", names(l), vapply(l, function(x) as.character(x), "")), collapse = ", "), "}")
P_main <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (spec in list(list("20 genes (current manuscript)", P_GENES20), list("17 genes (after exclusion)", P_GENES17))) {
    a <- P_analyze(spec[[2]], label = spec[[1]]); o <- a$out
    cat("=====", spec[[1]], "| genes", o$n_genes, "\n")
    for (k in names(o$counts)) cat("  ", k, .py_counts(o$counts[[k]]), "\n")
    cat("   pooled consensus (p<=0.05 in all 4 schemes, any sign):", length(o$consensus), "| sign-consistent:", length(o$consensus_signcons), "\n")
    for (nm in c("perm", "perm_signcons")) cat(sprintf("   permutation null (%s): mean %.1f SD %.1f, P(null>=%d)=%.3f, 5-95%%: %.0f-%.0f\n", nm, o[[nm]]$mean, o[[nm]]$sd, o[[nm]]$observed, o[[nm]]$p_ge_obs, o[[nm]]$q05, o[[nm]]$q95))
    cat("   treatment-adjusted q<0.05 per scheme:", .py_counts(o$adj_sig_q05), "| all four:", length(o$adj_all_four), "| consensus pairs sig in >=1 scheme:", length(o$cons_adj_any), "of", length(o$consensus), "\n")
    cat("   consensus pairs sig in >=3 schemes:", paste0("[", paste(sprintf("'%s'", names(o$cons_adj_counts)[as.integer(unlist(o$cons_adj_counts)) >= 3]), collapse = ", "), "]"), "\n")
    cat("   sign-inconsistent pairs (sig all 4 but sign differs):", length(o$sign_inconsistent), "\n")
    cat(sprintf("   |r| needed at df=8: p<0.05 %.3f ; best pair q<0.05 raw %.3f / norm %.3f\n", o$rcrit_df8_p05, o$rcrit_df8_q05_best$raw, o$rcrit_df8_q05_best$norm))
    write_text_bytes(py_json_typed(o, indent = NULL), file.path(out_dir, sprintf("pipeline_%dgenes.json", o$n_genes)))
  }
}
if (sys.nframe() == 0L) {
  HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
  source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R"))
  P_main(Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")))
}
