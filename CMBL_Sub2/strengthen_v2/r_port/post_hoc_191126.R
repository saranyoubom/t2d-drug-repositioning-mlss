## post_hoc_191126.R - R port of ars_review_191126/post_hoc/post_hoc_checks_191126.py, post_hoc_validation_191126.py and post_hoc_ins1_191126.py (19 Nov 2026). Two POST HOC checks of the benchmark
## (method in POSTHOC_METHODS_191126.md, fixed before any result):
##   check 1  share of the pooled covariance carried by group means, observed against its value when within-group coupling is removed (independent within-group permutation, selection of
##            pooled-significant pairs repeated in every permuted matrix);
##   check 2  the within-group coupling test T_r2 after removal of the first principal component of the residuals (the same removal in every permuted matrix; same permutations as the original test).
## Gates (before any summary is written): the run must reproduce the stored share_between and T_r2_p of dev_core_results.jsonl (WP2_DEV_CORE: default the R run results_R/wp2/compact/dev_core_results.jsonl,
## else the archived Python file; the archived file differs from an R run in the degenerate matrices with constant genes, see R_PORT_PLAN.md).
## Usage: Rscript post_hoc_191126.R dev | validation | ins1       (results into <POSTHOC_OUT>/results_<mode>; default results_R/post_hoc)
##   dev         development series (37), stored matrices (WP2_DEV_SETS: oracle/dev_sets_archive.json or the R extract), no download
##   validation  the 59 validation sets rebuilt from the cached GEO files (results_R/wp2_data; fetched by wp2_prefetch.R)
##   ins1        the INS-1 matrix of the manuscript (addendum POSTHOC_METHODS_ADDENDUM_191126.md), 10000 permutations; INS1_MATRIX overrides the matrix path
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- HERE
W2H <- file.path(R_PORT, "wp2_benchmark2"); source(file.path(W2H, "wp2_common.R")); source(file.path(R_PORT, "coexp_checks.R")); source(file.path(R_PORT, "cc_fast.R")); source(file.path(R_PORT, "bench_pipeline.R"))
source(file.path(W2H, "bench_data_io.R")); source(file.path(W2H, "bench2_core.R")); suppressMessages(library(readxl))
PH_ARC <- normalizePath(file.path(R_PORT, "..", "ars_review_191126", "post_hoc")); POUT <- Sys.getenv("POSTHOC_OUT", file.path(R_PORT, "results_R", "post_hoc")); NPERM <- 2000L
METHODS <- file.path(PH_ARC, "POSTHOC_METHODS_191126.md"); DES <- BC_DES
share_expected <- function(X, nperm, seed) {      # observed median share and the distribution of the same statistic without within-group coupling (cc_group_mean_share_expected, mean added)
  X <- as.matrix(X); obs <- cc_group_mean_share(X, DES); e <- cc_group_mean_share_expected(X, DES, nperm = nperm, seed = seed)
  list(obs_median = obs$median, obs_mean = obs$mean, obs_n_pairs = as.integer(obs$n_pairs), null_n_pairs_mean = e$null_n_pairs_mean, null_no_sig_frac = e$null_no_sig_frac, null_mean = e$null_mean, null_q025 = e$null_q025, null_q975 = e$null_q975,
       excess = if (is.na(obs$median)) NaN else obs$median - e$null_mean, above_q975 = if (is.na(obs$median)) NA else obs$median > e$null_q975)
}
tr2_after_pc1 <- function(Rb) {      # T_r2 of a stack (c, G, N) after subtracting the first singular component of every matrix; also the share of the residual variance in PC1
  d <- dim(Rb); cn <- d[1]; G <- d[2]; N <- d[3]; t2 <- numeric(cn); sh <- numeric(cn)
  for (b in seq_len(cn)) {
    M <- matrix(Rb[b, , ], G, N); s <- svd(M); R1 <- M - s$d[1] * s$u[, 1] %o% s$v[, 1]; r <- cc_corr(R1); r2 <- r * r
    t2[b] <- (sum(r2) - sum(diag(r2))) / (G * (G - 1)); sh[b] <- s$d[1]^2 / sum(s$d^2)
  }
  list(t = t2, share = sh)
}
coupling_pc1 <- function(X, nperm, seed) {
  X <- as.matrix(X); R <- DES$Q %*% log(X); G <- ncol(X); N <- nrow(X); st <- array(t(R), c(1, G, N))
  obs0 <- .coupling_stats(st, DES, NULL, FALSE)$T_r2[1]; o1 <- tr2_after_pc1(st); obs1 <- o1$t[1]; pc1_share <- o1$share[1]
  ch <- cc_perm_apply(R, DES, nperm, seed, NULL, fun = function(Rb) list(n0 = .coupling_stats(Rb, DES, NULL, FALSE)$T_r2, n1 = tr2_after_pc1(Rb)$t))
  n0 <- unlist(lapply(ch, function(x) x$n0)); n1 <- unlist(lapply(ch, function(x) x$n1)); p0 <- .pvals(n0, obs0)$p_ge_plus1; p1 <- .pvals(n1, obs1)$p_ge_plus1
  list(T_r2 = obs0, T_r2_p = p0, T_r2_pc1 = obs1, T_r2_pc1_p = p1, pc1_share = pc1_share, null_pc1_mean = mean(n1))
}
analyse <- function(matrices, label, stored = NULL) {      # matrices: list of list(gse, platform, tag, X); stored: named list "gse tag" -> flat
  rows <- list(); t0 <- Sys.time(); gate <- list(share_max_abs_diff = 0, tr2_p_mismatch = 0L, n_checked = 0L)
  for (m in matrices) {
    seed <- 2026 + as.integer(substring(m$gse, 4)); a <- share_expected(m$X, NPERM, seed); b <- coupling_pc1(m$X, NPERM, seed); rows[[length(rows) + 1]] <- c(list(gse = m$gse, platform = m$platform, set = m$tag), a, b)
    k <- paste(m$gse, m$tag)
    if (!is.null(stored) && !is.null(stored[[k]])) {
      f <- stored[[k]]; gate$n_checked <- gate$n_checked + 1L
      if (!(is.na(a$obs_median) && (is.null(f$share_between) || is.na(f$share_between)))) gate$share_max_abs_diff <- max(gate$share_max_abs_diff, abs(a$obs_median - f$share_between))
      if (abs(b$T_r2_p - f$T_r2_p) > 1e-12) gate$tr2_p_mismatch <- gate$tr2_p_mismatch + 1L
    }
    cat(sprintf("  %s %s %s: share %s, T_r2 P %s, after PC1 P %s (%.0f s)\n", m$gse, m$platform, m$tag, py_fmt_f(a$obs_median, 3), py_fmt_f(b$T_r2_p, 4), py_fmt_f(b$T_r2_pc1_p, 4), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  list(rows = rows, gate = gate)
}
qs <- function(v) if (length(v)) vapply(c(0.25, 0.5, 0.75), function(p) np_quantile(v, p), 1) else rep(NaN, 3)
summarize <- function(rows, label, gate, methods_sha) {
  S <- list(label = label, methods_sha256 = methods_sha, nperm = as.integer(NPERM), gate = gate, sets = list())
  for (tag in c("A", "B", "A+B")) {
    R <- Filter(function(r) tag == "A+B" || r$set == tag, rows); n <- length(R); defined <- Filter(function(r) !is.na(r$obs_median) && !is.na(r$null_mean), R)
    f0 <- as.integer(sum(vapply(R, function(r) r$T_r2_p < 0.05, TRUE))); f1 <- as.integer(sum(vapply(R, function(r) r$T_r2_pc1_p < 0.05, TRUE))); both <- as.integer(sum(vapply(R, function(r) r$T_r2_p < 0.05 && r$T_r2_pc1_p < 0.05, TRUE))); only0 <- f0 - both; only1 <- f1 - both
    S$sets[[tag]] <- list(n_matrices = as.integer(n), n_share_defined = as.integer(length(defined)), obs_share_quartiles = as.list(qs(vapply(defined, function(r) r$obs_median, 1))), null_share_quartiles = as.list(qs(vapply(defined, function(r) r$null_mean, 1))),
                          excess_quartiles = as.list(qs(vapply(defined, function(r) r$excess, 1))), n_above_null_q975 = as.integer(sum(vapply(defined, function(r) isTRUE(r$above_q975), TRUE))), pc1_share_quartiles = as.list(qs(vapply(R, function(r) r$pc1_share, 1))),
                          flagged_original = f0, flagged_original_ci = as.list(py_wilson(f0, n)), flagged_pc1 = f1, flagged_pc1_ci = as.list(py_wilson(f1, n)), paired = list(both = both, only_original = only0, only_after_removal = only1, neither = as.integer(n - f0 - only1)))
  }
  S
}
write_md_post <- function(S, path) {
  f3 <- function(x) py_fmt_f(unlist(x), 3); qf <- function(x) paste(f3(x), collapse = " / ")
  L <- c(sprintf("# Post hoc checks of the benchmark: results (%s)", S$label), "", sprintf("Method: `POSTHOC_METHODS_191126.md` (SHA-256 %s, self-recorded), %d permutations per matrix, seed 2026 + GSE number.", S$methods_sha256, S$nperm), "",
         sprintf("Gates: %d matrices compared with the stored results; largest difference in the observed share %s; %d mismatches of the stored T_r2 P.", S$gate$n_checked, py_fmt_g(S$gate$share_max_abs_diff, 3), S$gate$tr2_p_mismatch), "")
  for (tag in names(S$sets)) {
    v <- S$sets[[tag]]
    L <- c(L, sprintf("## Gene set %s (%d matrices)", tag, v$n_matrices), "",
           sprintf("- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (%d matrices with a defined value): %s", v$n_share_defined, qf(v$obs_share_quartiles)),
           sprintf("- Expected without within-group coupling (mean of the null medians), quartiles: %s", qf(v$null_share_quartiles)),
           sprintf("- Excess (observed minus expected), quartiles: %s; observed above the 97.5th null percentile in %d of %d", qf(v$excess_quartiles), v$n_above_null_q975, v$n_share_defined),
           sprintf("- Tr2 flagged (P < 0.05): original %d of %d (95%% CI %s to %s); after removal of PC1 %d of %d (%s to %s)", v$flagged_original, v$n_matrices, py_fmt_f(v$flagged_original_ci[[1]], 2), py_fmt_f(v$flagged_original_ci[[2]], 2), v$flagged_pc1, v$n_matrices, py_fmt_f(v$flagged_pc1_ci[[1]], 2), py_fmt_f(v$flagged_pc1_ci[[2]], 2)),
           sprintf("- Paired: flagged both times %d, only before removal %d, only after removal %d, neither %d", v$paired$both, v$paired$only_original, v$paired$only_after_removal, v$paired$neither),
           sprintf("- Share of the residual variance carried by PC1, quartiles: %s", qf(v$pc1_share_quartiles)), "")
  }
  write_md(L, path); L
}
stored_flat <- function(path) { s <- list(); for (j in read_jsonl(path)) s[[paste(j$gse, j$tag)]] <- j$flat; s }
write_rows <- function(rows, path) dictwriter_csv(rows, names(rows[[1]]), path, fmt = "repr")
mat12 <- function(rows) do.call(rbind, lapply(rows, function(r) as.numeric(unlist(r))))
main <- function(mode) {
  out <- file.path(POUT, paste0("results_", mode)); dir.create(out, showWarnings = FALSE, recursive = TRUE); msha <- sha256_file(METHODS); cat("methods file sha256", msha, "\n")
  fail <- function(g) { write_json_typed(g, file.path(out, "GATE_FAILED.json")); cat("GATE FAILED\n"); quit(status = 2) }
  if (mode == "dev") {
    sp <- Sys.getenv("WP2_DEV_SETS", file.path(R_PORT, "oracle", "dev_sets_archive.json")); Z <- fromJSON(sp, simplifyVector = FALSE); pp <- file.path(dirname(sp), "dev_prepared.json"); if (!file.exists(pp)) pp <- file.path(ARC, "compact", "dev_prepared.json")
    prep <- fromJSON(pp, simplifyVector = FALSE); sc <- Sys.getenv("WP2_DEV_CORE", file.path(R_PORT, "results_R", "wp2", "compact", "dev_core_results.jsonl")); if (!file.exists(sc)) sc <- file.path(ARC, "compact", "dev_core_results.jsonl"); stored <- stored_flat(sc)
    mats <- unlist(lapply(names(prep), function(g) lapply(c("A", "B"), function(tag) list(gse = g, platform = prep[[g]]$platform, tag = tag, X = mat12(Z[[sprintf("%s__%s", g, tag)]])))), recursive = FALSE); label <- "development series (37), stored matrices"
    only <- setdiff(strsplit(Sys.getenv("POSTHOC_ONLY", ""), ",", fixed = TRUE)[[1]], ""); if (length(only)) { mats <- Filter(function(m) m$gse %in% only, mats); label <- paste0(label, " [subset: ", paste(only, collapse = ","), "]") }      # POSTHOC_ONLY: a few series for a check
    r <- analyse(mats, label, stored); if (r$gate$share_max_abs_diff > 1e-9 || r$gate$tr2_p_mismatch != 0 || r$gate$n_checked != length(mats)) fail(r$gate)
    write_rows(r$rows, file.path(out, "per_matrix.csv")); S <- summarize(r$rows, label, r$gate, msha); write_json_typed(S, file.path(out, "summary.json")); cat(paste(write_md_post(S, file.path(out, "summary.md")), collapse = "\n"), "\n")
  } else if (mode == "ins1") {
    ADD <- file.path(PH_ARC, "POSTHOC_METHODS_ADDENDUM_191126.md"); G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r"); nperm <- 10000L; seed <- 2026
    mp <- Sys.getenv("INS1_MATRIX", ""); if (!nzchar(mp)) mp <- c(file.path(R_PORT, "..", "..", "..", "SBA", "Figures", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"), file.path(R_PORT, "..", "..", "R_17genes", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"))[c(file.exists(file.path(R_PORT, "..", "..", "..", "SBA", "Figures", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx")), TRUE)][1]
    cat("INS-1 matrix:", normalizePath(mp), "\n"); X <- as.matrix(as.data.frame(read_excel(mp))[, G17]); storage.mode(X) <- "double"; stopifnot(identical(dim(X), c(12L, 17L)))
    a <- share_expected(X, nperm, seed); b <- coupling_pc1(X, nperm, seed)
    gate <- list(share_ok = isTRUE(round(a$obs_median, 3) == 0.971), tr2_ok = isTRUE(round(b$T_r2, 3) == 0.133), p_ok = isTRUE(abs(b$T_r2_p - 0.061) <= 0.004)); cat("gates", paste(names(gate), unlist(gate), sep = " = ", collapse = ", "), "\n")
    res <- list(main_methods_sha256 = msha, addendum_sha256 = sha256_file(ADD), nperm = nperm, seed = as.integer(seed), gate = gate, share = a, coupling = b); if (!all(unlist(gate))) { write_json_typed(res, file.path(out, "GATE_FAILED.json")); quit(status = 2) }
    write_json_typed(res, file.path(out, "summary.json"))
    pyb <- function(x) if (isTRUE(x)) "True" else "False"; plus3 <- function(x) paste0(if (x < 0) "-" else "+", py_fmt_f(abs(x), 3))      # Python's True / False and '%+.3f'
    L <- c(sprintf("# Post hoc checks on the INS-1 matrix (addendum; %d permutations, seed %d)", nperm, as.integer(seed)), "",
           sprintf("Methods: `POSTHOC_METHODS_191126.md` (SHA-256 %s) and the addendum (SHA-256 %s), both self-recorded.", res$main_methods_sha256, res$addendum_sha256), "",
           sprintf("Gates: observed share 0.971, T_r2 0.133 and P within 0.004 of 0.061 (manuscript): {'share_ok': %s, 'tr2_ok': %s, 'p_ok': %s}.", pyb(gate$share_ok), pyb(gate$tr2_ok), pyb(gate$p_ok)), "",
           "## Check 1: share of the pooled covariance carried by group means",
           sprintf("- Observed (median over %d pooled-significant pairs): %s", a$obs_n_pairs, py_fmt_f(a$obs_median, 3)),
           sprintf("- Expected without within-group coupling (mean of the null medians, 10000 permutations): %s (2.5th to 97.5th percentile %s to %s); observed minus expected %s; above the 97.5th percentile: %s", py_fmt_f(a$null_mean, 3), py_fmt_f(a$null_q025, 3), py_fmt_f(a$null_q975, 3), plus3(a$excess), pyb(a$above_q975)),
           sprintf("- Pooled-significant pairs in the permuted matrices: mean %s (observed %d); permutations without a significant pair: %s", py_fmt_f(a$null_n_pairs_mean, 1), a$obs_n_pairs, py_fmt_f(a$null_no_sig_frac, 3)), "",
           "## Check 2: T_r2 before and after removal of the first principal component",
           sprintf("- Original: T_r2 %s, P = %s", py_fmt_f(b$T_r2, 4), py_fmt_f(b$T_r2_p, 4)), sprintf("- After PC1 removal: T_r2 %s, P = %s (null mean %s)", py_fmt_f(b$T_r2_pc1, 4), py_fmt_f(b$T_r2_pc1_p, 4), py_fmt_f(b$null_pc1_mean, 4)),
           sprintf("- Share of the residual variance carried by PC1: %s", py_fmt_f(b$pc1_share, 3)))
    write_md(L, file.path(out, "summary.md"), final_newline = TRUE); print(res)
  } else if (mode == "validation") {
    io_configure(file.path(DATA, "v2_data", "arrays_heldout", "matrices"), file.path(DATA, "bench_data", "platforms")); arr <- Filter(function(r) identical(r$status, "analysed"), read_jsonl(file.path(ARC, "heldout_arrays_results.jsonl"))); rna <- Filter(function(r) identical(r$status, "analysed"), read_jsonl(file.path(ARC, "heldout_rnaseq_results.jsonl")))
    fa <- read.csv(file.path(ARC, "frozen_heldout_arrays.csv"), stringsAsFactors = FALSE); fr <- read.csv(file.path(ARC, "frozen_heldout_rnaseq.csv"), stringsAsFactors = FALSE); mats <- list(); stored <- list(); arm_of <- list(); mism <- list()
    for (rec in arr) {
      r2 <- prepare_array(as.list(fa[fa$gse == rec$gse, ][1, ])); if (is.null(r2$sets)) stop(rec$gse, ": not analysable when rebuilt: ", r2$rec$status)
      for (tag in c("A", "B")) { if (!identical(as.character(r2$sets[[tag]]$genes), as.character(unlist(rec[[tag]]$genes)))) mism[[length(mism) + 1]] <- c(rec$gse, tag); mats[[length(mats) + 1]] <- list(gse = rec$gse, platform = rec$platform, tag = tag, X = r2$sets[[tag]]$X); stored[[paste(rec$gse, tag)]] <- rec[[tag]]$flat }
      arm_of[[rec$gse]] <- paste("arrays", rec$platform)
    }
    for (rec in rna) {
      row <- as.list(fr[fr$gse == rec$gse, ][1, ]); org <- row$organism; fc <- fetch_counts(DATA, rec$gse, org); if (is.null(fc$path)) stop(rec$gse, ": ", fc$why)
      cc <- read_counts_file(fc$path, strsplit(row$gsms, ";", fixed = TRUE)[[1]]); r0 <- list(gse = rec$gse, platform = org, organism = org, order = as.integer(row$order), how = row$how, groups = row$groups)
      r2 <- prepare_rnaseq(r0, cc$C, cc$ids, NULL, as.integer(substring(rec$gse, 4))); if (is.null(r2$sets)) stop(rec$gse, ": not analysable when rebuilt: ", r2$rec$status)
      for (tag in c("A", "B")) { if (!identical(as.character(r2$sets[[tag]]$genes), as.character(unlist(rec[[tag]]$genes)))) mism[[length(mism) + 1]] <- c(rec$gse, tag); mats[[length(mats) + 1]] <- list(gse = rec$gse, platform = paste("rnaseq", org), tag = tag, X = r2$sets[[tag]]$X); stored[[paste(rec$gse, tag)]] <- rec[[tag]]$flat }
      arm_of[[rec$gse]] <- paste("rnaseq", org)
    }
    cat(sprintf("rebuilt %d matrices; gene sets different from the stored ones: %d\n", length(mats), length(mism))); if (length(mism)) { write_json_typed(lapply(mism, as.list), file.path(out, "GATE_FAILED_genes.json")); stop("GATE FAILED: gene sets differ from the stored ones") }
    label <- "validation sets (59), rebuilt matrices"; r <- analyse(mats, label, stored); if (r$gate$share_max_abs_diff > 1e-9 || r$gate$tr2_p_mismatch != 0 || r$gate$n_checked != length(mats)) fail(r$gate)
    rows <- lapply(r$rows, function(x) { x$arm <- arm_of[[x$gse]]; x }); write_rows(rows, file.path(out, "per_matrix.csv")); S <- summarize(rows, label, r$gate, msha); write_json_typed(S, file.path(out, "summary.json")); write_md_post(S, file.path(out, "summary.md"))
    for (arm in py_sort(unique(unlist(arm_of)))) { sub <- Filter(function(x) x$arm == arm, rows); Sa <- summarize(sub, paste("validation sets,", arm), r$gate, msha); fn <- gsub(" ", "_", arm); write_md_post(Sa, file.path(out, sprintf("summary_%s.md", fn))); write_json_typed(Sa, file.path(out, sprintf("summary_%s.json", fn))) }
  } else stop("usage: Rscript post_hoc_191126.R dev | validation | ins1")
}
if (sys.nframe() == 0L) { a <- commandArgs(trailingOnly = TRUE); main(if (length(a)) a[1] else "dev") }
