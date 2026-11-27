## design_table.R - R port of wp2_benchmark2/design_table.py: protocol v2 section 3.7 (WP2 deliverable F): design guidance by simulation.
## G = 17 genes, 3 groups, n per group in 4, 6, 8, 12, 20, 30; one module of m = 6 genes with within-group correlation r in 0.3, 0.5, 0.7, 0.9 (cc_inject_factor), the other 11 genes independent;
## r = 0 is added as a type I error check (ADD-ON, same code). Data: independent log-normal values from cc_simulate (group effects N(0, 1), residual SD 0.3 on the log scale), seed 2026.
## Per cell R independent experiments (replicates here are independent experiments: wells from one preparation do not count). Outputs per cell:
##   power_T_r2       share of experiments with P < 0.05 for T_r2 (cc_coupling_test, 500 within-group permutations, P = (1 + k)/(500 + 1))
##   pair_power       mean over the 15 true module pairs of the share detected at BH q < 0.05 (partial correlation of log values after removing group means, residual df = N - 4, 136 pairs, BH over all pairs)
##   any_pair         share of experiments in which at least one true pair is detected
##   false_pairs      mean number of detected pairs outside the module (expected near 0)
##   r_best_pair      analytic |r| the single best of the 136 pairs needs for BH q < 0.05 (p <= 0.05/136) at the residual df; r_single is |r| for an uncorrected p < 0.05
## Usage: Rscript design_table.R run [worker=1 [nworkers=1 [R=300]]]  |  Rscript design_table.R summarize      (cells n = 30, 20, ... are written one file per n to <OUT>/compact/design_cells; a cell is skipped if its file exists)
if (!exists("HERE") || !file.exists(file.path(HERE, "wp2_common.R"))) { HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, "..")) }
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "coexp_checks.R")); source(file.path(R_PORT, "cc_fast.R"))
DIR <- file.path(OUT, "compact", "design_cells"); dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
NS <- c(4L, 6L, 8L, 12L, 20L, 30L); RS <- c(0.0, 0.3, 0.5, 0.7, 0.9); G <- 17L; M <- 6L; NPERM <- 500L
cell <- function(n, R) {
  rows <- list(); t0 <- Sys.time(); up <- upper.tri(matrix(0, G, G)); true <- (row(up)[up] <= M) & (col(up)[up] <= M)       # the 15 pairs inside the module (order of the pairs does not matter for BH and the means)
  for (ri in seq_along(RS)) {
    r <- RS[ri]
    for (rep in 0:(R - 1L)) {
      rng <- np_rng(c(2026, 37, n, ri - 1, rep)); sm <- cc_simulate(c(n, n, n), G, rng, effect_sd = 1.0, noise_sd = 0.3); X <- sm$X; codes <- sm$codes
      if (r > 0) X <- cc_inject_factor(X, codes, 1:M, r, rng)
      des <- cc_design(codes); ct <- cc_coupling_test(X, des, nperm = NPERM, seed = np_integers(rng, 0, 2^31 - 1, 1), eig = FALSE)
      Rm <- cc_residual_corr(X, des); pv <- .pval(Rm[up], des$df_resid); q <- .bh(pv); det <- q < 0.05
      rows[[length(rows) + 1]] <- list(n = n, r = r, rep = as.integer(rep), T_r2_p = ct$T_r2$p_ge_plus1, pair_power = sum(det[true]) / sum(true), any_pair = any(det[true]), false_pairs = as.integer(sum(det[!true])), n_detected = as.integer(sum(det)))
    }
  }
  saveRDS(rows, file.path(DIR, sprintf("n%d.rds", n))); cat(sprintf("done n = %d %d experiments %.1f s\n", n, length(rows), as.numeric(difftime(Sys.time(), t0, units = "secs")))); invisible(rows)
}
run <- function(worker, nw, R) {
  todo <- rev(sort(NS)); todo <- todo[!file.exists(file.path(DIR, sprintf("n%d.rds", todo)))]; todo <- todo[seq_along(todo) %% nw == (worker %% nw)]
  cat(sprintf("cells to run (n per group): %s, R = %d, worker %d of %d\n", paste(todo, collapse = " "), R, worker, nw)); for (n in todo) cell(n, R); cat("design simulation finished\n")
}
summarize <- function() {
  indir <- Sys.getenv("WP2_DESIGN_IN", DIR)           # WP2_DESIGN_IN = the archived Python folder compact/design_cells (n<N>.json) to test the summary alone; default: the RDS files of the R run
  rd <- function(n) { f <- file.path(indir, sprintf("n%d.rds", n)); if (file.exists(f)) rows_df(readRDS(f)) else rows_df(fromJSON(file.path(indir, sprintf("n%d.json", n)), simplifyVector = FALSE)) }
  R <- do.call(rbind, lapply(NS, rd)); out <- list()
  for (n in NS) for (r in RS) {
    g <- R[R$n == n & R$r == r, , drop = FALSE]; df <- 3L * n - 4L; N <- nrow(g); k <- as.integer(sum(g$T_r2_p < 0.05)); ci <- py_wilson(k, N); k2 <- as.integer(sum(g$any_pair)); ci2 <- py_wilson(k2, N)
    out[[length(out) + 1]] <- list(n_per_group = n, N = 3L * n, resid_df = df, r = r, replicates = as.integer(N), power_T_r2 = k / N, T_r2_lo = ci[1], T_r2_hi = ci[2], pair_power = np_mean(g$pair_power), any_pair = k2 / N, any_pair_lo = ci2[1], any_pair_hi = ci2[2],
                                   false_pairs = np_mean(as.numeric(g$false_pairs)), r_best_pair = cc_rcrit(df, 0.05 / 136), r_single = cc_rcrit(df, 0.05))
  }
  T <- rows_df(out); write_csv_crlf(T, file.path(OUT, "design_table.csv"))
  md <- c("# Design table (protocol 3.7): 17 genes, one module of 6 genes, 3 groups, simulated independent experiments", "",
          "Power of T_r2 (P < 0.05, 500 permutations) / share of true module pairs detected at BH q < 0.05 / experiments with at least one true pair detected. r = 0 rows give the type I error of T_r2 and the mean number of false pairs.", "",
          paste0("| n per group (N, residual df) | analytic abs(r) for the best pair | ", paste(sprintf("r = %s", py_fmt_f(RS, 1)), collapse = " | "), " |"), paste0("|---|---|", strrep("---|", length(RS))))
  for (n in NS) { g <- T[T$n_per_group == n, ]; g <- g[order(g$r), ]
    md <- c(md, paste0(sprintf("| %d (%d, %d) | %s | ", n, 3L * n, 3L * n - 4L, py_fmt_f(g$r_best_pair[1], 2)), paste(sprintf("%s / %s / %s", py_fmt_f(g$power_T_r2, 2), py_fmt_f(g$pair_power, 2), py_fmt_f(g$any_pair, 2)), collapse = " | "), " |")) }
  md <- c(md, "", sprintf("Replicates per cell: %d independent simulated experiments. Wells or replicates from one preparation are not independent experiments and do not count toward n.", T$replicates[1]), "")
  write_md(md, file.path(OUT, "design_table.md")); cat(paste(md, collapse = "\n"), "\n")
}
if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE); step <- if (length(a)) a[1] else "summarize"
  if (step == "run") run(if (length(a) >= 2) as.integer(a[2]) else 1L, if (length(a) >= 3) as.integer(a[3]) else 1L, if (length(a) >= 4) as.integer(a[4]) else 300L) else if (step == "summarize") summarize() else stop("run [worker nworkers R] | summarize")
}
