## calibrate_t2.R - R port of wp2_benchmark2/calibrate_t2.py: protocol v2 section 3.2: calibration of the within-group coupling tests on decoupled development data (WP2 deliverable B).
## For each of the 37 development series, gene set A (the 12 x 17 matrices saved by gate_part1 extract):
##   decouple = permute every gene's log values independently within each group (group means kept, coupling removed), then
##   (i)  no injection, 100 replicates per series: type I error of every test at P < 0.05;
##   (ii) inject ONE latent factor, centred within group, into m genes (drawn at random among the 14 non-denominator genes, as in the part 1 positive control), loading sqrt(r/(1 - r)) x the gene's pooled
##        within-group SD, so that two injected genes have within-group correlation r: r in 0.3, 0.5, 0.7, 0.9; m in 3, 6, 10; 20 replicates per series and condition.
## Tests (cc_coupling_test, 500 within-group permutations each, P = (1 + #{null >= obs})/(500 + 1)): T_r2 (primary), T_p, T_eig. The consensus-COUNT test of part 1 (cc_consensus_null, 500 permutations)
## is run at the REDUCED replicate number declared here: m = 6 (all four r) and no injection, the first 20 replicates of each series x condition (the protocol requires at least 5).
## Requirement fixed in the protocol: T_r2 type I error <= 0.06 and power >= 80% at r = 0.7, m = 6 (pooled over all series x replicates; the spread over series is reported as well).
## Seeds: numpy SeedSequence([2026, GSE number, condition, replicate]).
## Usage: Rscript calibrate_t2.R run [worker=1 of nworkers=1]  |  Rscript calibrate_t2.R summarize        Environment: WP2_DEV_SETS (JSON of the 12 x 17 matrices; default: <OUT>/compact/dev_sets.json made by gate_part1.R extract;
## the archived Python matrices are oracle/dev_sets_archive.json), WP2_DEV_CORE (dev_core_results.jsonl for the degenerate-series sensitivity; default the archived file).
if (!exists("HERE") || !file.exists(file.path(HERE, "wp2_common.R"))) { HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, "..")) }   # (HERE can be preset when the file is sourced)
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "coexp_checks.R")); source(file.path(R_PORT, "cc_fast.R")); source(file.path(HERE, "bench2_core.R")); source(file.path(HERE, "endpoints_lib.R"))
SERDIR <- file.path(OUT, "compact", "calib_series"); dir.create(SERDIR, showWarnings = FALSE, recursive = TRUE)
RS <- c(0.3, 0.5, 0.7, 0.9); MS <- c(3L, 6L, 10L); NREP_INJ <- 20L; NREP_NULL <- 100L; NREP_COUNT <- 20L; NPERM <- 500L
CONDITIONS <- c(list(list(0.0, 0L, NREP_NULL)), unlist(lapply(RS, function(r) lapply(MS, function(m) list(r, m, NREP_INJ))), recursive = FALSE))
dev_sets_path <- Sys.getenv("WP2_DEV_SETS", file.path(OUT, "compact", "dev_sets.json"))
load_dev_sets <- function() {
  Z <- fromJSON(dev_sets_path, simplifyVector = FALSE); pp <- file.path(dirname(dev_sets_path), "dev_prepared.json"); if (!file.exists(pp)) pp <- file.path(ARC, "compact", "dev_prepared.json")
  list(Z = Z, prep = fromJSON(pp, simplifyVector = FALSE))
}
mat12 <- function(rows) do.call(rbind, lapply(rows, function(r) as.numeric(unlist(r))))
one_series <- function(gse, X0) {
  n <- as.integer(substring(gse, 4)); rows <- list(); t0 <- Sys.time()
  for (ci in seq_along(CONDITIONS)) {
    r <- CONDITIONS[[ci]][[1]]; m <- CONDITIONS[[ci]][[2]]; nrep <- CONDITIONS[[ci]][[3]]
    for (rep in 0:(nrep - 1L)) {
      rng <- np_rng(c(2026, n, ci - 1, rep)); Xd <- cc_decouple(X0, BC_DES, rng)
      if (m > 0) { idx <- sort(np_choice_noreplace_idx(rng, 14, m) + 3); Xi <- cc_inject_factor(Xd, BC_DES, idx + 1, r, rng) } else { idx <- integer(0); Xi <- Xd }      # gene indices 3..16 (0-based as in Python)
      sd <- np_integers(rng, 0, 2^31 - 1, 1); ct <- cc_coupling_test(Xi, BC_DES, nperm = NPERM, seed = sd, eig = TRUE)
      row <- list(gse = gse, r = r, m = m, rep = as.integer(rep), injected = paste(idx, collapse = " "), T_r2_obs = ct$T_r2$observed, T_r2_p = ct$T_r2$p_ge_plus1, T_p_obs = ct$T_p$observed, T_p_p = ct$T_p$p_ge_plus1,
                  T_eig_obs = ct$T_eig$observed, T_eig_p = ct$T_eig$p_ge_plus1, count_obs = NA_real_, count_p = NA_real_, count_p_plus1 = NA_real_)
      if ((m == 0 || m == 6) && rep < NREP_COUNT) { cn <- cc_consensus_null(Xi, BC_DES, BC_DEN, nperm = NPERM, seed = sd); row$count_obs <- as.numeric(cn$observed); row$count_p <- cn$p_ge; row$count_p_plus1 <- cn$p_ge_plus1 }
      rows[[length(rows) + 1]] <- row
    }
  }
  saveRDS(rows, file.path(SERDIR, paste0(gse, ".rds"))); cat(sprintf("done %s %d rows %.1f s\n", gse, length(rows), as.numeric(difftime(Sys.time(), t0, units = "secs")))); invisible(rows)
}
run <- function(worker, nw) {
  D <- load_dev_sets(); todo <- names(D$prep); only <- setdiff(strsplit(Sys.getenv("WP2_CALIB_ONLY", ""), ",", fixed = TRUE)[[1]], ""); if (length(only)) todo <- c(only[only %in% todo], todo[!(todo %in% only)])      # WP2_CALIB_ONLY: these series first (e.g. the degenerate ones)
  todo <- todo[!file.exists(file.path(SERDIR, paste0(todo, ".rds")))]; todo <- todo[seq_along(todo) %% nw == (worker %% nw)]
  cat(sprintf("series to run: %d (worker %d of %d)\n", length(todo), worker, nw))
  for (gse in todo) one_series(gse, mat12(D$Z[[paste0(gse, "__A")]]))
  cat("calibration run finished\n")
}
wilson <- function(k, n) py_wilson(k, n)
summarize <- function() {
  D <- load_dev_sets(); prep <- D$prep; frames <- list()
  indir <- Sys.getenv("WP2_CALIB_IN", SERDIR)         # WP2_CALIB_IN = the archived Python folder compact/calib_series (JSON per series) to test the summary alone; default: the RDS files of the R run
  for (g in names(prep)) {
    f <- file.path(indir, paste0(g, ".rds")); fj <- file.path(indir, paste0(g, ".json"))
    if (file.exists(f)) frames[[length(frames) + 1]] <- rows_df(readRDS(f)) else if (file.exists(fj)) frames[[length(frames) + 1]] <- rows_df(lapply(fromJSON(fj, simplifyVector = FALSE), nullNA)) else stop("missing series file ", f)
  }
  R <- do.call(rbind, frames); rownames(R) <- NULL; R$platform <- vapply(R$gse, function(g) prep[[g]]$platform, ""); write_csv_crlf(R, file.path(OUT, "calibration_replicates.csv"))
  stats_ <- list(c("T_r2", "T_r2_p"), c("T_p", "T_p_p"), c("T_eig", "T_eig_p"), c("count", "count_p_plus1")); per <- list(); tab <- list()
  keys <- unique(R[, c("r", "m")]); keys <- keys[order(keys$r, keys$m), ]
  for (i in seq_len(nrow(keys))) {
    r <- keys$r[i]; m <- keys$m[i]; g <- R[R$r == r & R$m == m, , drop = FALSE]
    for (st in stats_) {
      sname <- st[1]; col <- st[2]; gg <- if (sname == "count") g[!is.na(g[[col]]), , drop = FALSE] else g
      if (!nrow(gg)) next
      for (gse in py_sort(unique(gg$gse))) { h <- gg[gg$gse == gse, , drop = FALSE]; per[[length(per) + 1]] <- list(gse = gse, platform = prep[[gse]]$platform, r = r, m = m, statistic = sname, n_rep = as.integer(nrow(h)), n_reject = as.integer(sum(h[[col]] < 0.05))) }
      k <- as.integer(sum(gg[[col]] < 0.05)); n <- as.integer(nrow(gg)); ci <- wilson(k, n)
      ser <- vapply(py_sort(unique(gg$gse)), function(gse) { x <- gg[[col]][gg$gse == gse]; sum(x < 0.05) / length(x) }, 1)
      tab[[length(tab) + 1]] <- list(r = r, m = m, condition = if (m == 0) "no injection (type I error)" else "injected", statistic = sname, n_series = as.integer(length(unique(gg$gse))), n_tests = n, n_reject = k, rate = k / n, ci_lo = ci[1], ci_hi = ci[2],
                                     series_min = min(ser), series_q1 = np_quantile(ser, 0.25), series_median = median(ser), series_q3 = np_quantile(ser, 0.75), series_max = max(ser),
                                     n_series_rate_ge_0_80 = as.integer(sum(ser >= 0.8)), n_series_rate_gt_0_06 = as.integer(sum(ser > 0.06)))
    }
  }
  write_csv_crlf(rows_df(per), file.path(OUT, "calibration_per_condition.csv")); T <- rows_df(tab); T <- T[order(T$m, T$r, T$statistic, method = "radix"), ]; rownames(T) <- NULL; write_csv_crlf(T, file.path(OUT, "calibration_tables.csv"))
  by_pl <- list()
  for (pl in py_sort(unique(R$platform))) { g <- R[R$platform == pl, ]; kk <- unique(g[, c("r", "m")]); kk <- kk[order(kk$r, kk$m), ]
    for (i in seq_len(nrow(kk))) { r <- kk$r[i]; m <- kk$m[i]; if ((r == 0.0 && m == 0) || (r == 0.7 && m == 6)) { h <- g[g$r == r & g$m == m, ]
      by_pl[[length(by_pl) + 1]] <- list(platform = pl, r = r, m = m, n_tests = as.integer(nrow(h)), T_r2_rate = sum(h$T_r2_p < 0.05) / nrow(h), T_p_rate = sum(h$T_p_p < 0.05) / nrow(h), T_eig_rate = sum(h$T_eig_p < 0.05) / nrow(h)) } } }
  write_csv_crlf(rows_df(by_pl), file.path(OUT, "calibration_by_platform.csv"))
  # ---- verdict (rules fixed in the protocol and in the header above)
  pick <- function(m, r = NULL, s) { d <- T[T$m == m & T$statistic == s & (if (is.null(r)) TRUE else T$r == r), , drop = FALSE]; if (nrow(d)) d[1, , drop = FALSE] else NULL }
  t1 <- pick(0, NULL, "T_r2"); pw <- pick(6, 0.7, "T_r2"); ok1 <- t1$rate <= 0.06; ok2 <- pw$rate >= 0.80; f4 <- function(x) py_fmt_f(x, 4); f2 <- function(x) py_fmt_f(x, 2)
  V <- c(sprintf("PROTOCOL 3.2 REQUIREMENT (T_r2 type I error <= 0.06 and power >= 80%% at r = 0.7, m = 6; pooled over %d series)", t1$n_series), "",
         sprintf("T_r2 type I error: %d/%d = %s (Wilson 95%% CI %s to %s) -> %s", t1$n_reject, t1$n_tests, f4(t1$rate), f4(t1$ci_lo), f4(t1$ci_hi), if (ok1) "PASS" else "FAIL"),
         sprintf("T_r2 power at r = 0.7, m = 6: %d/%d = %s (Wilson 95%% CI %s to %s) -> %s; per series: min %s, median %s, %d of %d series at or above 0.80", pw$n_reject, pw$n_tests, f4(pw$rate), f4(pw$ci_lo), f4(pw$ci_hi),
                 if (ok2) "PASS" else "FAIL", f2(pw$series_min), f2(pw$series_median), pw$n_series_rate_ge_0_80, pw$n_series),
         sprintf("OVERALL: %s", if (ok1 && ok2) "PASS" else "FAIL"), "")
  cand <- list(); for (s in c("T_r2", "T_p", "T_eig", "count")) { a <- pick(0, NULL, s); b <- pick(6, 0.7, s); if (!is.null(a) && !is.null(b)) cand[[length(cand) + 1]] <- list(s, a$rate, b$rate) }
  f3 <- function(x) py_fmt_f(x, 3)
  V <- c(V, paste0("All pre-listed statistics (type I error, power at r = 0.7, m = 6): ", paste(vapply(cand, function(c) sprintf("%s %s, %s", c[[1]], f3(c[[2]]), f3(c[[3]])), ""), collapse = "; ")))
  good <- Filter(function(c) c[[2]] <= 0.06 && c[[3]] >= 0.80, cand)
  if (!(ok1 && ok2)) {
    q <- Filter(function(c) c[[2]] <= 0.06, cand); best <- if (length(q)) q[[which.max(vapply(q, function(c) c[[3]], 1))]] else NULL
    V <- c(V, paste0("Secondary statistic named by the fixed rule (highest power at r = 0.7, m = 6 among statistics with type I error <= 0.06): ", if (!is.null(best)) best[[1]] else "none qualifies"))
  } else V <- c(V, paste0("Pre-listed statistics meeting both requirements: ", if (length(good)) paste(vapply(good, function(c) c[[1]], ""), collapse = ", ") else "none"))
  # ---- post hoc sensitivity (protocol section 8, WP2 entry of 5 Nov): without the series whose gene set A has genes with no within-group variation
  core <- read_jsonl(Sys.getenv("WP2_DEV_CORE", file.path(ARC, "compact", "dev_core_results.jsonl")))
  deg <- py_sort(unique(unlist(lapply(core, function(j) if (j$tag == "A" && (if (is.null(j$flat$n_zero_var_genes)) 0 else j$flat$n_zero_var_genes) > 0) j$gse))))
  Sx <- R[!(R$gse %in% deg), , drop = FALSE]; sens <- list()
  V <- c(V, "", sprintf("POST HOC SENSITIVITY (not the pre-specified analysis): without the %d series whose gene set A has genes with no within-group variation (%s); %d series remain", length(deg), paste(deg, collapse = ", "), length(unique(Sx$gse))))
  for (st in stats_) {
    sname <- st[1]; col <- st[2]; a <- Sx[Sx$m == 0, , drop = FALSE]; b <- Sx[Sx$m == 6 & Sx$r == 0.7, , drop = FALSE]
    if (sname == "count") { a <- a[!is.na(a[[col]]), , drop = FALSE]; b <- b[!is.na(b[[col]]), , drop = FALSE] }
    ka <- as.integer(sum(a[[col]] < 0.05)); na <- as.integer(nrow(a)); kb <- as.integer(sum(b[[col]] < 0.05)); nb <- as.integer(nrow(b)); wa <- wilson(ka, na); wb <- wilson(kb, nb)
    sens[[length(sens) + 1]] <- list(statistic = sname, type1_k = ka, type1_n = na, type1 = ka / na, type1_lo = wa[1], type1_hi = wa[2], power_k = kb, power_n = nb, power = kb / nb, power_lo = wb[1], power_hi = wb[2])
    V <- c(V, sprintf("  %-6s type I error %d/%d = %s (%s to %s); power at r = 0.7, m = 6: %d/%d = %s (%s to %s)", sname, ka, na, f4(ka / na), f4(wa[1]), f4(wa[2]), kb, nb, f4(kb / nb), f4(wb[1]), f4(wb[2])))
  }
  write_csv_crlf(rows_df(sens), file.path(OUT, "calibration_sensitivity_no_degenerate.csv")); write_md(V, file.path(OUT, "calibration_verdict.txt"), final_newline = TRUE); cat(paste(V, collapse = "\n"), "\n")
}
if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE); step <- if (length(a)) a[1] else "summarize"
  if (step == "run") run(if (length(a) >= 2) as.integer(a[2]) else 1L, if (length(a) >= 3) as.integer(a[3]) else 1L) else if (step == "summarize") summarize() else stop("run [worker nworkers] | summarize")
}
