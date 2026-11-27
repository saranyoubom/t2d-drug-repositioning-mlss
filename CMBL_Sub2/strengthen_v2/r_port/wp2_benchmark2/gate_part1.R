## gate_part1.R - R port of wp2_benchmark2/gate_part1.py: WP2 deliverable A gate: reload the development series from the downloaded part 1 data and reproduce the part 1 per-dataset numbers before
## anything new is done. Steps (Rscript gate_part1.R <step>):
##   extract      re-run the part 1 loop (frozen order, duplicate rule, 20 analysable per platform) with bench_data_io.R on <data>/bench_data/{matrices,platforms}; compare status, DE counts, scale and the two
##                17-gene draws with results_per_dataset.jsonl; save the 12 x 17 matrices (<OUT>/compact/dev_sets.json, dev_prepared.json, extract_check.json)
##   merge        puts the archived matrices of the series that GEO no longer serves (GSE31851) into dev_sets.json / dev_prepared.json of an extract run that did not (extract does it itself now)
##   compare      B_analyse (bench_pipeline.R, the part 1 function, 2000 permutations, seed 2026 + GSE number: the part 1 call) on gene sets A and B; compare observed consensus count, null mean, null SD, P with
##                per_dataset_results.csv (exact match required); resumes from <OUT>/compact/gate_bench_pipeline.jsonl; optional arguments: worker nworkers
##   consistency  coexp_checks (bench2_core.R analyse_matrix) on the same 74 matrices: observed values must be identical to part 1, null means and P agree within Monte Carlo error (not a gate);
##                results cached in <OUT>/compact/dev_core_results.jsonl and used by apply_dev.R; WP2_DEV_SETS can point to the archived matrices (oracle/dev_sets_archive.json)
##   report       writes gate_part1_reproduction.txt from the three result files
## Downloaded data are untrusted: read as numbers only.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "coexp_checks.R")); source(file.path(R_PORT, "cc_fast.R")); source(file.path(R_PORT, "bench_pipeline.R")); source(file.path(HERE, "bench_data_io.R")); source(file.path(HERE, "bench2_core.R"))
COMPACT <- file.path(OUT, "compact"); F_EXTRACT <- file.path(COMPACT, "extract_check.json"); F_SETS <- file.path(COMPACT, "dev_sets.json"); F_PREP <- file.path(COMPACT, "dev_prepared.json")
F_BENCH <- file.path(COMPACT, "gate_bench_pipeline.jsonl"); F_CORE <- file.path(COMPACT, "dev_core_results.jsonl"); F_REPORT <- file.path(OUT, "gate_part1_reproduction.txt")
part1_records <- function() read_jsonl(file.path(P1, "results_per_dataset.jsonl"))
mat_rows <- function(M) lapply(seq_len(nrow(M)), function(i) as.numeric(M[i, ]))        # 12 x 17 matrix -> list of rows (a list of 17 numbers each)
## series that GEO no longer serves (GSE31851, HTTP 403 since 2026): the archived 12 x 17 matrices (r_port/oracle/dev_sets_archive.json) and the archived prepared record stand in, in the frozen order of the rows,
## so that dev_sets.json holds all 74 matrices and the later steps (calibration, post hoc, Fig. 10 inputs) pool the same 37 series as the Python analysis; the series are listed in extract_check.json
fill_archived <- function(arrays, prep, unavailable, order) {
  if (!length(unavailable)) return(list(arrays = arrays, prep = prep))
  AZ <- fromJSON(file.path(R_PORT, "oracle", "dev_sets_archive.json"), simplifyVector = FALSE); AP <- fromJSON(file.path(ARC, "compact", "dev_prepared.json"), simplifyVector = FALSE)
  for (g in unavailable) { for (tag in c("A", "B")) arrays[[sprintf("%s__%s", g, tag)]] <- AZ[[sprintf("%s__%s", g, tag)]]; prep[[g]] <- AP[[g]] }
  have <- unique(sub("__[AB]$", "", names(arrays))); ord <- order[order %in% have]
  list(arrays = arrays[unlist(lapply(ord, function(g) paste0(g, c("__A", "__B"))))], prep = prep[ord])
}
merge_archived <- function() {      # step merge: apply fill_archived to the outputs of an extract run that did not do it
  ex <- fromJSON(F_EXTRACT, simplifyVector = FALSE); un <- as.character(unlist(ex$series_not_readable_from_GEO_archived_record_used))
  fa <- fill_archived(fromJSON(F_SETS, simplifyVector = FALSE), fromJSON(F_PREP, simplifyVector = FALSE), un, vapply(ex$checks, function(c) c$gse, ""))
  write_json_typed(fa$arrays, F_SETS, indent = NULL); write_json_typed(fa$prep, F_PREP); cat(sprintf("merged archived matrices of %s: %d matrices, %d series\n", paste(un, collapse = ", "), length(fa$arrays), length(fa$prep)))
}
extract <- function() {
  io_configure(file.path(DATA, "bench_data", "matrices"), file.path(DATA, "bench_data", "platforms"))
  rows <- read.csv(file.path(P1, "frozen_datasets.csv"), stringsAsFactors = FALSE, check.names = FALSE); recs <- part1_records(); ref <- setNames(recs, vapply(recs, function(r) r$gse, ""))
  done <- c(GPL570 = 0L, GPL1261 = 0L); used <- list(); checks <- list(); arrays <- list(); prep <- list(); unavailable <- character(0)
  for (i in seq_len(nrow(rows))) {
    row <- as.list(rows[i, ]); gpl <- row$platform; if (done[[gpl]] >= 20L) next
    t0 <- Sys.time(); gs <- strsplit(row$gsms, ";", fixed = TRUE)[[1]]; clash <- py_sort(unique(unlist(used[gs[gs %in% names(used)]]))); sets <- NULL
    if (length(clash)) rec <- list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = paste0("not analysable: duplicate of ", clash[1]))
    else {
      r <- tryCatch(prepare_array(row), error = function(e) {
        r0 <- ref[[row$gse]]; msg <- conditionMessage(e)
        if (grepl("HTTP 403", msg, fixed = TRUE) && !is.null(r0) && identical(r0$status, "analysed")) {   # GEO answers 403 for the whole directory of GSE31851 since 2026 (checked 25 Nov 2026): the series cannot be re-read; the archived record stands in
          unavailable <<- c(unavailable, row$gse)
          list(rec = c(r0[c("gse", "platform", "order", "how", "groups", "n_genes_mapped", "scale", "q99", "n_expressed", "n_de")], list(status = "prepared")), sets = NULL)
        } else list(rec = list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = sprintf("not analysable: %s %s", class(e)[1], substr(msg, 1, 120))), sets = NULL)
      })
      rec <- r$rec; sets <- r$sets
    }
    if (rec$status == "prepared") rec$status <- "analysed"
    if (rec$status %in% c("analysed", "analysable (dry run)")) { done[[gpl]] <- done[[gpl]] + 1L; for (g in gs) used[[g]] <- row$gse }
    r0 <- ref[[row$gse]]; ok <- list()
    if (is.null(r0)) ok$in_part1_log <- FALSE
    else {
      eqn <- function(a, b) !is.null(a) && !is.null(b) && !is.na(a) && !is.na(b) && a == b
      ok$status <- identical(rec$status, r0$status); ok$n_de <- identical(rec$n_de, r0$n_de) || eqn(rec$n_de, r0$n_de); ok$n_expressed <- eqn(rec$n_expressed, r0$n_expressed) || (is.null(rec$n_expressed) && is.null(r0$n_expressed))
      ok$n_genes_mapped <- eqn(rec$n_genes_mapped, r0$n_genes_mapped) || (is.null(rec$n_genes_mapped) && is.null(r0$n_genes_mapped)); ok$scale <- identical(rec$scale, r0$scale)
      ok$q99 <- if (!is.null(r0$q99)) isTRUE(abs(if (is.null(rec$q99)) 0 else rec$q99 - r0$q99) < 1e-9) else is.null(rec$q99)
      if (!is.null(sets)) for (tag in c("A", "B")) { ok[[paste0("genes_", tag)]] <- identical(as.character(sets[[tag]]$genes), as.character(unlist(r0[[tag]]$genes))); ok[[paste0("n_de_in_set_", tag)]] <- eqn(sets[[tag]]$n_de_in_set, r0[[tag]]$n_de_in_set) }
    }
    checks[[length(checks) + 1]] <- list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = rec$status, n_de = rec$n_de, ok = ok, all_ok = all(unlist(ok)), seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    if (!is.null(sets)) {
      for (tag in c("A", "B")) arrays[[sprintf("%s__%s", row$gse, tag)]] <- mat_rows(sets[[tag]]$X)
      prep[[row$gse]] <- list(platform = gpl, order = as.integer(row$order), n_de = rec$n_de, genes_A = as.list(sets$A$genes), genes_B = as.list(sets$B$genes), n_de_in_set_A = sets$A$n_de_in_set, n_de_in_set_B = sets$B$n_de_in_set, scale = rec$scale, groups = rec$groups)
    }
    cat(sprintf("%s %s (#%s): %s | DE %s | all matches %s | %ss\n", gpl, row$gse, row$order, rec$status, if (is.null(rec$n_de)) "-" else rec$n_de, if (checks[[length(checks)]]$all_ok) "True" else "False", checks[[length(checks)]]$seconds))
  }
  summary <- list(n_rows_processed = as.integer(length(checks)), n_in_part1_log = as.integer(sum(vapply(checks, function(c) is.null(c$ok$in_part1_log), TRUE))), n_part1_records = as.integer(length(ref)),
                  n_all_ok = as.integer(sum(vapply(checks, function(c) c$all_ok, TRUE))), n_analysed = as.integer(sum(vapply(checks, function(c) c$status == "analysed", TRUE))),
                  analysed_per_platform = as.list(done), series_not_readable_from_GEO_archived_record_used = as.list(unavailable), checks = checks)
  fa <- fill_archived(arrays, prep, unavailable, vapply(checks, function(c) c$gse, ""))
  write_json_typed(summary, F_EXTRACT); write_json_typed(fa$arrays, F_SETS, indent = NULL); write_json_typed(fa$prep, F_PREP)
  cat("extract finished:", paste(sprintf("%s=%s", names(summary)[names(summary) != "checks"], vapply(summary[names(summary) != "checks"], function(v) paste(unlist(v), collapse = "/"), "")), collapse = ", "), "\n")
}
load_sets <- function() {
  sp <- Sys.getenv("WP2_DEV_SETS", F_SETS); pp <- file.path(dirname(sp), "dev_prepared.json"); if (!file.exists(pp)) pp <- file.path(ARC, "compact", "dev_prepared.json")
  list(Z = fromJSON(sp, simplifyVector = FALSE), prep = fromJSON(pp, simplifyVector = FALSE))
}
mat12 <- function(rows) do.call(rbind, lapply(rows, function(r) as.numeric(unlist(r))))
compare <- function(worker = 1L, nw = 1L) {
  D <- load_sets(); tasks <- list(); have <- character(0)
  if (file.exists(F_BENCH)) have <- vapply(read_jsonl(F_BENCH), function(j) paste(j$gse, j$tag), "")
  for (gse in names(D$prep)) for (tag in c("A", "B")) if (!(paste(gse, tag) %in% have)) tasks[[length(tasks) + 1]] <- list(gse = gse, tag = tag)
  tasks <- tasks[seq_along(tasks) %% nw == (worker %% nw)]; cat(sprintf("tasks to run: %d already done: %d worker %d of %d\n", length(tasks), length(have), worker, nw))
  fn <- if (nw > 1) sub("[.]jsonl$", sprintf("_w%d.jsonl", worker), F_BENCH) else F_BENCH
  for (t in tasks) {
    t0 <- Sys.time(); X <- mat12(D$Z[[sprintf("%s__%s", t$gse, t$tag)]]); r <- B_analyse(X, 1:3, nperm = 2000L, seed = 2026 + as.integer(substring(t$gse, 4))); sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(py_json_typed(list(gse = t$gse, tag = t$tag, res = r, sec = sec), indent = NULL), "\n", file = fn, sep = "", append = TRUE)
    cat(sprintf("done %s %s obs %d P %s %ss\n", t$gse, t$tag, r$observed, py_fmt_g(r$perm_p_plain, 6), sec))
  }
  cat("compare finished\n")
}
consistency <- function() {
  D <- load_sets(); have <- character(0); if (file.exists(F_CORE)) have <- vapply(read_jsonl(F_CORE), function(j) paste(j$gse, j$tag), "")
  for (gse in names(D$prep)) for (tag in c("A", "B")) {
    if (paste(gse, tag) %in% have) next
    t0 <- Sys.time(); X <- mat12(D$Z[[sprintf("%s__%s", gse, tag)]]); res <- analyse_matrix(X, unlist(D$prep[[gse]][[paste0("genes_", tag)]]), as.integer(substring(gse, 4)))
    cat(py_json_typed(list(gse = gse, platform = D$prep[[gse]]$platform, tag = tag, flat = res$flat, detail = res$detail), indent = NULL), "\n", file = F_CORE, sep = "", append = TRUE)
    cat(sprintf("core %s %s obs %d T_r2 P %s %.1fs\n", gse, tag, res$flat$observed, py_fmt_f(res$flat$T_r2_p, 3), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  cat("consistency finished\n")
}
report <- function() {
  ex <- fromJSON(F_EXTRACT, simplifyVector = FALSE); T <- read.csv(file.path(P1, "per_dataset_results.csv"), stringsAsFactors = FALSE, check.names = FALSE); recs <- part1_records()
  ref <- list(); for (r in recs) if (identical(r$status, "analysed")) for (t in c("A", "B")) ref[[paste(r$gse, t)]] <- r[[t]]
  new <- list(); for (j in read_jsonl(F_BENCH)) new[[paste(j$gse, j$tag)]] <- j$res
  L <- c("GATE A: reproduction of the part 1 per-dataset numbers on the development series", "", sprintf("Run: %s; R %s", format(Sys.time(), "%Y-%m-%d %H:%M"), getRversion()), "",
         "1. Reload of the cached series (bench_data), same loop as run_benchmark (frozen order, duplicate rule, at most 20 analysable per platform)",
         sprintf("   rows processed %d (part 1 log: %d records); status, DE counts, expressed counts, scale, 99th percentile, gene lists of sets A and B and DE-in-set counts identical for %d of %d rows; analysed per platform {%s}",
                 ex$n_rows_processed, ex$n_part1_records, ex$n_all_ok, ex$n_rows_processed, paste(sprintf("'%s': %d", names(ex$analysed_per_platform), unlist(ex$analysed_per_platform)), collapse = ", ")))
  bad <- Filter(function(c) !c$all_ok, ex$checks)
  L <- c(L, paste0("   mismatching rows: ", if (!length(bad)) "none" else paste(vapply(bad, function(c) sprintf("%s [%s]", c$gse, paste(names(c$ok)[!unlist(c$ok)], collapse = ", ")), ""), collapse = "; ")))
  L <- c(L, "", "2. B_analyse (part 1 function; 2000 permutations; seed 2026 + GSE number) against per_dataset_results.csv, gene sets A and B", "   criterion: observed count identical; null mean within 1e-9; permutation P (plain fraction) identical", "",
         "series        set  observed(p1/now)  null_mean(p1/now)         P(p1/now)       null_sd(p1/now)       match")
  keys <- names(ref); keys <- keys[order(sub(" .*", "", keys), sub(".* ", "", keys), method = "radix")]; n_ok <- 0L; n_tot <- 0L; maxd <- 0
  for (k in keys) {
    gse <- sub(" .*", "", k); tag <- sub(".* ", "", k); n_tot <- n_tot + 1L; t <- T[T$gse == gse & T$set == tag, , drop = FALSE]
    if (is.null(new[[k]]) || nrow(t) != 1) { L <- c(L, sprintf("%-12s %s  MISSING (new %s, table rows %d)", gse, tag, if (!is.null(new[[k]])) "True" else "False", nrow(t))); next }
    r <- new[[k]]; ok <- (as.integer(t$observed) == r$observed) && abs(t$null_mean - r$null_mean) <= 1e-9 && abs(t$perm_p - r$perm_p_plain) <= 1e-12 && abs(t$null_sd - r$null_sd) <= 1e-9
    maxd <- max(maxd, abs(t$null_mean - r$null_mean)); n_ok <- n_ok + as.integer(ok)
    L <- c(L, sprintf("%-12s %s  %3d / %3d        %9s / %9s   %6s / %6s    %8s / %8s   %s", gse, tag, as.integer(t$observed), r$observed, py_fmt_f(t$null_mean, 4), py_fmt_f(r$null_mean, 4), py_fmt_f(t$perm_p, 3), py_fmt_f(r$perm_p_plain, 3), py_fmt_f(t$null_sd, 4), py_fmt_f(r$null_sd, 4), if (ok) "OK" else "MISMATCH"))
  }
  L <- c(L, "", sprintf("   matches: %d of %d (series x set); largest absolute difference in the null mean %s", n_ok, n_tot, sprintf("%.2e", maxd)))
  gate <- n_ok == n_tot && n_tot == 2L * ex$n_analysed && ex$n_all_ok == ex$n_rows_processed && ex$n_analysed == 37L
  L <- c(L, "", sprintf("GATE RESULT: %s", if (gate) sprintf("PASS (all %d analysed series reload, the draws are identical and the part 1 numbers are reproduced exactly)", ex$n_analysed) else "FAIL"), "")
  ## 3. consistency of analyse_matrix (coexp_checks) with part 1 (ported 25 Nov 2026 from the Python report; reads the R core results)
  if (file.exists(F_CORE)) {
    core <- list(); for (j in read_jsonl(F_CORE)) core[[paste(j$gse, j$tag)]] <- j; R <- list()
    for (k in names(ref)) {
      if (is.null(core[[k]])) next; f <- core[[k]]$flat; r0 <- ref[[k]]; u <- function(x) as.numeric(unlist(x)); s0 <- r0$share_between_raw; s0 <- if (is.null(s0)) NaN else as.numeric(s0)
      R[[length(R) + 1]] <- data.frame(obs_eq = f$observed == r0$observed, same_eq = f$observed_same_sign == r0$observed_same_sign,
        pooled_eq = isTRUE(all(u(list(f$pooled_raw, f$pooled_den1, f$pooled_den2, f$pooled_den3)) == u(r0$pooled_sig_per_scheme))),
        adj_eq = isTRUE(all(u(list(f$adjq_raw_obs, f$adjq_den1_obs, f$adjq_den2_obs, f$adjq_den3_obs)) == u(r0$adj_q05_per_scheme))), any_eq = f$adj_any == r0$cons_adj_any,
        all4_eq = f$adj_all4 == r0$cons_adj_all4, pairs4_eq = f$adj_pairs_all4 == r0$adj_pairs_all4, share_eq = if (is.nan(s0)) TRUE else isTRUE(abs(as.numeric(f$share_between) - s0) < 1e-9),
        dmean = f$null_mean - r0$null_mean, dsd = f$null_sd - r0$null_sd, dp = f$perm_p - r0$perm_p_plain, seflag = abs(f$null_mean - r0$null_mean) / max(r0$null_sd, 1e-9) * sqrt(2000))
    }
    R <- do.call(rbind, R); n <- nrow(R); se <- R$seflag / sqrt(2); se <- se[!is.na(se)]
    L <- c(L, sprintf("3. Consistency of the new module (coexp_checks via bench2_core, 2000 permutations) with part 1 on the same %d matrices (not a gate: random streams differ by design)", n),
           sprintf("   identical observed values: consensus count %d/%d, same-sign count %d/%d, pooled counts per scheme %d/%d, adjusted q<0.05 counts per scheme %d/%d, consensus pairs retained (any scheme) %d/%d, (all schemes) %d/%d,",
                   sum(R$obs_eq), n, sum(R$same_eq), n, sum(R$pooled_eq), n, sum(R$adj_eq), n, sum(R$any_eq), n, sum(R$all4_eq), n),
           sprintf("   pairs significant in all schemes after adjustment %d/%d, covariance share %d/%d", sum(R$pairs4_eq), n, sum(R$share_eq), n),
           sprintf("   null mean difference (new - part 1): mean %s, SD %s, largest absolute %s; null SD difference: mean %s; P difference: mean %s, largest absolute %s;", py_fmt_f(mean(R$dmean), 3), py_fmt_f(sd(R$dmean), 3),
                   py_fmt_f(max(abs(R$dmean)), 3), py_fmt_f(mean(R$dsd), 3), py_fmt_f(mean(R$dp), 4), py_fmt_f(max(abs(R$dp)), 3)),
           sprintf("   |difference of null means| in Monte Carlo SEs of the difference (sqrt(2/2000) x SD): median %s, 95th percentile %s, maximum %s", py_fmt_f(median(se), 2), py_fmt_f(unname(quantile(se, 0.95, type = 7)), 2), py_fmt_f(max(se), 2)), "")
  }
  write_md(L, F_REPORT, final_newline = TRUE); cat(paste(L, collapse = "\n"), "\n")
}
a <- commandArgs(trailingOnly = TRUE); step <- if (length(a)) a[1] else "report"
if (sys.nframe() == 0L) {
  if (step == "extract") extract() else if (step == "compare") compare(if (length(a) >= 2) as.integer(a[2]) else 1L, if (length(a) >= 3) as.integer(a[3]) else 1L)
  else if (step == "consistency") consistency() else if (step == "report") report() else if (step == "merge") merge_archived() else stop("step: extract | merge | compare [worker nworkers] | consistency | report")
}
