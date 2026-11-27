## summarize_benchmark.R - R port of benchmark_route3/summarize_benchmark.py (tables and summary; the figure is made with the figures of step 4).
## Summarizes results_per_dataset.jsonl into the pre-specified endpoints (PROTOCOL sections 4 and 5): processing_log.csv, per_dataset_results.csv and summary_endpoints.json. Applies the interpretive bands fixed
## in the protocol; nothing is tuned here. Usage: Rscript summarize_benchmark.R   (RB_RESULTS = the jsonl to read, default the archived Python results; RB_OUT = output folder)
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R"))
IN <- Sys.getenv("RB_RESULTS", normalizePath(file.path(R_PORT, "..", "..", "benchmark_route3", "results_per_dataset.jsonl"))); OUT <- Sys.getenv("RB_OUT", file.path(R_PORT, "results_R", "benchmark_route3"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
recs <- lapply(readLines(IN, encoding = "UTF-8", warn = FALSE), fromJSON, simplifyVector = FALSE); g <- function(r, k) if (is.null(r[[k]])) NA else r[[k]]
status <- data.frame(gse = vapply(recs, function(r) r$gse, ""), platform = vapply(recs, function(r) r$platform, ""), order = vapply(recs, function(r) as.integer(r$order), 1L), status = vapply(recs, function(r) r$status, ""),
                     n_de = vapply(recs, function(r) as.numeric(g(r, "n_de")), 1), stringsAsFactors = FALSE)
py_write_csv_repr(status, file.path(OUT, "processing_log.csv"), eol = .NL)   # pandas to_csv ends lines with os.linesep (the archived files are CRLF); before 26 Nov 2026 this wrote LF (audit C3)
ok <- Filter(function(r) r$status == "analysed", recs); rows <- list()
for (r in ok) for (tag in c("A", "B", "A_spike")) {
  x <- r[[tag]]; obs <- x$observed
  rows[[length(rows) + 1]] <- data.frame(gse = r$gse, platform = r$platform, set = tag, scale = r$scale, n_de = as.integer(r$n_de), n_de_in_set = as.numeric(g(x, "n_de_in_set")), observed = as.integer(obs), null_mean = as.numeric(x$null_mean), null_sd = as.numeric(x$null_sd),
                                         null_q95 = as.numeric(x$null_q95), z = as.numeric(x$z), perm_p = as.numeric(x$perm_p_plain), observed_same_sign = as.integer(x$observed_same_sign), adj_any = as.integer(x$cons_adj_any), adj_all4 = as.integer(x$cons_adj_all4),
                                         adj_pairs_all4 = as.integer(x$adj_pairs_all4), retained_any = if (obs) x$cons_adj_any / obs else NA_real_, retained_all4 = if (obs) x$cons_adj_all4 / obs else NA_real_, share_between = as.numeric(x$share_between_raw),
                                         pooled_raw = as.integer(x$pooled_sig_per_scheme[[1]]), adj_q05_raw = as.integer(x$adj_q05_per_scheme[[1]]), stringsAsFactors = FALSE)
}
T <- do.call(rbind, rows); py_write_csv_repr(T, file.path(OUT, "per_dataset_results.csv"), eol = .NL)
anal <- status[status$status == "analysed", ]; pp <- table(anal$platform); na_ <- table(gsub("[0-9]+", "N", status$status[status$status != "analysed"])); na_ <- na_[order(-as.integer(na_), method = "radix")]
S <- list(n_series_logged = length(recs), n_analysed = length(ok), analysed_per_platform = lapply(as.list(pp[py_sort(names(pp))]), as.integer), not_analysable = lapply(as.list(na_), as.integer))
wilson <- function(k, n, alpha = 0.05) { z <- qnorm(1 - alpha / 2); p <- k / n; den <- 1 + z^2 / n; centre <- (p + z^2 / 2 / n) / den; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den; c(centre - half, centre + half) }     # statsmodels proportion_confint(method='wilson')
prop <- function(mask) { k <- as.integer(sum(mask)); n <- as.integer(length(mask)); ci <- if (n) wilson(k, n) else c(NaN, NaN); list(k = k, n = n, prop = if (n) k / n else NaN, ci95 = list(ci[1], ci[2])) }
med <- function(x) { x <- x[!is.na(x)]; list(n = as.integer(length(x)), median = if (length(x)) median(x) else NaN, q1 = if (length(x)) unname(quantile(x, .25, type = 7)) else NaN, q3 = if (length(x)) unname(quantile(x, .75, type = 7)) else NaN) }
for (tag in c("A", "B", "A_spike")) {
  d <- T[T$set == tag, ]
  S[[tag]] <- list(E1_not_above_null_P_ge_0_05 = prop(d$perm_p >= 0.05), E1_flagged_P_lt_0_05 = prop(d$perm_p < 0.05), z = med(d$z), observed = med(d$observed), null_mean = med(d$null_mean), observed_minus_null = med(d$observed - d$null_mean),
                   E2_retained_any_scheme = med(d$retained_any), E2_retained_all_schemes = med(d$retained_all4), E2_pairs_sig_all4_after_adjustment = med(d$adj_pairs_all4), E3_share_between = med(d$share_between), datasets_with_zero_observed = as.integer(sum(d$observed == 0)))
  for (gpl in c("GPL570", "GPL1261")) { dd <- d[d$platform == gpl, ]; S[[tag]][[paste0("by_", gpl)]] <- list(n = as.integer(nrow(dd)), P_ge_0_05 = if (nrow(dd)) prop(dd$perm_p >= 0.05) else NULL, retained_any = med(dd$retained_any)) }
}
band <- S$A$E1_not_above_null_P_ge_0_05$prop
S$interpretive_band_gene_set_A <- if (is.nan(band)) "no analysable dataset" else if (band >= 0.8) "pooled consensus counts are generally explained by the design (>= 80%)" else if (band >= 0.5) "often explained (50-80%)" else "generally not explained: within-group coupling is common (< 50%)"
S$positive_control_detected_P_lt_0_05 <- S$A_spike$E1_flagged_P_lt_0_05; S$positive_control_ok_ge_80pct <- if (S$A_spike$E1_flagged_P_lt_0_05$n) S$A_spike$E1_flagged_P_lt_0_05$prop >= 0.8 else NULL
write_json_typed(S, file.path(OUT, "summary_endpoints.json")); cat(py_json_typed(S), "\n"); cat("figure: made with the figures (step 4)\n")
