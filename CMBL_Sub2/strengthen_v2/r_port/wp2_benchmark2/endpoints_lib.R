## endpoints_lib.R - R port of wp2_benchmark2/endpoints_lib.py: endpoints of protocol v2 sections 3.3, 3.4 and 3.6 from per-matrix results (a data frame, one row per series and gene set).
## Nothing is tuned here: P < 0.05 everywhere, Wilson 95% intervals (statsmodels operation order, z = scipy norm.isf), medians with quartiles (numpy linear interpolation).
## P conventions: the consensus-count endpoint E1 uses the plain fraction P(null >= observed) as in part 1 (perm_p); every other test uses (1 + k)/(n + 1) (columns *_p, adjq_*_p).
## R integer = Python int, R double = Python float (py_json_typed): counters are made with as.integer.
GEO_EXPECTED <- 0.6086          # (2/pi) arccos(0.576): share of pooled-significant raw pairs if group means carry the correlation (3 groups, df 10)
SCHEMES <- c("raw", "den1", "den2", "den3")
SCHEME_LABEL <- c(raw = "raw", den1 = "ratio to gene 1", den2 = "ratio to gene 2", den3 = "ratio to gene 3")

prop <- function(mask) {            # NaN comparisons are False in pandas: NA -> FALSE
  mask <- as.logical(mask); mask[is.na(mask)] <- FALSE; k <- as.integer(sum(mask)); n <- as.integer(length(mask)); ci <- py_wilson(k, n)
  list(k = k, n = n, prop = if (n) k / n else NaN, lo = ci[1], hi = ci[2])
}
med <- function(x) {
  x <- as.numeric(x); x <- x[!is.na(x)]
  if (!length(x)) return(list(n = 0L, median = NaN, q1 = NaN, q3 = NaN))
  list(n = as.integer(length(x)), median = median(x), q1 = np_quantile(x, 0.25), q3 = np_quantile(x, 0.75))
}
endpoints <- function(T) {           # T: data frame, one row per series (ONE gene set). Returns a nested list of all endpoints
  g <- function(nm) T[[nm]]
  out <- list(n_series = as.integer(nrow(T)))
  out$E1_count_not_above_null <- prop(g("perm_p") >= 0.05); out$E1_count_flagged <- prop(g("perm_p") < 0.05)
  out$coupling_T_r2 <- prop(g("T_r2_p") < 0.05)
  if (!is.null(g("T_p_p"))) out$coupling_T_p <- prop(g("T_p_p") < 0.05)
  if (!is.null(g("T_eig_p"))) out$coupling_T_eig <- prop(g("T_eig_p") < 0.05)
  out$T_r2_observed <- med(g("T_r2_obs")); out$T_r2_null_mean <- med(g("T_r2_null_mean"))
  ob <- g("observed") > 0; ob[is.na(ob)] <- FALSE
  out$E2_retained_any_scheme <- med(g("retained_any")[ob]); out$E2_retained_all_schemes <- med(g("retained_all4")[ob]); out$E2_retained_raw_scheme <- med(g("retained_raw")[ob])
  out$E2_n_series_with_consensus_pairs <- as.integer(sum(ob)); out$E2_pairs_sig_all_after_adjustment <- med(g("adj_pairs_all4"))
  out["E2_raw_any_consensus_pair_retained"] <- list(if (sum(ob)) prop(g("adj_cons_raw")[ob] > 0) else NULL)
  out$E3_share_between <- med(g("share_between"))
  out$raw_sig_share_group_effect <- med(g("share_sig_group_effect")); out$raw_sig_share_all_pairs <- med(g("share_sig_all"))
  ok <- !is.na(g("share_sig_group_effect")); out$raw_sig_share_above_expectation <- prop(g("share_sig_group_effect")[ok] > GEO_EXPECTED)
  out$observed_consensus <- med(g("observed")); out$null_mean_consensus <- med(g("null_mean")); out$consensus_minus_null <- med(g("observed") - g("null_mean"))
  for (s in SCHEMES) {
    o <- g(paste0("adjq_", s, "_obs")); nm <- g(paste0("adjq_", s, "_null_mean")); p <- g(paste0("adjq_", s, "_p"))
    ratio <- ifelse(!is.na(nm) & nm > 0, o / nm, NA_real_)
    out[[paste0("adj_", s)]] <- list(share_P_lt_0_05 = prop(p < 0.05), median_observed = med(o), median_null_mean = med(nm), median_ratio_obs_over_null = med(ratio),
                                     n_ratio_defined = as.integer(sum(!is.na(ratio))), median_excess = med(o - nm), n_series_observed_gt_0 = as.integer(sum(o > 0, na.rm = TRUE)))
  }
  out
}
.row <- function(arm, gset, key, scheme, d, kind = "proportion") {
  if (kind == "proportion") list(arm = arm, gene_set = gset, endpoint = key, scheme = scheme, k = as.numeric(d$k), n = d$n, estimate = d$prop, lo = d$lo, hi = d$hi)
  else list(arm = arm, gene_set = gset, endpoint = key, scheme = scheme, k = NA_real_, n = d$n, estimate = d$median, lo = d$q1, hi = d$q3)
}
long_table <- function(arm, gset, e) {     # flatten endpoints() into rows (proportions with Wilson limits, medians with quartiles as lo/hi); k is a float column (NaN in the median rows) like pandas
  rows <- list()
  for (key in names(e)) {
    d <- e[[key]]; if (is.null(d) || !is.list(d)) next
    if ("prop" %in% names(d)) rows[[length(rows) + 1]] <- .row(arm, gset, key, "", d)
    else if ("median" %in% names(d)) rows[[length(rows) + 1]] <- .row(arm, gset, key, "", d, "median")
    else if (startsWith(key, "adj_")) {
      s <- substring(key, 5); rows[[length(rows) + 1]] <- .row(arm, gset, "adj_share_P_lt_0_05", s, d$share_P_lt_0_05)
      for (sub in c("median_observed", "median_null_mean", "median_ratio_obs_over_null", "median_excess")) rows[[length(rows) + 1]] <- .row(arm, gset, paste0("adj_", sub), s, d[[sub]], "median")
    }
  }
  rows
}
fmt_prop <- function(d) if (!is.null(d) && d$n) sprintf("%d/%d = %s%% (%s-%s%%)", d$k, d$n, py_fmt_f(100 * d$prop, 0), py_fmt_f(100 * d$lo, 0), py_fmt_f(100 * d$hi, 0)) else "n/a"
fmt_med <- function(d, pct = FALSE, nd = 2) {      # py_fmt_f: Python's correctly rounded '%.Nf' (R's sprintf rounds near-ties differently on Windows)
  if (is.null(d) || d$n == 0) return("n/a")
  f <- if (pct) function(v) paste0(py_fmt_f(100 * v, 0), "%") else function(v) py_fmt_f(v, nd)
  sprintf("%s [%s, %s] (n=%d)", f(d$median), f(d$q1), f(d$q3), d$n)
}
