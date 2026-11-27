## summarize_validation.R - R port of wp2_benchmark2/summarize_validation.py: protocol v2 section 3.6 (WP2 deliverable F): endpoints on the validation sets, pooled and per set, next to the development series.
## Reads dev_application_per_series.csv (development, 37 series) and heldout_arrays_results.jsonl (held-out arrays) and, if the RNA-seq arm was run, heldout_rnaseq_results.jsonl.
## Writes validation_per_series.csv, validation_endpoints_long.csv, validation_summary.json, validation_summary.md and WP2_summary_table.csv / .md into OUT (env WP2_OUT).
## Inputs: WP2_DEV_CSV, WP2_HELDOUT, WP2_RNASEQ (defaults: the archived Python files; point them to the files of the R run).  Usage: Rscript summarize_validation.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp2_common.R")); source(file.path(HERE, "endpoints_lib.R")); suppressMessages(library(data.table))
load_jsonl <- function(path, arm_of) {
  if (!file.exists(path)) return(NULL); rows <- list()
  for (r in read_jsonl(path)) {
    if (!identical(r$status, "analysed")) next
    for (tag in c("A", "B")) rows[[length(rows) + 1]] <- c(list(gse = r$gse, platform = if (!is.null(r$platform)) r$platform else r$organism, gene_set = tag, arm = arm_of(r), n_de = r$n_de), nullNA(r[[tag]]$flat))
  }
  if (length(rows)) rows_df(rows) else NULL
}
dev <- read.csv(Sys.getenv("WP2_DEV_CSV", file.path(ARC, "dev_application_per_series.csv")), stringsAsFactors = FALSE, check.names = FALSE); dev$arm <- "development arrays"
ho <- load_jsonl(Sys.getenv("WP2_HELDOUT", file.path(ARC, "heldout_arrays_results.jsonl")), function(r) paste("held-out", r$platform))
rs <- load_jsonl(Sys.getenv("WP2_RNASEQ", file.path(ARC, "heldout_rnaseq_results.jsonl")), function(r) paste("RNA-seq", r$platform))
ALL <- as.data.frame(rbindlist(Filter(Negate(is.null), list(dev, ho, rs)), use.names = TRUE, fill = TRUE)); write_csv_crlf(ALL, file.path(OUT, "validation_per_series.csv"))
sw <- function(x, p) startsWith(x, p)
ARMS <- list(list("Development arrays (37)", function(d) d$arm == "development arrays"))
if (!is.null(ho)) {
  for (pl in sort(unique(ho$platform), method = "radix")) ARMS[[length(ARMS) + 1]] <- list(sprintf("Held-out %s", pl), local({ p <- pl; function(d) d$arm == paste("held-out", p) }))
  ARMS[[length(ARMS) + 1]] <- list("Held-out arrays pooled", function(d) sw(d$arm, "held-out"))
}
if (!is.null(rs)) {
  for (pl in sort(unique(rs$platform), method = "radix")) ARMS[[length(ARMS) + 1]] <- list(sprintf("RNA-seq %s", pl), local({ p <- pl; function(d) d$arm == paste("RNA-seq", p) }))
  ARMS[[length(ARMS) + 1]] <- list("RNA-seq pooled", function(d) sw(d$arm, "RNA-seq"))
}
if (!is.null(ho) || !is.null(rs)) ARMS[[length(ARMS) + 1]] <- list("All validation sets pooled", function(d) sw(d$arm, "held-out") | sw(d$arm, "RNA-seq"))
# post hoc sensitivity (protocol section 8, WP2 entry): leave out the matrices in which at least one of the 17 genes has no within-group variation
ARMS[[length(ARMS) + 1]] <- list("Development arrays without degenerate series (post hoc)", function(d) d$arm == "development arrays" & d$n_zero_var_genes == 0)
if (!is.null(ho) || !is.null(rs)) ARMS[[length(ARMS) + 1]] <- list("All validation sets pooled without degenerate series (post hoc)", function(d) (sw(d$arm, "held-out") | sw(d$arm, "RNA-seq")) & d$n_zero_var_genes == 0)
S <- list(); long_rows <- list()
for (gs in c("A", "B")) for (a in ARMS) {
  sel <- a[[2]](ALL); sel[is.na(sel)] <- FALSE; d <- ALL[ALL$gene_set == gs & sel, , drop = FALSE]
  if (!nrow(d)) next
  e <- endpoints(d); S[[sprintf("%s | set %s", a[[1]], gs)]] <- e; long_rows <- c(long_rows, long_table(a[[1]], gs, e))
}
write_csv_crlf(rows_df(long_rows), file.path(OUT, "validation_endpoints_long.csv")); write_json_typed(S, file.path(OUT, "validation_summary.json"))
binom_p <- function(d) if (d$n) pbinom(d$k - 1, d$n, 0.05, lower.tail = FALSE) else NaN                 # scipy binomtest(k, n, 0.05, alternative = 'greater')
fisher <- function(d1, d2) if (!d1$n || !d2$n) NaN else fisher.test(matrix(c(d1$k, d2$k, d1$n - d1$k, d2$n - d2$k), 2))$p.value
g2 <- function(x) sprintf("%.2g", x)                                                                    # '%.2g' of Python and C agree except that R keeps 'e-05' exponents like Python
md <- c("# Validation sets: endpoints of protocol 3.6 (WP2 deliverable F)", ""); rowdefs <- NULL
for (gs in c("A", "B")) {
  md <- c(md, sprintf("## Gene set %s (%s)", gs, if (gs == "A") "17 random differentially expressed genes" else "17 random expressed genes"), "")
  keys <- lapply(ARMS, function(a) list(a[[1]], S[[sprintf("%s | set %s", a[[1]], gs)]])); keys <- Filter(function(k) !is.null(k[[2]]), keys)
  md <- c(md, paste0("| Endpoint | ", paste(vapply(keys, function(k) sprintf("%s (n = %d)", k[[1]], k[[2]]$n_series), ""), collapse = " | "), " |"), paste0("|---|", strrep("---|", length(keys))))
  rowdefs <- list(list("Within-group coupling: T_r2 P < 0.05 (primary)", function(e) fmt_prop(e$coupling_T_r2)), list("T_p P < 0.05", function(e) fmt_prop(e$coupling_T_p)), list("T_eig P < 0.05", function(e) fmt_prop(e$coupling_T_eig)),
                  list("E1: consensus count not above null (P >= 0.05)", function(e) fmt_prop(e$E1_count_not_above_null)),
                  list("E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3]", function(e) fmt_med(e$E2_retained_raw_scheme, pct = TRUE)),
                  list("E2 any scheme (part 1 definition)", function(e) fmt_med(e$E2_retained_any_scheme, pct = TRUE)), list("E2 all four schemes", function(e) fmt_med(e$E2_retained_all_schemes, pct = TRUE)),
                  list("Pairs significant in all schemes after adjustment", function(e) fmt_med(e$E2_pairs_sig_all_after_adjustment, nd = 1)),
                  list("E3 covariance share carried by group means", function(e) fmt_med(e$E3_share_between, pct = TRUE)),
                  list("Raw pooled-significant share among group-effect genes (expectation 0.61)", function(e) fmt_med(e$raw_sig_share_group_effect)), list("Series with share above 0.61", function(e) fmt_prop(e$raw_sig_share_above_expectation)))
  for (s in SCHEMES) {
    rowdefs[[length(rowdefs) + 1]] <- list(sprintf("Adjusted q < 0.05 count above its null (P < 0.05), %s scheme", SCHEME_LABEL[[s]]), local({ s <- s; function(e) fmt_prop(e[[paste0("adj_", s)]]$share_P_lt_0_05) }))
    rowdefs[[length(rowdefs) + 1]] <- list(sprintf("Median observed / null mean adjusted count, %s scheme", SCHEME_LABEL[[s]]), local({ s <- s; function(e) fmt_med(e[[paste0("adj_", s)]]$median_ratio_obs_over_null, nd = 2) }))
  }
  for (rd in rowdefs) md <- c(md, paste0("| ", rd[[1]], " | ", paste(vapply(keys, function(k) rd[[2]](k[[2]]), ""), collapse = " | "), " |"))
  md <- c(md, "", paste0("T_r2 prevalence against the 5% expected under no coupling (exact binomial, one-sided): ", paste(vapply(keys, function(k) sprintf("%s P = %s", k[[1]], g2(binom_p(k[[2]]$coupling_T_r2))), ""), collapse = "; ")))
  d0 <- S[[sprintf("Development arrays (37) | set %s", gs)]]; dh <- S[[sprintf("Held-out arrays pooled | set %s", gs)]]
  if (!is.null(d0) && !is.null(dh)) md <- c(md, sprintf("Held-out arrays versus development (Fisher exact, two-sided): T_r2 prevalence P = %s; consensus count not above null P = %s.", g2(fisher(dh$coupling_T_r2, d0$coupling_T_r2)), g2(fisher(dh$E1_count_not_above_null, d0$E1_count_not_above_null))))
  md <- c(md, "")
}
write_md(md, file.path(OUT, "validation_summary.md")); cat(paste(md, collapse = "\n"), "\n")
# ---- one summary table of numbers (arms in columns), gene set A and B stacked
tab <- list()
for (gs in c("A", "B")) for (rd in rowdefs) {
  row <- list(gene_set = gs, endpoint = rd[[1]])
  for (a in ARMS) { k <- sprintf("%s | set %s", a[[1]], gs); if (!is.null(S[[k]])) row[[a[[1]]]] <- rd[[2]](S[[k]]) }
  tab[[length(tab) + 1]] <- row
}
Tt <- rows_df(tab); write_csv_crlf(Tt, file.path(OUT, "WP2_summary_table.csv"))
cols <- names(Tt); lines <- c(paste0("| ", paste(cols, collapse = " | "), " |"), paste0("|", strrep("---|", length(cols))), apply(Tt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
write_md(lines, file.path(OUT, "WP2_summary_table.md"), final_newline = TRUE)
