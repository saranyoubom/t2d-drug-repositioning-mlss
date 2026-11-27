## apply_dev.R - R port of wp2_benchmark2/apply_dev.py: protocol v2 sections 3.3 and 3.4 on the 37 development series (WP2 deliverable C).
## Reads compact/dev_core_results.jsonl (bench2_core analyse_matrix on the 74 matrices of gate_part1: gene sets A and B; coupling tests 2000 permutations, consensus null 2000, calibrated adjusted counts 500).
## Writes dev_application_per_series.csv, dev_endpoints_long.csv, dev_application_summary.json and dev_application_summary.md into OUT (env WP2_OUT).
## Input: WP2_DEV_CORE (default: the archived Python file compact/dev_core_results.jsonl; point it to <OUT>/compact/dev_core_results.jsonl for the R run).  Usage: Rscript apply_dev.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp2_common.R")); source(file.path(HERE, "endpoints_lib.R"))
recs <- read_jsonl(Sys.getenv("WP2_DEV_CORE", file.path(ARC, "compact", "dev_core_results.jsonl")))
T <- rows_df(lapply(recs, function(r) c(list(gse = r$gse, platform = r$platform, gene_set = r$tag), nullNA(r$flat)))); write_csv_crlf(T, file.path(OUT, "dev_application_per_series.csv"))
p1 <- read.csv(file.path(P1, "per_dataset_results.csv"), stringsAsFactors = FALSE); S <- list(); long_rows <- list(); md <- c("# Development series: protocol 3.3 and 3.4 (WP2 deliverable C)", "")
arms <- list(list("development, all 37 series", NULL), list("development, human GPL570", "GPL570"), list("development, mouse GPL1261", "GPL1261"),
             list("development, without series with degenerate genes (post hoc sensitivity)", "ND"))     # ND: no gene with zero within-group variance in the 17-gene set (protocol section 8, WP2 entry)
for (gs in c("A", "B")) {
  for (a in arms) {
    label <- a[[1]]; pl <- a[[2]]
    sel <- T$gene_set == gs & (if (is.null(pl)) TRUE else if (pl == "ND") T$n_zero_var_genes == 0 else T$platform == pl); d <- T[sel, , drop = FALSE]
    if (!nrow(d)) next
    e <- endpoints(d); S[[sprintf("%s | set %s", label, gs)]] <- e; long_rows <- c(long_rows, long_table(label, gs, e))
  }
  q <- p1[p1$set == gs, ]      # part 1 values for the same series (loop-based null, 2000 permutations) for comparison
  S[[sprintf("part 1 | set %s", gs)]] <- list(E1_count_not_above_null = prop(q$perm_p >= 0.05), retained_any = med(q$retained_any), retained_all4 = med(q$retained_all4), adj_pairs_all4 = med(q$adj_pairs_all4), share_between = med(q$share_between))
}
write_csv_crlf(rows_df(long_rows), file.path(OUT, "dev_endpoints_long.csv")); write_json_typed(S, file.path(OUT, "dev_application_summary.json"))
for (gs in c("A", "B")) {
  e <- S[[sprintf("development, all 37 series | set %s", gs)]]; q <- S[[sprintf("part 1 | set %s", gs)]]
  md <- c(md, sprintf("## Gene set %s (%s)", gs, if (gs == "A") "17 random differentially expressed genes, FDR < 0.05" else "17 random expressed genes"), "",
          "| Endpoint | Result (37 series) |", "|---|---|",
          sprintf("| Within-group coupling, T_r2 (primary), P < 0.05 | %s |", fmt_prop(e$coupling_T_r2)), sprintf("| T_p, P < 0.05 | %s |", fmt_prop(e$coupling_T_p)), sprintf("| T_eig, P < 0.05 | %s |", fmt_prop(e$coupling_T_eig)),
          sprintf("| Consensus count not above its permutation null (E1; P >= 0.05), this run | %s |", fmt_prop(e$E1_count_not_above_null)), sprintf("| same, part 1 (loop null) | %s |", fmt_prop(q$E1_count_not_above_null)),
          sprintf("| E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = %d) | median %s |", e$E2_n_series_with_consensus_pairs, fmt_med(e$E2_retained_raw_scheme, pct = TRUE)),
          sprintf("| E2 any scheme (part 1 definition), this run / part 1 | %s / %s |", fmt_med(e$E2_retained_any_scheme, pct = TRUE), fmt_med(q$retained_any, pct = TRUE)),
          sprintf("| E2 all four schemes, this run / part 1 | %s / %s |", fmt_med(e$E2_retained_all_schemes, pct = TRUE), fmt_med(q$retained_all4, pct = TRUE)),
          sprintf("| Pairs significant in all schemes after adjustment, this run / part 1 | %s / %s |", fmt_med(e$E2_pairs_sig_all_after_adjustment, nd = 1), fmt_med(q$adj_pairs_all4, nd = 1)),
          sprintf("| E3 share of pooled covariance carried by group means | %s |", fmt_med(e$E3_share_between, pct = TRUE)),
          sprintf("| Raw pooled-significant share among genes with a group effect (expectation 0.61) | %s; above 0.61 in %s |", fmt_med(e$raw_sig_share_group_effect), fmt_prop(e$raw_sig_share_above_expectation)), "")
  md <- c(md, "Ratio-scheme diagnostics (observed adjusted q < 0.05 count against its null from 500 within-group permutations):", "", "| Scheme | Series with P < 0.05 | Median observed | Median null mean | Median observed / null mean | Median excess |", "|---|---|---|---|---|---|")
  for (s in SCHEMES) {
    a <- e[[paste0("adj_", s)]]
    md <- c(md, sprintf("| %s | %s | %s | %s | %s (defined in %d) | %s |", SCHEME_LABEL[[s]], fmt_prop(a$share_P_lt_0_05), fmt_med(a$median_observed, nd = 1), fmt_med(a$median_null_mean, nd = 1), fmt_med(a$median_ratio_obs_over_null, nd = 2), a$n_ratio_defined, fmt_med(a$median_excess, nd = 1)))
  }
  md <- c(md, "")
  for (a in arms[2:3]) {
    label <- a[[1]]; pl <- a[[2]]; ee <- endpoints(T[T$gene_set == gs & T$platform == pl, , drop = FALSE])
    md <- c(md, sprintf("By platform, %s: T_r2 P < 0.05 in %s; consensus count not above null in %s.", strsplit(label, ", ", fixed = TRUE)[[1]][2], fmt_prop(ee$coupling_T_r2), fmt_prop(ee$E1_count_not_above_null)))
  }
  ee <- S[[sprintf("%s | set %s", arms[[4]][[1]], gs)]]
  md <- c(md, sprintf("Post hoc sensitivity, %s (n = %d): T_r2 P < 0.05 in %s; T_p in %s; consensus count not above null in %s; E2 raw scheme alone %s; covariance share carried by group means %s.",
                      strsplit(arms[[4]][[1]], ", ", fixed = TRUE)[[1]][2], ee$n_series, fmt_prop(ee$coupling_T_r2), fmt_prop(ee$coupling_T_p), fmt_prop(ee$E1_count_not_above_null), fmt_med(ee$E2_retained_raw_scheme, pct = TRUE), fmt_med(ee$E3_share_between, pct = TRUE)), "")
}
write_md(md, file.path(OUT, "dev_application_summary.md")); cat(paste(md, collapse = "\n"), "\n")
