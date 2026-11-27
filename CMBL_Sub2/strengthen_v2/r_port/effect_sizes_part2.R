## effect_sizes_part2.R - R port of wp5_reporting/effect_sizes.py, PART 2 (Python lines 174-302): rt_panel, raw_p, reproduce_reported, gene_rows.
## Group indices are 1-based here (Python 0-based): CTRL = 1, DAPT = 2, DKK-1 = 3.
rt_panel <- function(af3) {
  rows <- af4_rows(AF4, "Fig_4c_Raw_expression_data"); h <- which(vapply(rows, function(r) length(r) && is_str(r[[1]]) && r[[1]] == "Wnt5a", TRUE))[1]
  genes <- hdr_names(rows[[h]]); genes <- genes[!is.na(genes)]; X <- do.call(rbind, lapply(rows[(h + 1):(h + 12)], function(r) vapply(r[seq_along(genes)], function(c) as.numeric(c), 1)))
  stopifnot(identical(dim(X), c(12L, 17L)), setequal(genes, vapply(af3, function(r) r$gene, "")))
  if (file.exists(FIG4B)) { m <- as.matrix(as.data.frame(read_excel(FIG4B))[, genes]); storage.mode(m) <- "double"; GATE$matrix_AF4_vs_Fig4b_max_abs_rel_diff <<- max(abs(m - X) / abs(X)) } else GATE["matrix_AF4_vs_Fig4b_max_abs_rel_diff"] <<- list(NULL)
  res <- list(); over <- numeric(0)
  for (g in genes) {
    x <- X[, match(g, genes)]; fc <- x / mean(x[GI$CTRL]); res[[g]] <- list(x = x, fc = fc, groups = list(fc[GI$CTRL], fc[GI$DAPT], fc[GI[["DKK-1"]]]))
    cc <- x[GI$CTRL]; over <- c(over, mean(cc) / exp(mean(log(cc))))
  }
  GATE$ctrl_mean_if_normalized_to_mean_dCt_max <<- max(over); GATE$ctrl_mean_if_normalized_to_mean_dCt_median <<- median(over)
  list(genes = genes, res = res)
}
raw_p <- function(groups, i, j, use_kw) if (use_kw) pairwise_p_wilcox(groups[[j]], groups[[i]]) else pairwise_p_anova(groups, i, j)
## gate: Additional file 3 (means, SEMs, replicates, overall and pairwise p) against two recomputations: 'printed' = from the printed 3-decimal replicate values R1-R4,
## 'full' = from the 12 full-precision values of Additional file 4
reproduce_reported <- function(af3, res) {
  names_ <- c("CTRL", "DAPT", "DKK-1"); G <- list()
  for (basis in c("printed", "full")) {
    d_mean <- numeric(0); d_sem <- numeric(0); d_rep <- numeric(0); d_ctrl1 <- numeric(0); p_ok <- 0L; p_n <- 0L; ov_ok <- 0L; ov_n <- 0L; p_bad <- list(); ov_bad <- list()   # counters are integers (JSON ints)
    for (r in af3) {
      g <- r$gene; groups <- if (basis == "printed") lapply(names_, function(k) as.numeric(r$rep[[k]])) else res[[g]]$groups
      for (k in 1:3) { a <- groups[[k]]; nm <- names_[k]; d_mean <- c(d_mean, abs(mean(a) - r$mean[[nm]])); d_sem <- c(d_sem, abs(sem_(a) - r$sem[[nm]])); d_rep <- c(d_rep, abs(a - as.numeric(r$rep[[nm]]))) }
      d_ctrl1 <- c(d_ctrl1, abs(mean(groups[[1]]) - 1)); use_kw <- startsWith(r$test, "Kruskal")
      for (pr in list(list(c(1, 2), "DAPT"), list(c(1, 3), "DKK-1"), list(c(2, 3), "D-K"))) {
        praw <- raw_p(groups, pr[[1]][1], pr[[1]][2], use_kw); key <- pr[[2]]
        for (vs in list(list(praw, r$p_raw[[key]]), list(min(1.0, 3 * praw), r$p_adj[[key]]))) {
          p_n <- p_n + 1L
          if (printed_match(vs[[1]], vs[[2]])) p_ok <- p_ok + 1L else p_bad[[length(p_bad) + 1]] <- list(g, key, vs[[1]], as.character(vs[[2]]))
        }
      }
      ov <- if (use_kw) kruskal.test(groups)$p.value else oneway.test(unlist(groups) ~ factor(rep(1:3, vapply(groups, length, 1))), var.equal = TRUE)$p.value
      ov_n <- ov_n + 1L; if (printed_match(ov, r$overall_p)) ov_ok <- ov_ok + 1L else ov_bad[[length(ov_bad) + 1]] <- list(g, ov, as.character(r$overall_p))
    }
    G[[basis]] <- list(group_means_max_abs_diff = max(d_mean), group_sems_max_abs_diff = max(d_sem), n_group_means = length(d_mean), replicates_max_abs_diff = max(d_rep), n_replicates = length(d_rep),
                       ctrl_mean_is_1_max_abs_dev = max(d_ctrl1), pairwise_p_checked = p_n, pairwise_p_matching_printed = p_ok, pairwise_p_not_matching = p_bad, overall_p_checked = ov_n,
                       overall_p_matching_printed = ov_ok, overall_p_not_matching = ov_bad)
  }
  GATE$n_genes <<- length(af3); GATE$from_printed_replicates <<- G$printed; GATE$from_full_precision_values <<- G$full
}
## one record per gene and contrast (reported analysis)
gene_rows <- function(af3, genes, res) {
  names_ <- c("CTRL", "DAPT", "DKK-1"); ci_diff <- numeric(0)
  for (r in af3) {
    g <- r$gene; use_kw <- startsWith(r$test, "Kruskal"); printed <- lapply(names_, function(k) as.numeric(r$rep[[k]])); res[[g]]$praw <- list(); res[[g]]$praw_printed <- list()
    for (pr in list(list(c(1, 2), "DAPT"), list(c(1, 3), "DKK-1"))) {
      i <- pr[[1]][1]; j <- pr[[1]][2]; key <- pr[[2]]; res[[g]]$praw[[key]] <- raw_p(res[[g]]$groups, i, j, use_kw); res[[g]]$praw_printed[[key]] <- raw_p(printed, i, j, use_kw)
      a <- welch(log2(res[[g]]$groups[[j]]), log2(res[[g]]$groups[[i]])); b <- welch(log2(printed[[j]]), log2(printed[[i]])); ci_diff <- c(ci_diff, max(abs(a[1] - b[1]), abs(a[2] - b[2]), abs(a[3] - b[3])))
    }
  }
  GATE$ci_limits_full_vs_printed_replicates_max_abs_diff_log2 <<- max(ci_diff)
  rel <- unlist(lapply(af3, function(r) vapply(CONTRASTS[vapply(CONTRASTS, function(k) res[[r$gene]]$praw_printed[[k]] > 0, TRUE)], function(k) abs(res[[r$gene]]$praw[[k]] - res[[r$gene]]$praw_printed[[k]]) / res[[r$gene]]$praw_printed[[k]], 1)))
  GATE$p_raw_full_vs_printed_replicates_max_rel_diff <<- max(rel); GATE$p_raw_full_vs_printed_replicates_median_rel_diff <<- median(rel)
  q <- list(); gn <- vapply(af3, function(r) r$gene, "")
  for (key in CONTRASTS) {
    qq <- bh_(vapply(gn, function(g) res[[g]]$praw[[key]], 1)); q[[key]] <- setNames(as.list(unname(qq)), gn)
    qp <- bh_(vapply(gn, function(g) res[[g]]$praw_printed[[key]], 1)); qs <- bh_(vapply(af3, function(r) as.numeric(r$p_raw[[key]]), 1))
    GATE[[paste0("bh_q_max_abs_diff_full_vs_printed_replicates_", key)]] <<- max(abs(unlist(q[[key]]) - qp)); GATE[[paste0("bh_q_max_abs_diff_full_vs_printed_strings_", key)]] <<- max(abs(unlist(q[[key]]) - qs))
    GATE[[paste0("bh_q05_status_differs_", key)]] <<- as.list(gn[(unlist(q[[key]]) < 0.05) != (qp < 0.05)])
  }
  recs <- list()
  for (r in af3) {
    g <- r$gene; R <- res[[g]]; L <- lapply(R$groups, log2)
    for (kk in list(list("DAPT", 2), list("DKK-1", 3))) {
      key <- kk[[1]]; k <- kk[[2]]; w <- welch(L[[k]], L[[1]]); p3 <- pooled3(L, 1, k)
      recs[[length(recs) + 1]] <- list(readout = "RT-qPCR", feature = g, category = r$category, contrast = paste(key, "vs CTRL"), variant = "as reported", flag_26Aug_Gapdh = (key == "DAPT" && g %in% G5), test_in_manuscript = r$test,
        n_ctrl = 4L, n_trt = 4L, ctrl_mean_reported = r$mean[["CTRL"]], ctrl_sem_reported = r$sem[["CTRL"]], trt_mean_reported = r$mean[[key]], trt_sem_reported = r$sem[[key]], trt_mean_recomputed = mean(R$groups[[k]]),
        trt_sem_recomputed = sem_(R$groups[[k]]), fc_ratio_arith = mean(R$groups[[k]]) / mean(R$groups[[1]]), log2fc_geo = w[1], ci95_lo = w[2], ci95_hi = w[3], ci_method = "Welch t, log2 scale", df = w[4],
        fc_geo = 2^w[1], fc_geo_lo = 2^w[2], fc_geo_hi = 2^w[3], ci95_lo_pooled3 = p3[2], ci95_hi_pooled3 = p3[3], df_pooled3 = p3[4], ci_excludes_0 = (w[2] > 0 || w[3] < 0), ci_excludes_0_pooled3 = (p3[2] > 0 || p3[3] < 0),
        p_welch_log2 = w[5], p_raw_manuscript = as.character(r$p_raw[[key]]), p_adj_bonf_manuscript = as.character(r$p_adj[[key]]), p_raw_recomputed = R$praw[[key]], p_adj_bonf_recomputed = min(1.0, 3 * R$praw[[key]]),
        bh_q_17genes = q[[key]][[g]], interval_label = LABEL, note = "")
    }
  }
  list(recs = recs, q = q, res = res)       # res is returned because it now carries praw / praw_printed (Python mutates it in place)
}
