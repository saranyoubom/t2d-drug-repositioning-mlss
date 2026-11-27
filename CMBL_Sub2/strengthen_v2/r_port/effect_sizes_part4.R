## effect_sizes_part4.R - R port of wp5_reporting/effect_sizes.py, PART 4 (Python lines 359-725): other_readouts (the DNA block is dna_records() of part 3; 2-NBDG and C-peptide follow here),
## write_md and es_main. The csv.DictWriter-style CSV writer, the typed JSON writer (R integer = Python int, R double = Python float) and the line-end rule are in pycompat.R.
## Group indices are 1-based. The script name printed in the Markdown header can be set with ES_LABEL (default effect_sizes.R; "effect_sizes.py" reproduces the archived text).

## ---------------------------------------------------------------- 2-NBDG and C-peptide (after the DNA block of part 3)
other_readouts <- function() {
  d <- dna_records(); recs <- d$recs; extra <- d$extra; five <- extra$dna_five
  ## ---- 2-NBDG: 25,000 cells, second reading, as in Fig. 2b (S1), and the 11 other scenarios of Table S6
  nb <- af4_rows(AF4, "Fig_2b_2NBDG_raw_blocks"); h <- which(vapply(nb, function(r) length(r) > 0 && is_str(r[[1]]) && r[[1]] == "block", TRUE))[1]; blocks <- list()
  for (r in nb[(h + 1):length(nb)]) if (length(r) && !is.null(r[[1]]) && !identical(r[[1]], "")) {
    b <- r[[1]]; if (is.null(blocks[[b]])) blocks[[b]] <- list()
    blocks[[b]][[length(blocks[[b]]) + 1]] <- list(rep = as.integer(num(r[[2]])), blank = num(r[[3]]), CTRL = num(r[[4]]), DKK1 = num(r[[5]]), DAPT = num(r[[6]]))
  }
  sc <- af4_rows(AF4, "SuppTable6_2NBDG_scenarios"); h6 <- which(vapply(sc, function(r) length(r) > 0 && is_str(r[[1]]) && r[[1]] == "scenario", TRUE))[1]; sc_rep <- list()
  for (r in sc[(h6 + 1):length(sc)]) if (length(r) && is_str(r[[1]]) && nzchar(r[[1]])) sc_rep[[strsplit(r[[1]], " ", fixed = TRUE)[[1]][1]]] <- r
  sc_hdr <- hdr_names(sc[[h6]]); sc_col <- setNames(seq_along(sc_hdr), sc_hdr)
  vals <- function(block, grp_name, ctrl2 = NULL, drop_high_ctrl = FALSE, blank = "none") {
    b <- blocks[[block]]; v <- vapply(b, function(x) x[[grp_name]], 1); bl <- vapply(b, function(x) x$blank, 1)
    if (grp_name == "CTRL" && !is.null(ctrl2)) v[2] <- ctrl2
    if (blank == "row") v <- v - bl else if (blank == "mean") v <- v - mean(bl)
    if (grp_name == "CTRL" && drop_high_ctrl) v <- v[v != max(v)]
    v
  }
  spec <- list(S1 = list("25K_read2", 331.0, TRUE, "none"), S2 = list("25K_read2", NULL, TRUE, "none"), S3 = list("25K_read2", 331.0, FALSE, "none"), S4 = list("25K_read2", NULL, FALSE, "none"),
               S5 = list("25K_read2", NULL, TRUE, "row"), S6 = list("25K_read2", NULL, FALSE, "row"), S7 = list("25K_read2", NULL, FALSE, "mean"), S8 = list("25K_read1", NULL, FALSE, "none"),
               S9 = list("25K_read1", NULL, TRUE, "none"), S10 = list("50K_read1", NULL, FALSE, "none"), S11 = list("50K_read2", NULL, FALSE, "none"), S12 = list("50K_read2", NULL, FALSE, "row"))
  nb_rows <- list(); gate_nb <- numeric(0)
  for (s in names(spec)) {
    sp <- spec[[s]]; blk <- sp[[1]]; c2 <- sp[[2]]; drop <- sp[[3]]; bl <- sp[[4]]
    cc <- vals(blk, "CTRL", c2, drop, bl); t <- list(DAPT = vals(blk, "DAPT", NULL, FALSE, bl), `DKK-1` = vals(blk, "DKK1", NULL, FALSE, bl))
    norm <- mean(cc); r6 <- sc_rep[[s]]; g6 <- function(k) num(r6[[sc_col[[k]]]])
    gate_nb <- c(gate_nb, max(abs(mean(cc / norm) - g6("CTRL_mean")), abs(sem_(cc / norm) - g6("CTRL_sem")), abs(mean(t$DAPT / norm) - g6("DAPT_mean")), abs(sem_(t$DAPT / norm) - g6("DAPT_sem")),
                                abs(mean(t[["DKK-1"]] / norm) - g6("DKK1_mean")), abs(sem_(t[["DKK-1"]] / norm) - g6("DKK1_sem"))))
    if (min(min(cc), min(t$DAPT), min(t[["DKK-1"]])) <= 0) next
    for (nm in c("DAPT", "DKK-1")) {
      note <- paste0("CTRL values used: ", paste(sprintf("%g", cc), collapse = ", "), if (s == "S1") "; published analysis = S1" else "")
      rc <- ratio_rec("2-NBDG", paste0("2-NBDG fluorescence, scenario ", s), paste(nm, "vs CTRL"), t[[nm]], cc, "Kruskal-Wallis, Bonferroni Wilcoxon",
                      variant = if (s == "S1") "published (Fig. 2b)" else paste("Table S6 scenario", s), note = note)
      pre <- if (nm == "DAPT") "DAPT" else "DKK1"
      rc$trt_mean_recomputed <- mean(t[[nm]] / norm); rc$trt_sem_recomputed <- sem_(t[[nm]] / norm); rc$ctrl_mean_reported <- g6("CTRL_mean"); rc$ctrl_sem_reported <- g6("CTRL_sem")
      rc$trt_mean_reported <- g6(paste0(pre, "_mean")); rc$trt_sem_reported <- g6(paste0(pre, "_sem"))
      nb_rows[[length(nb_rows) + 1]] <- rc
    }
  }
  GATE[["NBDG_scenarios_max_abs_diff_mean_sem_vs_Table_S6"]] <<- max(gate_nb); GATE[["NBDG_scenarios_checked"]] <<- length(gate_nb)
  pw2 <- stat_pairwise(af4_rows(AF4, "Fig_2b_stats"))              # published S1 p values (Fig_2b_stats)
  for (i in seq_along(nb_rows)) if (startsWith(nb_rows[[i]]$variant, "published")) {
    nm <- strsplit(nb_rows[[i]]$contrast, " ", fixed = TRUE)[[1]][1]; nb_rows[[i]]$p_raw_manuscript <- pw2[[nm]][1]; nb_rows[[i]]$p_adj_bonf_manuscript <- pw2[[nm]][2]
  }
  recs <- c(recs, nb_rows)
  ## ---- C-peptide
  cp <- af4_rows(AF4, "Fig_3b_Cpeptide_all_conditions"); h <- which(vapply(cp, function(r) length(r) > 0 && is_str(r[[1]]) && r[[1]] == "group", TRUE))[1]
  cp_hdr <- hdr_names(cp[[h]]); cc3 <- setNames(seq_along(cp_hdr), cp_hdr); W <- list(); gn3 <- c("CTRL", "DAPT", "DKK-1")
  for (r in cp[(h + 1):length(cp)]) if (length(r) && is_str(r[[1]]) && r[[1]] %in% gn3) {
    key <- paste(r[[1]], sprintf("%g", num(r[[cc3[["condition_mM"]]]])))
    W[[key]] <- c(W[[key]], list(list(dod = num(r[[cc3[["dOD"]]]]), per = num(r[[cc3[["ng_per_ug_groupmeanDNA"]]]]), od590 = num(r[[cc3[["OD590"]]]]), od450 = num(r[[cc3[["OD450"]]]]), gdna = num(r[[cc3[["group_mean_DNA_ug"]]]]))))
  }
  Wk <- function(g, cond) W[[paste(g, sprintf("%g", cond))]]
  ## spread of the four wells per group at 5.5 mM (all wells above the lowest standard): coefficient of variation of dOD and the OD450 readings
  extra[["cpeptide_5.5mM_cv_dod"]] <- setNames(lapply(gn3, function(g) { v <- vapply(Wk(g, 5.5), function(w) w$dod, 1); sd(v) / mean(v) }), gn3)
  extra[["cpeptide_5.5mM_od450"]] <- setNames(lapply(gn3, function(g) vapply(Wk(g, 5.5), function(w) w$od450, 1)), gn3)
  stats_rows <- af4_rows(AF4, "Fig_3b_stats")
  pw3 <- list(per = stat_pairwise(stats_rows), sig = stat_pairwise(stats_rows, "Signal (dOD)")); blocks3 <- list(per = stat_means(stats_rows), sig = stat_means(stats_rows, "Signal (dOD)"))
  supported <- list(); below <- list()
  for (cond in c(0.0, 2.8, 5.5, 22.0, 44.0)) {
    b <- setNames(lapply(gn3, function(g) as.integer(sum(vapply(Wk(g, cond), function(w) w$dod < 0.044 - 1e-9, TRUE)))), gn3)
    below[[.json_num(cond)]] <- b; if (all(unlist(b) == 0L)) supported[[length(supported) + 1]] <- cond
  }
  extra[["cpeptide_below_floor"]] <- below; extra[["cpeptide_supported_conditions"]] <- supported
  gd <- five                                                       # denominators: five DNA values per group (Table S7); the group-mean DNA is the workbook's denominator
  for (cond in supported) {
    for (kl in list(list("per", sprintf("C-peptide per ug group-mean DNA (workbook units), %g mM glucose", cond)), list("dod", sprintf("ELISA signal dOD, %g mM glucose", cond)))) {
      key <- kl[[1]]; lab <- kl[[2]]; ctl <- vapply(Wk("CTRL", cond), function(w) w[[key]], 1)
      for (nm in c("DAPT", "DKK-1")) {
        trt <- vapply(Wk(nm, cond), function(w) w[[key]], 1); pp <- pw3[[if (key == "per") "per" else "sig"]][[nm]]
        rc <- ratio_rec("C-peptide", lab, paste(nm, "vs CTRL"), trt, ctl, "Kruskal-Wallis, Bonferroni Wilcoxon", p_raw = if (is.null(pp)) NULL else pp[1], p_adj = if (is.null(pp)) NULL else pp[2],
                        note = if (key == "per") "wells only; the DNA denominator is a constant" else "")
        rc$trt_mean_recomputed <- mean(trt); rc$trt_sem_recomputed <- sem_(trt); rc$ctrl_mean_reported <- mean(ctl); rc$ctrl_sem_reported <- sem_(ctl)
        recs[[length(recs) + 1]] <- rc
        if (key == "per") {
          ## interval that also carries the uncertainty of the DNA denominator (mean of five DNA values per group), delta method on the log scale
          lt <- log(trt); lc <- log(ctl); comp <- c(var(lt) / 4, var(lc) / 4, var(log(gd[[nm]])) / 5, var(log(gd[["CTRL"]])) / 5); dfc <- c(3, 3, 4, 4)
          se <- sqrt(sum(comp)); dfw <- se^4 / sum(comp^2 / dfc); dd <- (mean(lt) - mean(lc)) / log(2); t975 <- qt(0.975, dfw); lo <- dd - t975 * se / log(2); hi <- dd + t975 * se / log(2)
          rc2 <- rc; rc2$variant <- "with DNA-denominator uncertainty"; rc2$log2fc_geo <- dd; rc2$ci95_lo <- lo; rc2$ci95_hi <- hi; rc2$df <- dfw; rc2$fc_geo <- 2^dd; rc2$fc_geo_lo <- 2^lo; rc2$fc_geo_hi <- 2^hi
          rc2$ci_excludes_0 <- (lo > 0 || hi < 0); rc2$p_welch_log2 <- ""; rc2$p_raw_manuscript <- ""; rc2$p_adj_bonf_manuscript <- ""
          rc2$note <- paste0("post hoc: Welch-Satterthwaite combination of the log variances of the four wells per group and of the five DNA values per group (the group means used as denominators); ",
                             "the DNA values come from 5 different wells, one per glucose condition")
          recs[[length(recs) + 1]] <- rc2
          ## sensitivity: CTRL well with the optical artefact left out (OD590 0.105 against 0.036-0.040), post hoc
          if (cond == 5.5) {
            keep <- vapply(Wk("CTRL", cond), function(w) w$od590 < 0.08, TRUE)
            if (sum(!keep) == 1) {
              rc3 <- ratio_rec("C-peptide", lab, paste(nm, "vs CTRL"), trt, ctl[keep], "not tested by the authors", variant = "sensitivity: CTRL artefact well left out", note = "post hoc; the authors kept the well (Table S7 note); CTRL n = 3")
              rc3$ctrl_mean_reported <- ""; rc3$ctrl_sem_reported <- ""; recs[[length(recs) + 1]] <- rc3
            }
          }
        }
      }
    }
  }
  ## gate: 5.5 mM means and SEMs of the four wells against Fig_3b_stats (per ug group-mean DNA and dOD)
  dm <- numeric(0)
  for (sk in list(c("per", "per"), c("sig", "dod"))) for (g in gn3) { v <- vapply(Wk(g, 5.5), function(w) w[[sk[2]]], 1); dm <- c(dm, abs(mean(v) - blocks3[[sk[1]]][[g]]$mean), abs(sem_(v) - blocks3[[sk[1]]][[g]]$sem)) }
  GATE[["Cpeptide_5.5mM_mean_sem_max_abs_diff_vs_Fig_3b_stats"]] <<- max(dm)
  list(recs = recs, extra = extra, W = W, below = below, pw3 = pw3)
}

## ---------------------------------------------------------------- effect_sizes.md
write_md <- function(af3, genes, res, recs, var, oth, extra, below, path, script_label) {
  L <- character(0); A <- function(s) L <<- c(L, s); gi <- function(g) paste0("*", g, "*"); DAG <- "†"
  A("# Effect sizes with 95% intervals, across-gene BH q and Gapdh sensitivity (WP5a)")
  A("")
  A(sprintf("Produced by `%s` from Additional files 3 and 4 and `gapdh_exact.json`. Every interval below is **%s**.", script_label, LABEL))
  A("")
  A("## How to read the tables")
  A("")
  A("* **Values used.** For each gene, the 12 relative expression values (2^-dCt, one per well-level replicate; Additional file 4, sheet `Fig_4c_Raw_expression_data`, full precision) divided by the arithmetic mean of the four CTRL values of that gene. These are the replicate fold changes R1-R4 of Additional file 3 (reproduced to rounding, see README gate). Rows 1-4 are CTRL, 5-8 DAPT, 9-12 DKK-1. DAPT values of the five 26 August genes are as reported unless a variant is named.")
  A(paste0("* **log2 FC [95% CI].** Each replicate fold change is log2-transformed. Effect = mean log2 (treated) minus mean log2 (CTRL), that is log2 of the ratio of geometric means. Interval = Welch t-interval (n = 4 vs 4; Welch-Satterthwaite df, range in the summary below), so no equal-variance assumption is made. Not adjusted for the three pairwise contrasts per gene or for the 17 genes. ",
           "The same values rounded to the 3 decimals printed in Additional file 3 move the interval limits by at most ", sprintf("%.4f", GATE$ci_limits_full_vs_printed_replicates_max_abs_diff_log2), " log2 units."))
  A("* **Mean +/- SEM.** As printed in Additional file 3 (fold change of the arithmetic mean relative to the CTRL arithmetic mean; SEM of the four replicate fold changes). The manuscript fold change is the ratio of arithmetic means, so it differs slightly from 2^(log2 FC).")
  A("* **p (Bonf.).** The manuscript's Bonferroni-corrected pairwise p (Additional file 3): pooled-variance t-tests after one-way ANOVA, or Wilcoxon rank-sum tests after Kruskal-Wallis for Wnt5a and Lef1 (marked KW; the smallest attainable adjusted p is 0.086).")
  G1 <- GATE$from_printed_replicates
  A(sprintf("* **BH q.** Benjamini-Hochberg adjustment over the 17 raw pairwise p of the same contrast (CTRL vs DAPT, or CTRL vs DKK-1). The raw p are those of the manuscript's tests (pooled-variance t-test, or Wilcoxon rank-sum for the two KW genes) recomputed from the 12 full-precision values. The p printed in Additional file 3 were computed from its 3-decimal replicate values (all %d of %d printed pairwise p are reproduced from those values); the full-precision p differ from them by a median of %.1f%% and at most %.0f%% (relative), and q changes by at most %.4f with no change in whether q < 0.05. q is reported beside the manuscript's tests and does not replace them; it answers a different question (17 genes) from the Bonferroni p (3 contrasts per gene). q < 0.05 is marked in bold.",
            G1$pairwise_p_matching_printed, G1$pairwise_p_checked, 100 * GATE$p_raw_full_vs_printed_replicates_median_rel_diff, 100 * GATE$p_raw_full_vs_printed_replicates_max_rel_diff,
            max(GATE$bh_q_max_abs_diff_full_vs_printed_replicates_DAPT, GATE[["bh_q_max_abs_diff_full_vs_printed_replicates_DKK-1"]])))
  A("* **Excludes 0 (pooled).** Sensitivity: whether the interval excludes 0 when the within-group variance is pooled over the three groups (df = 9), the model behind the manuscript's pairwise tests.")
  A(paste0("* **", DAG, "** The DAPT value depends on the Gapdh wells of the 26 August run (Tcf7l2, Kcnj11, Cacna1d, Cacna1c, Ptbp1); see the section on Gapdh variants."))
  A("")
  ## summary counts
  A("## Summary of the 17-gene panel")
  A("")
  dfs <- vapply(recs, function(x) x$df, 1); widths <- vapply(recs, function(x) x$ci95_hi - x$ci95_lo, 1)
  A(sprintf("Welch-Satterthwaite df range %.1f to %.1f; median width of the 95%% interval %.2f log2 units (a %.1f-fold range). Counts below are out of 17 genes per contrast.", min(dfs), max(dfs), median(widths), 2^median(widths)))
  A("")
  A("| Contrast | Manuscript Bonferroni p < 0.05 | BH q < 0.05 | Welch interval excludes 0 | Pooled interval excludes 0 | Genes where Bonferroni and BH disagree | Genes where Welch interval and Bonferroni disagree |")
  A("|---|---|---|---|---|---|---|")
  join_d <- function(d) if (!length(d)) "none" else paste(vapply(d, function(x) paste0("*", sub(" (", "* (", x, fixed = TRUE)), ""), collapse = ", ")
  for (key in CONTRASTS) {
    rr <- Filter(function(x) x$contrast == paste(key, "vs CTRL"), recs); fe <- vapply(rr, function(x) x$feature, "")
    bon <- setNames(vapply(rr, function(x) as.numeric(x$p_adj_bonf_manuscript) < 0.05, TRUE), fe); qq <- setNames(vapply(rr, function(x) x$bh_q_17genes < 0.05, TRUE), fe)
    cw <- setNames(vapply(rr, function(x) isTRUE(x$ci_excludes_0), TRUE), fe); cpp <- setNames(vapply(rr, function(x) isTRUE(x$ci_excludes_0_pooled3), TRUE), fe)
    d1 <- vapply(fe[bon != qq], function(g) sprintf("%s (%s)", g, if (qq[[g]]) "q only" else "Bonferroni only"), ""); d2 <- vapply(fe[bon != cw], function(g) sprintf("%s (%s)", g, if (cw[[g]]) "interval only" else "Bonferroni only"), "")
    A(sprintf("| %s vs CTRL | %d | %d | %d | %d | %s | %s |", key, sum(bon), sum(qq), sum(cw), sum(cpp), join_d(unname(d1)), join_d(unname(d2))))
  }
  A("")
  A("Bonferroni p < 0.05 is the manuscript's own criterion. The intervals are unadjusted, so a few genes whose Bonferroni p is just above 0.05 have intervals that exclude 0.")
  A("")
  for (key in CONTRASTS) {
    A(sprintf("## RT-qPCR: %s vs CTRL (17 genes)", key)); A(""); A(sprintf("Interval label: %s.", LABEL)); A("")
    A(paste0("| Gene | Test | CTRL mean +/- SEM | ", key, " mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.), manuscript | BH q | Excludes 0 (Welch / pooled) |"))
    A("|---|---|---|---|---|---|---|---|---|")
    for (r in Filter(function(x) x$contrast == paste(key, "vs CTRL"), recs)) {
      g <- r$feature; dag <- if (r$flag_26Aug_Gapdh) DAG else ""; qv <- r$bh_q_17genes; qs <- if (qv < 0.05) sprintf("**%s**", fmtp(qv)) else fmtp(qv)
      A(sprintf("| %s%s | %s | %s +/- %s | %s +/- %s | %s [%s, %s] | %s [%s, %s] | %s | %s | %s / %s |", gi(g), dag, if (startsWith(r$test_in_manuscript, "Kruskal")) "KW" else "ANOVA", sprintf("%g", r$ctrl_mean_reported),
                sprintf("%g", r$ctrl_sem_reported), sprintf("%g", r$trt_mean_reported), sprintf("%g", r$trt_sem_reported), fmt(r$log2fc_geo, 2), fmt(r$ci95_lo, 2), fmt(r$ci95_hi, 2), fmt(r$fc_geo, 2), fmt(r$fc_geo_lo, 2), fmt(r$fc_geo_hi, 2),
                r$p_adj_bonf_manuscript, qs, if (r$ci_excludes_0) "yes" else "no", if (r$ci_excludes_0_pooled3) "yes" else "no"))
    }
    A("")
  }
  ## Gapdh variants
  A("## The five genes whose DAPT values depend on the 26 August Gapdh wells")
  A("")
  A(sprintf("DAPT vs CTRL; CTRL and DKK-1 values are identical in all three variants. Reported = Gapdh values as analyzed in the authors' calculation sheet (C10 and D10 entered twice, A10 and B10 set aside). C0 = reference is the mean of the two valid wells (Cq 22.68) for all four DAPT replicates. C1 = reference is the mean of all four wells (Cq 23.62), failed wells included (worst case). Same interval method as above, label: %s. BH q is recomputed over the 17 genes of the DAPT contrast with the five p replaced by the variant values.", LABEL))
  A("")
  A("| Gene | Variant | Mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.) | BH q | Sign of log2 FC |")
  A("|---|---|---|---|---|---|---|---|")
  vlab <- c(reported = "reported", C0 = "C0 (two valid wells)", C1 = "C1 (all four wells)")
  for (v in var) {
    qs <- if (v$bh_q < 0.05) sprintf("**%s**", fmtp(v$bh_q)) else fmtp(v$bh_q); sgn <- if (v$log2fc < 0) "down" else "up"
    A(sprintf("| %s | %s | %.3f +/- %.3f | %s [%s, %s] | %s [%s, %s] | %s | %s | %s |", gi(v$gene), vlab[[v$variant]], v$mean, v$sem, fmt(v$log2fc, 2), fmt(v$lo, 2), fmt(v$hi, 2), fmt(2^v$log2fc, 2), fmt(2^v$lo, 2), fmt(2^v$hi, 2), fmtp(v$p_bonf), qs, sgn))
  }
  A("")
  byvar <- function(nm) { o <- Filter(function(v) v$variant == nm, var); setNames(o, vapply(o, function(v) v$gene, "")) }
  vr <- byvar("reported"); v0 <- byvar("C0"); v1 <- byvar("C1")
  flip <- G5[vapply(G5, function(g) (v1[[g]]$log2fc > 0) != (vr[[g]]$log2fc > 0), TRUE)]
  inc0 <- function(vv) G5[vapply(G5, function(g) vv[[g]]$lo <= 0 && 0 <= vv[[g]]$hi, TRUE)]
  jn <- function(x) if (length(x)) paste0("*", x, "*", collapse = ", ") else "none"
  A(sprintf("The \"reported\" rows are the main-table values; all p and q in this table are recomputed from the full-precision values (pooled-variance t-tests, df = 9), so their third digit can differ from the printed p of Additional file 3. Fold changes and Bonferroni p of all three variants match `gapdh_exact.json` (max absolute difference of fold changes %s, max relative difference of p %s). With C1 the sign of the log2 FC reverses for %s. The 95%% interval includes 0 under C1 for %s, under C0 for %s and under the reported values for %s.",
            .py_e2(GATE$gapdh_variants_max_abs_diff_fc_vs_gapdh_exact_json), .py_e2(GATE$gapdh_variants_max_rel_diff_bonf_p_vs_gapdh_exact_json), jn(flip), jn(inc0(v1)), jn(inc0(v0)), jn(inc0(vr))))
  A("")
  ## other readouts
  A("## DNA content, 2-NBDG and C-peptide per ug DNA")
  A("")
  A(sprintf("Interval label: %s. Intervals are Welch t-intervals on the log2 scale of the replicate values, as above.", LABEL))
  A("")
  A("| Readout | Contrast | Variant | n (CTRL / treated) | Mean ratio (arith.) | log2 ratio [95% CI] | 2^log2 ratio [95% CI] | p (Bonf.), manuscript |")
  A("|---|---|---|---|---|---|---|---|")
  for (r in oth) {
    if (r$readout == "2-NBDG" && !startsWith(r$variant, "published")) next
    pj <- r$p_adj_bonf_manuscript; pj <- if (identical(pj, "")) "" else if (is.numeric(pj)) fmtp(pj) else as.character(pj)
    A(sprintf("| %s | %s | %s | %d / %d | %s | %s [%s, %s] | %s [%s, %s] | %s |", r$feature, r$contrast, r$variant, r$n_ctrl, r$n_trt, fmt(r$fc_ratio_arith, 3), fmt(r$log2fc_geo, 2), fmt(r$ci95_lo, 2), fmt(r$ci95_hi, 2),
              fmt(r$fc_geo, 2), fmt(r$fc_geo_lo, 2), fmt(r$fc_geo_hi, 2), pj))
  }
  A("")
  pick <- function(readout, feat, contrast, variant) { for (x in oth) if (x$readout == readout && grepl(feat, x$feature, fixed = TRUE) && x$contrast == contrast && x$variant == variant) return(x); stop("pick: no record") }
  dna_d <- pick("DNA content", "wells at 2.8", "DAPT vs CTRL", "as reported"); dna_k <- pick("DNA content", "wells at 2.8", "DKK-1 vs CTRL", "as reported")
  dna5_d <- pick("DNA content", "five wells", "DAPT vs CTRL", "sensitivity: five wells"); dna5_k <- pick("DNA content", "five wells", "DKK-1 vs CTRL", "sensitivity: five wells")
  nbs <- list(); for (x in oth) if (x$readout == "2-NBDG") { k <- strsplit(x$contrast, " ", fixed = TRUE)[[1]][1]; nbs[[k]] <- c(nbs[[k]], list(x)) }
  n_ex <- sum(vapply(unlist(nbs, recursive = FALSE), function(x) isTRUE(x$ci_excludes_0), TRUE)); n_all <- sum(vapply(nbs, length, 1L))
  cp_d <- pick("C-peptide", "per ug", "DAPT vs CTRL", "as reported"); cp_k <- pick("C-peptide", "per ug", "DKK-1 vs CTRL", "as reported")
  cs_d <- pick("C-peptide", "dOD", "DAPT vs CTRL", "as reported"); cs_k <- pick("C-peptide", "dOD", "DKK-1 vs CTRL", "as reported")
  l2 <- function(k) vapply(nbs[[k]], function(x) x$log2fc_geo, 1)
  A("Notes on this table.")
  A("")
  A(sprintf("* **DNA content.** The four values per group are wells that differed in their 20-minute glucose exposure (Methods), not four replicates of one condition, and the CTRL wells range from %s to %s ug. The lower interval limit is %s (DAPT) and %s (DKK-1) log2 units above 0, that is, barely above no change. Adding the fifth (basal) well, which the authors left out for an unrecorded reason, gives %s and %s log2 units with lower limits %s and %s (post hoc sensitivity, not a replacement).",
            sprintf("%.3f", min(extra$dna_four$CTRL)), sprintf("%.3f", max(extra$dna_four$CTRL)), fmt(dna_d$ci95_lo, 2), fmt(dna_k$ci95_lo, 2), fmt(dna5_d$log2fc_geo, 2), fmt(dna5_k$log2fc_geo, 2), fmt(dna5_d$ci95_lo, 2), fmt(dna5_k$ci95_lo, 2)))
  A(sprintf("* **2-NBDG: not supported.** Across the 12 scenarios of Table S6, %d of %d intervals exclude 0, the log2 ratio of DAPT runs from %s to %s and that of DKK-1 from %s to %s, so even the direction depends on the analysis choice. The published analysis (S1) uses the CTRL values 195, 331 and 331 (the value 331 appears twice because the analysis block of the workbook holds 331 for CTRL replicate 2, where the reading block holds 184; the reason is not recorded). Only S1 is shown in the table above; all twelve are in the scenario table below and in the CSV.",
            n_ex, n_all, fmt(min(l2("DAPT")), 2), fmt(max(l2("DAPT")), 2), fmt(min(l2("DKK-1")), 2), fmt(max(l2("DKK-1")), 2)))
  A(sprintf("* **C-peptide at 5.5 mM glucose.** The raw signal (dOD) is higher than CTRL in both treated groups (mean ratio %s and %s) whereas the value per ug of group-mean DNA is lower (%s and %s): the sign of the difference is set by the DNA denominator, which is 2.3 and 2.6 times larger in the treated groups (Table S7). Three of the four CTRL wells read 0.047 to 0.048, only 0.003 to 0.004 above the lowest standard (0.044), and the fourth (0.083) carries the optical artefact noted in Table S7. The \"wells only\" interval treats the DNA denominator as an exact constant; the row \"with DNA-denominator uncertainty\" adds the variance of the five DNA values per group (one well per glucose condition) and is much wider. The narrow intervals of the row \"CTRL artefact well left out\" come from near-identical readings (the readings are 0.001 apart in dOD) and describe the spread of the readings, not the uncertainty of the biology.",
            fmt(cs_d$fc_ratio_arith, 2), fmt(cs_k$fc_ratio_arith, 2), fmt(cp_d$fc_ratio_arith, 2), fmt(cp_k$fc_ratio_arith, 2)))
  A("")
  A("### 2-NBDG: range over the twelve analysis scenarios of Table S6")
  A("")
  A("| Scenario | CTRL values used | DAPT mean ratio | DAPT log2 ratio [95% CI] | DKK-1 mean ratio | DKK-1 log2 ratio [95% CI] |")
  A("|---|---|---|---|---|---|")
  scn <- list(); for (r in oth) if (r$readout == "2-NBDG") { s <- strsplit(r$feature, "scenario ", fixed = TRUE)[[1]][2]; k <- strsplit(r$contrast, " ", fixed = TRUE)[[1]][1]; scn[[s]][[k]] <- r }
  for (s in names(scn)) {
    dd <- scn[[s]]$DAPT; kk <- scn[[s]][["DKK-1"]]
    A(sprintf("| %s | %s | %s | %s [%s, %s] | %s | %s [%s, %s] |", s, strsplit(strsplit(dd$note, "CTRL values used: ", fixed = TRUE)[[1]][2], ";", fixed = TRUE)[[1]][1], fmt(dd$fc_ratio_arith, 2), fmt(dd$log2fc_geo, 2), fmt(dd$ci95_lo, 2), fmt(dd$ci95_hi, 2),
              fmt(kk$fc_ratio_arith, 2), fmt(kk$log2fc_geo, 2), fmt(kk$ci95_lo, 2), fmt(kk$ci95_hi, 2)))
  }
  A("")
  A("### C-peptide: which glucose conditions the data support")
  A("")
  A("A comparison is estimated only where no well of any of the three groups lies below the lowest standard (dOD < 0.044). Wells below the lowest standard per group (CTRL / DAPT / DKK-1):")
  A("")
  A("| Glucose (mM) | Below lowest standard | Estimated here |")
  A("|---|---|---|")
  for (cond in names(below)) { d <- below[[cond]]; A(sprintf("| %s | %d / %d / %d | %s |", cond, d$CTRL, d$DAPT, d[["DKK-1"]], if (all(unlist(d) == 0L)) "yes" else "no: at least one group has wells below the lowest standard")) }
  A("")
  write_text_bytes(paste(L, collapse = .NL), path)             # Python wrote '\n'.join(L) in text mode on Windows, so the archived file has CRLF line ends
}
.py_e2 <- function(x) sprintf("%.1e", x)                        # Python '%.1e' (two-digit exponent, as the UCRT sprintf of R on Windows)

## ---------------------------------------------------------------- main
es_main <- function(out_dir = file.path(ES_HERE, "results_R", "wp5"), script_label = Sys.getenv("ES_LABEL", "effect_sizes.R")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  af3 <- load_af3(); rp <- rt_panel(af3); genes <- rp$genes; res <- rp$res
  reproduce_reported(af3, res)
  gr <- gene_rows(af3, genes, res); recs <- gr$recs; q <- gr$q; res <- gr$res
  var <- variants_g5(af3, res, q)
  ot <- other_readouts(); oth <- ot$recs; extra <- ot$extra; below <- ot$below
  allrecs <- c(recs, oth)
  cols <- c("readout", "feature", "category", "contrast", "variant", "flag_26Aug_Gapdh", "test_in_manuscript", "n_ctrl", "n_trt", "ctrl_mean_reported", "ctrl_sem_reported", "trt_mean_reported", "trt_sem_reported",
            "trt_mean_recomputed", "trt_sem_recomputed", "fc_ratio_arith", "log2fc_geo", "ci95_lo", "ci95_hi", "ci_method", "df", "fc_geo", "fc_geo_lo", "fc_geo_hi", "ci95_lo_pooled3", "ci95_hi_pooled3", "df_pooled3",
            "ci_excludes_0", "ci_excludes_0_pooled3", "p_welch_log2", "p_raw_manuscript", "p_adj_bonf_manuscript", "p_raw_recomputed", "p_adj_bonf_recomputed", "bh_q_17genes", "interval_label", "note")
  vnote <- c(reported = "Gapdh reference of the 26 Aug DAPT group as analyzed (C10 and D10 entered twice)", C0 = "Gapdh = mean of the two valid wells C10 and D10 (22.68) for all four DAPT replicates",
             C1 = "Gapdh = mean of all four wells A10-D10 (23.62), failed wells included (worst case)")
  vrecs <- lapply(var, function(v) {                           # 26 August variants as additional rows
    r0 <- af3[[which(vapply(af3, function(r) r$gene, "") == v$gene)[1]]]
    list(readout = "RT-qPCR", feature = v$gene, category = r0$category, contrast = "DAPT vs CTRL", variant = if (v$variant != "reported") paste("Gapdh", v$variant) else "as reported (26 Aug Gapdh)",
         flag_26Aug_Gapdh = TRUE, test_in_manuscript = "One-way ANOVA", n_ctrl = 4L, n_trt = 4L, trt_mean_recomputed = v$mean, trt_sem_recomputed = v$sem, fc_ratio_arith = v$fc_arith, log2fc_geo = v$log2fc, ci95_lo = v$lo, ci95_hi = v$hi,
         ci_method = "Welch t, log2 scale", df = v$df, fc_geo = 2^v$log2fc, fc_geo_lo = 2^v$lo, fc_geo_hi = 2^v$hi, ci95_lo_pooled3 = v$lo3, ci95_hi_pooled3 = v$hi3, df_pooled3 = 9L,
         ci_excludes_0 = (v$lo > 0 || v$hi < 0), ci_excludes_0_pooled3 = (v$lo3 > 0 || v$hi3 < 0), p_raw_recomputed = v$p_raw, p_adj_bonf_recomputed = v$p_bonf, bh_q_17genes = v$bh_q, interval_label = LABEL, note = vnote[[v$variant]])
  })
  dictwriter_csv(c(allrecs, Filter(function(x) !startsWith(x$variant, "as reported"), vrecs)), cols, file.path(out_dir, "effect_sizes.csv"))
  write_text_bytes(gsub("\n", .NL, py_json_typed(list(gate = GATE, extra = extra, below_floor = below)), fixed = TRUE), file.path(out_dir, "gate_reproduction.json"))
  write_md(af3, genes, res, recs, var, oth, extra, below, file.path(out_dir, "effect_sizes.md"), script_label)
  cat(py_json_typed(GATE), "\n")
  invisible(list(recs = allrecs, vrecs = vrecs, gate = GATE, extra = extra, below = below))
}
