## effect_sizes_part3.R - R port of wp5_reporting/effect_sizes.py, PART 3 (Python lines 303-395): variants_g5, ratio_rec and the DNA-content block of other_readouts (dna_records).
## The 2-NBDG, C-peptide, summary and output parts (Python lines 396-725) follow in part 4. Group indices 1-based.
library(jsonlite)
## DAPT values of the five 26 August genes under the reported, C0 and C1 Gapdh references
variants_g5 <- function(af3, res, q_reported) {
  J <- fromJSON(GAPDH_JSON, simplifyVector = FALSE); HK <- as.numeric(unlist(J$hk_used)); REF <- list(reported = NULL, C0 = as.numeric(J$ref$C0), C1 = as.numeric(J$ref$C1)); out <- list(); chk <- list()
  gn <- vapply(af3, function(r) r$gene, "")
  for (var in c("reported", "C0", "C1")) {
    praw_all <- setNames(lapply(gn, function(g) res[[g]]$praw[["DAPT"]]), gn); recs <- list()
    for (g in G5) {
      R <- res[[g]]; x <- R$x
      if (!is.null(REF[[var]])) x[GI$DAPT] <- x[GI$DAPT] * 2.0^(REF[[var]] - HK)
      fc <- x / mean(x[GI$CTRL]); groups <- list(fc[GI$CTRL], fc[GI$DAPT], fc[GI[["DKK-1"]]]); L <- lapply(groups, log2)
      w <- welch(L[[2]], L[[1]]); p3 <- pooled3(L, 1, 2); praw <- pairwise_p_anova(groups, 1, 2); praw_all[[g]] <- praw
      recs[[length(recs) + 1]] <- list(gene = g, variant = var, mean = mean(groups[[2]]), sem = sem_(groups[[2]]), fc_arith = mean(groups[[2]]) / mean(groups[[1]]), log2fc = w[1], lo = w[2], hi = w[3], df = w[4], lo3 = p3[2], hi3 = p3[3], p_raw = praw, p_bonf = min(1.0, 3 * praw))
      ref <- J[[var]]$per_gene[[g]]; rc <- recs[[length(recs)]]; chk[[length(chk) + 1]] <- c(abs(rc$fc_arith - ref$fc_dapt), abs(rc$p_bonf - ref[["CTRL-DAPT"]]) / max(ref[["CTRL-DAPT"]], 1e-300))
    }
    qv <- setNames(bh_(unlist(praw_all)), names(praw_all))
    for (i in seq_along(recs)) { recs[[i]]$bh_q <- unname(qv[[recs[[i]]$gene]]); recs[[i]]["sign_vs_reported"] <- list(NULL) }
    out <- c(out, recs)
  }
  s1c <- list(Kcnj11 = c(0.40, 0.40, 0.76), Cacna1d = c(0.76, 0.76, 1.46), Cacna1c = c(0.71, 0.72, 1.37), Ptbp1 = c(0.62, 0.62, 1.18), Tcf7l2 = c(0.78, 0.78, 1.49))   # Table S1c, DAPT fold change as reported / C0 / C1, as printed
  byv <- list(); for (r in out) byv[[paste(r$gene, r$variant)]] <- r$fc_arith
  GATE$gapdh_variants_match_Table_S1c_at_2_decimals <<- all(vapply(G5, function(g) all(vapply(1:3, function(k) abs(round(byv[[paste(g, c("reported", "C0", "C1")[k])]], 2) - s1c[[g]][k]) < 1e-9, TRUE)), TRUE))
  GATE$gapdh_variants_max_abs_diff_fc_vs_gapdh_exact_json <<- max(vapply(chk, function(c) c[1], 1)); GATE$gapdh_variants_max_rel_diff_bonf_p_vs_gapdh_exact_json <<- max(vapply(chk, function(c) c[2], 1))
  GATE$gapdh_refs <<- list(C0 = REF$C0, C1 = REF$C1, hk_used = HK)
  out
}
ratio_rec <- function(readout, feature, contrast, a_trt, a_ctl, test, p_raw = NULL, p_adj = NULL, note = "", variant = "as reported") {
  w <- welch(log2(a_trt), log2(a_ctl))
  list(readout = readout, feature = feature, category = "", contrast = contrast, variant = variant, flag_26Aug_Gapdh = FALSE, test_in_manuscript = test, n_ctrl = length(a_ctl), n_trt = length(a_trt),
       ctrl_mean_reported = "", ctrl_sem_reported = "", trt_mean_reported = "", trt_sem_reported = "", trt_mean_recomputed = "", trt_sem_recomputed = "", fc_ratio_arith = mean(a_trt) / mean(a_ctl),
       log2fc_geo = w[1], ci95_lo = w[2], ci95_hi = w[3], ci_method = "Welch t, log2 scale", df = w[4], fc_geo = 2^w[1], fc_geo_lo = 2^w[2], fc_geo_hi = 2^w[3], ci95_lo_pooled3 = "", ci95_hi_pooled3 = "",
       df_pooled3 = "", ci_excludes_0 = (w[2] > 0 || w[3] < 0), ci_excludes_0_pooled3 = "", p_welch_log2 = w[5], p_raw_manuscript = if (is.null(p_raw)) "" else p_raw, p_adj_bonf_manuscript = if (is.null(p_adj)) "" else p_adj,
       p_raw_recomputed = "", bh_q_17genes = "", interval_label = LABEL, note = note)
}
## DNA content: the four values per group listed in Table S7 (wells at 2.8, 5.5, 22 and 44 mM); basal well = sensitivity
dna_records <- function() {
  recs <- list(); extra <- list(); s7 <- af4_rows(AF4, "SuppTable7_Cpeptide_all_gluc"); h <- which(vapply(s7, function(r) length(r) && is_str(r[[1]]) && r[[1]] == "group", TRUE))[1]
  hdr <- hdr_names(s7[[h]]); col <- setNames(seq_along(hdr), hdr); gnames <- c("CTRL", "DAPT", "DKK-1"); dna <- setNames(lapply(gnames, function(.) numeric(0)), gnames)
  for (r in s7[(h + 1):length(s7)]) if (length(r) && is_str(r[[1]]) && r[[1]] %in% gnames) dna[[r[[1]]]][sprintf("%g", as.numeric(r[[col[["condition_mM"]]]]))] <- as.numeric(r[[col[["well_DNA_ug"]]]])
  four <- lapply(dna, function(v) unname(v[sprintf("%g", c(2.8, 5.5, 22.0, 44.0))])); five <- lapply(dna, function(v) unname(v[sprintf("%g", c(0.0, 2.8, 5.5, 22.0, 44.0))]))
  extra$dna_four <- four; extra$dna_five <- five
  rows <- af4_rows(AF4, "Fig_1c_stats"); rep <- stat_means(rows); pw <- stat_pairwise(rows)
  GATE$DNA_mean_max_abs_diff_vs_Fig_1c_stats <<- max(vapply(gnames, function(g) abs(mean(four[[g]]) - rep[[g]]$mean), 1)); GATE$DNA_sem_max_abs_diff_vs_Fig_1c_stats <<- max(vapply(gnames, function(g) abs(sem_(four[[g]]) - rep[[g]]$sem), 1))
  grp <- lapply(gnames, function(g) four[[g]])
  GATE$DNA_bonf_p_max_abs_diff_vs_Fig_1c_stats <<- max(vapply(list(list(2, "DAPT"), list(3, "DKK-1")), function(kn) abs(min(1.0, 3 * pairwise_p_anova(grp, 1, kn[[1]])) - pw[[kn[[2]]]][2]), 1))
  for (kn in list(list("DAPT", 2), list("DKK-1", 3))) {
    nm <- kn[[1]]; k <- kn[[2]]
    rc <- ratio_rec("DNA content", "Total DNA (ug), wells at 2.8, 5.5, 22, 44 mM glucose", paste(nm, "vs CTRL"), grp[[k]], grp[[1]], "One-way ANOVA, Bonferroni t", p_raw = pw[[nm]][1], p_adj = pw[[nm]][2],
                    note = "four wells that differed in their 20-minute glucose exposure; the fifth (basal) well was left out by the authors for an unrecorded reason")
    L <- lapply(grp, log2); p3 <- pooled3(L, 1, k)
    rc$ctrl_mean_reported <- rep$CTRL$mean; rc$ctrl_sem_reported <- rep$CTRL$sem; rc$trt_mean_reported <- rep[[nm]]$mean; rc$trt_sem_reported <- rep[[nm]]$sem; rc$trt_mean_recomputed <- mean(grp[[k]]); rc$trt_sem_recomputed <- sem_(grp[[k]])
    rc$ci95_lo_pooled3 <- p3[2]; rc$ci95_hi_pooled3 <- p3[3]; rc$df_pooled3 <- p3[4]; rc$ci_excludes_0_pooled3 <- (p3[2] > 0 || p3[3] < 0); recs[[length(recs) + 1]] <- rc
    five_g <- lapply(gnames, function(g) five[[g]])
    recs[[length(recs) + 1]] <- ratio_rec("DNA content", "Total DNA (ug), all five wells incl. basal (sensitivity)", paste(nm, "vs CTRL"), five_g[[k]], five_g[[1]], "not tested by the authors", variant = "sensitivity: five wells", note = "post hoc sensitivity; the authors analyzed four wells")
  }
  list(recs = recs, extra = extra)
}
