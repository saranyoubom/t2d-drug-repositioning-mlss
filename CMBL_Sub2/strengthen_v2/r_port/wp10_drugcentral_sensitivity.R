## wp10_drugcentral_sensitivity.R - R port of wp10_versions/drugcentral_sensitivity.py: WP10 part B: the MLSS of the reported run (DrugCentral file of the reported analysis) versus the same MLSS computed
## with the current DrugCentral drug-target file (database dump of 25 Sep 2026, Database/drug.target.interaction.tsv.gz). The port mlss_port.R is used unchanged (rules: exact GENE match to the reported hubs,
## human target, pIC50 = largest human IC50/Ki/EC50 with relation '=', same Hes1-normalized correlation lookup, same weights); only the interaction file differs.
## Gate first: the port on the OLD file reproduces the reported scores (R_17genes/Fig_6c_MLSS_ALL_Combinations.csv, the same gate as wp3_mlss.R). Pre-specification: PRESPEC_WP10_database_versions.md (section 4).
## Usage: Rscript wp10_drugcentral_sensitivity.R   (reads only; writes <OUT>/mlss_new_all_combinations.csv, mlss_comparison.csv, pic50_changes.csv, drugcentral_summary.json)
## Environment: WP10_OUT (default results_R/wp10), MLSS_RDIR (R_17genes folder), WP_PROJECT (project root holding Database/).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- HERE
source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "mlss_port.R")); suppressMessages({ library(jsonlite); library(data.table) })
RDIR <- Sys.getenv("MLSS_RDIR", normalizePath(file.path(R_PORT, "..", "..", "R_17genes"))); PROJECT <- Sys.getenv("WP_PROJECT", normalizePath(file.path(R_PORT, "..", "..", "..")))
NEW <- file.path(PROJECT, "Database", "drug.target.interaction.tsv.gz"); RES <- Sys.getenv("WP10_OUT", file.path(R_PORT, "results_R", "wp10")); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
read_new <- function(path) {            # dump of 25 Sep 2026: tab separated, text as is, ACT_VALUE numeric
  dc <- as.data.frame(fread(text = readLines(gzfile(path, "rt"), encoding = "UTF-8", warn = FALSE), sep = "\t", quote = "\"", colClasses = "character", na.strings = "", data.table = FALSE, header = TRUE))
  names(dc)[1] <- sub("^﻿", "", names(dc)[1]); dc$ACT_VALUE <- suppressWarnings(as.numeric(dc$ACT_VALUE)); dc$G <- toupper(dc$GENE); dc
}
pk <- function(a, b) mlss_pkey(a, b)
main <- function() {
  hb <- mlss_hubs_from_centrality(file.path(RDIR, "Fig_6b_Gene_Centrality_Summary.csv"), 5); hubs <- hb$hubs; corr <- mlss_corr_from_matrix_csv(file.path(RDIR, "Fig_4e_Hes1_significant_correlations_only.csv"))
  old <- mlss_read_drugcentral(file.path(RDIR, "drug.target.interaction.csv")); new <- read_new(NEW)
  R_CSV <- read.csv(file.path(RDIR, "Fig_6c_MLSS_ALL_Combinations.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  P0 <- mlss_run_reported(old, hubs, corr); gate_ok <- nrow(P0) == nrow(R_CSV) && max(abs(P0[["MLSS_v4.0"]] - R_CSV[order(R_CSV$Rank), "MLSS_v4.0"])) < 1e-9
  cat("GATE on the old file:", if (gate_ok) "PASS" else "FAIL", "\n"); if (!gate_ok) quit(status = 1)
  P1 <- mlss_run_reported(new, hubs, corr)
  out <- list(hubs = as.list(hubs), rows_old = nrow(old), rows_new = nrow(new), distinct_drugs_old = length(unique(old$DRUG_NAME)), distinct_drugs_new = length(unique(new$DRUG_NAME)))
  t0 <- mlss_targets_from(old, hubs); t1 <- mlss_targets_from(new, hubs); out$mlss_drugs_old <- length(t0); out$mlss_drugs_new <- length(t1)
  out$drugs_added <- as.list(py_sort(setdiff(names(t1), names(t0)))); out$drugs_removed <- as.list(py_sort(setdiff(names(t0), names(t1))))
  per_hub <- function(t) { r <- lapply(hubs, function(g) sum(vapply(t, function(gs) g %in% gs, TRUE))); names(r) <- hubs; r }
  out$drugs_per_hub_old <- lapply(per_hub(t0), as.integer); out$drugs_per_hub_new <- lapply(per_hub(t1), as.integer)
  out$combinations_old <- nrow(P0); out$combinations_new <- nrow(P1); out$cross_target_old <- as.integer(sum(P0$Complementarity_C == 1)); out$cross_target_new <- as.integer(sum(P1$Complementarity_C == 1))
  P0$key <- pk2 <- mlss_pkeys(P0$Drug_A, P0$Drug_B); P1$key <- mlss_pkeys(P1$Drug_A, P1$Drug_B); r0 <- P0; rownames(r0) <- r0$key; r1 <- P1; rownames(r1) <- r1$key
  shared <- r0$key[r0$key %in% r1$key]; out$combinations_shared <- as.integer(length(shared))
  s0 <- r0[shared, "MLSS_v4.0"]; s1 <- r1[shared, "MLSS_v4.0"]; k0 <- r0[shared, "Rank"]; k1 <- r1[shared, "Rank"]
  out$spearman_scores_shared <- round(cor(s0, s1, method = "spearman"), 4); out$spearman_ranks_shared <- round(cor(k0, k1, method = "spearman"), 4)
  out$score_changed_shared <- as.integer(sum(abs(s0 - s1) > 1e-9)); out$max_abs_score_change_shared <- round(max(abs(s0 - s1)), 4)
  top0 <- P0$key[order(P0$Rank)][1:15]; top1 <- P1$key[order(P1$Rank)][1:15]; out$top15_old_still_in_top15_new <- as.integer(sum(top0 %in% top1))
  out$top15_old_ranks_in_new <- lapply(top0, function(k) list(r0[k, "Combination"], as.integer(r0[k, "Rank"]), if (k %in% r1$key) as.integer(r1[k, "Rank"]) else NULL))
  out$top15_new <- lapply(top1, function(k) list(as.integer(r1[k, "Rank"]), r1[k, "Combination"], round(r1[k, "MLSS_v4.0"], 4)))
  sc1 <- round(P1[order(P1$Rank), "MLSS_v4.0"], 9); sc0 <- round(P0[order(P0$Rank), "MLSS_v4.0"], 9); out$tie_straddles_rank15_new <- sc1[15] == sc1[16]; out$tie_straddles_rank15_old <- sc0[15] == sc0[16]
  nd <- unlist(out$drugs_added); mask <- P1$Drug_A %in% nd | P1$Drug_B %in% nd; rows_new <- if (length(nd)) P1[mask, , drop = FALSE] else P1[0, , drop = FALSE]; rows_new <- rows_new[order(rows_new$Rank), , drop = FALSE]
  out$combinations_with_new_drugs <- as.integer(nrow(rows_new)); out$new_drug_combination_ranks <- lapply(head(seq_len(nrow(rows_new)), 25), function(i) list(as.integer(rows_new$Rank[i]), rows_new$Combination[i], round(rows_new[["MLSS_v4.0"]][i], 4)))
  out$new_drug_in_top15 <- as.integer(sum(!(top1 %in% r0$key)))
  # potency changes for drugs and hub targets of the reported run
  pot0 <- mlss_potency_from(old); pot1 <- mlss_potency_from(new); hubU <- toupper(hubs); ks <- sort(union(names(pot0), names(pot1)), method = "radix")
  dg <- do.call(rbind, strsplit(ks, "\t", fixed = TRUE)); ord <- order(dg[, 1], dg[, 2], method = "radix"); ks <- ks[ord]; dg <- dg[ord, , drop = FALSE]; chg <- list()
  for (i in seq_along(ks)) { d <- dg[i, 1]; g <- dg[i, 2]
    if (g %in% hubU && (d %in% names(t0) || d %in% names(t1))) { a <- if (ks[i] %in% names(pot0)) pot0[[ks[i]]] else NA_real_; b <- if (ks[i] %in% names(pot1)) pot1[[ks[i]]] else NA_real_
      if (is.na(a) || is.na(b) || abs(a - b) > 1e-9) chg[[length(chg) + 1]] <- list(drug = d, gene = g, pIC50_old = a, pIC50_new = b) } }
  out$pic50_changes <- as.integer(length(chg))
  py_write_csv_repr(if (length(chg)) do.call(rbind, lapply(chg, as.data.frame, stringsAsFactors = FALSE)) else data.frame(drug = character(0), gene = character(0), pIC50_old = numeric(0), pIC50_new = numeric(0)), file.path(RES, "pic50_changes.csv"), eol = "\r\n")
  cols <- intersect(c("Rank", "Combination", "MLSS_v4.0", "Complementarity_C", "Balance_B", "Potency_P", "Mean_pIC50", "Targets_A", "Targets_B", "Drug_A", "Drug_B"), names(P1)); py_write_csv_repr(P1[, cols], file.path(RES, "mlss_new_all_combinations.csv"), eol = "\r\n")
  union_k <- union(r0$key, r1$key); rk1 <- ifelse(union_k %in% r1$key, r1[union_k, "Rank"], 1e6); rk0 <- ifelse(union_k %in% r0$key, r0[union_k, "Rank"], 1e6); union_k <- union_k[order(rk1, rk0, method = "radix")]
  cmp <- do.call(rbind, lapply(union_k, function(k) { a <- if (k %in% r0$key) r0[k, ] else NULL; b <- if (k %in% r1$key) r1[k, ] else NULL
    data.frame(combination = (if (!is.null(b)) b else a)$Combination, rank_old = if (!is.null(a)) as.character(as.integer(a$Rank)) else "", MLSS_old = if (!is.null(a)) py_fmt_f(a[["MLSS_v4.0"]], 6) else "",
               rank_new = if (!is.null(b)) as.character(as.integer(b$Rank)) else "", MLSS_new = if (!is.null(b)) py_fmt_f(b[["MLSS_v4.0"]], 6) else "", status = if (!is.null(a) && !is.null(b)) "both" else if (!is.null(b)) "new only" else "old only", stringsAsFactors = FALSE) }))
  py_write_csv_repr(cmp, file.path(RES, "mlss_comparison.csv"), eol = "\r\n"); write_json_typed(out, file.path(RES, "drugcentral_summary.json")); cat("written to", RES, "\n")
}
if (sys.nframe() == 0L) main()
