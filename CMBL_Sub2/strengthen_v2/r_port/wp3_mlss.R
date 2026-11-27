## wp3_mlss.R - R port of wp3_mlss/mlss_port.py (gate, part 1 and 2) and wp3_mlss/mlss_scenarios.py (protocol section 5 scenarios and extras).
## Run: Rscript wp3_mlss.R   (writes gate_scores.csv, mlss_scenarios.csv, mlss_scenario_scores.csv, mlss_permutation_draws.csv, mlss_scenarios_log.txt into R_OUT_DIR = results_R/wp3)
options(warn = 1)
args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "rng", "numpy_rng.R")); source(file.path(here, "pycompat.R")); source(file.path(here, "mlss_port.R"))
suppressMessages({ library(readxl); library(dplyr); library(jsonlite) })
RDIR <- Sys.getenv("MLSS_RDIR", file.path(here, "..", "..", "R_17genes")); OUTDIR <- Sys.getenv("R_OUT_DIR", file.path(here, "results_R", "wp3"))
GAPDH_JSON <- file.path(here, "..", "..", "reanalysis_17genes", "gapdh_exact.json")
SEED <- 2026; NPERM <- 1000; LOG <- character(0); log <- function(s = "") { LOG <<- c(LOG, as.character(s)); cat(s, "\n", sep = "") }
t_start <- Sys.time()

## ------------------------------------------------------------------ inputs and the reported reference
dc <- mlss_read_drugcentral(file.path(RDIR, "drug.target.interaction.csv")); POT <- mlss_potency_from(dc)
hb <- mlss_hubs_from_centrality(file.path(RDIR, "Fig_6b_Gene_Centrality_Summary.csv"), 5); HUBS5 <- hb$hubs
CORR_REP <- mlss_corr_from_matrix_csv(file.path(RDIR, "Fig_4e_Hes1_significant_correlations_only.csv")); TG_REP <- mlss_targets_from(dc, HUBS5)
## as mlss_port.py: the authors' original matrix (SBA/Figures, 21 genes, the archived results were made with it) if present, else the 17-gene copy shipped in R_17genes
## (same values rounded to 16 digits: identical except for rank ties, which can move tie-sensitive nulls in the third decimal)
load_matrix <- function() {
  cand <- c(Sys.getenv("MLSS_INS1_MATRIX", ""), file.path(here, "..", "..", "..", "..", "SBA", "Figures", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"), file.path(RDIR, "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"))
  p <- cand[nzchar(cand) & file.exists(cand)][1]; cat("INS-1 matrix:", normalizePath(p), "\n")
  X <- as.matrix(as.data.frame(read_excel(p))[, MLSS_GENES17]); storage.mode(X) <- "double"; X
}
X <- load_matrix(); NET_REP <- mlss_network_from_matrix(X)
R_CSV <- read.csv(file.path(RDIR, "Fig_6c_MLSS_ALL_Combinations.csv"), stringsAsFactors = FALSE, check.names = FALSE)
R_KEYS <- mlss_pkeys(R_CSV$Drug_A, R_CSV$Drug_B); REP_SCORE <- setNames(R_CSV[["MLSS_v4.0"]], R_KEYS)
R15 <- R_CSV[order(R_CSV$Rank), ][1:15, ]; REP_TOP15 <- mlss_pkeys(R15$Drug_A, R15$Drug_B)
REP <- mlss_table(TG_REP, POT, CORR_REP); stopifnot(max(abs(REP[["MLSS_v4.0"]] - R_CSV[order(R_CSV$Rank), "MLSS_v4.0"])) < 1e-9)      # the gate, repeated

## ------------------------------------------------------------------ gate part 1 (the protocol gate) and part 2 (upstream network step)
log("MLSS port gate: R port vs existing R output Fig_6c_MLSS_ALL_Combinations.csv")
log(sprintf("Inputs (read only, R_17genes): hubs %s, %d correlated pairs, drug.target.interaction.csv (%d rows)", paste(HUBS5, collapse = ", "), length(CORR_REP), nrow(dc)))
P <- REP; P$key <- mlss_pkeys(P$Drug_A, P$Drug_B); Rr <- R_CSV; Rr$key <- R_KEYS
Mg <- merge(Rr, P, by = "key", suffixes = c("_R", "_py")); d <- abs(Mg[["MLSS_v4.0_R"]] - Mg[["MLSS_v4.0_py"]])
log(sprintf("Port: %d combinations; R output: %d rows; matched by drug pair: %d", nrow(P), nrow(Rr), nrow(Mg)))
log(sprintf("MLSS score: max |difference| = %.3e (criterion < 0.005); mean |difference| = %.3e", max(d), mean(d)))
for (cn in c("Complementarity_C", "Balance_B", "Network_Coverage_V", "Potency_P", "Mean_pIC50", "Mean_Abs_r", "Mean_Pos_r", "Mean_Neg_r", "Base_Score", "Synergy_Bonus", "Antagonism_Penalty"))
  log(sprintf("  %-20s max |difference| = %.3e", cn, max(abs(Mg[[paste0(cn, "_R")]] - Mg[[paste0(cn, "_py")]]))))
rk <- abs(Mg$Rank_R - Mg$Rank_py); gate_ok <- nrow(Mg) == 153 && max(d) < 0.005 && max(rk) == 0
log(sprintf("Absolute rank difference per combination: max %d, number differing %d", max(rk), sum(rk > 0)))
log(sprintf("GATE (153 scores within 0.005 and rank order reproduced): %s", if (gate_ok) "PASS" else "FAIL"))
gs <- data.frame(Rank_R = Mg$Rank_R, Rank_py = Mg$Rank_py, Combination_R = Mg$Combination_R, MLSS_v4.0_R = Mg[["MLSS_v4.0_R"]], MLSS_v4.0_py = Mg[["MLSS_v4.0_py"]], abs_diff = d, stringsAsFactors = FALSE)
gs <- gs[order(gs$Rank_R), ]; py_write_csv(gs, file.path(OUTDIR, "gate_scores.csv"), nd = 14)
log(""); log("Part 2: upstream network step (INS-1 matrix -> four schemes -> consensus -> degree) against stored results")
cons_R <- py_sort(vapply(strsplit(read.csv(file.path(RDIR, "Fig_4d_consensus_pairs.csv"), stringsAsFactors = FALSE)$Consensus_Pairs, "_", fixed = TRUE), function(p) mlss_pkey(p[1], p[2]), ""))
log(sprintf("  consensus pairs equal to Fig_4d_consensus_pairs.csv: %s (%d vs %d)", identical(py_sort(NET_REP$consensus), cons_R), length(NET_REP$consensus), length(cons_R)))
hub_df <- hb$table; log(sprintf("  degree table equal to Fig_6b_Gene_Centrality_Summary.csv: %s", identical(sort(unlist(NET_REP$degree)[hub_df$Gene]) , sort(setNames(hub_df$Degree, hub_df$Gene)))))
J <- fromJSON(GAPDH_JSON, simplifyVector = FALSE)
for (w in c("reported", "C0", "C1")) { n <- mlss_network_from_matrix(mlss_gapdh_variant(X, w)); S <- J[[w]]; sc <- vapply(S$consensus, function(p) mlss_pkey(p[[1]], p[[2]]), "")
  log(sprintf("  %-8s consensus %2d (gapdh_exact.json %2d, identical %s); degree identical %s; significant pairs %s", w, length(n$consensus), length(S$consensus), identical(n$consensus, py_sort(sc)),
              identical(as.list(n$degree[order(names(n$degree))]), S$degree[order(names(S$degree))]), paste(sprintf("%s %d", names(n$sig_counts), n$sig_counts), collapse = ", "))) }

## ------------------------------------------------------------------ scenario machinery
safe_rho <- function(a, b) if (length(a) > 2 && diff(range(a)) > 0 && diff(range(b)) > 0) cor(a, b, method = "spearman") else NaN
rank_r2 <- function(y, cols) {
  yr <- rank(y); A <- cbind(1, do.call(cbind, lapply(cols, as.numeric))); res <- lm.fit(A, yr)$residuals; tot <- sum((yr - mean(yr))^2)
  if (tot > 0) 1 - sum(res^2) / tot else NaN
}
endpoints <- function(T) {
  if (nrow(T) == 0) return(list())
  keys <- mlss_pkeys(T$Drug_A, T$Drug_B); sc <- T[["MLSS_v4.0"]]; s9 <- round(sc, 9); com <- which(keys %in% names(REP_SCORE)); ref <- unname(REP_SCORE[keys[com]])
  top15 <- head(keys, 15); genes <- unique(unlist(strsplit(c(T$Targets_A, T$Targets_B), "|", fixed = TRUE))); genes <- py_sort(genes[genes != ""]); c1 <- T$Complementarity_C >= 1 - 1e-12
  list(n_drugs = length(unique(c(T$Drug_A, T$Drug_B))), n_combinations = nrow(T), n_targets_involved = length(genes), targets_involved = paste(genes, collapse = "|"), n_common = length(com),
       spearman_vs_reported = safe_rho(sc[com], ref), top15_overlap = length(intersect(top15, REP_TOP15)), top15_tie_at_cutoff = length(sc) > 15 && s9[15] == s9[16], top_pair = T$Combination[1],
       n_tied_at_top = sum(s9 == s9[1]), top_score = sc[1], top_pair_in_reported_top15 = keys[1] %in% REP_TOP15, max_abs_score_change = if (length(com)) max(abs(sc[com] - ref)) else NaN,
       r2_potency_all = safe_rho(sc, T$Mean_pIC50)^2, r2_potency_within_C1 = if (sum(c1) > 2) safe_rho(sc[c1], T$Mean_pIC50[c1])^2 else NaN,
       r2_C_plus_potency = rank_r2(sc, list(T$Complementarity_C, rank(T$Mean_pIC50))))
}
ROWS <- list(); SCORES <- list(); DRAWS <- list()
COLS <- c("scenario_id", "family", "description", "n_hubs", "hubs", "n_drugs", "n_combinations", "n_targets_involved", "targets_involved", "r_key_pair", "n_common", "spearman_vs_reported",
          "top15_overlap", "top15_tie_at_cutoff", "top_pair", "n_tied_at_top", "top_score", "top_pair_in_reported_top15", "max_abs_score_change", "r2_potency_all", "r2_potency_within_C1", "r2_C_plus_potency",
          "n_draws", "spearman_p2_5", "spearman_p97_5", "spearman_min", "spearman_max", "top15_overlap_p2_5", "top15_overlap_p97_5", "top15_overlap_max", "top_pair_in_reported_top15_frac", "r2_potency_all_p2_5", "r2_potency_all_p97_5")
add <- function(sid, family, desc, T, hubs, corr, extra = NULL) {
  e <- endpoints(T)
  row <- c(list(scenario_id = sid, family = family, description = desc, n_hubs = length(hubs), hubs = paste(py_sort(hubs), collapse = "|"), r_key_pair = if (!is.null(corr) && !is.na(corr[MLSS_KEY])) unname(corr[MLSS_KEY]) else NA_real_), e, extra)
  ROWS[[length(ROWS) + 1]] <<- row; SCORES[[sid]] <<- setNames(T[["MLSS_v4.0"]], mlss_pkeys(T$Drug_A, T$Drug_B))
  log(sprintf("%-34s rho %s  top15 overlap %s  top %s (%s tied)  targets %s  drugs %s", sid, if (!is.na(e$spearman_vs_reported)) sprintf("%.4f", e$spearman_vs_reported) else "nan", e$top15_overlap, e$top_pair, e$n_tied_at_top, e$n_targets_involved, e$n_drugs))
}
write_outputs <- function() {
  out <- bind_rows(lapply(ROWS, function(r) as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE))); for (cn in setdiff(COLS, names(out))) out[[cn]] <- NA_real_; out <- out[, COLS]
  py_write_csv(out, file.path(OUTDIR, "mlss_scenarios.csv"), nd = 4)
  if (length(DRAWS)) py_write_csv(bind_rows(DRAWS), file.path(OUTDIR, "mlss_permutation_draws.csv"), nd = 4)
  W <- data.frame(Rank = R_CSV$Rank, Combination = R_CSV$Combination, reported_MLSS = R_CSV[["MLSS_v4.0"]], stringsAsFactors = FALSE)
  for (sid in names(SCORES)) W[[sid]] <- unname(SCORES[[sid]][R_KEYS])
  py_write_csv(W, file.path(OUTDIR, "mlss_scenario_scores.csv"), nd = 7); writeLines(LOG, file.path(OUTDIR, "mlss_scenarios_log.txt")); nrow(out)
}
npq <- function(x, p) if (anyNA(x)) NaN else np_percentile(x, p)
fmtg <- function(x) sprintf("%+.2f", x)

## ------------------------------------------------------------------ 0 reported
add("S00_reported_port", "gate", "Port with the reported inputs (reproduces the R output; gate in gate_output.txt)", REP, HUBS5, CORR_REP)

## ------------------------------------------------------------------ 1 Gapdh variants carried through the pipeline
for (w in c("C0", "C1")) {
  net <- mlss_network_from_matrix(mlss_gapdh_variant(X, w)); dg <- unlist(net$degree); hubs <- toupper(names(dg)[dg >= 5]); corr <- mlss_hes1_lookup(net); tg <- mlss_targets_from(dc, hubs)
  add(paste0("S01_gapdh_", w), "gapdh", sprintf("Gapdh %s: 26 Aug DAPT reference %s; consensus %d pairs, hubs (degree >= 5) %s; Hes1-normalized r lookup from the recalculated matrix", w,
      if (w == "C0") "mean of the two valid wells" else "mean of all four wells", length(net$consensus), paste0("[", paste0("'", py_sort(hubs), "'", collapse = ", "), "]")), mlss_table(tg, POT, corr), hubs, corr)
}

## ------------------------------------------------------------------ 2 source of r in the Synergy Bonus / Antagonism Penalty
sources <- list(list("kcnj11", "Kcnj11-normalized pooled r", mlss_hes1_lookup(NET_REP, source = "Kcnj11-normalized")), list("ins1", "Ins1-normalized pooled r", mlss_hes1_lookup(NET_REP, source = "Ins1-normalized")),
  list("raw_pooled", "raw pooled r", mlss_hes1_lookup(NET_REP, source = "Raw")),
  list("raw_adj_linear", "raw treatment-adjusted r (linear scale, df 8), kept only if p < 0.05 as in the R lookup", mlss_adjusted_lookup(X, log = FALSE, signif_only = TRUE)),
  list("raw_adj_linear_nofilter", "raw treatment-adjusted r (linear scale, df 8), value used whatever its p", mlss_adjusted_lookup(X, log = FALSE, signif_only = FALSE)),
  list("raw_adj_log", "raw treatment-adjusted r (log scale, df 8; protocol 2.2 scale), kept only if p < 0.05", mlss_adjusted_lookup(X, log = TRUE, signif_only = TRUE)),
  list("raw_adj_log_nofilter", "raw treatment-adjusted r (log scale, df 8), value used whatever its p", mlss_adjusted_lookup(X, log = TRUE, signif_only = FALSE)),
  list("none", "no correlation in the Synergy Bonus and Antagonism Penalty (V kept as reported)", setNames(numeric(0), character(0))))
for (s in sources) {
  name <- s[[1]]; desc <- s[[2]]; lk <- s[[3]]; T <- mlss_table(TG_REP, POT, lk, corr_v = CORR_REP)
  add(paste0("S02_r_", name), "r_source", paste0("r in Synergy/Antagonism: ", desc, sprintf("; r(Cacna1d,Glp1r) = %s", if (is.na(lk[MLSS_KEY])) "absent" else sprintf("%.4f", lk[MLSS_KEY]))), T, HUBS5, lk)
  T2 <- mlss_table(TG_REP, POT, lk); k1 <- paste(T$Drug_A, T$Drug_B, sep = "\t"); k2 <- paste(T2$Drug_A, T2$Drug_B, sep = "\t")
  if (max(abs(T2[["MLSS_v4.0"]][order(k2)] - T[["MLSS_v4.0"]][order(k1)])) > 1e-12)
    add(paste0("S02_r_", name, "_all"), "r_source", "as above but the same lookup also feeds V (as in the R script; an empty lookup makes V = 0 for every pair)", T2, HUBS5, lk)
}

## ------------------------------------------------------------------ 3 hub rule
for (k in c(4, 5, 6)) {
  hh <- mlss_hubs_from_centrality(file.path(RDIR, "Fig_6b_Gene_Centrality_Summary.csv"), k)$hubs; tg <- mlss_targets_from(dc, hh)
  add(paste0("S03_hub_ge", k), "hub_rule", sprintf("Hub rule degree >= %d on the reported network: hubs %s", k, paste0("[", paste0("'", py_sort(hh), "'", collapse = ", "), "]")), mlss_table(tg, POT, CORR_REP), hh, CORR_REP)
}
write_outputs()          # checkpoint: deterministic scenarios done

## ------------------------------------------------------------------ 4 permutations
DRUGS <- names(TG_REP)[mlss_order_key(names(TG_REP))]; GENE_OF <- setNames(vapply(DRUGS, function(d) TG_REP[[d]][1], ""), DRUGS)
perm_potency <- function(rng) {
  slots <- unname(POT[paste(DRUGS, GENE_OF[DRUGS], sep = "\t")]); idx <- np_permutation(rng, length(DRUGS)) + 1; p <- numeric(0)
  for (i in seq_along(DRUGS)) if (!is.na(slots[idx[i]])) p[paste(DRUGS[i], GENE_OF[[DRUGS[i]]], sep = "\t")] <- slots[idx[i]]
  list(TG_REP, p)
}
perm_target <- function(rng) {
  genes <- unname(GENE_OF[DRUGS]); idx <- np_permutation(rng, length(DRUGS)) + 1; tg <- list(); p <- numeric(0)
  for (i in seq_along(DRUGS)) { d <- DRUGS[i]; tg[[d]] <- genes[idx[i]]; kk <- paste(d, genes[i], sep = "\t"); if (!is.na(POT[kk])) p[paste(d, genes[idx[i]], sep = "\t")] <- POT[[kk]] }
  list(tg, p)
}
run_perm <- function(name, fn) {
  rng <- np_rng(SEED); out <- list()
  for (b in seq_len(NPERM) - 1) { tp <- fn(rng); T <- mlss_table(tp[[1]], tp[[2]], CORR_REP); e <- endpoints(T)
    out[[b + 1]] <- c(list(family = name, draw = as.integer(b)), e[c("spearman_vs_reported", "top15_overlap", "top_pair", "n_tied_at_top", "top_pair_in_reported_top15", "r2_potency_all", "r2_potency_within_C1", "r2_C_plus_potency", "n_targets_involved")]) }
  D <- bind_rows(lapply(out, as.data.frame, stringsAsFactors = FALSE)); DRAWS[[length(DRAWS) + 1]] <<- D; q <- function(cn, p_) npq(D[[cn]], p_)
  u <- unique(D$top_pair); cnt <- vapply(u, function(v) sum(D$top_pair == v), 1); o <- order(-cnt, seq_along(u)); top_name <- u[o[1]]; top_n <- cnt[o[1]]
  extra <- list(n_draws = NPERM, spearman_p2_5 = q("spearman_vs_reported", 2.5), spearman_p97_5 = q("spearman_vs_reported", 97.5), spearman_min = min(D$spearman_vs_reported, na.rm = TRUE), spearman_max = max(D$spearman_vs_reported, na.rm = TRUE),
                top15_overlap_p2_5 = q("top15_overlap", 2.5), top15_overlap_p97_5 = q("top15_overlap", 97.5), top15_overlap_max = max(D$top15_overlap), top_pair_in_reported_top15_frac = mean(D$top_pair_in_reported_top15),
                r2_potency_all_p2_5 = q("r2_potency_all", 2.5), r2_potency_all_p97_5 = q("r2_potency_all", 97.5))
  row <- c(list(scenario_id = paste0("S04_perm_", name), family = "permutation",
     description = sprintf("%s permuted across the 18 drugs, %d draws, seed %d; spearman/overlap/r2 columns are means over draws, spread in the p2_5/p97_5 columns", if (name == "potency") "potency values (missing values travel with the permutation)" else "drug-to-target mapping (each drug keeps its potency)", NPERM, SEED),
     n_hubs = length(HUBS5), hubs = paste(py_sort(HUBS5), collapse = "|"), r_key_pair = unname(CORR_REP[MLSS_KEY]), n_drugs = 18, n_combinations = 153, n_targets_involved = trunc(mean(D$n_targets_involved)), targets_involved = "CACNA1D|GLP1R", n_common = 153,
     spearman_vs_reported = mean(D$spearman_vs_reported, na.rm = TRUE), top15_overlap = mean(D$top15_overlap), top_pair = sprintf("most frequent: %s (%d of %d draws)", top_name, as.integer(top_n), NPERM), n_tied_at_top = mean(D$n_tied_at_top),
     top_pair_in_reported_top15 = mean(D$top_pair_in_reported_top15) > 0.5, r2_potency_all = mean(D$r2_potency_all, na.rm = TRUE), r2_potency_within_C1 = mean(D$r2_potency_within_C1, na.rm = TRUE), r2_C_plus_potency = mean(D$r2_C_plus_potency, na.rm = TRUE)), extra)
  ROWS[[length(ROWS) + 1]] <<- row
  log(sprintf("%-34s rho mean %.3f (2.5-97.5%%: %.3f to %.3f; max %.3f)  top15 overlap mean %.2f (%.0f to %.0f; max %d)  top pair in reported top 15: %.1f%%  R2 potency all %.3f, within C=1 %.3f, C+potency %.3f", row$scenario_id, row$spearman_vs_reported,
              extra$spearman_p2_5, extra$spearman_p97_5, extra$spearman_max, row$top15_overlap, extra$top15_overlap_p2_5, extra$top15_overlap_p97_5, as.integer(extra$top15_overlap_max), 100 * extra$top_pair_in_reported_top15_frac, row$r2_potency_all, row$r2_potency_within_C1, row$r2_C_plus_potency))
}
run_perm("potency", perm_potency); write_outputs()
run_perm("target", perm_target); write_outputs()

## ------------------------------------------------------------------ extras (not in the protocol)
for (r in c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)) {
  lk <- CORR_REP; lk[MLSS_KEY] <- r
  add(paste0("X1_r_sweep_", fmtg(r)), "r_sweep", sprintf("EXTRA: r(Cacna1d,Glp1r) set to %s in the single lookup used for V, Synergy and Antagonism (r = 0 makes V = 0 for every pair)", fmtg(r)), mlss_table(TG_REP, POT, lk), HUBS5, lk)
}
net_h <- NET_REP$res[["Hes1-normalized"]]; out <- list()
for (i in seq_along(net_h$r)) { lk <- CORR_REP; lk[MLSS_KEY] <- net_h$r[i]; e <- endpoints(mlss_table(TG_REP, POT, lk)); out[[i]] <- list(r = net_h$r[i], rho = e$spearman_vs_reported, overlap = e$top15_overlap, top_in = e$top_pair_in_reported_top15) }
O <- bind_rows(lapply(out, as.data.frame))
ROWS[[length(ROWS) + 1]] <- list(scenario_id = "X2_r_any_pair_120", family = "r_any_pair", description = sprintf("EXTRA: r(Cacna1d,Glp1r) replaced, in turn, by each of the 120 Hes1-normalized pooled r (range %.3f to %.3f); columns are means over the 120 replacements, spread in min/max", min(O$r), max(O$r)),
  n_hubs = length(HUBS5), hubs = paste(py_sort(HUBS5), collapse = "|"), r_key_pair = mean(O$r), n_drugs = 18, n_combinations = 153, n_targets_involved = 2, targets_involved = "CACNA1D|GLP1R", n_common = 153, spearman_vs_reported = mean(O$rho), top15_overlap = mean(O$overlap),
  n_draws = 120, spearman_min = min(O$rho), spearman_max = max(O$rho), spearman_p2_5 = np_percentile(O$rho, 2.5), spearman_p97_5 = np_percentile(O$rho, 97.5), top15_overlap_p2_5 = np_percentile(O$overlap, 2.5), top15_overlap_p97_5 = np_percentile(O$overlap, 97.5),
  top15_overlap_max = max(O$overlap), top_pair_in_reported_top15_frac = mean(O$top_in))
log(sprintf("%-34s rho mean %.4f (min %.4f, max %.4f), top15 overlap mean %.1f (min %d, max %d) over 120 replacements of r(Cacna1d,Glp1r) [%.3f to %.3f]", "X2_r_any_pair_120", mean(O$rho), min(O$rho), max(O$rho), mean(O$overlap), as.integer(min(O$overlap)), as.integer(max(O$overlap)), min(O$r), max(O$r)))
WS <- list(list("MLSS (Standard)", .45, .30, .10, .15), list("C-Heavy (60:25:15)", .60, .25, .10, .05), list("C-Heavy (55:20:25)", .55, .20, .15, .10), list("B-Heavy (35:50:15)", .35, .50, .10, .05),
           list("V-Heavy (35:30:35)", .35, .30, .30, .05), list("P-Heavy (30:20:10:40)", .30, .20, .10, .40), list("Balanced (40:40:10:10)", .40, .40, .10, .10), list("Network-Focus (30:20:40:10)", .30, .20, .40, .10))
for (i in seq_along(WS)) { w <- WS[[i]]
  add(paste0("X3_weights_", i - 1), "weights_all153", sprintf("EXTRA: Fig_6d weight scenario \"%s\" (C:B:V:P) on all 153 combinations", w[[1]]), mlss_table(TG_REP, POT, CORR_REP, weights = c(w[[2]], w[[3]], w[[4]], w[[5]])), HUBS5, CORR_REP) }
AB <- list(list("Without_C", c(0, .30 / .55, .10 / .55, .15 / .55)), list("Without_B", c(.45 / .70, 0, .10 / .70, .15 / .70)), list("Without_V", c(.45 / .90, .30 / .90, 0, .15 / .90)),
           list("Without_P", c(.45 / .85, .30 / .85, .10 / .85, 0)), list("Base_Only", MLSS_W_STD))
for (a in AB) { nm <- a[[1]]; w <- a[[2]]
  if (nm == "Base_Only") { T <- mlss_table(TG_REP, POT, setNames(numeric(0), character(0)), weights = w, corr_v = CORR_REP); lk <- setNames(numeric(0), character(0)) } else { T <- mlss_table(TG_REP, POT, CORR_REP, weights = w); lk <- CORR_REP }
  add(paste0("X4_ablation_", nm), "ablation_all153", sprintf("EXTRA: Fig_6d leave-one-out \"%s\" (remaining weights renormalized; Base_Only = no Synergy/Antagonism) on all 153 combinations", nm), T, HUBS5, lk) }
log(""); n_rows <- write_outputs()
log(sprintf("written: mlss_scenarios.csv (%d rows), mlss_permutation_draws.csv (%d rows), mlss_scenario_scores.csv, mlss_scenarios_log.txt", n_rows, sum(vapply(DRAWS, nrow, 1L))))
write_outputs(); cat(sprintf("done in %.0f s\n", as.numeric(Sys.time() - t_start, units = "secs")))
