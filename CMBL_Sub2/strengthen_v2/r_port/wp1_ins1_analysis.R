## wp1_ins1_analysis.R - R port of wp1_ins1_structure/wp1_ins1_analysis.py (protocol v2 section 2, items 2.2-2.8 on the INS-1 17-gene matrix).
## POST HOC and exploratory: one cell preparation, n = 12 (3 groups x 4 wells). Needs gene_run_map.csv and ins1_sheet_cq.csv (inputs of the Python script, kept in
## wp1_ins1_structure/), numpy_rng.R, coexp_checks.R and pycompat.R. Seed 2026 throughout, drawing the same numpy streams as the Python script.
## Writes results_ins1.json and tbl_*.csv into R_OUT_DIR (default results_R/wp1 next to this script); the archived Python files are not touched.
options(warn = 1)
args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "rng", "numpy_rng.R")); source(file.path(here, "coexp_checks.R")); source(file.path(here, "pycompat.R"))
suppressMessages({ library(readxl); library(dplyr) })
WP <- Sys.getenv("WP1_DIR", file.path(here, "..", "wp1_ins1_structure"))
MATRIX <- Sys.getenv("INS1_MATRIX", file.path(here, "..", "..", "R_17genes", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"))
BENCH <- Sys.getenv("BENCH_DIR", file.path(here, "..", "..", "benchmark_route3"))
OUTDIR <- Sys.getenv("R_OUT_DIR", file.path(here, "results_R", "wp1"))
G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r"); DEN <- c("Ins1", "Kcnj11", "Hes1"); SEED <- 2026
PROVISIONAL <- c(Glp1r = "Jun03", Glut2 = "Jun03", Hes1 = "Jun03", Ins1 = "Jun03", Lef1 = "Jun24", Tcf7 = "Jun24", Wnt5a = "Jun24", Wnt5b = "Jun24", Wnt9a = "Jun24", Cacna1c = "Aug26", Cacna1d = "Aug26",
                 Kcnj11 = "Aug26", Ptbp1 = "Aug26", Tcf7l2 = "Aug26", Itpr1 = "Nov09", Ryr2 = "Nov09", Ryr3 = "Nov09")        # planning guess from Table S1c, SUPERSEDED by the verified map
OUT <- list(); t_start <- Sys.time()
save_tbl <- function(df, name, nd = 4) { py_write_csv(df, file.path(OUTDIR, name), nd); df }
rows_df <- function(rows) bind_rows(lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE)))
upper_mean <- function(M) mean(M[upper.tri(M)])
corr_cols <- function(D) { Z <- D - matrix(colMeans(D), nrow(D), ncol(D), byrow = TRUE); Z <- Z / matrix(sqrt(colSums(Z^2)), nrow(D), ncol(D), byrow = TRUE); crossprod(Z) }
arr3 <- function(z, d) aperm(array(z, dim = rev(d)), 3:1)         # numpy C-order (d1, d2, d3) array from a flat vector of draws
kruskal_list <- function(gs) {                                    # scipy.stats.kruskal(*samples).statistic (tie-corrected)
  x <- unlist(gs); N <- length(x); r <- rank(x); idx <- rep(seq_along(gs), lengths(gs)); s <- tapply(r, idx, sum); n <- lengths(gs)
  h <- 12 / (N * (N + 1)) * sum(s^2 / n) - 3 * (N + 1); tc <- table(x); h / (1 - sum(tc^3 - tc) / (N^3 - N))
}
wilcoxon_py <- function(x, y) {                                   # scipy.stats.wilcoxon(x, y).pvalue: zeros dropped, exact for n <= 50 without ties, else normal approximation without continuity correction
  d <- x - y; d <- d[d != 0]; n <- length(d); if (n == 0) return(NaN)
  r <- rank(abs(d)); rp <- sum(r[d > 0]); rm <- sum(r[d < 0]); tt <- min(rp, rm); ties <- table(r); has_ties <- any(ties > 1)
  if (n <= 50 && !has_ties) return(min(1, 2 * psignrank(tt, n)))
  mn <- n * (n + 1) / 4; se <- sqrt(n * (n + 1) * (2 * n + 1) / 24 - sum(ties^3 - ties) / 48); 2 * pnorm(-abs((tt - mn) / se))
}

## ------------------------------------------------------------------------------------------------------------------ load
X <- as.matrix(as.data.frame(read_excel(MATRIX))[, G17]); storage.mode(X) <- "double"; grp <- rep(c("CTRL", "DAPT", "DKK-1"), each = 4); des <- cc_design(grp)
mp <- read.csv(file.path(WP, "gene_run_map.csv"), stringsAsFactors = FALSE, check.names = FALSE); rownames(mp) <- mp$gene; cq <- read.csv(file.path(WP, "ins1_sheet_cq.csv"), stringsAsFactors = FALSE)
stopifnot(all(mp[G17, "status"] == "resolved"))
runs <- mp[G17, "run"]; refs <- mp[G17, "reference_set"]
getcol <- function(col) vapply(G17, function(g) { d <- cq[cq$gene == g, ]; d[order(d$replicate), col] }, numeric(12))
Tq <- getcol("target_cq"); RF <- getcol("reference_cq")             # 12 x 17
stopifnot(all(abs(2^(-(Tq - RF)) - X) <= 1e-12 * abs(X)))           # the sheet Cq values reproduce the matrix
base <- cc_observed_summary(X, grp, DEN, G17); iu17 <- triu_rowmajor(17); NAMES <- base$scheme_names
OUT$status <- "POST HOC exploratory analysis of one cell preparation (3 groups x 4 wells); protocol v2 section 2; not confirmatory"
OUT$map <- setNames(lapply(G17, function(g) list(run = mp[g, "run"], reference_set = mp[g, "reference_set"], sheet = mp[g, "sheet"], rows = mp[g, "sheet_rows"])), G17)
OUT$run_sizes_17genes <- as.list(table(runs)); OUT$reference_set_sizes_17genes <- as.list(table(refs))
OUT$baseline <- list(consensus_any = base$consensus_any, consensus_same_sign = base$consensus_same_sign, pooled_sig_per_scheme = as.list(base$pooled_sig_per_scheme),
                     adj_q05_per_scheme = as.list(base$adj_q_lt_alpha_per_scheme), cons_adj_any = base$cons_adj_any, degree = as.list(base$degree), consensus_pairs = base$consensus_pairs, scheme_names = as.list(NAMES))
cat(sprintf("loaded; baseline consensus %d, pooled %s, adjusted %s\n", base$consensus_any, paste(base$pooled_sig_per_scheme, collapse = ","), paste(base$adj_q_lt_alpha_per_scheme, collapse = ",")))

## ------------------------------------------------------------------------------------------------------------------ 2.2 run-structure test
cat("2.2 run structure\n"); rows22 <- list()
flatten <- function(res, label) for (tag in c("r", "abs_r")) {
  t <- res[[tag]]
  rows22[[length(rows22) + 1]] <<- list(labelling = label, statistic = tag, scope = "all runs", n_genes = res$n_genes, n_inside_pairs = res$n_same_pairs, n_cross_pairs = res$n_between_pairs,
                                        mean_inside = t$mean_same, mean_cross = t$mean_between, D = t$D, null_sd = t$null_sd, p_two_sided = t$p_two_sided)
  for (run in names(t$per_run)) { v <- t$per_run[[run]]
    rows22[[length(rows22) + 1]] <<- list(labelling = label, statistic = tag, scope = paste("run", run), n_genes = v$n_genes, n_inside_pairs = v$n_pairs_inside, mean_inside = v$mean_inside, mean_cross = v$mean_cross, D = v$D, p_two_sided = v$p_two_sided) }
}
rs_primary <- cc_run_structure_test(X, grp, runs, 20000, SEED, G17); flatten(rs_primary, "verified map (run = sheet): primary")
rs_refset <- cc_run_structure_test(X, grp, refs, 20000, SEED, G17); flatten(rs_refset, "verified map, run = reference set (ADD-ON secondary: Jun24 splits in two Gapdh vectors)")
rs_prov <- cc_run_structure_test(X, grp, unname(PROVISIONAL[G17]), 20000, SEED, G17); flatten(rs_prov, "provisional map of the planning reconnaissance (superseded)")
t22 <- save_tbl(rows_df(rows22), "tbl_2_2_run_structure.csv"); OUT[["2_2_run_structure"]] <- list(primary = rs_primary, reference_set = rs_refset, provisional_map = rs_prov)
r_log <- cc_residual_corr(X, grp); r_lin <- cc_residual_corr(X, grp, log = FALSE)

## ------------------------------------------------------------------------------------------------------------------ 2.3 consensus pairs and hubs by run
cat("2.3 consensus by run\n")
elig <- which(!(G17 %in% DEN)); ge <- G17[elig]; ne <- length(elig); iue <- triu_rowmajor(ne)
lab_e <- runs[elig]; uniq_e <- unique(lab_e); code_e <- match(lab_e, uniq_e)
consv <- as.numeric(base$.cons_any[elig, elig][iue]); same_e <- code_e[iue[, 1]] == code_e[iue[, 2]]
a <- sum(same_e & consv == 1); b <- sum(same_e & consv == 0); c_ <- sum(!same_e & consv == 1); d <- sum(!same_e & consv == 0)
m22 <- matrix(c(a, b, c_, d), 2, byrow = TRUE); fo <- (a * d) / (b * c_); fp <- fisher.test(m22)$p.value; fpg <- fisher.test(m22, alternative = "greater")$p.value
rng <- np_rng(SEED); u <- matrix(np_random(rng, 20000 * ne), 20000, ne, byrow = TRUE)
cp <- t(vapply(seq_len(20000), function(r) code_e[order(u[r, ])], numeric(ne)))
null_same <- drop(((cp[, iue[, 1]] == cp[, iue[, 2]]) + 0) %*% consv); obs_same <- sum(same_e * consv)
p_ge <- (1 + sum(null_same >= obs_same)) / 20001; p_le <- (1 + sum(null_same <= obs_same)) / 20001
p_two <- (1 + sum(abs(null_same - mean(null_same)) >= abs(obs_same - mean(null_same)) - 1e-12)) / 20001
deg <- vapply(ge, function(g) { v <- base$degree[g]; if (is.na(v)) 0 else as.numeric(v) }, 1); nrun <- table(lab_e)
hub_rows <- lapply(seq_len(ne), function(i) {
  k <- lab_e[i]; s_possible <- as.integer(nrun[[k]]) - 1L
  partners <- Filter(function(j) j != i && base$.cons_any[elig[min(i, j)], elig[max(i, j)]], seq_len(ne)); same_edges <- sum(lab_e[partners] == k)
  list(gene = ge[i], run = k, reference_set = refs[elig[i]], genes_in_run_17 = sum(runs == k), eligible_genes_in_run = as.integer(nrun[[k]]), degree = deg[i], is_hub_deg_ge_5 = deg[i] >= 5,
       same_run_edges = same_edges, between_run_edges = deg[i] - same_edges, possible_same_run_partners = s_possible, expected_same_run_edges = deg[i] * s_possible / (ne - 1))
})
hubs_df <- rows_df(hub_rows); hubs_df <- hubs_df[order(-hubs_df$degree, hubs_df$gene, method = "radix"), ]; save_tbl(hubs_df, "tbl_2_3_hubs_vs_run.csv")
runrows <- list(); n_all_hub <- sum(deg >= 5)
for (k in uniq_e) {
  m <- lab_e == k; n_k <- sum(m); n_hub <- sum(deg[m] >= 5); mean_deg <- mean(deg[m])
  null_mean_deg <- vapply(1:5000, function(.) mean(deg[np_permutation(rng, ne)[1:n_k] + 1]), 1)
  runrows[[length(runrows) + 1]] <- list(run = k, eligible_genes = n_k, mean_degree = mean_deg, hubs_deg_ge_5 = n_hub, hubs_expected = n_all_hub * n_k / ne,
                                         hypergeometric_p_ge = phyper(n_hub - 1, n_all_hub, ne - n_all_hub, n_k, lower.tail = FALSE), p_mean_degree_ge = (1 + sum(null_mean_deg >= mean_deg - 1e-12)) / 5001)
}
runtab <- save_tbl(rows_df(runrows), "tbl_2_3_degree_by_run.csv")
Hs <- kruskal_list(lapply(seq_along(uniq_e), function(k) deg[code_e == k])); Hp <- pchisq(Hs, length(uniq_e) - 1, lower.tail = FALSE)
## as in the Python script, the permutation is drawn anew for every group inside the comprehension (rng.permutation(ne) is called K times per null draw)
nullH <- vapply(1:3000, function(.) kruskal_list(lapply(seq_along(uniq_e), function(k) deg[code_e[np_permutation(rng, ne) + 1] == k])), 1)
pH <- (1 + sum(nullH >= Hs - 1e-12)) / 3001
inv <- lapply(uniq_e, function(k) {
  m_inv <- lab_e[iue[, 1]] == k | lab_e[iue[, 2]] == k; n_inv <- sum(m_inv); n_c <- sum(m_inv & consv == 1)
  list(run = k, eligible_pairs_involving_run = n_inv, consensus_pairs_involving_run = n_c, consensus_expected = sum(consv) * n_inv / length(consv),
       hypergeometric_p_ge = phyper(n_c - 1, n_inv, length(consv) - n_inv, sum(consv), lower.tail = FALSE))
})
save_tbl(rows_df(inv), "tbl_2_3_pairs_involving_run.csv")
pairs_e <- data.frame(gene_a = ge[iue[, 1]], gene_b = ge[iue[, 2]], run_a = lab_e[iue[, 1]], run_b = lab_e[iue[, 2]], same_run = same_e, consensus = consv == 1, stringsAsFactors = FALSE)
save_tbl(pairs_e, "tbl_2_3_eligible_pairs.csv")
OUT[["2_3_consensus_by_run"]] <- list(eligible_pairs = length(consv), same_run_pairs = sum(same_e), between_run_pairs = sum(!same_e), consensus_same_run = a, same_run_not_consensus = b, consensus_between_run = c_,
  between_run_not_consensus = d, share_consensus_same_run = a / max(a + b, 1), share_consensus_between_run = c_ / max(c_ + d, 1), fisher_odds_ratio = fo, fisher_p_two_sided = fp, fisher_p_greater = fpg,
  label_permutation_null_mean = mean(null_same), label_permutation_null_sd = np_std(null_same), label_permutation_p_ge = p_ge, label_permutation_p_le = p_le, label_permutation_p_two_sided = p_two,
  pairs_involving_run = inv, hubs = as.list(ge[deg >= 5]), hub_table = hub_rows, degree_by_run = runrows, kruskal_H_degree_by_run = Hs, kruskal_p_asymptotic = Hp, kruskal_p_label_permutation = pH)
cat(sprintf("  consensus same-run %d/%d, between-run %d/%d; Fisher OR %.2f P %.3f; label permutation P(>=) %.4f (null mean %.1f)\n", a, a + b, c_, c_ + d, fo, fp, p_ge, mean(null_same)))

## ------------------------------------------------------------------------------------------------------------------ 2.4 shared-reference contribution
cat("2.4 shared reference\n")
r_cq <- cc_residual_corr(Tq, grp, log = FALSE); r_dcq <- cc_residual_corr(X, grp); r_ref <- cc_residual_corr(RF, grp, log = FALSE); stopifnot(all(abs(r_dcq - r_log) <= 1e-8 + 1e-5 * abs(r_log)))
Qm <- des$Q; E <- Qm %*% RF; Tr <- Qm %*% Tq; svar <- colSums(Tr^2); r_exp <- crossprod(E) / sqrt(outer(svar + colSums(E^2), svar + colSums(E^2)))
P <- data.frame(gene_a = G17[iu17[, 1]], gene_b = G17[iu17[, 2]], run_a = runs[iu17[, 1]], run_b = runs[iu17[, 2]], ref_a = refs[iu17[, 1]], ref_b = refs[iu17[, 2]], stringsAsFactors = FALSE)
P$same_run <- P$run_a == P$run_b; P$same_reference_set <- P$ref_a == P$ref_b; P$r_target_cq <- r_cq[iu17]; P$r_dcq <- r_dcq[iu17]; P$r_reference_vectors <- r_ref[iu17]; P$r_expected_shared_reference_only <- r_exp[iu17]
save_tbl(P, "tbl_2_4_pairs.csv"); rc8 <- des$rc_resid
summ <- function(mask, label) {
  s <- P[mask, ]; if (nrow(s) == 0) return(list(scope = label, n_pairs = 0L))
  list(scope = label, n_pairs = nrow(s), mean_r_target_cq = mean(s$r_target_cq), mean_r_dcq = mean(s$r_dcq), mean_r_expected_shared_reference_only = mean(s$r_expected_shared_reference_only),
       mean_abs_r_target_cq = mean(abs(s$r_target_cq)), mean_abs_r_dcq = mean(abs(s$r_dcq)), n_abs_r_ge_0632_target_cq = sum(abs(s$r_target_cq) >= rc8), n_abs_r_ge_0632_dcq = sum(abs(s$r_dcq) >= rc8),
       wilcoxon_p_dcq_vs_target_cq = if (nrow(s) >= 6) wilcoxon_py(s$r_dcq, s$r_target_cq) else NaN)
}
srows <- c(list(summ(P$same_run, "all same-run pairs"), summ(!P$same_run, "all between-run pairs")), lapply(py_sort(unique(runs)), function(k) summ(P$run_a == k & P$run_b == k, paste("same run:", k))),
           lapply(py_sort(unique(refs)), function(k) summ(P$ref_a == k & P$ref_b == k, paste("same reference set:", k))), list(summ(P$same_run & !P$same_reference_set, "same run, different reference sets (Jun24 refA x refC)")))
t24 <- save_tbl(rows_df(srows), "tbl_2_4_shared_reference_summary.csv")
## ADD-ON: layout control for the 26 Aug run (rows in ascending Cq order within group) against a random pairing, simulated from the observed group means and SDs
aug <- which(runs == "Aug26"); naug <- length(aug); B <- 5000; rngs <- np_rng(SEED); sims <- list()
for (layout in c("ascending rows (as in the sheet)", "random row order")) {
  Tm <- array(0, dim = c(B, 12, naug)); Rm <- matrix(0, B, 12)
  for (k in seq_along(des$rows)) {
    rows_ <- des$rows[[k]]; Ta <- Tq[rows_, aug, drop = FALSE]; m <- colMeans(Ta); s <- apply(Ta, 2, sd)
    z <- arr3(np_standard_normal(rngs, B * 4 * naug), c(B, 4, naug)); draw <- z * array(rep(s, each = B * 4), dim = c(B, 4, naug)) + array(rep(m, each = B * 4), dim = c(B, 4, naug))
    if (startsWith(layout, "ascending")) draw <- aperm(apply(draw, c(1, 3), sort), c(2, 1, 3))
    Tm[, rows_, ] <- draw; ref_obs <- RF[rows_, aug[1]]
    if (k == 2) Rm[, rows_] <- matrix(ref_obs, B, 4, byrow = TRUE)
    else { rd <- matrix(np_normal(rngs, mean(ref_obs), sd(ref_obs), B * 4), B, 4, byrow = TRUE); Rm[, rows_] <- if (startsWith(layout, "ascending")) t(apply(rd, 1, sort)) else rd }
  }
  mt <- numeric(B); md <- numeric(B)
  for (bb in seq_len(B)) { Tb <- matrix(Tm[bb, , ], 12, naug); mt[bb] <- upper_mean(corr_cols(Qm %*% Tb)); md[bb] <- upper_mean(corr_cols(Qm %*% (Tb - Rm[bb, ]))) }
  sims[[layout]] <- list(mean_r_target_cq = mt, mean_r_dcq = md)
}
aug_obs <- P[P$run_a == "Aug26" & P$run_b == "Aug26", ]
OUT[["2_4_addon_layout_control_Aug26"]] <- list(note = "ADD-ON (not in the protocol): simulated Aug26 data with independent genes, observed group means/SDs; compare with the observed mean within-group r of the 10 Aug26 pairs",
  observed = list(mean_r_target_cq = mean(aug_obs$r_target_cq), mean_r_dcq = mean(aug_obs$r_dcq)),
  simulated = lapply(sims, function(v) list(mean_r_target_cq = list(mean = mean(v$mean_r_target_cq), q025 = np_percentile(v$mean_r_target_cq, 2.5), q975 = np_percentile(v$mean_r_target_cq, 97.5)),
                                            mean_r_dcq = list(mean = mean(v$mean_r_dcq), q025 = np_percentile(v$mean_r_dcq, 2.5), q975 = np_percentile(v$mean_r_dcq, 97.5), p_ge_observed = mean(v$mean_r_dcq >= mean(aug_obs$r_dcq))))))
## sensitivity versions of the pooled analysis
base_set <- vapply(base$consensus_pairs, function(p) paste(py_sort(p), collapse = "|"), ""); var_rows <- list()
variant <- function(kind) { R2 <- RF
  if (kind == "a") for (rr in des$rows) R2[rr, ] <- matrix(colMeans(RF[rr, , drop = FALSE]), length(rr), 17, byrow = TRUE)
  else if (kind == "b") R2[] <- matrix(colMeans(RF), 12, 17, byrow = TRUE)
  2^(-(Tq - R2)) }
for (vv in list(list("reported (primary)", X), list("(a) reference = run-and-group mean", variant("a")), list("(b) reference = run-wide mean (no normalization)", variant("b")))) {
  name <- vv[[1]]; Xv <- vv[[2]]; ov <- cc_observed_summary(Xv, grp, DEN, G17); cs <- vapply(ov$consensus_pairs, function(p) paste(py_sort(p), collapse = "|"), "")
  cnv <- cc_consensus_null(Xv, grp, DEN, 10000, SEED, genes = G17); ctv <- cc_coupling_test(Xv, grp, 10000, SEED); rsv <- cc_run_structure_test(Xv, grp, runs, 20000, SEED, G17)
  hubs_v <- names(ov$degree)[ov$degree >= 5]; sh <- cc_group_mean_share(Xv, grp, ov$.S[[1]] & upper.tri(ov$.S[[1]]))
  rn <- function(p) { g <- strsplit(p, "|", fixed = TRUE)[[1]]; paste0(paste(g, collapse = "-"), " (", paste(py_sort(unique(runs[match(g, G17)])), collapse = "/"), ")") }
  Q0 <- ov$.Q[[1]]; extra <- list()
  for (k in names(rsv$r$per_run)) for (sn in c("mean_inside", "D", "p")) extra[[paste0(sn, "_", k)]] <- rsv$r$per_run[[k]][[c(mean_inside = "mean_inside", D = "D", p = "p_two_sided")[[sn]]]]
  qpairs <- which(Q0[iu17] < 0.05)
  var_rows[[length(var_rows) + 1]] <- c(list(version = name, consensus_pairs = ov$consensus_any, kept = length(intersect(cs, base_set)), lost = length(setdiff(base_set, cs)), new = length(setdiff(cs, base_set)),
    null_mean = cnv$null_mean, p_null_ge = cnv$p_ge, pooled_raw = ov$pooled_sig_per_scheme[1], pooled_Ins1 = ov$pooled_sig_per_scheme[2], pooled_Kcnj11 = ov$pooled_sig_per_scheme[3], pooled_Hes1 = ov$pooled_sig_per_scheme[4],
    adj_raw = ov$adj_q_lt_alpha_per_scheme[1], adj_Ins1 = ov$adj_q_lt_alpha_per_scheme[2], adj_Kcnj11 = ov$adj_q_lt_alpha_per_scheme[3], adj_Hes1 = ov$adj_q_lt_alpha_per_scheme[4],
    hubs_deg_ge_5 = paste(sprintf("%s %d", hubs_v, as.integer(ov$degree[hubs_v])), collapse = "; "), max_degree = if (length(ov$degree)) max(ov$degree) else 0,
    T_r2 = ctv$T_r2$observed, T_r2_null_mean = ctv$T_r2$null_mean, T_r2_p = ctv$T_r2$p_ge_plus1, T_p = ctv$T_p$observed, T_p_null_mean = ctv$T_p$null_mean, T_p_p = ctv$T_p$p_ge_plus1,
    share_median = sh$median, D_same_minus_between = rsv$r$D, D_p = rsv$r$p_two_sided, mean_r_same_run = rsv$r$mean_same, mean_r_between_run = rsv$r$mean_between), extra,
    list(raw_scheme_pairs_q_lt_005 = paste(vapply(qpairs, function(q) { i <- iu17[q, 1]; j <- iu17[q, 2]; sprintf("%s-%s (%s) q=%.3f", G17[i], G17[j], paste(py_sort(unique(c(runs[i], runs[j]))), collapse = "/"), Q0[i, j]) }, ""), collapse = "; "),
         consensus_pairs_lost = paste(py_sort(vapply(setdiff(base_set, cs), rn, "")), collapse = "; "), consensus_pairs_new = paste(py_sort(vapply(setdiff(cs, base_set), rn, "")), collapse = "; ")))
}
t24v <- save_tbl(rows_df(var_rows), "tbl_2_4_sensitivity_versions.csv"); OUT[["2_4_shared_reference"]] <- list(summary = srows, versions = var_rows)

## ------------------------------------------------------------------------------------------------------------------ 2.5 calibrated adjusted counts
cat("2.5 calibrated adjusted counts\n")
an <- cc_adjusted_null(X, grp, DEN, 2000, SEED, genes = G17); ac <- list()
for (n in an$scheme_names) { v <- an$schemes[[n]]
  ac[[length(ac) + 1]] <- list(scheme = n, observed_q_lt_005 = v$observed, null_mean = v$null_mean, null_sd = v$null_sd, null_q95 = v$null_q95, excess = v$excess, p_null_ge_obs = v$p_ge, p_plus1 = v$p_ge_plus1, calibrated = (n == "raw"),
    reading = if (n == "raw") "calibrated: no pair has q < 0.05 (null mean ~0)" else if (v$p_ge >= 0.05) "artefact-prone: shared denominator, observed count not above its null" else "artefact-prone: shared denominator; observed above null at P < 0.05 (still not calibrated)") }
t25 <- save_tbl(rows_df(ac), "tbl_2_5_adjusted_calibrated.csv"); Qs <- base$.Q
rq <- P[, c("gene_a", "gene_b", "run_a", "run_b", "same_run")]; rq$r_pooled_raw_linear <- cor(X)[iu17]; rq$r_adjusted_partial_linear <- r_lin[iu17]
rq$p_adjusted <- 2 * pt(-abs(rq$r_adjusted_partial_linear * sqrt(des$df_resid / (1 - rq$r_adjusted_partial_linear^2))), des$df_resid)
for (s in seq_along(NAMES)) rq[[paste0("q_", NAMES[s])]] <- Qs[[s]][iu17]
rq$consensus_pair <- base$.cons_any[iu17]; rq <- rq[order(rq$q_raw), ]; save_tbl(rq, "tbl_2_5_q_per_pair.csv")
OUT[["2_5_adjusted_calibrated"]] <- list(schemes = ac, raw_scheme_min_q = min(Qs[[1]], na.rm = TRUE), raw_scheme_pairs_q_lt_005 = sum(Qs[[1]][iu17] < 0.05),
  smallest_raw_q_pairs = lapply(seq_len(5), function(i) as.list(rq[i, c("gene_a", "gene_b", "r_adjusted_partial_linear", "p_adjusted", "q_raw")])))

## ------------------------------------------------------------------------------------------------------------------ 2.6 within-group coupling and run-block null
cat("2.6 coupling and run-block null\n")
smk <- cc_same_run_mask(runs); masks <- list(same_run = smk$same, between_run = smk$between)
ct_i <- cc_coupling_test(X, grp, 10000, SEED, pair_masks = masks); ct_b <- cc_coupling_test(X, grp, 10000, SEED, blocks = runs, pair_masks = masks)
cn_i <- cc_consensus_null(X, grp, DEN, 10000, SEED, genes = G17); cn_b <- cc_consensus_null(X, grp, DEN, 10000, SEED, blocks = runs, genes = G17); ctrows <- list()
for (nm in list(list("independent genes (protocol null)", ct_i), list("run-block (replicate labels permuted jointly within each run)", ct_b))) for (k in names(nm[[2]])) { v <- nm[[2]][[k]]
  ctrows[[length(ctrows) + 1]] <- list(null = nm[[1]], statistic = k, observed = v$observed, null_mean = v$null_mean, null_sd = v$null_sd, null_q95 = v$null_q95, p_ge_plus1 = v$p_ge_plus1) }
for (nm in list(list("independent genes (protocol null)", cn_i), list("run-block (replicate labels permuted jointly within each run)", cn_b))) { cn <- nm[[2]]
  ctrows[[length(ctrows) + 1]] <- list(null = nm[[1]], statistic = "consensus_count", observed = cn$observed, null_mean = cn$null_mean, null_sd = cn$null_sd, null_q95 = cn$null_q95, p_ge_plus1 = cn$p_ge_plus1, p_ge = cn$p_ge)
  if (startsWith(nm[[1]], "run-block")) ctrows[[length(ctrows) + 1]] <- list(null = nm[[1]], statistic = "consensus_count_between_run_only", observed = cn$observed - a, null_mean = cn$null_mean - a, null_sd = cn$null_sd, null_q95 = cn$null_q95 - a, p_ge_plus1 = cn$p_ge_plus1, p_ge = cn$p_ge) }
t26 <- save_tbl(rows_df(ctrows), "tbl_2_6_coupling.csv"); strip <- function(l) l[!startsWith(names(l), ".")]
OUT[["2_6_coupling"]] <- list(independent = ct_i, run_block = ct_b, consensus_independent = strip(cn_i), consensus_run_block = strip(cn_b), consensus_same_run_pairs_observed = a, consensus_between_run_pairs_observed = c_)

## ------------------------------------------------------------------------------------------------------------------ 2.7 geometric expectation
cat("2.7 geometry\n")
gz <- list(`ANOVA on log scale` = cc_geometry_check(X, grp, log = TRUE), `ANOVA on linear scale` = cc_geometry_check(X, grp, log = FALSE))
bm_ <- read.csv(file.path(BENCH, "per_dataset_results.csv"), stringsAsFactors = FALSE); A <- bm_[bm_$set == "A", ]; A$share_raw <- A$pooled_raw / 136
rcB <- cc_rcrit(10); expB <- cc_geometric_expectation(rcB, 3)
bsum <- list(n_series = nrow(A), median_share = median(A$share_raw), q25 = np_percentile(A$share_raw, 25), q75 = np_percentile(A$share_raw, 75), mean_share = mean(A$share_raw), expected = expB, n_above_expected = sum(A$share_raw > expB),
             by_platform = lapply(split(A, A$platform), function(g) list(n = nrow(g), median = median(g$share_raw))))
L1 <- gz[["ANOVA on log scale"]]; L2 <- gz[["ANOVA on linear scale"]]
geo <- data.frame(data = c(sprintf("INS-1, genes with ANOVA group effect (log scale, BH q<0.05): %d genes", L1$n_group_effect_genes), "INS-1, all 17 genes", sprintf("INS-1, genes with ANOVA group effect (linear scale): %d genes", L2$n_group_effect_genes),
                           sprintf("benchmark part 1, gene set A (17 DE genes), median of %d series", nrow(A))), n_pairs = c(L1$n_pairs_group_effect, 136, L2$n_pairs_group_effect, 136),
                  n_sig = c(L1$n_sig_group_effect, L1$n_sig_all, L2$n_sig_group_effect, median(A$pooled_raw)), share = c(L1$share_sig_group_effect, L1$share_sig_all, L2$share_sig_group_effect, bsum$median_share),
                  expected = c(L1$expected, L1$expected, L2$expected, expB), stringsAsFactors = FALSE)
save_tbl(geo, "tbl_2_7_geometry.csv"); OUT[["2_7_geometry"]] <- list(ins1 = gz, benchmark_set_A = bsum, note = "anchor, not a test; assumes the pooled correlation is carried by the group means with uniform orientation of the centred group-mean vectors")

## ------------------------------------------------------------------------------------------------------------------ 2.8 scale sensitivity
cat("2.8 scale sensitivity\n"); srows8 <- list()
for (kd in list(c("pearson", "linear-scale Pearson (the manuscript)"), c("log", "log-scale Pearson"), c("spearman", "Spearman"))) {
  kind <- kd[1]; ov <- cc_observed_summary(X, grp, DEN, G17, kind = kind); cnk <- cc_consensus_null(X, grp, DEN, 10000, SEED, kind = kind, genes = G17); sh <- cc_group_mean_share(X, grp, ov$.S[[1]] & upper.tri(ov$.S[[1]]), kind = kind)
  srows8[[length(srows8) + 1]] <- list(scale = kd[2], pooled_raw = ov$pooled_sig_per_scheme[1], pooled_Ins1 = ov$pooled_sig_per_scheme[2], pooled_Kcnj11 = ov$pooled_sig_per_scheme[3], pooled_Hes1 = ov$pooled_sig_per_scheme[4], consensus_pairs = ov$consensus_any,
    same_sign = ov$consensus_same_sign, null_mean = cnk$null_mean, null_sd = cnk$null_sd, p_null_ge = cnk$p_ge, covariance_share_median = sh$median, covariance_share_n_pairs = sh$n_pairs, adj_q_raw = ov$adj_q_lt_alpha_per_scheme[1],
    adj_q_Ins1 = ov$adj_q_lt_alpha_per_scheme[2], adj_q_Kcnj11 = ov$adj_q_lt_alpha_per_scheme[3], adj_q_Hes1 = ov$adj_q_lt_alpha_per_scheme[4], consensus_with_adj_q_any = ov$cons_adj_any,
    hubs_deg_ge_5 = paste(sprintf("%s %d", names(ov$degree)[ov$degree >= 5], as.integer(ov$degree[ov$degree >= 5])), collapse = "; "))
}
t28 <- save_tbl(rows_df(srows8), "tbl_2_8_scale_sensitivity.csv"); OUT[["2_8_scale_sensitivity"]] <- srows8
OUT$runtime_seconds <- round(as.numeric(Sys.time() - t_start, units = "secs"), 1); py_write_json(OUT, file.path(OUTDIR, "results_ins1.json"))
cat(sprintf("done in %.0fs; wrote results_ins1.json and tbl_*.csv to %s\n", as.numeric(Sys.time() - t_start, units = "secs"), OUTDIR))
