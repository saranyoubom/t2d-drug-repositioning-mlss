## check_refA_vs_plate1.R - R port of wp11_plate_provenance/check_refA_vs_plate1.py: WP11 (8 Nov 2026), author item A12: the Wnt-family genes of the 24 June sheet (Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7 and their Gapdh vector 'refA') lie about 2.6 cycles
## above the 2025 re-export of 24 June plate 1 (file 2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 without ' 2', threshold 45.70 RFU). reconstruct_threshold.R finds that a single threshold of about 210 RFU applied to the same curves reproduces
## the sheet (mean 0.04 cycles). This script asks the question that matters for the analysis: how much do the quantities that the analysis uses (dCq = target - Gapdh, and ddCq between groups) differ (a) between the sheet and the 2025
## re-export (45.70 RFU) and (b) between the sheet and the curves read at the fitted threshold? Group means of the four wells are used (the sheet carries no well labels for this run), so the comparison does not depend on how replicates
## are paired. Read only. Output: results_R/wp11/refA_vs_plate1.txt.  Usage: Rscript check_refA_vs_plate1.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R")); LG <- plog_new(); T_FIT <- 209.91      # fitted by reconstruct_threshold.R (results/threshold_fit.txt)
pl <- read_plate(file.path(SBA, "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1")); cyc <- pl$cyc; curves <- pl$curves; tab <- pl$tab
sheet <- list(); S <- read_csv_chr(file.path(W1, "ins1_sheet_cq.csv"))
for (i in seq_len(nrow(S))) { k <- key_of(S$gene[i], S$group[i]); if (is.null(sheet[[k]])) sheet[[k]] <- list(); sheet[[k]][[length(sheet[[k]]) + 1L]] <- c(as.numeric(S$target_cq[i]), as.numeric(S$reference_cq[i])) }
G <- c("CTRL", "DAPT", "DKK-1"); genesA <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7"); wells <- split(tab$well, key_of(tab$gene, tab$group)); cqs <- split(tab$cq, key_of(tab$gene, tab$group))
mean_export <- function(g, k) np_mean(cqs[[key_of(g, k)]]); mean_T <- function(g, k) np_mean(vapply(wells[[key_of(g, k)]], function(w) cq_at(cyc, curves[[w]], T_FIT), 1))
sh <- function(g, k, j) vapply(sheet[[key_of(g, k)]], function(x) x[j], 1)
gs <- vapply(G, function(k) np_mean(sh("Wnt5a", k, 2)), 1)      # refA is the reference vector of Wnt5a (WP1 gene-to-run map)
plog(LG, py_sprintf("Gapdh group means (4 wells): sheet refA %s | 2025 re-export (45.70 RFU) %s | curves at %.2f RFU %s", paste(vapply(G, function(k) py_sprintf("%s %.2f", k, gs[[k]]), ""), collapse = " "), paste(vapply(G, function(k) py_sprintf("%s %.2f", k, mean_export("gapdh", k)), ""), collapse = " "), T_FIT, paste(vapply(G, function(k) py_sprintf("%s %.2f", k, mean_T("gapdh", k)), ""), collapse = " ")))
D <- list(); plog(LG)
plog(LG, py_sprintf("%-6s %-6s %9s | %9s %9s %8s | %9s %9s %8s", "gene", "group", "Cq sheet", "dCq sheet", "dCq 45.70", "diff", "Cq sheet", py_sprintf("dCq %.0f", T_FIT), "diff"))
for (g in genesA) for (k in G) {
  ts <- np_mean(sh(g, k, 1)); ds <- ts - gs[[k]]; de <- mean_export(norm(g), k) - mean_export("gapdh", k); dt <- mean_T(norm(g), k) - mean_T("gapdh", k); D[[key_of(g, k)]] <- c(ds, de, dt)
  plog(LG, py_sprintf("%-6s %-6s %9.2f | %9.2f %9.2f %+8.2f | %9.2f %9.2f %+8.2f", g, k, ts, ds, de, ds - de, ts, dt, ds - dt))
}
d1 <- vapply(D, function(v) v[1] - v[2], 1); d2 <- vapply(D, function(v) v[1] - v[3], 1)
plog(LG); plog(LG, py_sprintf("dCq (target minus Gapdh), sheet minus re-export at 45.70 RFU, 15 gene-group means: mean %+.2f, SD %.2f, range %+.2f to %+.2f", np_mean(d1), np_std1(d1), min(d1), max(d1)))
plog(LG, py_sprintf("dCq (target minus Gapdh), sheet minus curves at %.2f RFU, 15 gene-group means: mean %+.2f, SD %.2f, range %+.2f to %+.2f", T_FIT, np_mean(d2), np_std1(d2), min(d2), max(d2)))
plog(LG); plog(LG, py_sprintf("ddCq (group dCq minus CTRL dCq): sheet | re-export 45.70 | difference || curves at %.0f | difference", T_FIT)); dd1 <- numeric(0); dd2 <- numeric(0)
for (g in genesA) { cells <- character(0)
  for (k in c("DAPT", "DKK-1")) { s <- D[[key_of(g, k)]][1] - D[[key_of(g, "CTRL")]][1]; e <- D[[key_of(g, k)]][2] - D[[key_of(g, "CTRL")]][2]; t <- D[[key_of(g, k)]][3] - D[[key_of(g, "CTRL")]][3]
    dd1 <- c(dd1, s - e); dd2 <- c(dd2, s - t); cells <- c(cells, py_sprintf("%s-CTRL %+.2f | %+.2f | %+.2f || %+.2f | %+.2f", k, s, e, s - e, t, s - t)) }
  plog(LG, py_sprintf("  %-6s %s", g, paste(cells, collapse = "      "))) }
plog(LG, py_sprintf("difference of the 10 ddCq values, sheet minus re-export (45.70 RFU): mean %+.2f, SD %.2f, largest absolute %.2f", np_mean(dd1), np_std1(dd1), max(abs(dd1))))
plog(LG, py_sprintf("difference of the 10 ddCq values, sheet minus curves at %.2f RFU: mean %+.3f, SD %.3f, largest absolute %.3f", T_FIT, np_mean(dd2), np_std1(dd2), max(abs(dd2))))
sd_s <- np_mean(unlist(lapply(genesA, function(g) vapply(G, function(k) np_std1(sh(g, k, 1)), 1)))); sd_e <- np_mean(unlist(lapply(genesA, function(g) vapply(G, function(k) np_std1(cqs[[key_of(norm(g), k)]]), 1))))
plog(LG, py_sprintf("mean SD of the four target Cq values per block: sheet %.3f, re-export %.3f", sd_s, sd_e))
write_txt(LG$lines, file.path(RES, "refA_vs_plate1.txt"))
