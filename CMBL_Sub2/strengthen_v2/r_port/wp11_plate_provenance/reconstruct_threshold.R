## reconstruct_threshold.R - R port of wp11_plate_provenance/reconstruct_threshold.py: WP11 (8 Nov 2026): which threshold reproduces the values that the calculation sheet holds for the seven Wnt-family genes and their Gapdh (the 'refA' genes)?
## Method: Cq is recomputed from the baseline-subtracted amplification curves ('Quantification Amplification Results') of a plate by linear interpolation of the first crossing of a threshold T.
## Validation: for plates with a known threshold (the user's exports of 8 Nov: 'Before' as found and 'After' = 150) the fitted T must come back as that threshold, and the recomputed Cq must equal the exported Cq.
## Then T is fitted to the sheet values (4 wells per gene and group on the 24 June plate 1, file 2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 without ' 2', 2025 re-export): the T with the smallest mean absolute difference of the sorted blocks.
## Blocks in which a well never reaches T are skipped and listed (one well, G04 of the Wnt2b DAPT block, peaks at 89.6 RFU). The sheet values are matched without any well mapping: sorted block against sorted block.
## Read only. Usage: Rscript reconstruct_threshold.R   (writes results_R/wp11/threshold_fit.txt; WP11_OUT overrides the folder)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R")); LG <- plog_new()
fit_T <- function(cyc, curves, tab, grid) {      # the T of the grid with the smallest mean |recomputed - exported Cq| over all wells (a T at which some well does not cross is skipped); c(mean, T, max)
  best <- NULL; use <- which(tab$well %in% names(curves))
  for (T in grid) { d <- vapply(use, function(i) abs(cq_at(cyc, curves[[tab$well[i]]], T) - tab$cq[i]), 1); if (anyNA(d)) next; s <- np_mean(d); if (is.null(best) || s < best[1]) best <- c(s, T, max(d)) }
  best
}
## 1. validation on the user's exports of 8 Nov (known thresholds)
plog(LG, "=== validation: the fitted threshold must come back as the known one ===")
for (v in list(c("After", "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2", 150.0), c("Before", "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2", 138.061), c("Before", "2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1", 36.072), c("After", "2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1", 150.0))) {
  pl <- read_plate(file.path(EXD, v[1], v[2])); coarse <- fit_T(pl$cyc, pl$curves, pl$tab, np_arange(10, 400, 1)); fine <- fit_T(pl$cyc, pl$curves, pl$tab, np_arange(coarse[2] - 1.5, coarse[2] + 1.5, 0.01))
  plog(LG, py_sprintf("  %-6s %-46s known T %.3f  fitted T %.3f  mean |recomputed - exported Cq| %.4f (max %.4f)", v[1], v[2], as.numeric(v[3]), fine[2], fine[1], fine[3]))
}
## 2. the 24 June plate 1 (2025 re-export): threshold of the export, then the threshold of the sheet
pl <- read_plate(file.path(SBA, "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1")); cyc <- pl$cyc; curves <- pl$curves; tab <- pl$tab
coarse <- fit_T(cyc, curves, tab, np_arange(10, 400, 1)); fine <- fit_T(cyc, curves, tab, np_arange(coarse[2] - 1.5, coarse[2] + 1.5, 0.01))
plog(LG); plog(LG, py_sprintf("=== 24 June plate 1 (file without \" 2\", 2025 re-export): threshold of the export = %.3f (mean |diff| %.4f, max %.4f) ===", fine[2], fine[1], fine[3]))
sheet <- list(); S <- read_csv_chr(file.path(W1, "ins1_sheet_cq.csv"))
for (i in seq_len(nrow(S))) { k <- key_of(S$gene[i], S$group[i]); if (is.null(sheet[[k]])) sheet[[k]] <- list(); sheet[[k]][[length(sheet[[k]]) + 1L]] <- c(as.numeric(S$target_cq[i]), as.numeric(S$reference_cq[i])) }
G <- c("CTRL", "DAPT", "DKK-1"); genesA <- c("Wnt2", "Wnt2b", "Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7")      # Wnt2 and Wnt2b are not among the 17 analyzed genes; they share the plate and the vector refA (Wnt2b has own DAPT values on the sheet)
wells <- split(tab$well, key_of(tab$gene, tab$group)); peaks <- vapply(curves, function(y) max(y, na.rm = TRUE), 1); weak <- tab[which(peaks[tab$well] < 150), ]
plog(LG, "wells per gene and group on this plate: " %+% paste(vapply(Filter(function(g) key_of(norm(g), "CTRL") %in% names(wells), c(genesA, "Gapdh")), function(g) py_sprintf("%s:%d", g, length(wells[[key_of(norm(g), "CTRL")]])), ""), collapse = ", "))
sel7 <- tab$gene %in% c(norm(genesA), "gapdh"); pk <- unname(peaks[tab$well[sel7]])
plog(LG, py_sprintf("peak baseline-subtracted RFU over the 96 wells of the seven assays and Gapdh: median %.0f, minimum %.1f; wells with a peak below 150 RFU: %s", np_median(pk), min(pk), paste(vapply(seq_len(nrow(weak)), function(i) py_sprintf("%s (%s %s, %.1f RFU)", weak$well[i], weak$gene[i], weak$group[i], peaks[[weak$well[i]]]), ""), collapse = ", ")))
BL <- list(); for (g in genesA) for (k in G) BL[[length(BL) + 1L]] <- list(key = key_of(norm(g), k), k = k, name = g, s = sort(vapply(sheet[[key_of(g, k)]], function(x) x[1], 1)))
for (k in G) BL[[length(BL) + 1L]] <- list(key = key_of("gapdh", k), k = k, name = "Gapdh(refA)", s = sort(vapply(sheet[[key_of("Wnt5a", k)]], function(x) x[2], 1)))
block_cq <- function(b, T) vapply(wells[[b$key]], function(w) cq_at(cyc, curves[[w]], T), 1, USE.NAMES = FALSE)
mean_cost <- function(T, sel = function(n) TRUE) {      # c(mean |sorted recomputed - sorted sheet| over the blocks that every well crosses T, number of skipped blocks)
  d <- list(); nsk <- 0L
  for (b in BL) { if (!sel(b$name)) next; e <- block_cq(b, T); if (anyNA(e)) { nsk <- nsk + 1L; next }; d[[length(d) + 1L]] <- abs(sort(e) - b$s) }
  c(if (length(d)) np_mean(unlist(d)) else Inf, nsk)
}
best_T <- function(sel = function(n) TRUE, lo = 60, hi = 420) {      # a threshold that skips many blocks is not a valid fit: at most 2 of the 24 blocks may be skipped
  grid <- np_arange(lo, hi, 1); cs <- t(vapply(grid, function(T) mean_cost(T, sel), c(0, 0))); ok <- is.finite(cs[, 1]) & cs[, 2] <= 2; T0 <- min(grid[ok][cs[ok, 1] == min(cs[ok, 1])])
  fine <- np_arange(T0 - 1.5, T0 + 1.5, 0.01); cf <- t(vapply(fine, function(T) mean_cost(T, sel), c(0, 0))); idx <- which(cf[, 2] <= 2); i <- idx[which.min(cf[idx, 1])]
  list(T = fine[i], v = cf[i, 1], n = cf[i, 2])
}
b1 <- best_T(); T1 <- b1$T; v1 <- b1$v; n1 <- as.integer(b1$n)
plog(LG, py_sprintf("threshold that best reproduces the sheet values (seven assays and Gapdh refA, blocks of four wells): T = %.2f RFU, mean |recomputed - sheet| = %.4f cycles over %d values (blocks skipped because a well never reaches T: %d)", T1, v1, 4L * (24L - n1), n1))
plog(LG, "mean |recomputed - sheet| at other thresholds: " %+% paste(vapply(c(45.7, 100, 150, 180, 200, T1, 220, 250), function(T) py_sprintf("T %.0f: %.3f", T, mean_cost(T)[1]), ""), collapse = "; "))
bg <- best_T(function(n) n == "Gapdh(refA)"); bt <- best_T(function(n) n != "Gapdh(refA)"); b5 <- best_T(function(n) n %in% c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7"))
plog(LG, py_sprintf("thresholds fitted separately: Gapdh blocks T = %.2f (mean %.4f); target blocks T = %.2f (mean %.4f); the five analyzed genes only (Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7) T = %.2f (mean %.4f)", bg$T, bg$v, bt$T, bt$v, b5$T, b5$v))
plog(LG, py_sprintf("per block at T = %.2f (largest absolute difference of the four sorted values, in cycles; CTRL / DAPT / DKK-1):", T1))
for (g in c(genesA, "Gapdh(refA)")) {
  row <- vapply(G, function(k) { b <- Filter(function(b) b$name == g && b$k == k, BL)[[1]]; e <- block_cq(b, T1); if (anyNA(e)) "skipped" else py_sprintf("%.3f", max(abs(sort(e) - b$s))) }, "")
  plog(LG, py_sprintf("  %-12s %s", g, paste(row, collapse = " / ")))
}
write_txt(LG$lines, file.path(RES, "threshold_fit.txt"))
