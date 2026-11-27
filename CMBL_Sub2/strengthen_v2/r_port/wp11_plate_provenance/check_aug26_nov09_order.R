## check_aug26_nov09_order.R - R port of wp11_plate_provenance/check_aug26_nov09_order.py: WP11 (8 Nov 2026), author item A13: who put the 26 August rows in ascending Cq order, the instrument export or the workbook?
## Inputs (read only): the instrument exports of the .pcrd files made in 2025 ('Quantification Cq Results', wells in plate order), the user's 26 August 'Quantification Summary' workbook (sheet '0': CFX columns B-H plus the user's
## calculation columns J-O), the matrix rows of WP1 (ins1_sheet_cq.csv, gene_run_map.csv). Outputs: results_R/wp11/aug26_nov09_order_check.txt and .json.  Usage: Rscript check_aug26_nov09_order.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R")); LG <- plog_new(); G <- c("CTRL", "DAPT", "DKK-1"); TOL <- 0.0051
F_EXP <- file.path(RT, "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx"); F_SUM <- file.path(RT, "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx"); F_N9 <- file.path(RT, "2020-11-09_ITPR, RYR1-3 -  Quantification Cq Results.xlsx")
wkey_order <- function(w) order(as.integer(substring(w, 2L)), substr(w, 1L, 1L))      # plate order: column by column, row A to H
meta <- list(); M <- read_csv_chr(file.path(W1, "gene_run_map.csv")); for (i in seq_len(nrow(M))) meta[[M$gene[i]]] <- M$run[i]
sheet <- list(); S <- read_csv_chr(file.path(W1, "ins1_sheet_cq.csv")); for (i in seq_len(nrow(S))) { k <- key_of(S$gene[i], S$group[i]); if (is.null(sheet[[k]])) sheet[[k]] <- list(); sheet[[k]][[length(sheet[[k]]) + 1L]] <- list(rep = as.integer(S$replicate_in_group[i]), t = as.numeric(S$target_cq[i]), r = as.numeric(S$reference_cq[i])) }
ss <- function(g, k) { x <- sheet[[key_of(g, k)]]; x[order(vapply(x, function(z) z$rep, 1L))] }
genes26 <- names(meta)[vapply(meta, function(r) r == "Aug26", TRUE)]; genes9 <- names(meta)[vapply(meta, function(r) r == "Nov09", TRUE)]
plog(LG, "genes of the 26 August run in the matrix: " %+% paste(genes26, collapse = ", ") %+% "; of the 9 November run: " %+% paste(genes9, collapse = ", "))
nondecr <- function(v) all(v[-length(v)] <= v[-1])      # Python: all(v[i] <= v[i+1]) (TRUE for fewer than 2 values)
## A. instrument export, 26 August
ex <- cq_table(F_EXP); blocks <- list()
for (g in c(genes26, "Gapdh")) for (k in G) { it <- ex[[key_of(norm(g), k)]]; o <- wkey_order(vapply(it, function(z) z$well, "")); blocks[[key_of(g, k)]] <- it[o] }
bw <- function(g, k) vapply(blocks[[key_of(g, k)]], function(z) z$well, ""); bq <- function(g, k) vapply(blocks[[key_of(g, k)]], function(z) z$cq, 1)
asc <- vapply(unlist(lapply(genes26, function(g) lapply(G, function(k) c(g, k))), recursive = FALSE), function(gk) nondecr(bq(gk[1], gk[2])), TRUE)
plog(LG); plog(LG, "A. Instrument export of the 26 August plate (re-export of the .pcrd, wells in plate order)")
plog(LG, py_sprintf("   target blocks (5 genes x 3 groups) whose four wells are in ascending Cq order in plate order: %d of %d (chance expectation about 1 in 24 per block, i.e. 0.6 of 15)", sum(asc), length(asc)))
## B. the user's Quantification Summary workbook
wsum <- sheet0(F_SUM); rows <- list()
for (r in 2:length(wsum)) { cells <- wsum[[r]]; get <- function(c) if (c <= length(cells)) cells[[c]] else NA; w <- get(2); tg <- get(4); sm <- get(6); cq <- get(7); if (is.na(w) || !isnum(cq)) next
  rows[[length(rows) + 1L]] <- list(r = r, w = as.character(w), g = norm(tg), k = grp(sm), c = as.numeric(cq)) }
inst <- list(); for (key in names(ex)) { g <- split_key(key)[1]; for (z in ex[[key]]) inst[[paste(z$well, g, sep = "\r")]] <- z$cq }
ikey <- function(x) paste(x$w, x$g, sep = "\r"); tot <- length(rows); same <- sum(vapply(rows, function(x) !is.null(inst[[ikey(x)]]) && abs(inst[[ikey(x)]] - x$c) <= TOL, TRUE))
tgt_rows <- Filter(function(x) x$g != "gapdh", rows); ref_rows <- Filter(function(x) x$g == "gapdh", rows)
same_t <- sum(vapply(tgt_rows, function(x) abs(inst[[ikey(x)]] - x$c) <= TOL, TRUE)); same_r <- sum(vapply(ref_rows, function(x) !is.null(inst[[ikey(x)]]) && abs(inst[[ikey(x)]] - x$c) <= TOL, TRUE))
plog(LG); plog(LG, "B. The user's 26 August 'Quantification Summary' workbook (sheet '0', Well and Cq columns) against the instrument value of the same well")
plog(LG, py_sprintf("   rows compared: %d (60 target rows, 12 Gapdh rows); Cq equal to the instrument value of the named well: %d of %d (target rows %d of %d, Gapdh rows %d of %d)", tot, same, tot, same_t, length(tgt_rows), same_r, length(ref_rows)))
wbv <- function(g, k) vapply(Filter(function(x) x$g == norm(g) && identical(x$k, k), rows), function(x) x$c, 1)      # workbook values of a block, in sheet order
ordered <- 0L; nblk <- 0L
for (g in c(genes26, "Gapdh")) for (k in G) { v <- wbv(g, k); if (length(v) == 4) { nblk <- nblk + 1L; ordered <- ordered + as.integer(nondecr(v)) } }
plog(LG, py_sprintf("   blocks (gene x group) whose four Cq values run in ascending order down the sheet: %d of %d (Gapdh DAPT is the exception: its four cells hold two values twice)", ordered, nblk))
plog(LG, "   block by block, sorted workbook values against sorted instrument values (tolerance 0.005):"); setres <- list()
for (g in c(genes26, "Gapdh")) { cells <- character(0)
  for (k in G) { v <- sort(wbv(g, k)); e <- sort(bq(g, k)); ok <- length(v) == length(e) && all(abs(v - e) <= TOL); setres[[key_of(g, k)]] <- ok; cells <- c(cells, if (ok) "same set" else "DIFFERENT") }
  plog(LG, py_sprintf("     %-8s %s", g, paste(cells, collapse = " / "))) }
out <- list(); diffs <- sort(vapply(names(setres)[!unlist(setres)], function(key) paste(split_key(key), collapse = " "), "", USE.NAMES = FALSE), method = "radix")
out$summary_workbook <- list(rows = as.integer(tot), cq_equal_to_instrument_value_of_named_well = as.integer(same), target_rows = list(as.integer(same_t), length(tgt_rows)), gapdh_rows = list(as.integer(same_r), length(ref_rows)),
                             blocks_same_value_set_as_instrument = as.integer(sum(unlist(setres))), blocks_total = length(setres), blocks_different = as.list(diffs))
dap <- vapply(Filter(function(x) x$g == "gapdh" && identical(x$k, "DAPT"), rows), function(x) x$c, 1)
plog(LG, py_sprintf("   Gapdh DAPT: workbook cells %s; instrument wells %s", paste(py_fmt_f(dap, 3L), collapse = ", "), paste(vapply(seq_along(bw("Gapdh", "DAPT")), function(i) py_sprintf("%s %.3f", bw("Gapdh", "DAPT")[i], bq("Gapdh", "DAPT")[i]), ""), collapse = ", ")))
## C. the matrix rows of WP1
plog(LG); plog(LG, "C. The matrix rows (WP1, row order of the sheet) against the instrument export"); asc_m <- 0L; set_m <- 0L; nb <- 0L
for (g in genes26) for (k in G) { v <- vapply(ss(g, k), function(z) z$t, 1); e <- bq(g, k); nb <- nb + 1L; asc_m <- asc_m + as.integer(nondecr(v)); set_m <- set_m + as.integer(all(abs(sort(v) - sort(e)) <= TOL)) }
plog(LG, py_sprintf("   26 August target blocks ascending in the sheet row order: %d of %d; same value set as the instrument: %d of %d", asc_m, nb, set_m, nb))
for (k in G) { v <- vapply(ss(genes26[1], k), function(z) z$r, 1); e <- bq("Gapdh", k); plog(LG, py_sprintf("   Gapdh %-5s matrix reference vector %s | instrument wells %s", k, paste(py_fmt_f(v, 3L), collapse = " "), paste(py_fmt_f(e, 3L), collapse = " "))) }
out$matrix_aug26 <- list(target_blocks_ascending = list(asc_m, nb), target_blocks_same_set_as_instrument = list(set_m, nb), plate_order_ascending_instrument = list(as.integer(sum(asc)), length(asc)))
## D. pairing of target and Gapdh wells (descriptive, post hoc)
plog(LG); plog(LG, "D. Pairing of target and Gapdh wells on the 26 August plate (descriptive; the reported analysis keeps the sheet pairing)")
dcq_matrix <- function(variant) {
  t(vapply(genes26, function(g) unlist(lapply(G, function(k) {
    if (variant == "sheet") vapply(ss(g, k), function(z) z$t - z$r, 1)      # the matrix as analyzed: sheet vectors, k-th smallest with k-th smallest
    else { tt <- bq(g, k); rr <- bq("Gapdh", k); if (variant == "rank_instrument") { tt <- sort(tt); rr <- sort(rr) }; tt - rr } })), numeric(12)))
}
desc <- list()
for (cs in list(c("as analyzed (sheet vectors, rank pairing, DAPT Gapdh wells A10 and B10 replaced)", "sheet"), c("instrument values, rank pairing (all four Gapdh wells)", "rank_instrument"), c("instrument values, plate-order pairing (all four Gapdh wells)", "plate_instrument"))) {
  Mx <- dcq_matrix(cs[2]); ng <- length(genes26); Gm <- array(0, c(ng, 3, 4)); for (i in seq_len(ng)) for (j in 1:3) Gm[i, j, ] <- Mx[i, ((j - 1) * 4 + 1):(j * 4)]
  means <- apply(Gm, c(1, 2), np_mean); res <- t(vapply(seq_len(ng), function(i) unlist(lapply(1:3, function(j) Gm[i, j, ] - means[i, j])), numeric(12))); Cm <- cor(t(res)); iu <- Cm[t(upper.tri(Cm))]      # the upper triangle in row-major order (numpy triu_indices)
  sds <- np_mean(as.numeric(t(apply(Gm, c(1, 2), np_std1)))); dd <- np_mean(means[, 2] - means[, 1])      # the 15 SDs flattened in C order
  desc[[cs[2]]] <- list(within_group_sd_of_dCq = py_round(sds, 3L), mean_residual_r_between_genes = py_round(np_mean(iu), 3L), mean_DAPT_minus_CTRL_dCq = py_round(dd, 3L))
  plog(LG, py_sprintf("   %-80s within-group SD of dCq %.3f | mean residual r between the 5 genes %+.3f | mean (DAPT - CTRL) dCq %+.3f", cs[1], desc[[cs[2]]]$within_group_sd_of_dCq, desc[[cs[2]]]$mean_residual_r_between_genes, desc[[cs[2]]]$mean_DAPT_minus_CTRL_dCq))
}
out$pairing_descriptive <- desc
## E. 9 November
ex9 <- cq_table(F_N9); plog(LG); plog(LG, "E. 9 November run: matrix rows against the instrument export (exported order)"); same_seq <- 0L; same_set9 <- 0L; nb9 <- 0L; asc9_sheet <- 0L; asc9_exp <- 0L
for (g in genes9) for (k in G) { v <- vapply(ss(g, k), function(z) z$t, 1); it <- ex9[[key_of(norm(g), k)]]; e <- if (is.null(it)) numeric(0) else vapply(it, function(z) z$cq, 1); if (length(e) != 4) next
  nb9 <- nb9 + 1L; same_seq <- same_seq + as.integer(all(abs(v - e) <= TOL)); same_set9 <- same_set9 + as.integer(all(abs(sort(v) - sort(e)) <= TOL)); asc9_sheet <- asc9_sheet + as.integer(nondecr(v)); asc9_exp <- asc9_exp + as.integer(nondecr(e)) }
plog(LG, py_sprintf("   blocks compared: %d; same four values in the same row order as the export: %d; same set of values: %d; ascending in the sheet: %d; ascending in the export: %d", nb9, same_seq, same_set9, asc9_sheet, asc9_exp))
out$nov09 <- list(blocks = nb9, same_sequence_as_export = same_seq, same_set_as_export = same_set9, ascending_in_sheet = asc9_sheet, ascending_in_export = asc9_exp)
write_txt(LG$lines, file.path(RES, "aug26_nov09_order_check.txt")); write_json_typed(out, file.path(RES, "aug26_nov09_order_check.json"))
