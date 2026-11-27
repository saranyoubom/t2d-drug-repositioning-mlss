## export_tables.R - R port of wp11_plate_provenance/export_tables.py: WP11 (9 Nov 2026): structured copies of WP11 results for Additional file 6 (so that no sheet has to be parsed from prose).
## Writes, in results_R/wp11/:  aug26_workbook_vs_instrument.csv (72 rows), threshold_validation.csv, refA_block_residuals.csv, refA_dcq_comparison.csv.
## Inputs: the 26 August instrument export and the user's 'Quantification Summary' workbook (as in check_aug26_nov09_order.R), results_R/wp11/threshold_fit.txt, results_R/wp11/refA_vs_plate1.txt (made by the other WP11 scripts).
## Read only; usage: Rscript export_tables.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R")); TOL <- 0.0051
F_EXP <- file.path(RT, "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx"); F_SUM <- file.path(RT, "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx")
ex <- cq_table(F_EXP); inst <- list()      # (well, gene) -> (group, Cq)
for (key in names(ex)) { sk <- split_key(key); for (z in ex[[key]]) inst[[paste(z$well, sk[1], sep = "\r")]] <- list(k = sk[2], c = z$cq) }
wsum <- sheet0(F_SUM); rows <- list()
for (r in 2:length(wsum)) { cells <- wsum[[r]]; get <- function(c) if (c <= length(cells)) cells[[c]] else NA; w <- get(2); tg <- get(4); sm <- get(6); cq <- get(7); if (is.na(w) || !isnum(cq)) next
  rows[[length(rows) + 1L]] <- list(r = r, w = as.character(w), tg = as.character(tg), k = grp(sm), cq = as.numeric(cq)) }
blocks <- list(); for (x in rows) { key <- key_of(norm(x$tg), x$k); blocks[[key]] <- c(blocks[[key]], x$cq) }
iblocks <- lapply(ex, function(it) sort(vapply(it, function(z) z$cq, 1))); recs <- list()
for (x in rows) {
  g <- norm(x$tg); ii <- inst[[paste(x$w, g, sep = "\r")]]; key <- key_of(g, x$k); bv <- sort(blocks[[key]]); ib <- iblocks[[key]]; same_set <- length(bv) == length(ib) && all(abs(bv - ib) <= TOL)
  b <- blocks[[key]]; asc <- all(b[-length(b)] <= b[-1]); note <- if (g == "gapdh" && identical(x$k, "DAPT")) "DAPT Gapdh: workbook holds the values of C10 and D10 twice; instrument wells A10 and B10 (24.987 and 24.132) are not in the workbook" else ""
  recs[[length(recs) + 1L]] <- list(workbook_row = as.integer(x$r), well = x$w, target = x$tg, group_in_workbook = x$k, workbook_Cq = py_round(x$cq, 6L), instrument_group_of_named_well = ii$k, instrument_Cq_of_named_well = py_round(ii$c, 6L), abs_difference = py_round(abs(x$cq - ii$c), 6L),
                                    workbook_Cq_equals_instrument_Cq_of_named_well = as.integer(abs(x$cq - ii$c) <= TOL), block_value_set_equals_instrument = as.integer(same_set), block_ascending_down_the_workbook = as.integer(asc), note = note)
}
cols <- names(recs[[1]]); dictwriter_csv(recs, cols, file.path(RES, "aug26_workbook_vs_instrument.csv"), fmt = "repr")
## ---- parse threshold_fit.txt
tf <- readLines(file.path(RES, "threshold_fit.txt"), warn = FALSE, encoding = "UTF-8"); mt <- function(pat, ln) { m <- regmatches(ln, regexec(pat, ln, perl = TRUE))[[1]]; if (length(m)) m[-1] else NULL }
val <- list(); for (ln in tf) { m <- mt("^\\s+(After|Before)\\s+(.+?)\\s+known T ([\\d.]+)\\s+fitted T ([\\d.]+)\\s+mean \\|recomputed - exported Cq\\| ([\\d.]+) \\(max ([\\d.]+)\\)", ln)
  if (!is.null(m)) val[[length(val) + 1L]] <- list(export_version = m[1], plate = m[2], known_threshold_RFU = m[3], fitted_threshold_RFU = m[4], mean_abs_difference_cycles = m[5], max_abs_difference_cycles = m[6]) }
dictwriter_csv(val, names(val[[1]]), file.path(RES, "threshold_validation.csv"))
res <- list(); for (ln in tf) { m <- mt("^\\s{2}(Wnt\\w+|Lef1|Tcf7|Gapdh\\(refA\\))\\s+(skipped|[\\d.]+) / (skipped|[\\d.]+) / (skipped|[\\d.]+)\\s*$", ln)
  if (!is.null(m)) for (i in 1:3) res[[length(res) + 1L]] <- list(assay = m[1], group = c("CTRL", "DAPT", "DKK-1")[i], largest_abs_difference_of_4_sorted_values_cycles_at_fitted_threshold = m[1 + i]) }
dictwriter_csv(res, names(res[[1]]), file.path(RES, "refA_block_residuals.csv"))
rv <- readLines(file.path(RES, "refA_vs_plate1.txt"), warn = FALSE, encoding = "UTF-8"); dq <- list()
for (ln in rv) { m <- mt("^(\\w+)\\s+(CTRL|DAPT|DKK-1)\\s+([\\d.]+) \\|\\s+([\\d.]+)\\s+([\\d.]+)\\s+([+-][\\d.]+) \\|\\s+([\\d.]+)\\s+([\\d.]+)\\s+([+-][\\d.]+)\\s*$", ln)
  if (!is.null(m)) dq[[length(dq) + 1L]] <- list(gene = m[1], group = m[2], mean_Cq_sheet = m[3], dCq_sheet = m[4], "dCq_reexport_45.70RFU" = m[5], dCq_sheet_minus_reexport = m[6], mean_Cq_sheet_again = m[7], "dCq_curves_at_209.91RFU" = m[8], dCq_sheet_minus_curves = m[9]) }
dictwriter_csv(dq, names(dq[[1]]), file.path(RES, "refA_dcq_comparison.csv"))
cat(py_sprintf("rows: aug26 %d, validation %d, residuals %d, dcq %d\n", length(recs), length(val), length(res), length(dq)))
cat(py_sprintf("equal named-well Cq: %d of %d; blocks with the same value set: %d of %d rows\n", sum(vapply(recs, function(r) r$workbook_Cq_equals_instrument_Cq_of_named_well, 1L)), length(recs), sum(vapply(recs, function(r) r$block_value_set_equals_instrument, 1L)), length(recs)))
