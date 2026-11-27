## check_wp11_numbers.R - R port of wp11_plate_provenance/check_wp11_numbers.py: WP11 (8 Nov 2026): independent check that the numbers of README.md and WORDING_provenance.md are the numbers of the result files, and that three key facts hold when
## recomputed from the raw exports. Part 1 (transcription): every number is looked up in the saved result files (results_R/wp11/*.txt, *.json; the folder WP11_OUT, which can also be the archived results/ to test the check itself) and
## must appear in the README / wording text (the archived texts of wp11_plate_provenance/). Part 2 (recomputation from raw files): 26 August DAPT Gapdh wells, well/Cq agreement of the workbook, and the equality of the plate-2 blocks with the sheet.
## Read only. Usage: Rscript check_wp11_numbers.R   (writes check_wp11_last_run.txt into WP11_CHECK_OUT, default results_R/wp11)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R")); suppressMessages(library(jsonlite)); CHK_OUT <- Sys.getenv("WP11_CHECK_OUT", file.path(R_PORT, "results_R", "wp11")); dir.create(CHK_OUT, showWarnings = FALSE, recursive = TRUE)
rd_text <- function(p) { b <- readBin(p, "raw", file.size(p)); s <- rawToChar(b); Encoding(s) <- "UTF-8"; gsub("\r\n", "\n", s, fixed = TRUE) }
README <- rd_text(file.path(ARC, "README.md")); WORDING <- if (file.exists(file.path(ARC, "WORDING_provenance.md"))) rd_text(file.path(ARC, "WORDING_provenance.md")) else ""; TXT <- paste0(README, "\n", WORDING)
out <- character(0); fails <- 0L
check <- function(name, ok, detail = "") { ok <- isTRUE(ok); out <<- c(out, py_sprintf("%s  %s %s", if (ok) "PASS" else "FAIL", name, detail)); if (!ok) fails <<- fails + 1L }
has <- function(s) grepl(enc2utf8(s), TXT, fixed = TRUE, useBytes = TRUE); read <- function(f) rd_text(file.path(RES, f)); inn <- function(s, x) grepl(s, x, fixed = TRUE)
mt <- function(pat, s) { m <- regmatches(s, regexec(pat, s, perl = TRUE))[[1]]; if (length(m)) m else NULL }; tup <- function(g) paste0("('", paste(g, collapse = "', '"), "')")
eqn <- function(x, y) identical(unlist(x), y)      # a list of numbers equal to a named / plain vector
## ---------------------------------------------------------------- part 1: transcription
tf <- read("threshold_fit.txt"); rv <- read("refA_vs_plate1.txt"); ao <- read("aug26_nov09_order_check.txt"); aj <- fromJSON(file.path(RES, "aug26_nov09_order_check.json"), simplifyVector = FALSE); t150 <- read("threshold150_summary.txt")
sj <- fromJSON(file.path(Sys.getenv("WP10_RES", file.path(R_PORT, "results_R", "wp10")), "string_web_export_check.json"), simplifyVector = FALSE)      # made by wp10_string_web_export_check.R (byte-identical to the archived file)
m <- mt("T = ([\\d.]+) RFU, mean \\|recomputed - sheet\\| = ([\\d.]+) cycles over (\\d+) values", tf)
check("fit: best T 209.91, mean 0.0398 (0.040), 92 values", !is.null(m) && abs(as.numeric(m[2]) - 209.91) < 0.006 && abs(as.numeric(m[3]) - 0.0398) < 0.0006 && m[4] == "92" && has("209.91") && has("0.040 cycles over 92 values"), if (!is.null(m)) m[1] else "")
check("fit: plate-1 export threshold 45.700", inn("threshold of the export = 45.700", tf) && has("45.70"))
for (kf in list(c("150.000", "150.000"), c("138.061", "138.060"), c("36.072", "36.070"))) check(sprintf("fit validation: known %s returns %s", kf[1], kf[2]), grepl(sprintf("known T %s\\s+fitted T %s", gsub(".", "\\.", kf[1], fixed = TRUE), gsub(".", "\\.", kf[2], fixed = TRUE)), tf, perl = TRUE))
check("fit: Gapdh-only T 211.59, target-only 209.53, five genes 209.91 (0.0293)", inn("Gapdh blocks T = 211.59", tf) && inn("target blocks T = 209.53", tf) && inn("T = 209.91 (mean 0.0293)", tf) && has("211.6") && has("209.5") && has("0.029"))
check("fit: failed well G04 of Wnt2b DAPT, peak 89.6", inn("G04 (wnt2b DAPT, 89.6 RFU)", tf) && has("G04") && has("89.6"))
dm <- mt("sheet minus re-export at 45\\.70 RFU, 15 gene-group means: mean ([+-][\\d.]+), SD ([\\d.]+), range ([+-][\\d.]+) to ([+-][\\d.]+)", rv)
check("refA: dCq re-export diff +0.22 (SD 0.11; +0.06 to +0.51)", !is.null(dm) && identical(dm[-1], c("+0.22", "0.11", "+0.06", "+0.51")) && has("+0.22") && has("+0.06 to +0.51"), if (!is.null(dm)) tup(dm[-1]) else "")
dt <- mt("sheet minus curves at 209\\.91 RFU, 15 gene-group means: mean ([+-][\\d.]+), SD ([\\d.]+), range ([+-][\\d.]+) to ([+-][\\d.]+)", rv)
check("refA: dCq diff at 209.91 -0.03 (SD 0.04; -0.10 to +0.02)", !is.null(dt) && identical(dt[-1], c("-0.03", "0.04", "-0.10", "+0.02")) && has("-0.03") && has("-0.10 to +0.02"), if (!is.null(dt)) tup(dt[-1]) else "")
d1 <- mt("sheet minus re-export \\(45\\.70 RFU\\): mean ([+-][\\d.]+), SD ([\\d.]+), largest absolute ([\\d.]+)", rv); d2 <- mt("sheet minus curves at 209\\.91 RFU: mean ([+-][\\d.]+), SD ([\\d.]+), largest absolute ([\\d.]+)", rv)
check("refA: ddCq diff re-export -0.14 (SD 0.09; 0.30)", !is.null(d1) && identical(d1[-1], c("-0.14", "0.09", "0.30")), if (!is.null(d1)) tup(d1[-1]) else "")
check("refA: ddCq diff at 209.91 +0.007 (SD 0.025; 0.042)", !is.null(d2) && identical(d2[-1], c("+0.007", "0.025", "0.042")) && has("+0.007") && has("0.025") && has("0.042"), if (!is.null(d2)) tup(d2[-1]) else "")
sw <- aj$summary_workbook
check("26 Aug: 32 of 72 rows equal the instrument value of the named well (target 26 of 60, Gapdh 6 of 12)", sw$cq_equal_to_instrument_value_of_named_well == 32 && sw$rows == 72 && eqn(sw$target_rows, c(26L, 60L)) && eqn(sw$gapdh_rows, c(6L, 12L)) && has("32 of 72") && has("26 of 60") && has("6 of 12"))
check("26 Aug: 40 of 72 rows no longer match", 72 - sw$cq_equal_to_instrument_value_of_named_well == 40 && has("40 of 72"))
check("26 Aug: 17 of 18 blocks ascending in the workbook; plate order ascending 1 of 15", inn("17 of 18", ao) && eqn(aj$matrix_aug26$plate_order_ascending_instrument, c(1L, 15L)) && has("17 of 18") && has("1 of the 15"))
check("26 Aug: matrix target blocks ascending 15 of 15 and same set 15 of 15", eqn(aj$matrix_aug26$target_blocks_ascending, c(15L, 15L)) && eqn(aj$matrix_aug26$target_blocks_same_set_as_instrument, c(15L, 15L)))
check("26 Aug: only the Gapdh DAPT block differs in value set", identical(unlist(sw$blocks_different), "Gapdh DAPT") && sw$blocks_same_value_set_as_instrument == 17 && sw$blocks_total == 18)
pd <- aj$pairing_descriptive; pk <- function(a, b, c) c(within_group_sd_of_dCq = a, mean_residual_r_between_genes = b, mean_DAPT_minus_CTRL_dCq = c)
check("pairing table: as analyzed 0.132 / +0.536 / +0.662", eqn(pd$sheet, pk(0.132, 0.536, 0.662)) && has("| 0.132 | +0.536 | +0.662 |"))
check("pairing table: rank instrument 0.397 / +0.975 / -0.276", eqn(pd$rank_instrument, pk(0.397, 0.975, -0.276)) && has("| 0.397 | +0.975 | -0.276 |"))
check("pairing table: plate-order instrument 0.481 / +0.938 / -0.276", eqn(pd$plate_instrument, pk(0.481, 0.938, -0.276)) && has("| 0.481 | +0.938 | -0.276 |"))
check("9 Nov: 11 blocks, same sequence 11, ascending 0", eqn(aj$nov09, c(blocks = 11L, same_sequence_as_export = 11L, same_set_as_export = 11L, ascending_in_sheet = 0L, ascending_in_export = 0L)) && has("11 gene-by-group blocks"))
check("threshold-150 shifts: Gr1 +2.23, Gr2 +1.97, Round2 +0.16", inn("Gr1     n 78  median +2.23", t150) && inn("Gr2     n 78  median +1.97", t150) && inn("Round2  n 72  median +0.16", t150) && has("+0.16"))
before <- strsplit(strsplit(t150, "=== exports with the threshold as found (Before) ===", fixed = TRUE)[[1]][2], "=== shift", fixed = TRUE)[[1]][1]; g5 <- c("Hes1", "Hey1", "Glut2", "Ins1", "Glp1r")
ok5 <- all(vapply(g5, function(g) grepl(sprintf("%s\\s+target Cq blocks equal to the export \\(CTRL / DAPT / DKK-1\\): Round2 / Round2 / Round2", g), before, perl = TRUE), TRUE))
check("Before: Hes1, Hey1, Glut2, Ins1, Glp1r target and Gapdh blocks = Round2 (plate 2) in all three groups", ok5 && all(vapply(g5, function(g) grepl(sprintf("%s\\s+reference blocks equal to Gapdh of: Round2 / Round2 / Round2", g), before, perl = TRUE), TRUE)))
tl <- strsplit(t150, "\n", fixed = TRUE)[[1]]; sel <- tl[grepl("target Cq blocks equal", tl, fixed = TRUE) & vapply(tl, function(l) any(startsWith(trimws(l), c("Wnt2 ", "Wnt2b", "Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7"))), TRUE)]
check("Wnt-family genes match no block of Gr1, Gr2 or plate 2 at either threshold", all(grepl("- / - / -", sel, fixed = TRUE)) && inn("Hes1   target Cq blocks equal to the export (CTRL / DAPT / DKK-1): - / - / -", t150))
alias <- sj$classification_from_web_export_12_5$alias; ea <- list("Known Validated" = 5, "Novel Discovery" = 16, "STRING Only" = 34, "Not Validated" = 65, in_STRING = 39, coverage_percent = 32.5)
check("STRING: web export = API 12.5 (59 pairs, 0 score difference, same classification)", sj$distinct_pairs$web_12_5 == 59 && sj$web_vs_api_12_5$max_abs_score_difference_in_common_pairs == 0 && isTRUE(sj$classification_web_equals_api) && identical(sort(names(alias)), sort(names(ea))) && all(unlist(alias)[sort(names(alias))] == unlist(ea)[sort(names(alias))]) && has("5 / 16 / 34 / 65"))
cc <- sj$web_12_5_vs_old_12_0
check("STRING: old 12.0 export, 54 of 55 common pairs changed score, 4 gained, 4 lost", cc$pairs_in_both == 55 && cc$score_change_in_common_pairs$n_changed == 54 && length(cc$pairs_gained_in_12_5) == 4 && length(cc$pairs_lost_in_12_5) == 4 && has("54 of the 55") && has("4 pairs are gained and 4 lost"))
## ---------------------------------------------------------------- part 2: recomputation from the raw files (the reader of wp11_lib.R, not the table() of the other scripts)
ex <- sheet0(file.path(RT, "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx")); h <- vapply(ex[[1]], function(x) if (is.na(x)[1]) "None" else as.character(x), ""); iw <- match("Well", h); it <- match("Target", h); ic <- match("Cq", h); inst <- list()
for (r in ex[-1]) if (!blank(r[[iw]]) && isnum(r[[ic]])) inst[[paste(as.character(r[[iw]]), tolower(strsplit(as.character(r[[it]]), "_", fixed = TRUE)[[1]][1]), sep = "\r")]] <- r[[ic]]
gw <- vapply(c("A10", "B10", "C10", "D10"), function(w) inst[[paste(w, "gapdh", sep = "\r")]], 1)
check("raw: 26 Aug Gapdh DAPT wells A10 24.987, B10 24.132, C10 22.569, D10 22.798", all(abs(gw - c(24.987, 24.132, 22.569, 22.798)) < 0.0006), "[" %+% paste(vapply(py_round(gw, 3L), function(x) py_sprintf("%s", x), ""), collapse = ", ") %+% "]")
sm <- sheet0(file.path(RT, "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx")); n <- 0L; okn <- 0L; nt <- 0L; okt <- 0L
for (r in sm[-1]) { w <- r[[2]]; tg <- r[[4]]; cq <- r[[7]]; if (is.na(w) || !isnum(cq)) next; n <- n + 1L; same <- abs(inst[[paste(as.character(w), tolower(as.character(tg)), sep = "\r")]] - cq) < 0.0051; okn <- okn + as.integer(same); if (tolower(as.character(tg)) != "gapdh") { nt <- nt + 1L; okt <- okt + as.integer(same) } }
check("raw: workbook Cq equals the instrument Cq of the named well in 32 of 72 rows (target 26 of 60)", identical(c(n, okn, nt, okt), c(72L, 32L, 60L, 26L)), sprintf("(%d, %d, %d, %d)", n, okn, nt, okt))
dap <- unlist(lapply(sm[-1], function(r) if (identical(r[[4]], "Gapdh") && identical(r[[6]], "DAPT")) r[[7]] else NULL))
check("raw: workbook Gapdh DAPT cells 22.569, 22.798, 22.569, 22.798", identical(py_round(dap, 3L), c(22.569, 22.798, 22.569, 22.798)), "[" %+% paste(vapply(py_round(dap, 3L), function(x) py_sprintf("%s", x), ""), collapse = ", ") %+% "]")
b <- sheet0(file.path(EXD, "Before", "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2 -  Quantification Cq Results.xlsx")); hh <- vapply(b[[1]], function(x) if (is.na(x)[1]) "None" else as.character(x), ""); iw <- match("Well", hh); it <- match("Target", hh); isamp <- match("Sample", hh); ic <- match("Cq", hh); blk <- list()
for (r in b[-1]) { if (is.na(r[[it]]) || !isnum(r[[ic]])) next; ts <- .target_sample(r[[it]], r[[isamp]]); g <- norm(sub("^[0-9]+[.]\\s*", "", ts$t)); sl <- tolower(if (is.na(ts$s)[1]) "none" else as.character(ts$s)); k <- if (grepl("ctrl", sl, fixed = TRUE)) "CTRL" else if (grepl("dapt", sl, fixed = TRUE)) "DAPT" else "DKK-1"
  blk[[key_of(g, k)]] <- c(blk[[key_of(g, k)]], r[[ic]]) }
sh <- list(); S <- read_csv_chr(file.path(W1, "ins1_sheet_cq.csv")); for (i in seq_len(nrow(S))) { k <- key_of(norm(S$gene[i]), S$group[i]); sh[[k]] <- c(sh[[k]], list(c(as.numeric(S$target_cq[i]), as.numeric(S$reference_cq[i])))) }
allclose <- function(a, b, atol) { a <- sort(a); b <- sort(b); length(a) == length(b) && all(abs(a - b) <= atol + 1e-5 * abs(b)) }      # numpy.allclose with the default rtol
eq <- 0L; eqg <- 0L; for (g in c("hes1", "hey1", "glut2", "ins1", "glp1r")) for (k in c("CTRL", "DAPT", "DKK-1")) { x <- sh[[key_of(g, k)]]; eq <- eq + as.integer(allclose(vapply(x, function(z) z[1], 1), blk[[key_of(g, k)]], 0.0051)); eqg <- eqg + as.integer(allclose(vapply(x, function(z) z[2], 1), blk[[key_of("gapdh", k)]], 0.0051)) }
check("raw: plate-2 export (threshold as found) equals the sheet in 15 of 15 target blocks and 15 of 15 Gapdh vectors", identical(c(eq, eqg), c(15L, 15L)), sprintf("(%d, %d)", eq, eqg))
tabs <- read_csv_chr(file.path(RES, "aug26_workbook_vs_instrument.csv")); bad <- tabs$block_value_set_equals_instrument == "0"
check("export_tables: aug26 csv has 72 rows, 32 with the Cq of the named well, 4 rows outside a matching value set (Gapdh DAPT)", nrow(tabs) == 72 && sum(as.integer(tabs$workbook_Cq_equals_instrument_Cq_of_named_well)) == 32 && sum(bad) == 4 && all(tabs$target[bad] == "Gapdh" & tabs$group_in_workbook[bad] == "DAPT"))
vt <- read_csv_chr(file.path(RES, "threshold_validation.csv")); dq <- read_csv_chr(file.path(RES, "refA_dcq_comparison.csv")); rs <- read_csv_chr(file.path(RES, "refA_block_residuals.csv"))
check("export_tables: 4 validation rows, 15 dCq rows (mean sheet-minus-curves -0.03), 24 residual rows", nrow(vt) == 4 && nrow(dq) == 15 && py_round(np_sum(as.numeric(dq$dCq_sheet_minus_curves)) / 15, 2L) == -0.03 && nrow(rs) == 24)
out <- c(out, py_sprintf("%d checks, %d failures", length(out), fails)); write_txt(out, file.path(CHK_OUT, "check_wp11_last_run.txt")); cat(paste(out, collapse = "\n"), "\n", sep = "")
