## check_threshold150.R - R port of wp11_plate_provenance/check_threshold150.py: WP11 (8 Nov 2026): plate provenance of the 24 June / 3 June sheet values with the exports that the user re-made on 8 Nov 2026 from the native .pcrd files
## (Exported_data/Before = threshold as found: auto 36.072 / auto 48.532 / manual 138.061; Exported_data/After = manual threshold 150 RFU for every run).
## Question (author items A12): do the target and reference (Gapdh) Cq values that the calculation sheet 'Mix_3+24 June 2020' holds for the nine June genes of the analysis come from the plates 2020-Jun-03 Gr1 and Gr2
## (two wells per group each) and 2020-Jun-24 Round2, once the threshold of the export is 150? Read only; writes results_R/wp11/threshold150_matches.csv and threshold150_summary.txt.
## Sheet values: the WP1 matrix rows ins1_sheet_cq.csv and gene_run_map.csv. Usage: Rscript check_threshold150.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R")); LG <- plog_new(); G <- c("CTRL", "DAPT", "DKK-1")
RUNS <- c(Gr1 = "2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx", Gr2 = "2020-Jun-03_Gr2_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx", Round2 = "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2 -  Quantification Cq Results.xlsx")
EX <- list(); for (v in c("Before", "After")) { EX[[v]] <- list(); for (r in names(RUNS)) EX[[v]][[r]] <- cq_table(file.path(EXD, v, RUNS[[r]])) }
sheet <- list(); meta <- list(); M <- read_csv_chr(file.path(W1, "gene_run_map.csv")); for (i in seq_len(nrow(M))) meta[[M$gene[i]]] <- M$run[i]
S <- read_csv_chr(file.path(W1, "ins1_sheet_cq.csv")); for (i in seq_len(nrow(S))) { k <- key_of(S$gene[i], S$group[i]); if (is.null(sheet[[k]])) sheet[[k]] <- list(); sheet[[k]][[length(sheet[[k]]) + 1L]] <- list(rep = as.integer(S$replicate_in_group[i]), t = as.numeric(S$target_cq[i]), r = as.numeric(S$reference_cq[i])) }
sheet_sorted <- function(g, k) { x <- sheet[[key_of(g, k)]]; x[order(vapply(x, function(z) z$rep, 1L))] }      # sorted(tuples): by the replicate number
genes <- names(meta)[vapply(meta, function(r) r == "Jun24", TRUE)]; out <- list()
eq <- function(a, b, tol = 0.0051) { a <- sort(a); b <- sort(b); length(a) == length(b) && all(abs(a - b) <= tol) }
SOURCES <- list(Round2 = "Round2", "Gr1+Gr2" = c("Gr1", "Gr2"), Gr1 = "Gr1", Gr2 = "Gr2")
export_cq <- function(ver, runs, gene, k) unlist(lapply(runs, function(rn) { x <- EX[[ver]][[rn]][[key_of(gene, k)]]; if (is.null(x)) numeric(0) else vapply(x, function(z) z$cq, 1) }))
for (ver in c("After", "Before")) {
  plog(LG, "=== exports with the threshold " %+% (if (ver == "After") "150 (After)" else "as found (Before)") %+% " ===")
  for (g in genes) { res <- character(0)
    for (k in G) { tv <- vapply(sheet_sorted(g, k), function(z) z$t, 1); hit <- character(0)
      for (nm in names(SOURCES)) { ex <- export_cq(ver, SOURCES[[nm]], norm(g), k); if (length(ex) == 4 && eq(tv, ex)) hit <- c(hit, nm) }
      res <- c(res, if (length(hit)) paste(hit, collapse = "+") else "-"); out[[length(out) + 1L]] <- list(ver, g, k, paste(py_fmt_f(tv, 2L), collapse = ","), paste(hit, collapse = ";")) }
    plog(LG, py_sprintf("  %-6s target Cq blocks equal to the export (CTRL / DAPT / DKK-1): %s", g, paste(res, collapse = " / "))) }
  plog(LG, "  reference (Gapdh) vector of each gene against the Gapdh wells of the exports:")
  for (g in genes) { res <- character(0)
    for (k in G) { rv <- vapply(sheet_sorted(g, k), function(z) z$r, 1); hit <- character(0)
      for (nm in names(SOURCES)) { ex <- export_cq(ver, SOURCES[[nm]], "gapdh", k); if (length(ex) == 4 && eq(rv, ex)) hit <- c(hit, nm) }
      res <- c(res, if (length(hit)) paste(hit, collapse = "+") else "-") }
    plog(LG, py_sprintf("  %-6s reference blocks equal to Gapdh of: %s", g, paste(res, collapse = " / "))) }
  plog(LG)
}
## shift caused by the threshold: After minus Before on the same wells
plog(LG, "=== shift of the Cq values by the threshold (After minus Before, same wells, all targets) ===")
for (rn in names(RUNS)) { d <- numeric(0)
  for (key in names(EX$After[[rn]])) { b <- EX$Before[[rn]][[key]]; bw <- if (is.null(b)) list() else setNames(lapply(b, function(z) z$cq), vapply(b, function(z) z$well, "")); for (z in EX$After[[rn]][[key]]) if (!is.null(bw[[z$well]])) d <- c(d, z$cq - bw[[z$well]]) }
  plog(LG, py_sprintf("  %-7s n %d  median %+.2f  range %+.2f to %+.2f", rn, length(d), np_median(d), min(d), max(d))) }
## replicate -> plate and well
plog(LG, "=== which plate and well carries each sheet replicate (After exports; value match within 0.005) ===")
for (g in genes) { parts <- character(0)
  for (k in G) { cell <- character(0)
    for (z in sheet_sorted(g, k)) { m <- character(0); for (rn in names(RUNS)) for (w in EX$After[[rn]][[key_of(norm(g), k)]]) if (abs(w$cq - z$t) <= 0.0051) m <- c(m, paste0(rn, ":", w$well)); cell <- c(cell, sprintf("%d=%s", z$rep, if (length(m) == 1) m[1] else if (!length(m)) "?" else paste(m, collapse = "/"))) }
    parts <- c(parts, paste(k, paste(cell, collapse = " "))) }
  plog(LG, py_sprintf("  %-6s %s", g, paste(parts, collapse = " | "))) }
recs <- lapply(out, function(o) setNames(o, c("export_version", "gene", "group", "sheet_target_cq_sorted", "export_source_with_equal_values")))
dictwriter_csv(recs, c("export_version", "gene", "group", "sheet_target_cq_sorted", "export_source_with_equal_values"), file.path(RES, "threshold150_matches.csv")); write_txt(LG$lines, file.path(RES, "threshold150_summary.txt"))
