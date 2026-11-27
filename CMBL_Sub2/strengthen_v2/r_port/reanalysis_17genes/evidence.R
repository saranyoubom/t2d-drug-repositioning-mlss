## evidence.R - R port of reanalysis_17genes/evidence.py: primer evidence table (reference sequence, species, Tm of both primers, mismatches, amplicon length against the table, melt temperature
## measured vs predicted from the amplicon, corrections, inclusion in the analysis) from Table S1 of Additional file 1.docx, NCBI sequences (primer_lib.R, cached) and cfx_summary.json (melt.R).
## Writes primer_evidence.csv (UTF-8 with BOM) and evidence.json into RA_OUT. Needs the NCBI cache or network access for the sequences.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "docx_table.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE); source(file.path(HERE, "primer_lib.R"))
AF1 <- Sys.getenv("AF1_DOCX", normalizePath(file.path(R_PORT, "..", "..", "..", "Supplementary", "Additional file 1.docx")))
summ <- fromJSON(file.path(OUT, "cfx_summary.json"), simplifyVector = FALSE)
rows <- list()
for (c in p_table_rows(AF1)) { s <- ws_split(c[4]); if (length(s) >= 2) rows[[length(rows) + 1]] <- list(gene = ws_split(c[1])[1], acc = c[2], F = s[1], R = s[2], L = as.integer(.allm("\\d+", c[5])[1]), eff = c[7]) }
## corrections (reasons in the evidence file)
CORR <- list(Hes1 = list(acc = "NM_024360.4", F = "CACCGGACAAACCAAAGACAGC", R = "GAGGTGCTTCACTGTCATTTCC", note = "r_Hes1_Fw_NEW / r_Hes1_Rv_NEW (13 Feb 2020 primer sheet); printed pair was the canine K9_Hes1 design"),
             Tcf7 = list(acc = "NM_001427430.1", F = "GAGCACACTTCGCAGAGACTTTTC", R = "GTAGTGCATGAGCAGATTGAAGGC", note = "r_Tcf7_Fw_NEW / r_Tcf7_Rv_NEW (13 Feb 2020 primer sheet)"),
             Ryr3 = list(acc = "XM_087613279.1", note = "accession updated; primers unchanged (previous accession XM_039106806.1 no longer matches)"))
EXCL <- list(Hey1 = "identity unresolved (printed pair is a canine design; melt identical to Hes1)", Ryr1 = "non-specific melt (75.1 °C); 12 of 12 wells flagged by rule", Wnt2 = "17% of wells at Cq ≥ 35 or without melt peak",
             Wnt2b = "46% of wells at Cq ≥ 35 or without melt peak")
ampseq <- function(acc, F, R) { fz <- p_fetch(acc); f <- p_find(fz$seq, F); r <- p_find(fz$seq, p_rc(R)); if (is.na(f$mm) || is.na(r$mm)) return(list(a = NULL, mf = NULL, mr = NULL, h = fz$hdr)); list(a = substr(fz$seq, f$pos + 1, r$pos + nchar(R)), mf = f$mm, mr = r$mm, h = fz$hdr) }
base <- list()
for (r in rows) {
  g <- r$gene; cc <- CORR[[g]]; acc <- if (!is.null(cc$acc)) cc$acc else r$acc; F <- if (!is.null(cc$F)) cc$F else r$F; R <- if (!is.null(cc$R)) cc$R else r$R
  z <- ampseq(acc, F, R); if (g == "Hey1") z <- ampseq("NM_001191845.2", F, R)
  base[[g]] <- list(gene = g, acc = acc, F = F, R = R, amp = if (!is.null(z$a)) nchar(z$a) else NULL, Fmm = z$mf, Rmm = z$mr, hdr = substr(z$h, 2, 90), pred = if (!is.null(z$a)) p_pred(z$a) else NULL, tmF = round(p_tm(F), 2), tmR = round(p_tm(R), 2), eff = r$eff, old_acc = r$acc, old_L = r$L)
}
## calibration from well-matched genes
cal <- Filter(Negate(is.null), lapply(names(base), function(g) { v <- base[[g]]; if (!is.null(summ[[g]]) && !is.null(v$pred) && v$pred != 0 && !(g %in% c("Hes1", "Hey1", "Ryr3", "Ryr1", "Tcf7")) && !is.null(summ[[g]]$tm)) c(v$pred, summ[[g]]$tm) else NULL }))
X <- vapply(cal, function(x) x[1], 1); Y <- vapply(cal, function(x) x[2], 1); fit <- lm(Y ~ X); a0 <- unname(coef(fit)[1]); b <- unname(coef(fit)[2]); sd_ <- sqrt(sum(residuals(fit)^2) / (length(X) - 2))
cat("calibration genes", length(cal), "SD", round(sd_, 2), "\n")
out <- list()
for (g in names(base)) {
  v <- base[[g]]; s <- summ[[g]]; p <- if (!is.null(v$pred) && v$pred != 0) a0 + b * v$pred else NULL; has_tm <- !is.null(s) && !is.null(s$tm) && s$tm != 0
  out[[length(out) + 1]] <- list(gene = g, reference = v$acc, species = if (grepl("Rattus", v$hdr, fixed = TRUE)) "Rattus norvegicus" else if (grepl("Canis", v$hdr, fixed = TRUE)) "Canis lupus" else "?",
    forward = v$F, reverse = v$R, Tm_F = v$tmF, Tm_R = v$tmR, F_mismatch = v$Fmm, R_mismatch = v$Rmm, amplicon_bp = v$amp, table_bp_before = v$old_L,
    melt_mean = if (has_tm) round(s$tm, 2) else NULL, melt_sd = if (!is.null(s) && !is.null(s$sd)) round(s$sd, 2) else NULL, wells = if (!is.null(s)) s$n else NULL,
    wells_late_or_no_melt = if (!is.null(s)) max(s$late, s$nomelt) else NULL, predicted_melt = if (!is.null(p) && p != 0) round(p, 1) else NULL,
    melt_diff = if (has_tm && !is.null(p) && p != 0) round(s$tm - p, 1) else NULL, melt_diff_SD = if (has_tm && !is.null(p) && p != 0) round((s$tm - p) / sd_, 1) else NULL,
    correction = if (!is.null(CORR[[g]]$note)) CORR[[g]]$note else "", analysis = if (!is.null(EXCL[[g]])) paste0("excluded: ", EXCL[[g]]) else "included")
}
cols <- names(out[[1]]); dictwriter_csv(out, cols, file.path(OUT, "primer_evidence.csv"), fmt = "%.12g", bom = TRUE)
write_text_bytes(py_json_typed(list(rows = out, cal = list(a = a0, b = b, sd = sd_, n = length(cal))), indent = NULL), file.path(OUT, "evidence.json"))
nn <- function(x) if (is.null(x)) "None" else as.character(x)
for (o in out) cat(sprintf("%-8s %-15s %-6s mm=%s/%s amp=%s (was %s) melt=%s pred=%s diff=%s (%s SD) late/nomelt=%s/%s %s\n", o$gene, o$reference, substr(o$species, 1, 6), nn(o$F_mismatch), nn(o$R_mismatch), nn(o$amplicon_bp), nn(o$table_bp_before),
                          nn(o$melt_mean), nn(o$predicted_melt), nn(o$melt_diff), nn(o$melt_diff_SD), nn(o$wells_late_or_no_melt), nn(o$wells), substr(o$analysis, 1, 30)))
