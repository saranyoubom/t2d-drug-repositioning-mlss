## perturbation_decisions.R - R port of wp4_islets/perturbation_decisions.py: protocol 4.6, step 2: the keep / drop decision with its reason for EVERY candidate of perturbation/candidates.csv, made from title, organism, assay type and summary text
## only (no expression value was read). Rules (PRESPEC_WP4_choices.md section D): keep = expression profiling of islets, beta-cell lines or beta-cell-derived cells with an INHIBITION of the Notch (gamma-secretase) or Wnt (DKK1, Wnt inhibitor,
## beta-catenin or TCF7L2 loss) pathway against a control, at least 3 samples per arm, on a platform that covers the 17 genes. Everything else is dropped with the reason below. Series read in full before deciding are marked 'read'.
## Input: PERT_FROZEN/candidates.csv (the archived list of the two fixed queries of perturbation_scout.R); output WP4_PERT_DEC_OUT (default results_R/wp4/perturbation_decisions): candidates_decisions.csv.   Usage: Rscript perturbation_decisions.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R"))
PF <- Sys.getenv("PERT_FROZEN", normalizePath(file.path(R_PORT, "..", "wp4_islets", "perturbation"))); OUT <- Sys.getenv("WP4_PERT_DEC_OUT", file.path(R_PORT, "results_R", "wp4", "perturbation_decisions")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
D <- read.csv(file.path(PF, "candidates.csv"), stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
MANUAL <- list(
  GSE60803 = c("keep", "read: adult human beta-cell-derived cells, beta-catenin shRNA against control (WNT inhibition), several donors; expanded adult cells, not a stem-cell differentiation design (the regex flag fires on the word redifferentiation); sample arms to be checked in the series matrix before analysis"),
  GSE44045 = c("drop", "CORRECTED 6 Nov (addendum H): the sample records show that the arrays compare islets of rats fed a standard diet with a cafeteria diet (the series title refers to Sfrp5); no Notch or Wnt perturbation (K2). The first version of this reason described the series as Sfrp5 silencing (Wnt activation) from the title and summary; that description is withdrawn"),
  GSE44046 = c("drop", "CORRECTED 6 Nov (addendum H): same study as GSE44045 (30-day arm), a diet comparison in the sample records; no Notch or Wnt perturbation (K2)"),
  GSE30526 = c("drop", "read: COUP-TFII loss or gain of function; beta-catenin is a downstream target, the perturbation is of a transcription factor, not of the Wnt pathway"),
  GSE12802 = c("drop", "read: small-molecule screen for beta-cell proliferation in an immortalised mouse line; compounds ACTIVATE Wnt (thiophene-pyrimidines) or L-type calcium channels; no pathway inhibition arm"),
  GSE193888 = c("drop", "read: doxycycline-inducible beta-cell Notch GAIN of function in mouse islets; direction opposite to DAPT"), GSE249170 = c("drop", "read: beta-cell Jagged1 gain and loss of function in mice (ligand level, in vivo); not a gamma-secretase or pathway inhibition"),
  GSE139147 = c("drop", "read: HDAC inhibitor MS-275 in BRIN-BD11 cells; not a Notch or Wnt perturbation"), GSE224490 = c("drop", "read: BCAS2 splicing factor in mouse islets; not a Notch or Wnt perturbation"),
  GSE136887 = c("drop", "read: single-cell gain-of-function screen of T2D regulators in human islets; not a Notch or Wnt inhibition"), GSE195647 = c("drop", "read: mammary stem cell ageing; not islets"),
  GSE41284 = c("drop", "title: Tcf7l2 in liver (perinatal hepatic response); not islets or beta-cells"), GSE53352 = c("drop", "title: haematopoietic stem cells deficient of glycogen synthase kinase; not islets"),
  GSE50103 = c("drop", "title: bi-potent pancreas progenitor expansion (Lgr5 / R-spondin); progenitor design"), GSE14 = c("drop", "title: Cancer Genome Anatomy Project SAGE libraries; not a perturbation design"),
  GSE118725 = c("drop", "title: transcription-factor binding to non-coding variants (reporter and binding assays); not an expression perturbation"), GSE126550 = c("drop", "title: saturation mutagenesis of regulatory elements; reporter assay"),
  GSE308943 = c("drop", "title: microRNA profiling of rat islet cells; non-coding RNA assay, no Notch or Wnt perturbation"), GSE323218 = c("drop", "title: MPRA control from MIN6; reporter assay"), GSE323299 = c("drop", "title: MPRA from MIN6; reporter assay"),
  GSE228665 = c("drop", "title: isogenic hESC knock-outs of T2D genes and their differentiation; stem-cell design"))
rows <- list()
for (i in seq_len(nrow(D))) {
  r <- D[i, ]; if (r$accession %in% names(MANUAL)) { dec <- MANUAL[[r$accession]][1]; why <- MANUAL[[r$accession]][2] }
  else if (!grepl("Expression profiling", as.character(r$assay), fixed = TRUE)) { dec <- "drop"; why <- sprintf("not an expression profile (assay: %s)", r$assay) }
  else if (isTRUE(as.logical(r$stem_cell_flag))) { dec <- "drop"; why <- "title / summary describe a stem-cell, progenitor, reprogramming or lineage-differentiation design (flag from the regex, checked against the title)" }
  else { dec <- "drop"; why <- "no Notch or Wnt inhibition arm in the title or summary" }
  rows[[i]] <- data.frame(accession = r$accession, queries = r$queries, organism = r$organism, assay = r$assay, n_samples = r$n_samples, platform = r$platform, decision = dec, reason = why, title = r$title, stringsAsFactors = FALSE)
}
O <- do.call(rbind, rows); py_write_csv_repr(O, file.path(OUT, "candidates_decisions.csv"), eol = "\r\n"); cat(sprintf("%s %d\n", names(table(O$decision)), as.integer(table(O$decision))), "\n")
