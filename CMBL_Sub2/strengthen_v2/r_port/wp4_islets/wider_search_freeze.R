## wider_search_freeze.R - R port of wp4_islets/wider_search_freeze.py: B2/B3 of PRESPEC_WP4b_positive_control_and_wider_search.md: the keep / drop decision with its reason for EVERY one of the candidates of wider_search/candidates_all.csv, made from
## title, summary, design text and sample records only (no expression value was read), and the contrast of every kept series. Generic rules (organism, assay, sample number, stem-cell flag) give the mechanical reasons; the series read in full
## carry a specific reason. The file is fingerprinted (sha256 in ../PROTOCOL_v2_SHA256.txt) BEFORE any expression value of a candidate is read. Input: WS_FROZEN/candidates_all.csv (the archived list of the three queries of
## wider_search.R, frozen on 6 Nov 2026), output WP4_WS_FREEZE_OUT (default results_R/wp4/wider_search_freeze): frozen_candidates.csv.   Usage: Rscript wider_search_freeze.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R"))
WS <- Sys.getenv("WS_FROZEN", normalizePath(file.path(R_PORT, "..", "wp4_islets", "wider_search"))); OUT <- Sys.getenv("WP4_WS_FREEZE_OUT", file.path(R_PORT, "results_R", "wp4", "wider_search_freeze")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
D <- read.csv(file.path(WS, "candidates_all.csv"), stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
STEM <- "differentiat|iPSC|\\bESC\\b|embryonic stem|pluripoten|stem cell|stem-cell|progenitor|reprogram"
K <- function(sign, mirror, arm, weak, contrast, pairing, route, k4, why) list(sign = sign, mirror = mirror, arm = arm, weak = weak, contrast = contrast, pairing = pairing, route = route, k4 = k4, why = why)      # kept series: expected sign (+1 inhibition, -1 mirror), arm, contrast, pairing, data route, K4 status
KEEP <- list(
  GSE60803 = K(1L, 0L, "DKK-1", "confounded: controls are not expanded cells (addendum G)", "as fixed in addendum G of PRESPEC_WP4_choices.md: per donor, control = untreated or nontarget array, perturbed = mean of the beta-catenin shRNA arrays", "by donor", "series matrix (GPL11532); already analysed", "obtainable",
               "read: adult human beta-cell-derived cells, beta-catenin shRNA against control (Wnt inhibition), 4 donors (K1, K2, K3 met)"),
  GSE167073 = K(1L, 0L, "DKK-1", "receptor-level loss (Fzd9) in isolated islets of an insulinoma model under MycER activation", "perturbed = samples titled \"Fzd9 ko\" (4); control = samples titled \"Fzd9 wt\" (4)", "none (unpaired)", "series matrix values table (41346 rows, GPL17400); platform annotation from the GEO platform text export", "obtainable (annotation to be fetched at the analysis step)",
                "read: Fzd9 knock-out against wild type in islets of MycER-BclXL mice, 4 vs 4 (K1 islets of an insulinoma model, K2 Wnt receptor loss, K3 met)"),
  GSE99503 = K(-1L, 1L, "DKK-1", "GSK3 inhibitor 1-azakenpaullone is not specific to the Wnt pathway", "perturbed = \"Neg siRNA + 1-AKP #1-3\" (scramble siRNA plus GSK3 inhibitor); control = \"Neg siRNA #1-3\" (scramble siRNA, DMSO); the TUNAR siRNA arms are not used", "none (passages are not stated as pairs)", "NCBI-generated raw counts (GRCh38.p13, Entrez ids; HEAD 200, 566 KB)", "obtainable",
               "read: EndoC-betaH1 cells, GSK3 inhibition (Wnt activation, a MIRROR design) against vehicle in scramble-siRNA cells, 3 vs 3 (K1, K2, K3 met)"),
  GSE193888 = K(-1L, 1L, "DAPT", "", "perturbed = \"NICD-1..4\" (tetO_NICD+, doxycycline 8 weeks); control = \"Control-1..4\" (tetO_NICD-)", "none (unpaired)", "none by plain download: no supplementary file list (HTTP 404) and no NCBI-generated counts for mouse (HTTP 404)", "NOT obtainable",
                "read: beta-cell Notch gain of function in mouse islets (a MIRROR design for the DAPT arm), 4 vs 4 (K1, K2, K3 met); data not obtainable"),
  GSE249170 = K(-1L, 1L, "DAPT", "", "perturbed = \"Jag1-1..3\" (tetO_Jag1+, doxycycline 8 weeks); control = \"Control-1..3\" (tetO_Jag1-)", "none (unpaired)", "none by plain download: no supplementary file list (HTTP 404) and no NCBI-generated counts for mouse (HTTP 404)", "NOT obtainable",
                "read: beta-cell Jagged1 gain of function in mouse islets (a MIRROR design for the DAPT arm), 3 vs 3 (K1, K2, K3 met); data not obtainable"),
  `E-MTAB-2034` = K(1L, 0L, "DKK-1", "", "perturbed = \"TCF7L2 knockdown_Replicate 1-4\"; control = \"Control sample-Replicate 1-4\"; paired by replicate (passage) as named in the sample records", "by replicate (passages 62, 63, 81, 82)", "none: the ArrayExpress record holds raw reads only (ENA); no processed expression table; reads are not aligned here", "NOT obtainable",
                    "read: TCF7L2 knockdown against control in rat INS-1 832/13 cells, 4 vs 4 (K1, K2, K3 met; the closest model to the INS-1 experiment); no processed data"))
DROP <- c(GSE44045 = "read, sample records: the arrays compare islets of rats fed a standard diet with a cafeteria diet (the series title refers to Sfrp5, the array experiment is a diet comparison); no Notch or Wnt perturbation (K2)",
          GSE44046 = "read, sample records: same study as GSE44045 (30-day arm), a diet comparison; no Notch or Wnt perturbation (K2)", GSE44047 = "read: SuperSeries of GSE44045 and GSE44046 (diet comparison); no Notch or Wnt perturbation (K2)",
          GSE12802 = "read, sample records: compounds (BayK8644, exendin-4, their combination) against DMSO over a time course; no Wnt or Notch arm (K2)", GSE30526 = "read: COUP-TFII loss or gain of function; beta-catenin is a downstream target, the perturbation is of a transcription factor (K2)",
          GSE64444 = "read: MIN6 culture conditions that change Hes3; Hes3 is not a Notch pathway perturbation (K2)", GSE64449 = "read: Hes3 knockdown in MIN6 cells; Hes3 is not a Notch pathway perturbation (K2)",
          GSE7430 = "read: conditional beta-catenin knockout in the developing exocrine pancreas; embryonic pancreas tissue, not islets or beta-cells (K1)", GSE76189 = "read: Tspan2 and high glucose in beta-cells; no Notch or Wnt perturbation (K2)",
          GSE186944 = "read: Cdh13 null islets; no Notch or Wnt perturbation (K2)", `E-MEXP-3836` = "read: insulin receptor overexpression in INS-1 cells; Wnt activation is a consequence, not the perturbation (K2)",
          `E-MTAB-5325` = "title: pericytes of Tcf7l2-mutated mice; not islets or beta-cells (K1)", GSE24296 = "title: anti-DLL4 treatment, immune cell populations; not islets or beta-cells (K1)",
          GSE41284 = "title: Tcf7l2 in liver (K1)", GSE193888_dummy = "", GSE139147 = "read: HDAC inhibitor in BRIN-BD11 cells; no Notch or Wnt perturbation (K2)", GSE224490 = "read: BCAS2 splicing factor; no Notch or Wnt perturbation (K2)",
          GSE195647 = "read: mammary stem cell ageing (K1)", GSE136887 = "read: single-cell gain-of-function screen of T2D regulators in islets; no Notch or Wnt perturbation (K2)",
          GSE125151 = "title: serpin B13 in islets; no Notch or Wnt perturbation (K2)", GSE134272 = "title: glucolipotoxicity and a protective compound in rat islets; no Notch or Wnt perturbation (K2)",
          GSE17647 = "title: 4E-BP1 and HDL in beta-cells (K2)", GSE19877 = "title: paternal obesity and offspring islets (K2)", GSE44373 = "title: obese and lean rat islets, diet (K2)", GSE142572 = "title: hPSCs treated with CHIR99021; stem cells, not beta-cells (K1)")
rows <- list(); na0 <- function(x) if (is.na(x)) "" else as.character(x)
for (i in seq_len(nrow(D))) {
  r <- D[i, ]; acc <- r$accession; n <- suppressWarnings(as.numeric(r$n_samples)); org <- as.character(r$organism); assay <- as.character(r$assay); ae <- identical(as.character(r$entrytype), "ArrayExpress")
  stem <- grepl(STEM, paste(as.character(r$title), as.character(r$summary)), ignore.case = TRUE, perl = TRUE); base <- list(accession = acc, sources = r$sources, pmids = r$pmids, organism = org, assay = assay, n_samples = r$n_samples, title = r$title)
  if (acc %in% names(KEEP)) { k <- KEEP[[acc]]; rows[[length(rows) + 1]] <- c(base, list(decision = "keep", reason = k$why, expected_sign = as.character(k$sign), mirror = as.character(k$mirror), arm = k$arm, weak_match_note = k$weak, contrast = k$contrast, pairing = k$pairing, data_route = k$route, k4 = k$k4)); next }
  if (acc %in% names(DROP) && nzchar(DROP[[acc]])) { dec <- "drop"; why <- DROP[[acc]] }
  else if (!grepl("Homo sapiens|Mus musculus|Rattus norvegicus", org)) { dec <- "drop"; why <- "organism is not human, mouse or rat (K1)" }
  else if (!(ae || grepl("Expression profiling", assay, fixed = TRUE))) { dec <- "drop"; why <- sprintf("not an expression profile (assay: %s) (K1)", assay) }
  else if (!ae && !is.na(n) && n < 6) { dec <- "drop"; why <- "fewer than 6 samples, so fewer than 3 per arm (K3)" }
  else if (stem) { dec <- "drop"; why <- "stem-cell, progenitor or reprogramming design (K1; flag from title or summary, checked against the title)" }
  else { dec <- "drop"; why <- "no Notch or Wnt perturbation of islets or beta-cells in the title or summary (K2)" }
  rows[[length(rows) + 1]] <- c(base, list(decision = dec, reason = why, expected_sign = "", mirror = "", arm = "", weak_match_note = "", contrast = "", pairing = "", data_route = "", k4 = ""))
}
F <- do.call(rbind, lapply(rows, function(r) as.data.frame(lapply(r, function(v) if (is.null(v)) NA else v), stringsAsFactors = FALSE)))
for (cn in c("expected_sign", "mirror", "arm", "weak_match_note", "contrast", "pairing", "data_route", "k4")) F[[cn]] <- as.character(F[[cn]])
py_write_csv_repr(F, file.path(OUT, "frozen_candidates.csv"), eol = "\r\n"); cat("candidates", nrow(F), "| keep", sum(F$decision == "keep"), "| drop", sum(F$decision == "drop"), "\n")
