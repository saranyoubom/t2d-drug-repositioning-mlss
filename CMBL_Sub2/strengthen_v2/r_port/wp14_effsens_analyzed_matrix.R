## wp14_effsens_analyzed_matrix.R - R port of wp14_table_s8/wp14_effsens_analyzed_matrix.py: Table S8 (efficiency sensitivity at group level) recomputed from the Cq values of the
## ANALYZED matrix (wp1_ins1_structure/ins1_sheet_cq.csv) with the mean amplicon efficiencies read from Table S1 of Additional file 1 (a Word file, read through its XML).
## Writes results/efficiency_sensitivity_analyzed_matrix.csv into R_OUT_DIR (default results_R/wp14). Read only otherwise.
args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "pycompat.R"))
STR <- normalizePath(file.path(here, "..")); PROJ <- normalizePath(file.path(STR, "..", ".."))
AF1 <- Sys.getenv("AF1_DOCX", file.path(PROJ, "Supplementary", "Additional file 1_JTM_working.docx")); OUTDIR <- Sys.getenv("R_OUT_DIR", file.path(here, "results_R", "wp14"))
G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
html_unescape <- function(s) {
  s <- gsub("&lt;", "<", s, fixed = TRUE); s <- gsub("&gt;", ">", s, fixed = TRUE); s <- gsub("&quot;", "\"", s, fixed = TRUE); s <- gsub("&apos;", "'", s, fixed = TRUE)
  m <- gregexpr("&#[0-9]+;", s); for (x in unique(unlist(regmatches(s, m)))) s <- gsub(x, intToUtf8(as.integer(sub("&#([0-9]+);", "\\1", x))), s, fixed = TRUE)
  gsub("&amp;", "&", s, fixed = TRUE)
}
allm <- function(p, x) { r <- regmatches(x, gregexpr(p, x, perl = TRUE))[[1]]; if (length(r)) r else character(0) }
table0_rows <- function(path) {
  tmp <- tempfile(); dir.create(tmp); unzip(path, files = "word/document.xml", exdir = tmp)
  x <- readChar(file.path(tmp, "word", "document.xml"), file.size(file.path(tmp, "word", "document.xml")), useBytes = TRUE); Encoding(x) <- "UTF-8"
  tbl <- regmatches(x, regexpr("(?s)<w:tbl>.*?</w:tbl>", x, perl = TRUE)); lapply(allm("(?s)<w:tr[ >].*?</w:tr>", tbl), function(tr)
    vapply(allm("(?s)<w:tc>.*?</w:tc>", tr), function(tc) paste(vapply(allm("<w:t(?: [^>]*)?>[^<]*</w:t>", tc), function(t) html_unescape(sub("</w:t>$", "", sub("^<w:t(?: [^>]*)?>", "", t, perl = TRUE))), ""), collapse = ""), ""))
}
EFF <- list(); rows <- table0_rows(AF1)
for (c in rows[3:length(rows)]) {
  if (length(c) < 7 || c[1] == c[2]) next
  m <- allm("\\d+\\.\\d+|\\d+", c[7]); sym <- regmatches(c[1], regexpr("^[A-Za-z0-9]+", c[1]))
  if (length(m) && length(sym)) EFF[[sym]] <- as.numeric(m[1]) / 100
}
missing <- setdiff(c(G17, "Gapdh"), names(EFF))
cq <- read.csv(file.path(STR, "wp1_ins1_structure", "ins1_sheet_cq.csv"), stringsAsFactors = FALSE)
out <- list()
for (g in G17) {
  rec <- list(gene = g, E_gene = EFF[[g]], E_Gapdh = EFF[["Gapdh"]]); ct <- cq[cq$gene == g & cq$group == "CTRL", ]
  for (tg in list(c("DAPT", "DAPT"), c("DKK-1", "DKK1"))) {
    t <- tg[1]; tag <- tg[2]; tt <- cq[cq$gene == g & cq$group == t, ]
    dT <- mean(tt$target_cq) - mean(ct$target_cq); dG <- mean(tt$reference_cq) - mean(ct$reference_cq)
    fc100 <- 2^(-(dT - dG)); fcE <- (1 + EFF[[g]])^(-dT) / (1 + EFF[["Gapdh"]])^(-dG)
    rep <- mean(2^(-(tt$target_cq - tt$reference_cq))) / mean(2^(-(ct$target_cq - ct$reference_cq)))
    rec[[paste0(tag, "_FC_reported")]] <- rep; rec[[paste0(tag, "_FC_E100")]] <- fc100; rec[[paste0(tag, "_FC_Eamplicon")]] <- fcE; rec[[paste0(tag, "_ratio")]] <- fcE / fc100
  }
  flags <- character(0)
  for (tag in c("DAPT", "DKK1")) { A <- rec[[paste0(tag, "_FC_E100")]]; B <- rec[[paste0(tag, "_FC_Eamplicon")]]
    if ((A - 1) * (B - 1) < 0) flags <- c(flags, paste0(tag, ": direction flips")) else if (abs(log2(B / A)) > 0.585) flags <- c(flags, paste0(tag, ": >1.5x shift")) else if (abs(log2(B / A)) > 0.32) flags <- c(flags, paste0(tag, ": >1.25x shift")) }
  rec$flag <- paste(flags, collapse = "; "); out[[g]] <- rec
}
rd <- function(x) round(x, 3)
df <- do.call(rbind, lapply(out, function(r) data.frame(gene = r$gene, mean_amplicon_efficiency = rd(r$E_gene), DAPT_FC_reported = rd(r$DAPT_FC_reported), DAPT_FC_E100 = rd(r$DAPT_FC_E100), DAPT_FC_Eamplicon = rd(r$DAPT_FC_Eamplicon), DAPT_ratio = rd(r$DAPT_ratio),
  DKK1_FC_reported = rd(r$DKK1_FC_reported), DKK1_FC_E100 = rd(r$DKK1_FC_E100), DKK1_FC_Eamplicon = rd(r$DKK1_FC_Eamplicon), DKK1_ratio = rd(r$DKK1_ratio), flag = r$flag, stringsAsFactors = FALSE)))
py_write_csv(df, file.path(OUTDIR, "efficiency_sensitivity_analyzed_matrix.csv"), nd = 4)
flips <- vapply(out, function(r) grepl("flips", r$flag), TRUE)
cat("efficiencies missing for:", missing, "| E_Gapdh", EFF[["Gapdh"]], "\n"); cat("direction flips:", names(out)[flips], "| any shift > 1.25x:", sum(vapply(out, function(r) nzchar(r$flag), TRUE)), "of 17 genes | max |log2 ratio|:",
  sprintf("%.3f", max(abs(log2(unlist(lapply(out, function(r) c(r$DAPT_ratio, r$DKK1_ratio))))))), "\n")
