## effsens.R - R port of reanalysis_17genes/effsens.py: efficiency sensitivity at group level from the CFX Cq values of all plates (cfx_long.json of melt.R) with the mean amplicon efficiencies of Table S1
## (Additional file 1_CMBL_working.docx, read through its XML; AF1_DOCX overrides): fold change of DAPT and DKK-1 against CTRL per plate with E = 100% and with the amplicon efficiencies, geometric mean over plates,
## ratio, flag. Writes effsens.json and efficiency_sensitivity_group_level.csv (UTF-8 with BOM, CRLF) into RA_OUT; prints the table.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "docx_table.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
AF1 <- Sys.getenv("AF1_DOCX", normalizePath(file.path(R_PORT, "..", "..", "..", "Supplementary", "Additional file 1_CMBL_working.docx")))
rows <- fromJSON(file.path(OUT, "cfx_long.json"), simplifyVector = FALSE)
EFF <- list()
for (c in docx_table0_rows(AF1)[-(1:2)]) {
  if (length(c) < 7 || c[1] == c[2]) next
  m <- .allm("\\d+\\.\\d+|\\d+", c[7]); if (length(m)) EFF[[strsplit(trimws(c[1]), "\\s+")[[1]][1]]] <- as.numeric(m[1]) / 100
}
EXCL <- c("Hey1", "Wnt2", "Wnt2b", "Ryr1")
grp <- function(s) { s <- tolower(if (is.null(s)) "None" else as.character(s)); if (grepl("dapt", s, fixed = TRUE)) "DAPT" else if (grepl("dkk", s, fixed = TRUE)) "DKK" else if (grepl("ctrl", s, fixed = TRUE) || grepl("control", s, fixed = TRUE)) "CTRL" else NA_character_ }
Cq <- list()
for (r in rows) { g <- grp(r$sample); if (!is.na(g) && is.numeric(r$cq)) { k <- paste(substr(r$run, 1, 22), r$gene, g, sep = "\t"); Cq[[k]] <- c(Cq[[k]], r$cq) } }
keys <- do.call(rbind, strsplit(names(Cq), "\t", fixed = TRUE)); plates <- py_sort(unique(keys[, 1])); genes <- py_sort(setdiff(unique(keys[, 2]), c("Gapdh", EXCL)))
fc <- function(plate, gene, treat, E_t, E_g) {
  need <- c(paste(plate, gene, treat, sep = "\t"), paste(plate, gene, "CTRL", sep = "\t"), paste(plate, "Gapdh", treat, sep = "\t"), paste(plate, "Gapdh", "CTRL", sep = "\t"))
  if (any(vapply(need, function(k) length(Cq[[k]]) == 0, TRUE))) return(NULL)
  dt <- mean(Cq[[need[1]]]) - mean(Cq[[need[2]]]); dg <- mean(Cq[[need[3]]]) - mean(Cq[[need[4]]])
  (1 + E_t)^(-dt) / ((1 + E_g)^(-dg))
}
cat(sprintf("Gapdh efficiency used: %.3f\n", EFF[["Gapdh"]]))
cat(sprintf("%-8s %5s %6s | DAPT: FC(E=1) FC(E_amp) ratio | DKK-1: FC(E=1) FC(E_amp) ratio | flag\n", "gene", "E_t", "plates"))
out <- list()
for (g in genes) {
  rec <- list(gene = g, E = EFF[[g]])
  for (t in c("DAPT", "DKK")) {
    a <- Filter(function(x) !is.null(x) && x != 0, lapply(plates, function(p) fc(p, g, t, 1.0, 1.0))); b <- Filter(function(x) !is.null(x) && x != 0, lapply(plates, function(p) fc(p, g, t, EFF[[g]], EFF[["Gapdh"]])))
    a <- unlist(a); b <- unlist(b)
    rec[[t]] <- list(exp(mean(log(a))), exp(mean(log(b))), length(a), as.list(round(a, 3)), as.list(round(b, 3)))
  }
  flags <- character(0)
  for (t in c("DAPT", "DKK")) { A <- rec[[t]][[1]]; B <- rec[[t]][[2]]
    if ((A - 1) * (B - 1) < 0) flags <- c(flags, paste0(t, ": direction flips")) else if (abs(log2(B / A)) > 0.585) flags <- c(flags, paste0(t, ": >1.5x shift")) else if (abs(log2(B / A)) > 0.32) flags <- c(flags, paste0(t, ": >1.25x shift")) }
  rec$flag <- paste(flags, collapse = "; "); out[[length(out) + 1]] <- rec
  cat(sprintf("%-8s %5.2f %6d | %7.2f %8.2f %6.2f | %7.2f %8.2f %6.2f | %s\n", g, if (is.null(EFF[[g]])) 0 else EFF[[g]], rec$DAPT[[3]], rec$DAPT[[1]], rec$DAPT[[2]], rec$DAPT[[2]] / rec$DAPT[[1]], rec$DKK[[1]], rec$DKK[[2]], rec$DKK[[2]] / rec$DKK[[1]], rec$flag))
}
write_text_bytes(py_json_typed(out, indent = NULL), file.path(OUT, "effsens.json"))
rd <- function(x) round(x, 3)
recs <- lapply(out, function(r) list(gene = r$gene, mean_efficiency = r$E, plates = r$DAPT[[3]], DAPT_FC_E100 = rd(r$DAPT[[1]]), DAPT_FC_Eamplicon = rd(r$DAPT[[2]]), DAPT_ratio = rd(r$DAPT[[2]] / r$DAPT[[1]]),
                                    DKK1_FC_E100 = rd(r$DKK[[1]]), DKK1_FC_Eamplicon = rd(r$DKK[[2]]), DKK1_ratio = rd(r$DKK[[2]] / r$DKK[[1]]), flag = r$flag))
dictwriter_csv(recs, c("gene", "mean_efficiency", "plates", "DAPT_FC_E100", "DAPT_FC_Eamplicon", "DAPT_ratio", "DKK1_FC_E100", "DKK1_FC_Eamplicon", "DKK1_ratio", "flag"), file.path(OUT, "efficiency_sensitivity_group_level.csv"), fmt = "%.12g", bom = TRUE)
