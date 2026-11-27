## fingerprint_inputs.R - R port of wp11_plate_provenance/fingerprint_inputs.py: WP11 (8 Nov 2026): sha256 of every raw input that the WP11 scripts read, so the results can be tied to the exact files.
## Paths are relative to the project folder (CMBL_Sub2_291026). Usage: Rscript fingerprint_inputs.R   (writes results_R/wp11/INPUT_FINGERPRINTS.txt)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp11_lib.R"))
relpath <- function(path, start) {      # os.path.relpath
  sp <- function(x) strsplit(normalizePath(x, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]; p <- sp(path); s <- sp(start); n <- 0L
  while (n < min(length(p), length(s)) && tolower(p[n + 1L]) == tolower(s[n + 1L])) n <- n + 1L
  paste(c(rep("..", length(s) - n), p[seq_len(length(p) - n) + n]), collapse = "/")
}
PLATES <- c("2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1", "2020-Jun-03_Gr2_13genes-Ctrl-DAPT-DKK1", "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2"); files <- list(); add <- function(p, d) files[[length(files) + 1L]] <<- c(p, d)
for (ver in c("Before", "After")) for (pl in PLATES) for (kind in c("Quantification Cq Results.xlsx", "Quantification Amplification Results.xlsx")) add(file.path(P, "Exported_data", ver, paste0(pl, " -  ", kind)), sprintf("user export of 8 Nov 2026, threshold %s", if (ver == "Before") "as found" else "150"))
for (pl in PLATES) add(file.path(P, "Exported_data", "After", paste0(gsub(" 2", "", pl, fixed = TRUE), "_Threshold150.pcrd")), "user file saved with threshold 150 (8 Nov 2026)")      # the saved .pcrd of the second 24 June plate lost the ' 2' of its name (its exports keep it)
add(file.path(P, "Database", "STRING12.5", "string_interactions_21_nodes.tsv"), "user's STRING 12.5 web export, 8 Nov 2026"); add(file.path(P, "Database", "STRING12.5", "STRING12.5.txt"), "network URL of that export")
for (f in list(c("2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx", "26 Aug instrument export (2025 re-export of the .pcrd, wells in plate order)"),
               c("2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx", "26 Aug workbook with the user's calculation columns"), c("2020-11-09_ITPR, RYR1-3 -  Quantification Cq Results.xlsx", "9 Nov instrument export"))) add(file.path(RT, f[1]), f[2])
for (f in list(c("2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx", "24 June plate 1, 2025 re-export"), c("2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 -  Quantification Amplification Results.xlsx", "24 June plate 1, 2025 re-export (curves)"))) add(file.path(SBA, f[1]), f[2])
for (f in c("ins1_sheet_cq.csv", "gene_run_map.csv")) add(file.path(P, "CMBL_Sub2", "strengthen_v2", "wp1_ins1_structure", f), "WP1 matrix rows")
out <- format(Sys.time(), "Inputs read by the WP11 scripts; sha256, then description; fingerprinted %Y-%m-%d %H:%M:%S %z")
for (f in files) out <- c(out, sprintf("%s *%s", sha256_file(f[1]), relpath(f[1], P)), f[2])
write_txt(out, file.path(RES, "INPUT_FINGERPRINTS.txt")); cat(length(files), "files fingerprinted\n")
