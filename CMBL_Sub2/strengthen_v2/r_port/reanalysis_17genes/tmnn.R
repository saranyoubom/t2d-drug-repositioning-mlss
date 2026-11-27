## tmnn.R - R port of reanalysis_17genes/tmnn.py: nearest-neighbour melting temperatures of the primers of Table S1 against the Tm printed there (error mean, SD, maximum) and of the four new Hes1 / Tcf7 primers.
## The Tm function itself is p_tm() of primer_lib.R. Prints only.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "docx_table.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); source(file.path(HERE, "primer_lib.R"))
AF1 <- Sys.getenv("AF1_DOCX", normalizePath(file.path(R_PORT, "..", "..", "..", "Supplementary", "Additional file 1.docx")))
errs <- list()
for (c in p_table_rows(AF1)) {
  s <- ws_split(c[4]); tms <- as.numeric(.allm("\\d+\\.\\d+", c[6])); if (length(s) < 2 || length(tms) < 2) next
  for (k in 1:2) errs[[length(errs) + 1]] <- list(ws_split(c[1])[1], round(p_tm(s[k]), 2), tms[k], round(p_tm(s[k]) - tms[k], 2))
}
e <- vapply(errs, function(x) x[[4]], 1); cat("n", length(e), "mean err", round(mean(e), 2), "SD", round(np_std(e), 2), "max abs", max(abs(e)), "\n")
for (x in errs[1:8]) cat(sprintf("('%s', %s, %s, %s)\n", x[[1]], format(x[[2]], nsmall = 1), format(x[[3]], nsmall = 1), format(x[[4]], nsmall = 1)))
for (nm in list(c("r_Hes1_Fw_NEW", "CACCGGACAAACCAAAGACAGC"), c("r_Hes1_Rv_NEW", "GAGGTGCTTCACTGTCATTTCC"), c("r_Tcf7_Fw_NEW", "GAGCACACTTCGCAGAGACTTTTC"), c("r_Tcf7_Rv_NEW", "GTAGTGCATGAGCAGATTGAAGGC"))) cat(nm[1], round(p_tm(nm[2]), 2), nchar(nm[2]), "\n")
