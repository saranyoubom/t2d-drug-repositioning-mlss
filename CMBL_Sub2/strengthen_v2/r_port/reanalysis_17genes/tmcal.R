## tmcal.R - R port of reanalysis_17genes/tmcal.py: calibration of the melt temperature measured on the CFX (cfx_summary.json of melt.R) against the one predicted from the amplicon
## (81.5 + 0.41 GC% - 675/length - 21.6), over the well-matched genes, and the predicted melt temperature of seven candidate primer pairs on their reference sequences. Prints only.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "docx_table.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); source(file.path(HERE, "primer_lib.R"))
AF1 <- Sys.getenv("AF1_DOCX", normalizePath(file.path(R_PORT, "..", "..", "..", "Supplementary", "Additional file 1.docx")))
summ <- fromJSON(file.path(OUT, "cfx_summary.json"), simplifyVector = FALSE)
prim <- list(); for (c in p_table_rows(AF1)) { s <- ws_split(c[4]); if (length(s) >= 2) prim[[ws_split(c[1])[1]]] <- list(c[2], s[1], s[2]) }
amp <- function(acc, F, R) { fz <- p_fetch(acc); f <- p_find(fz$seq, F); r <- p_find(fz$seq, p_rc(R)); if (is.na(f$mm) || is.na(r$mm)) return(NULL); substr(fz$seq, f$pos + 1, r$pos + nchar(R)) }
pts <- list()
for (g in names(prim)) {
  if (g %in% c("Hes1", "Hey1", "Ryr3")) next
  a <- amp(prim[[g]][[1]], prim[[g]][[2]], prim[[g]][[3]]); if (is.null(a) || is.null(summ[[g]]) || is.null(summ[[g]]$tm) || g == "Ryr1") next
  pts[[length(pts) + 1]] <- list(g, nchar(a), round(p_pred(a), 2), summ[[g]]$tm)
}
X <- vapply(pts, function(p) p[[3]], 1); Y <- vapply(pts, function(p) p[[4]], 1); fit <- lm(Y ~ X); a0 <- unname(coef(fit)[1]); b <- unname(coef(fit)[2]); res <- unname(residuals(fit)); sd_ <- sqrt(sum(res^2) / (length(X) - 2))
cat(sprintf("calibration over %d genes: Tm_meas = %.2f + %.2f * Tm_pred ; residual SD %.2f C ; r = %.2f\n", length(pts), a0, b, sd_, cor(X, Y)))
for (i in seq_along(pts)) cat(sprintf("   %-8s N=%3d pred=%6.2f meas=%6.2f resid=%+.2f\n", pts[[i]][[1]], pts[[i]][[2]], pts[[i]][[3]], pts[[i]][[4]], res[i]))
cand <- list(`Hes1 printed(canine pair) on rat NM_024360` = list("NM_024360", "GAGAAGGCGGACATTCTGGA", "ACCTCGTTCATACACTCGCTG", summ$Hes1$tm), `Hes1 r_NEW on rat NM_024360` = list("NM_024360", "CACCGGACAAACCAAAGACAGC", "GAGGTGCTTCACTGTCATTTCC", summ$Hes1$tm),
             `Hes1 r_old (r_Hes1_Fw/Rv)` = list("NM_024360", "CTGGAGAAGGCAGACATTCTGG", "GTCACCTCGTTCATGCACTCG", summ$Hes1$tm), `Tcf7 printed` = list("NM_001427430", "CGTCCTTGATGCTGGGATCTG", "GATGATTGGCTTCTTCGCCTCC", summ$Tcf7$tm),
             `Tcf7 r_NEW` = list("NM_001427430", "GAGCACACTTCGCAGAGACTTTTC", "GTAGTGCATGAGCAGATTGAAGGC", summ$Tcf7$tm), `Ryr3 printed on XM_087613279` = list("XM_087613279", "CAACCTACTACCCCGAAGGAAC", "GAGTGGGCTCTGCTGTTAGG", summ$Ryr3$tm),
             `Hey1 printed(canine) on rat NM_001191845` = list("NM_001191845", "ACCTGAAAATGCTGCACACG", "GCTGGGAGGCGTAGTTGTTA", summ$Hey1$tm))
cat("\n")
for (k in names(cand)) {
  v <- cand[[k]]; fz <- p_fetch(v[[1]]); f <- p_find(fz$seq, v[[2]]); r <- p_find(fz$seq, p_rc(v[[3]])); a <- substr(fz$seq, f$pos + 1, r$pos + nchar(v[[3]])); p <- a0 + b * p_pred(a)
  cat(sprintf("%-44s N=%3d predicted Tm=%5.1f  measured=%5.1f  diff=%+.1f  (%+.1f SD)\n", k, nchar(a), p, v[[4]], v[[4]] - p, (v[[4]] - p) / sd_))
}
