## wp13_gapdh_versions.R - R port of wp13_gapdh_corrected/wp13_gapdh_versions.py (WP13, 9 Nov 2026, at the user's request: "I do not know why Gapdh wells A10 and B10 were replaced; can you calculate the corrected version"):
## the 26 August DAPT reference of the INS-1 matrix in three versions, carried through the pooled analysis and the Fig. 5 classification. No new rule: the versions were defined on 2 Nov 2026 (gapdh_exact) and are used unchanged.
##   reported : the DAPT replicates k = 1..4 of the five 26 August genes use the Gapdh values [22.5685, 22.7983, 22.5685, 22.7983] (wells C10 and D10, each entered twice; A10 and B10 not used)
##   C0       : A10 and B10 excluded WITHOUT duplicating: the mean of the two valid wells (22.6834) for all four DAPT replicates   ('corrected by exclusion')
##   C1       : all four wells (mean 23.6215) for all four DAPT replicates                                                        ('corrected, failed wells included')
## Reads only; writes gapdh_versions.json and gapdh_versions_summary.csv into R_OUT_DIR (default results_R/wp13). Usage: Rscript wp13_gapdh_versions.R
## Needs reanalysis_17genes/pipeline.R, wp10_string_sensitivity.R (sourced) and results_R/reanalysis_17genes/qc26_wells.csv (RA_OUT) from gapdh_qc_amplification.R.
suppressMessages({ library(jsonlite); library(readxl) })
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- HERE
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "reanalysis_17genes", "pipeline.R")); source(file.path(R_PORT, "wp10_string_sensitivity.R"))
RES13 <- Sys.getenv("R_OUT_DIR", file.path(R_PORT, "results_R", "wp13")); dir.create(RES13, showWarnings = FALSE, recursive = TRUE); RA <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes"))
G5 <- c("Tcf7l2", "Kcnj11", "Cacna1d", "Cacna1c", "Ptbp1")
HK <- c(22.5685093010354, 22.7982726003435, 22.5685093010354, 22.7982726003435)
RAW <- list(A10 = 24.9873996418614, B10 = 24.1318270358876, C10 = 22.5685093010354, D10 = 22.7982726003435)
REF <- list(C0 = mean(c(RAW$C10, RAW$D10)), C1 = mean(unlist(RAW)))
df0 <- P_load(P_GENES17)
variant <- function(ref) { d <- df0; for (g in G5) d[5:8, g] <- df0[5:8, g] * 2.0^(ref - HK); d }         # 2^-(Cq - ref) = 2^-(Cq - HK) * 2^(ref - HK)
x <- load_inputs(); oa <- x$old$node1; ob <- x$old$node2; sset_alias <- string_set(oa, ob, TRUE); sset_lit <- string_set(oa, ob, FALSE)
keyset <- function(l) vapply(l, function(p) paste(unlist(p), collapse = "-"), "")
out <- list()
for (name in c("reported", "C0", "C1")) {
  d <- if (name == "reported") df0 else variant(REF[[name]]); sch <- P_schemes(d); res <- P_sig_sets(sch); cons <- P_consensus(res, use_sign = FALSE); cn <- P_counts(res); adj <- P_adjusted(sch)
  deg <- list(); for (i in seq_len(nrow(cons))) for (a in cons[i, ]) deg[[a]] <- (if (is.null(deg[[a]])) 0L else deg[[a]]) + 1L
  h <- res[["Hes1-normalized"]]; w <- which(upper.tri(h$p) & !is.na(h$p) & h$p < 0.05 & abs(h$r) >= 0.7, arr.ind = TRUE)
  p07 <- setNames(h$r[w], pname(h$genes[w[, 1]], h$genes[w[, 2]])); consset <- if (nrow(cons)) pname(cons[, 1], cons[, 2]) else character(0)
  cl_a <- classify(as.list(p07), consset, sset_alias); cl_l <- classify(as.list(p07), consset, sset_lit)
  fc <- setNames(lapply(G5, function(g) mean(d[5:8, g]) / mean(d[1:4, g])), G5); dn <- names(deg); dv <- unlist(deg); dord <- order(-dv, dn, method = "radix")
  out[[name]] <- list(consensus = .pairs_list(cons), n_consensus = nrow(cons), sig = lapply(cn, function(v) v$sig), adj_q05 = lapply(adj, function(a) as.integer(sum(a$q[upper.tri(a$q)] < 0.05))), degree = as.list(dv[dord]),
                      hubs = as.list(dn[dv >= 5]), n_pairs_abs_r_ge_0_7 = length(p07), fig5_alias = counts(cl_a), fig5_literal = counts(cl_l), in_STRING_alias = as.integer(sum(names(p07) %in% sset_alias)), fc_dapt = fc,
                      r_cacna1d_glp1r = h$r[match("Cacna1d", h$genes), match("Glp1r", h$genes)])
}
base <- keyset(out$reported$consensus)
for (name in names(out)) { cs <- keyset(out[[name]]$consensus); out[[name]]$kept <- as.integer(length(intersect(cs, base))); out[[name]]$lost <- as.list(py_sort(setdiff(base, cs))); out[[name]]$new <- as.list(py_sort(setdiff(cs, base))) }
qc <- read.csv(file.path(RA, "qc26_wells.csv"), stringsAsFactors = FALSE); qc <- qc[qc$gene == "Gapdh", ]
out$qc_gapdh_wells <- lapply(seq_len(nrow(qc)), function(i) list(well = qc$well[i], sample = qc$sample[i], cq = round(qc$cq[i], 3), plateau_RFU = round(qc$plateau[i], 0)))
out$ref <- REF; out$hk_used <- HK; out$raw <- RAW
write_json_typed(out, file.path(RES13, "gapdh_versions.json"))
hdr <- c("version", "DAPT_reference_Gapdh_Cq", "consensus_pairs", "kept_of_21", "lost", "new", "hubs_degree_ge_5", "pooled_sig_raw/Ins1/Kcnj11/Hes1", "adjusted_q05_raw/Ins1/Kcnj11/Hes1", "r_Cacna1d_Glp1r_Hes1norm",
         "fig5_alias_KnownValidated/Novel/STRINGOnly/NotValidated", "in_STRING_alias_of_120", "fc_DAPT_Tcf7l2/Kcnj11/Cacna1d/Cacna1c/Ptbp1")
rows <- lapply(list(c("reported", "rows: 22.5685, 22.7983, 22.5685, 22.7983"), c("C0", sprintf("%.4f", REF$C0)), c("C1", sprintf("%.4f", REF$C1))), function(nr) {
  o <- out[[nr[1]]]; c5 <- o$fig5_alias
  c(nr[1], nr[2], as.character(o$n_consensus), as.character(o$kept), paste(unlist(o$lost), collapse = "; "), paste(unlist(o$new), collapse = "; "), paste(vapply(unlist(o$hubs), function(g) sprintf("%s %d", g, o$degree[[g]]), ""), collapse = "; "),
    paste(unlist(o$sig), collapse = "/"), paste(unlist(o$adj_q05), collapse = "/"), .json_num(round(o$r_cacna1d_glp1r, 3)), sprintf("%d/%d/%d/%d", c5$`Known Validated`, c5$`Novel Discovery`, c5$`STRING Only`, c5$`Not Validated`),
    as.character(o$in_STRING_alias), paste(sprintf("%.2f", unlist(o$fc_dapt[G5])), collapse = "/"))
})
write_text_bytes(paste0(paste(c(paste(vapply(hdr, .dw_quote, ""), collapse = ","), vapply(rows, function(r) paste(vapply(r, .dw_quote, ""), collapse = ","), "")), collapse = "\r\n"), "\r\n"), file.path(RES13, "gapdh_versions_summary.csv"))
for (n in c("reported", "C0", "C1")) {
  o <- out[[n]]; cat(sprintf("[%s] consensus %d (kept %d of 21; lost %s; new %s) | hubs %s | pooled sig %s | adj q<0.05 %s | |r|>=0.7 pairs %d | in STRING (alias) %d | r(Cacna1d,Glp1r) %.3f\n", n, o$n_consensus, o$kept,
              paste0("[", paste(sprintf("'%s'", unlist(o$lost)), collapse = ", "), "]"), paste0("[", paste(sprintf("'%s'", unlist(o$new)), collapse = ", "), "]"), paste0("{", paste(sprintf("'%s': %d", unlist(o$hubs), unlist(o$degree[unlist(o$hubs)])), collapse = ", "), "}"),
              paste0("[", paste(unlist(o$sig), collapse = ", "), "]"), paste0("[", paste(unlist(o$adj_q05), collapse = ", "), "]"), o$n_pairs_abs_r_ge_0_7, o$in_STRING_alias, o$r_cacna1d_glp1r))
  cat("     DAPT fold change vs CTRL:", paste(sprintf("%s: %.3f", G5, unlist(o$fc_dapt[G5])), collapse = ", "), "\n")
}
