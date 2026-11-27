## gapdh_exact.R - R port of reanalysis_17genes/gapdh_exact.py.
## Exact re-computation of the 26 August DAPT values from the author's original Quantification Summary sheet (Undergrad_Glucose uptake signaling/RT-PCR).
## Original analysis: DAPT replicate k of every target gene was paired with the Gapdh value HK_k = [22.5685, 22.7983, 22.5685, 22.7983] (raw wells C10 and D10, each used twice;
## raw A10 24.9874 and B10 24.1318 were parked in column V and not used). Relative expression = 2^-(Cq_target - Gapdh).
##  C0: reference = mean of the two valid wells (22.6834) for all four DAPT replicates (exclusion done without duplicating values)
##  C1: reference = mean of all four raw wells (23.6225)   (worst case: failed wells included)
## Writes gapdh_exact.json (indent 1) and gapdh_exact_pairs.csv. Run from any folder: Rscript gapdh_exact.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "pipeline.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
G5 <- c("Tcf7l2", "Kcnj11", "Cacna1d", "Cacna1c", "Ptbp1")
HK <- c(22.5685093010354, 22.7982726003435, 22.5685093010354, 22.7982726003435)      # as used (sheet rows 30-33)
RAW <- list(A10 = 24.9873996418614, B10 = 24.1318270358876, C10 = 22.5685093010354, D10 = 22.7982726003435)
REF <- list(C0 = mean(c(RAW$C10, RAW$D10)), C1 = mean(unlist(RAW)))
df0 <- P_load(P_GENES17)
variant <- function(ref) { d <- df0; for (g in G5) d[5:8, g] <- df0[5:8, g] * 2.0^(ref - HK); d }      # 2^-(Cq - ref) = 2^-(Cq - HK) * 2^(ref - HK)
anova_bonf <- function(x) {
  g <- list(x[1:4], x[5:8], x[9:12]); ss_w <- sum(vapply(g, function(a) sum((a - mean(a))^2), 1)); mse <- ss_w / 9
  gm <- mean(x); ss_b <- sum(vapply(g, function(a) length(a) * (mean(a) - gm)^2, 1)); p <- pf((ss_b / 2) / mse, 2, 9, lower.tail = FALSE); out <- list()
  for (pr in list(list(c(1, 2), "CTRL-DAPT"), list(c(1, 3), "CTRL-DKK1"), list(c(2, 3), "DAPT-DKK1"))) {
    t <- (mean(g[[pr[[1]][1]]]) - mean(g[[pr[[1]][2]]])) / sqrt(mse * 2 / 4); out[[pr[[2]]]] <- min(1.0, 2 * pt(abs(t), 9, lower.tail = FALSE) * 3)
  }
  c(list(omnibus_p = p, fc_dapt = mean(g[[2]]) / mean(g[[1]]), fc_dkk = mean(g[[3]]) / mean(g[[1]])), out)
}
summary_ <- function(d) {
  sch <- P_schemes(d); res <- P_sig_sets(sch); cons <- P_consensus(res, use_sign = FALSE); adj <- P_adjusted(sch); cn <- P_counts(res); deg <- list()
  for (i in seq_len(nrow(cons))) for (a in cons[i, ]) deg[[a]] <- (if (is.null(deg[[a]])) 0L else deg[[a]]) + 1L
  dn <- names(deg); dv <- unlist(deg); o <- order(-dv, dn, method = "radix"); deg <- as.list(dv[o])             # sorted by (-degree, name), C-locale names
  pget <- function(s, field) { i <- match("Cacna1d", s$genes); j <- match("Glp1r", s$genes); s[[field]][i, j] }
  list(r_cacna1d_glp1r = lapply(res, function(s) pget(s, "r")), p_cacna1d_glp1r = lapply(res, function(s) pget(s, "p")), sig = lapply(cn, function(v) v$sig), consensus = .pairs_list(cons),
       adj_q05 = lapply(adj, function(a) as.integer(sum(a$q[upper.tri(a$q)] < 0.05))), degree = deg)
}
out <- list()
for (nm in c("reported", "C0", "C1")) {
  d <- if (nm == "reported") df0 else variant(REF[[nm]])
  out[[nm]] <- c(list(per_gene = setNames(lapply(G5, function(g) anova_bonf(d[, g])), G5)), summary_(d))
}
keyset <- function(l) vapply(l, function(p) paste(unlist(p), collapse = "-"), "")
c0 <- keyset(out$reported$consensus)
cat("reference Cq: reported pairs HK", paste(sprintf("%.3f", HK), collapse = " "), sprintf("| C0 %.4f | C1 %.4f\n", REF$C0, REF$C1))
for (g in G5) cat(sprintf("%-8s DAPT fold change / padj(CTRL-DAPT):  reported %.3f / %.2g | C0 %.3f / %.2g | C1 %.3f / %.2g\n", g,
                          out$reported$per_gene[[g]]$fc_dapt, out$reported$per_gene[[g]][["CTRL-DAPT"]], out$C0$per_gene[[g]]$fc_dapt, out$C0$per_gene[[g]][["CTRL-DAPT"]], out$C1$per_gene[[g]]$fc_dapt, out$C1$per_gene[[g]][["CTRL-DAPT"]]))
for (n in c("reported", "C0", "C1")) {
  S <- out[[n]]; cs <- keyset(S$consensus); deg <- S$degree
  cat(sprintf("[%s] sig %s | consensus %d (kept %d, lost %d, new %d) | adj q<0.05 %s | hubs %s | Glp1r %d Cacna1d %d Ptbp1 %d Tcf7l2 %d\n", n, .py_counts(S$sig), length(cs), length(intersect(cs, c0)), length(setdiff(c0, cs)), length(setdiff(cs, c0)), .py_counts(S$adj_q05),
              paste0("[", paste(sprintf("'%s'", names(deg)[unlist(deg) >= 5]), collapse = ", "), "]"), if (is.null(deg$Glp1r)) 0L else deg$Glp1r, if (is.null(deg$Cacna1d)) 0L else deg$Cacna1d, if (is.null(deg$Ptbp1)) 0L else deg$Ptbp1, if (is.null(deg$Tcf7l2)) 0L else deg$Tcf7l2))
}
sets <- lapply(c("C0", "C1"), function(n) keyset(out[[n]]$consensus)); allp <- py_sort(unique(c(c0, unlist(sets))))        # sorted tuples = sorted 'a-b' strings (gene names contain no '-' or characters below it)
tab <- data.frame(pair = allp, involves_G5 = vapply(strsplit(allp, "-", fixed = TRUE), function(p) any(p %in% G5), TRUE), reported = allp %in% c0, C0 = allp %in% sets[[1]], C1 = allp %in% sets[[2]], stringsAsFactors = FALSE)
py_write_csv(tab, file.path(OUT, "gapdh_exact_pairs.csv"))
write_json_typed(c(list(ref = list(C0 = REF$C0, C1 = REF$C1), hk_used = HK, raw = RAW), out), file.path(OUT, "gapdh_exact.json"))
