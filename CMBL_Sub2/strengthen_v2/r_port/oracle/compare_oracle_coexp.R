## compare_oracle_coexp.R - compares the R port with the raw null distributions dumped by dump_oracle_coexp.py (counts must be equal; floats to 1e-12).
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1]))); port <- file.path(here, "..")
source(file.path(port, "rng", "numpy_rng.R")); source(file.path(port, "coexp_checks.R"))
suppressMessages(library(readxl))
G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
SRC <- Sys.getenv("INS1_MATRIX", file.path(here, "..", "..", "..", "R_17genes", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx"))
X <- as.matrix(as.data.frame(read_excel(SRC))[, G17]); storage.mode(X) <- "double"; grp <- rep(c("CTRL", "DAPT", "DKK-1"), each = 4); den <- c("Ins1", "Kcnj11", "Hes1")
rd <- function(f) readBin(file.path(here, f), "double", n = file.size(file.path(here, f)) / 8, size = 8)
cmp <- function(name, a, b, tol = 0) {
  d <- abs(a - b); nb <- sum(d > tol)
  cat(sprintf("%-34s n=%-6d differing entries: %d   max abs difference %.3g\n", name, length(a), nb, max(d)))
  if (nb > 0 && nb <= 8) { i <- which(d > tol); cat("     first differences at", head(i, 8), ": R", head(a[i], 8), " Python", head(b[i], 8), "\n") }
  nb == 0
}
t0 <- Sys.time()
an <- cc_adjusted_null(X, grp, den, nperm = 2000, seed = 2026, genes = G17); py <- matrix(rd("adjusted_null_2000.f64"), 2000, 4, byrow = TRUE)
ok1 <- cmp("adjusted null (2000 x 4)", as.vector(an$.null), as.vector(py))
cn <- cc_consensus_null(X, grp, den, nperm = 10000, seed = 2026, genes = G17)
ok2 <- cmp("consensus null, any sign", cn$.null, rd("consensus_null_10000.f64")); ok3 <- cmp("consensus null, same sign", cn$.null_same, rd("consensus_null_same_10000.f64"))
des <- cc_design(grp); R <- des$Q %*% log(X)
ch <- cc_perm_apply(R, des, 10000, 2026, NULL, fun = function(Rb) .coupling_stats(Rb, des))
ok4 <- cmp("coupling null T_p", unlist(lapply(ch, function(x) x$T_p)), rd("coupling_T_p_10000.f64"))
ok5 <- cmp("coupling null T_r2", unlist(lapply(ch, function(x) x$T_r2)), rd("coupling_T_r2_10000.f64"), 1e-12)
ok6 <- cmp("coupling null T_eig", unlist(lapply(ch, function(x) x$T_eig)), rd("coupling_T_eig_10000.f64"), 1e-9)
cat(sprintf("\nORACLE %s (%.0f s)\n", if (all(ok1, ok2, ok3, ok4, ok5, ok6)) "MATCH" else "DIFFERENCES", as.numeric(Sys.time() - t0, units = "secs")))
cat("adjusted null p_ge (R):", paste(sprintf("%s %.4f", names(an$schemes), vapply(an$schemes, function(s) s$p_ge, 1)), collapse = " | "), "\n")
