args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1]))); port <- file.path(here, "..")
source(file.path(port, "rng", "numpy_rng.R")); source(file.path(port, "coexp_checks.R")); suppressMessages(library(readxl))
G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
X <- as.matrix(as.data.frame(read_excel(file.path(here, "..", "..", "..", "R_17genes", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx")))[, G17]); storage.mode(X) <- "double"
grp <- rep(c("CTRL", "DAPT", "DKK-1"), each = 4); den <- c("Ins1", "Kcnj11", "Hes1"); des <- cc_design(grp)
py <- readBin(file.path(here, "spearman_null.f64"), "double", n = 20000, size = 8)
cn <- cc_consensus_null(X, grp, den, nperm = 10000, seed = 2026, kind = "spearman", genes = G17)
bad <- which(cn$.null != py[1:10000]); cat("differing permutations (any sign):", length(bad), "\n"); print(head(bad, 20))
cat("R:", cn$.null[head(bad, 10)], "\nPython:", py[head(bad, 10)], "\n")
if (length(bad)) {
  b <- bad[1]; stk <- cc_permuted_all(X, grp, 10000, 2026)           # (10000, 17, 12)
  Xb <- matrix(stk[b, , ], 17, 12); denx <- match(den, G17); elig <- setdiff(1:17, denx); rc <- des$rc_pooled
  cat(sprintf("\npermutation %d (1-based), rc = %.17g\n", b, rc))
  for (d in c(list(NULL), as.list(denx))) {
    M <- if (is.null(d)) Xb[elig, ] else Xb[elig, ] / matrix(Xb[d, ], length(elig), 12, byrow = TRUE); r <- cc_corr(.rank_rows(M)); a <- abs(r)[upper.tri(r)]
    cat(sprintf("scheme %-8s: pairs with |r| within 1e-9 of rc: %d; closest |r| - rc: %.3g; n ties in rows: %d\n", if (is.null(d)) "raw" else G17[d], sum(abs(a - rc) < 1e-9), a[which.min(abs(a - rc))] - rc, sum(apply(M, 1, function(v) anyDuplicated(v) > 0))))
  }
}
