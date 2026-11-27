args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1]))); port <- file.path(here, "..")
source(file.path(port, "rng", "numpy_rng.R")); source(file.path(port, "coexp_checks.R"))
v <- readBin(file.path(here, "runstruct_oracle.f64"), "double", n = 60 * 148, size = 8); runs <- c(rep("A", 5), rep("B", 4), rep("C", 3)); g <- rep(1:3, each = 4)
bad <- 0; pr <- numeric(60); pp <- numeric(60)
for (i in 0:59) {
  blk <- v[i * 148 + 1:148]; X <- matrix(blk[1:144], 12, 12, byrow = TRUE); py <- blk[145:148]
  res <- cc_run_structure_test(X, g, runs, nperm = 999, seed = i); r <- c(res$r$p_two_sided, res$r$D, res$r$null_sd, res$abs_r$p_two_sided); pr[i + 1] <- r[1]; pp[i + 1] <- py[1]
  if (any(abs(r - py) > 1e-9)) { bad <- bad + 1; cat(sprintf("data set %2d: R p %.6f D %.8f sd %.8f p_abs %.6f | Python p %.6f D %.8f sd %.8f p_abs %.6f\n", i, r[1], r[2], r[3], r[4], py[1], py[2], py[3], py[4])) }
}
cat(sprintf("%d of 60 data sets differ; mean p(r): R %.4f, Python %.4f\n", bad, mean(pr), mean(pp)))
