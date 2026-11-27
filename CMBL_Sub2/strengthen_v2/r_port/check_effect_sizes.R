## check_effect_sizes.R - WP5a: independent check of effect_sizes.csv (R port of wp5_reporting/check_effect_sizes.py).
## Recomputes, with other code paths, (1) the Welch intervals of the 34 gene contrasts with t.test(), (2) the pooled-variance intervals with lm() over the three groups,
## (3) the Benjamini-Hochberg q with a hand-written step-up routine. Prints the maximum absolute differences (expected below 1e-5, the rounding of the CSV to 6 significant digits).
## Usage: Rscript check_effect_sizes.R   (CMBL_ROOT overrides the project root = three folders above this script; ES_OUT = folder of effect_sizes.csv, default results_R/wp5)
suppressMessages(library(readxl))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
ROOT <- Sys.getenv("CMBL_ROOT", normalizePath(file.path(HERE, "..", "..", "..")))
AF4 <- file.path(ROOT, "Submission_package_JTM", "Additional file 4.xlsx"); OUT <- Sys.getenv("ES_OUT", file.path(HERE, "results_R", "wp5"))
sh <- as.matrix(read_excel(AF4, sheet = "Fig_4c_Raw_expression_data", col_names = FALSE, col_types = "text", .name_repair = "minimal"))
rows <- lapply(seq_len(nrow(sh)), function(i) { v <- sh[i, ]; v[!is.na(v)] })                       # empty cells dropped, as in the Python check
h <- which(vapply(rows, function(r) length(r) > 0 && r[1] == "Wnt5a", TRUE))[1]
genes <- rows[[h]]; X <- do.call(rbind, lapply(rows[(h + 1):(h + 12)], as.numeric)); colnames(X) <- genes
main <- read.csv(file.path(OUT, "effect_sizes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
main <- main[main$readout == "RT-qPCR" & main$variant == "as reported", ]; stopifnot(nrow(main) == 34)
design <- data.frame(dapt = c(rep(0, 4), rep(1, 4), rep(0, 4)), dkk = c(rep(0, 8), rep(1, 4)))
d_welch <- d_pool <- 0
for (i in seq_len(nrow(main))) {
  r <- main[i, ]; x <- X[, r$feature]; y <- log2(x / mean(x[1:4])); is_dapt <- startsWith(r$contrast, "DAPT"); k <- if (is_dapt) 5:8 else 9:12
  ci <- t.test(y[k], y[1:4])$conf.int                                                              # Welch (R default), difference treated minus CTRL
  d_welch <- max(d_welch, abs(ci[1] - r$ci95_lo), abs(ci[2] - r$ci95_hi), abs(mean(y[k]) - mean(y[1:4]) - r$log2fc_geo))
  cf <- confint(lm(y ~ dapt + dkk, data = cbind(y = y, design)), level = 0.95)[if (is_dapt) "dapt" else "dkk", ]    # pooled variance over the three groups, df = 9
  d_pool <- max(d_pool, abs(cf[1] - r$ci95_lo_pooled3), abs(cf[2] - r$ci95_hi_pooled3))
}
bh_manual <- function(p) { m <- length(p); o <- order(p); q <- numeric(m); prev <- 1; for (rank in m:1) { i <- o[rank]; prev <- min(prev, p[i] * m / rank); q[i] <- prev }; q }
d_bh <- 0
for (con in c("DAPT vs CTRL", "DKK-1 vs CTRL")) { s <- main[main$contrast == con, ]; d_bh <- max(d_bh, max(abs(bh_manual(s$p_raw_recomputed) - s$bh_q_17genes))) }
cat(sprintf("max abs difference, Welch interval and effect vs t.test():  %.2e\n", d_welch))
cat(sprintf("max abs difference, pooled interval vs lm() / confint():    %.2e\n", d_pool))
cat(sprintf("max abs difference, BH q vs hand-written step-up:           %.2e\n", d_bh))
