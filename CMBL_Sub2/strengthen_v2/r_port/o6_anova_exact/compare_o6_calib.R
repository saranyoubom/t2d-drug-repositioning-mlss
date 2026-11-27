## compare_o6_calib.R - 27 Nov 2026, impact check O6, last step: the WP2 numbers under the exact rule (results_R/o6_anova_exact/wp2_exact) next to the R run with numpy's order
## (calibration: results_R/wp2; development application and validation summary: results_R/o6_anova_exact/wp2_numpy_R, made by run_o6_calib.sh step 8 from the R run's inputs).
## Writes results_R/o6_anova_exact/o6_calibration_summary.md. Read-only on everything else.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
O <- file.path(R_PORT, "results_R", "o6_anova_exact"); X <- file.path(O, "wp2_exact"); N <- file.path(O, "wp2_numpy_R"); RUN <- file.path(R_PORT, "results_R", "wp2")
rd <- function(p) { if (!file.exists(p)) stop("missing ", p); x <- readLines(p, warn = FALSE, encoding = "UTF-8"); x[!grepl("^Warning message|built under R version", x)] }
md <- c("# Impact check O6, part 2: WP2 numbers with the exact ANOVA rule", "", sprintf("Written %s by o6_anova_exact/compare_o6_calib.R. Exact rule: development set of 36 series (GSE20677 and GSE49079 out, GSE30164 in; set A redrawn in GSE15892, GSE28477, GSE53183); held-out GSE68424 set A redrawn; everything else reused bit for bit (prep_exact_calib.log, merge_heldout.log).", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
## 1 calibration verdict
md <- c(md, "## 1. Calibration verdict (protocol 3.2)", "", "### R run, numpy order (37 series; equals the archived Python verdict except T_p and T_eig, BLAS noise)", "", "```", rd(file.path(RUN, "calibration_verdict.txt")), "```", "",
        "### Exact rule (36 series)", "", "```", rd(file.path(X, "calibration_verdict.txt")), "```", "")
## 2 per-series rates of the series that changed, and every series with a type I error of T_r2 above 0.06
pc <- function(d) { x <- read.csv(file.path(d, "calibration_per_condition.csv"), stringsAsFactors = FALSE); x$rate <- x$n_reject / x$n_rep; x }
a <- pc(RUN); b <- pc(X)
sel <- function(x, m, r) x[x$statistic == "T_r2" & x$m == m & abs(x$r - r) < 1e-9, c("gse", "platform", "n_rep", "n_reject")]
t1 <- merge(sel(a, 0, 0), sel(b, 0, 0), by = c("gse", "platform"), all = TRUE, suffixes = c("_numpy", "_exact")); pw <- merge(sel(a, 6, 0.7), sel(b, 6, 0.7), by = c("gse", "platform"), all = TRUE, suffixes = c("_numpy", "_exact"))
f <- function(k, n) ifelse(is.na(n), "-", sprintf("%d/%d", k, n))
show <- unique(c("GSE20677", "GSE49079", "GSE15892", "GSE28477", "GSE53183", "GSE30164", t1$gse[(!is.na(t1$n_reject_numpy) & t1$n_reject_numpy / t1$n_rep_numpy > 0.06) | (!is.na(t1$n_reject_exact) & t1$n_reject_exact / t1$n_rep_exact > 0.06)]))
md <- c(md, "## 2. T_r2 per series: the changed series and every series with a type I error above 0.06", "", "| series | platform | type I error, numpy | type I error, exact | power r = 0.7, m = 6, numpy | exact |", "|---|---|---|---|---|---|")
for (g in show) { u <- t1[t1$gse == g, ]; v <- pw[pw$gse == g, ]
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |", g, u$platform, f(u$n_reject_numpy, u$n_rep_numpy), f(u$n_reject_exact, u$n_rep_exact), f(v$n_reject_numpy, v$n_rep_numpy), f(v$n_reject_exact, v$n_rep_exact))) }
md <- c(md, "")
## 3 calibration tables: every row whose rate moved
ta <- read.csv(file.path(RUN, "calibration_tables.csv"), stringsAsFactors = FALSE); tb <- read.csv(file.path(X, "calibration_tables.csv"), stringsAsFactors = FALSE)
tt <- merge(ta[, c("r", "m", "statistic", "n_tests", "n_reject", "rate")], tb[, c("r", "m", "statistic", "n_tests", "n_reject", "rate")], by = c("r", "m", "statistic"), suffixes = c("_numpy", "_exact")); tt <- tt[order(tt$m, tt$r, tt$statistic), ]
md <- c(md, "## 3. Pooled calibration table (all conditions)", "", "| r | m | statistic | numpy | exact |", "|---|---|---|---|---|",
        sprintf("| %s | %d | %s | %d/%d = %.3f | %d/%d = %.3f |", format(tt$r), tt$m, tt$statistic, tt$n_reject_numpy, tt$n_tests_numpy, tt$rate_numpy, tt$n_reject_exact, tt$n_tests_exact, tt$rate_exact), "")
## 4 and 5: Markdown summaries, line by line
difflines <- function(title, p1, p2) {
  A <- rd(p1); B <- rd(p2); n <- max(length(A), length(B)); length(A) <- n; length(B) <- n; A[is.na(A)] <- "(none)"; B[is.na(B)] <- "(none)"; d <- which(A != B)
  c(paste("##", title), "", sprintf("%d of %d lines differ (left: numpy order, R run inputs; right: exact rule).", length(d), n), "", if (length(d)) unlist(lapply(d, function(i) c(sprintf("- line %d", i), sprintf("  - numpy: %s", A[i]), sprintf("  - exact: %s", B[i])))) else "- none", "")
}
md <- c(md, difflines("4. Development application (dev_application_summary.md)", file.path(N, "dev_application_summary.md"), file.path(X, "dev_application_summary.md")))
md <- c(md, difflines("5. Validation summary (validation_summary.md)", file.path(N, "validation_summary.md"), file.path(X, "validation_summary.md")))
md <- c(md, difflines("6. WP2 summary table (WP2_summary_table.md)", file.path(N, "WP2_summary_table.md"), file.path(X, "WP2_summary_table.md")))
md <- c(md, "Note: summarize_validation.R labels the development arm \"Development arrays (37)\" whatever the number of series; under the exact rule it holds 36 series.", "")
writeLines(md, file.path(O, "o6_calibration_summary.md"), useBytes = TRUE); cat(md, sep = "\n")
