## test_degenerate_inputs.R - checks of the R port on degenerate inputs (25 Nov 2026, after the GSE20677 finding; labels and three checks revised 26 Nov 2026 after the independent audit, audit_261126/AUDIT_REPORT.md C5d).
## Each check states what exact arithmetic gives and compares the port with it.
## PASS = the port gives the exact answer. KNOWN = the port deliberately reproduces a floating-point artefact of the Python analysis, so the archived numbers are reproduced (documented in R_PORT_PLAN.md; waits
## for the user's decision). NOISE = a rounding-noise artefact that both implementations have but the port does not reproduce bit for bit (the archived numbers differ for such inputs). FAIL = an error nobody knew about.
## REF lines show what exact arithmetic gives; they are not tests of the port.
## Synthetic data only (12 samples, 3 groups of 4, as the study design). Light: runs in seconds. Output: results_R/test_degenerate_inputs.txt
## Run: Rscript test_degenerate_inputs.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(HERE, "rng", "numpy_rng.R")); source(file.path(HERE, "pycompat.R")); source(file.path(HERE, "coexp_checks.R")); source(file.path(HERE, "bench_pipeline.R"))
defs <- function(file, names) { ex <- parse(file); for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% names) eval(e, globalenv()) }
defs(file.path(HERE, "benchmark_route3", "run_benchmark.R"), c("ANOVA_EXACT", "seq_sum_rows", "anova_fdr", "py_json_ascii"))     # only these definitions (sourcing the runner would create folders)
ANOVA_EXACT <- FALSE   # the archived behaviour whatever the environment says; check 4 switches it on for one extra line
OUTF <- file.path(HERE, "results_R", "test_degenerate_inputs.txt"); LINES <- character(0); NS <- c(PASS = 0L, KNOWN = 0L, NOISE = 0L, FAIL = 0L)
say <- function(...) { s <- paste0(...); LINES <<- c(LINES, s); cat(s, "\n", sep = "") }
check <- function(name, exact_ok, known = NULL, noise = NULL, detail = "") {
  st <- if (isTRUE(exact_ok)) "PASS" else if (!is.null(known)) "KNOWN" else if (!is.null(noise)) "NOISE" else "FAIL"; NS[[st]] <<- NS[[st]] + 1L
  say(sprintf("%-5s %s%s%s", st, name, if (nzchar(detail)) paste0(" | ", detail) else "", if (st == "KNOWN") paste0(" | known: ", known) else if (st == "NOISE") paste0(" | noise: ", noise) else ""))
}
seqmean <- function(v) Reduce(`+`, v) / length(v)

## data: 12 x 8, genes 1-4 lognormal, genes 5 and 6 constant at values whose numpy-style (sequential) mean is inexact, gene 7 = 2 x gene 1, gene 8 = gene 4 x small noise
set.seed(2026); N <- 12; grp <- rep(1:3, each = 4); des <- cc_design(grp)
cands <- seq(0.1, 9.9, by = 0.0137); inexact <- cands[vapply(cands, function(c) seqmean(rep(c, N)) != c, TRUE)]
stopifnot(length(inexact) >= 2); c5 <- inexact[1]; c6 <- inexact[2]
X <- matrix(exp(rnorm(N * 4, 0, 0.5)), N, 4); X <- cbind(X, c5, c6, 2 * X[, 1], X[, 4] * exp(rnorm(N, 0, 0.05))); colnames(X) <- paste0("g", 1:8)   # g8 tracks g4 closely: a real consensus pair
say("test_degenerate_inputs.R | ", R.version.string, " | ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say(sprintf("constant genes g5 = %.17g and g6 = %.17g: their mean as numpy computes it (sequential sum / 12) differs from the value by %.1e and %.1e", c5, c6, seqmean(rep(c5, N)) - c5, seqmean(rep(c6, N)) - c6))
say("")

say("== 1. pooled correlation (cc_corr, coexp_checks.R)")
r <- cc_corr(t(X))
check("a constant gene correlates 0 with every gene", all(r[5:6, -(5:6)] == 0) && r[5, 6] == 0, detail = sprintf("largest |r| %.1e", max(abs(r[5:6, ]))))
check("a gene and 2 x that gene correlate 1 (to within 1e-12)", abs(r[1, 7] - 1) < 1e-12, detail = sprintf("r = %.17g", r[1, 7]))
r_log <- cc_corr(t(log(X)))
check("constant gene on the log scale also correlates 0", all(r_log[5:6, -(5:6)] == 0) && r_log[5, 6] == 0, detail = sprintf("largest |r| %.1e", max(abs(r_log[5:6, ]))))

say(""); say("== 2. residual correlation after removing group means (cc_residual_corr: the within-group coupling statistics)")
rr <- cc_residual_corr(X, grp); resid <- t(log(X)) %*% des$Q; nonc <- setdiff(1:8, 5:6)
check("a constant gene has residual correlation 0 with every gene", all(rr[5:6, -(5:6)] == 0) && rr[5, 6] == 0,
      noise = "the residuals of a constant gene are rounding noise (largest |residual| above); cc_corr scales the noise to unit length, so a constant gene can correlate with real genes and with another constant gene. This is the case of GSE20677 and GSE28477 (calibration, analyse_matrix, route 3 adjusted endpoints), where R and Python give different numbers because the noise depends on the order of floating-point operations; neither is the exact value",
      detail = sprintf("largest |residual| of g5, g6 %.1e; largest |r| of a constant gene with a non-constant gene %.3f; r(g5, g6) = %.3f", max(abs(resid[5:6, ])), max(abs(rr[5:6, nonc])), rr[5, 6]))
const_in_groups <- X; const_in_groups[, 5] <- rep(c(c5, 2 * c5, 3 * c5), each = 4)
rr2 <- cc_residual_corr(const_in_groups, grp)
check("a gene constant within each group (different group levels) has residual correlation 0", all(rr2[5, -5] == 0),
      noise = "same rounding-noise case: the within-group residuals are exactly 0 in exact arithmetic", detail = sprintf("largest |r| %.3f", max(abs(rr2[5, -5]))))

say(""); say("== 3. part 1 / route 3 pipeline (B_pooled, bench_pipeline.R)")
p <- B_pooled(X, 1:3); cs <- B_consensus_from(p$S, p$SG)
check("two constant genes are not significant together in the raw scheme", !p$S[[1]][5, 6], detail = "Python's np.corrcoef gives r = +-1 here when numpy's mean of both columns is inexact (the GSE20677 A artefact pair)")
ratio_r <- vapply(2:4, function(k) { s <- B_scheme_list(8, 1:3)[[k]]; M <- B_scheme_matrix(X, s$keep, s$d); cor(M[, match(5, s$keep)], M[, match(6, s$keep)]) }, 0)
check("in each ratio scheme two constant genes divided by the same denominator correlate 1 to within 1e-12 (exactly 1 in exact arithmetic: both are a constant times 1/denominator)", all(abs(ratio_r - 1) < 1e-12),
      detail = sprintf("r = %s; this is why the consensus rule also needs the raw scheme (first version of this file expected 0 here: the expectation was wrong, not the port)", paste(format(ratio_r, digits = 17), collapse = ", ")))
check("a constant gene is never in a consensus pair", !any(cs$any[5:6, ]) && !any(cs$any[, 5:6]))
perm <- c(1:3, 8, 6, 5, 7, 4); inv <- order(perm); sym <- function(m) unname(m | t(m))
cons_any <- function(M) { q <- B_pooled(M, 1:3); B_consensus_from(q$S, q$SG)$any }
c1 <- sym(cons_any(X)); c2 <- sym(cons_any(X[, perm])[inv, inv])
check("the consensus pairs do not depend on the column order of the non-denominator genes", sum(c1[upper.tri(c1)]) >= 1 && identical(c1, c2),
      detail = sprintf("whole 8 x 8 consensus matrix compared after undoing the permutation; %d consensus pair(s)", sum(c1[upper.tri(c1)])))
f <- exp(rnorm(N, 0, 0.4)); Xf <- X * f; pf <- B_pooled(Xf, 1:3)
dil <- vapply(2:4, function(k) { s <- B_scheme_list(8, 1:3)[[k]]; a <- suppressWarnings(cor(B_scheme_matrix(X, s$keep, s$d))); b <- suppressWarnings(cor(B_scheme_matrix(Xf, s$keep, s$d))); max(abs(a - b), na.rm = TRUE) }, 0)
sameS <- all(vapply(2:4, function(k) identical(p$S[[k]], pf$S[[k]]), TRUE))
check("the port's ratio schemes are unchanged when every sample is multiplied by its own factor (dilution): correlations and significant pairs", max(dil) < 1e-12 && sameS,
      detail = sprintf("largest correlation difference %.1e; significant pairs identical in the 3 ratio schemes: %s", max(dil), sameS))

say(""); say("== 4. route 3 ANOVA (anova_fdr, run_benchmark.R) on a probe that is constant in all 12 samples")
L <- log2(cbind(X[, 1:4], c5, c6)); af <- anova_fdr(L)
check("a constant probe is not called differentially expressed (F = 0/0 is undefined; its p should not be 0)", all(af$p[5:6] > 0.05),
      known = "the port reproduces numpy's sequential grand mean on purpose: the mean differs from the value by an ulp, the between-group sum of squares is about 1e-31 against a within-group sum of exactly 0, so F = Inf and p = 0; on GPL1261 thousands of floor-level probes get p = 0 and make the BH step lenient (memory: route 3 ANOVA artefact, decision pending)",
      detail = sprintf("p of the constant probes: %s", paste(format(af$p[5:6], digits = 3), collapse = ", ")))
gm_exact <- t(vapply(1:3, function(g) colMeans(L[grp == g, , drop = FALSE]), numeric(ncol(L)))); grand_exact <- colMeans(L)
ssb_exact <- 4 * colSums((gm_exact - matrix(grand_exact, 3, ncol(L), byrow = TRUE))^2)
say(sprintf("REF   with long-double means (R's colMeans, not the port's anova_fdr) the between-group sum of squares of the two constant probes is %s: F = 0/0 in exact arithmetic, so they would not be called",
            paste(format(ssb_exact[5:6], digits = 3), collapse = " and ")))
ANOVA_EXACT <- TRUE; af_x <- anova_fdr(L); ANOVA_EXACT <- FALSE
check("with ANOVA_EXACT=1 (impact check O6, 26 Nov 2026) the constant probes get p = 1 and every other p is unchanged", all(af_x$p[5:6] == 1) && identical(af_x$p[1:4], af$p[1:4]),
      detail = sprintf("p of the constant probes with the switch: %s", paste(format(af_x$p[5:6], digits = 3), collapse = ", ")))

say(""); say("== 5. ranks with ties (Spearman version of the INS-1 analysis)")
check("R's rank(), which the port uses for the Spearman version, gives scipy.stats.rankdata's average ranks for ties", identical(rank(c(3, 1, 2, 2, 5)), c(4, 1, 2.5, 2.5, 5)))

say(""); say("== 6. output format (route 3 results_per_dataset.jsonl)")
check("non-ASCII characters are written as Python's json.dumps writes them (\\u escape, lowercase hex)", identical(py_json_ascii(py_json_typed(list(t = "17β-estradiol"), indent = NULL)), "{\"t\": \"17\\u03b2-estradiol\"}"),
      detail = "the beta of GSE29881 (17beta-estradiol); the 25 Nov run wrote it raw, so that line differed from the archive in its bytes only")

say(""); say("== 7. group-effect ANOVA of the coexpression tool with a constant gene (cc_anova_group_effect, cc_geometry_check; fix F3 of the audit of 27 Nov 2026)")
Xg <- X[, 1:4]; Xg[grp == 2, 1] <- Xg[grp == 2, 1] * 4; Xg[grp == 3, 1] <- Xg[grp == 3, 1] * 8; Xg[grp == 2, 2] <- Xg[grp == 2, 2] / 4; Xg[grp == 3, 2] <- Xg[grp == 3, 2] * 3   # genes 1 and 2: strong group effects
p_const <- function(c) cc_anova_group_effect(cbind(Xg, c), des)$p[5]; n_ge <- function(M) cc_geometry_check(M, des)$n_group_effect_genes
cands7 <- 1 + (1:3000) / 997; pv7 <- vapply(cands7, p_const, 0)                     # archived behaviour: which constant values give p = 0 (inexact grand mean) or NaN (exact)
ci <- cands7[which(pv7 == 0)[1]]; ce <- cands7[which(is.nan(pv7))[1]]; stopifnot(!is.na(ci), !is.na(ce)); ref <- n_ge(Xg)
old <- c(n_ge(cbind(Xg, ci)), n_ge(cbind(Xg, ce))); CC_CONSTANT_FIX <- TRUE; new <- c(n_ge(cbind(Xg, ci)), n_ge(cbind(Xg, ce))); ag <- cc_anova_group_effect(cbind(Xg, ce), des); CC_CONSTANT_FIX <- FALSE
check("a gene with the same value in every sample has no group effect and leaves the count of the other genes unchanged", identical(old, c(ref, ref)),
      known = "the port reproduces the archived module (coexp_checks.py 1.0, tool 1.1.0) on purpose: an inexact grand mean gives the constant gene p = 0 (counted), an exact one gives p = NaN and every q NaN (no gene counted); fixed in tool 1.1.1",
      detail = sprintf("group-effect genes without the constant gene %d; with it (inexact mean %.6g, exact mean %.6g): %d, %d", ref, ci, ce, old[1], old[2]))
check("with CC_CONSTANT_FIX=1 (coexp-checks 1.1.1) the constant gene gets p = 1, no q is NaN and the count of the other genes is unchanged", identical(new, c(ref, ref)) && ag$p[5] == 1 && !anyNA(ag$q),
      detail = sprintf("group-effect genes with the switch (inexact mean, exact mean): %d, %d", new[1], new[2]))

say(""); say(sprintf("SUMMARY: %d PASS, %d KNOWN (Python artefacts reproduced on purpose), %d NOISE (rounding noise of both implementations, not reproduced bit for bit), %d FAIL", NS[["PASS"]], NS[["KNOWN"]], NS[["NOISE"]], NS[["FAIL"]]))
writeLines(LINES, OUTF, useBytes = TRUE)
