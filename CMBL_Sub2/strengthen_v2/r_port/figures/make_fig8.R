## make_fig8.R - R port of wp6a_assets/make_fig8_final.py: FINAL Fig. 8 (public islet meta-analysis and perturbation concordance) for the JTM manuscript. Plots saved results of wp4_islets/meta and
## wp4_islets/wider_search only; no new analysis.
## (a) forest plot of the excess co-expression E per human dataset (square, +/-1.96 SE) and the pooled estimates (diamond, REML + Hartung-Knapp 95% CI) with the +/-0.10 margin shaded (no label: the legend
##     names the band). (b) pair-level pooled depth-adjusted rho of the 11 consensus pairs evaluable in at least 3 human datasets (dots = datasets, diamond = pooled estimate, gray line = 95% CI; number in
##     brackets = datasets). (c) log2 fold changes of the 17 genes in the three analyzable perturbation series (GSE60803, GSE167073, GSE99503; the last is a mirror design, fold changes multiplied by -1, said in
##     the legend) against the INS-1 DKK-1 arm, one small panel per series on shared axes; filled dots = same sign in both, open dots = opposite sign; shaded quadrants = same sign. The in-image texts that are
##     not an axis label, a panel letter, a row label, a number or a dataset title are removed (10 Nov 2026, finding F41 accepted by the user); the legend carries them.
## No key strip; 170 mm wide, 300 dpi PNG and PDF. Usage: Rscript make_fig8.R   Environment: WP4_META (default the archived wp4_islets/meta), WP4_WS (wp4_islets/wider_search), FIG_OUT (default results_R/figures).
## Output: Figure_8_islet_meta_final.png / .pdf
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R")); suppressMessages(library(jsonlite))
W4 <- normalizePath(file.path(R_PORT, "..", "wp4_islets")); MET <- Sys.getenv("WP4_META", file.path(W4, "meta")); WS <- Sys.getenv("WP4_WS", file.path(W4, "wider_search")); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures"))
ms <- fromJSON(file.path(MET, "meta_summary.json"), simplifyVector = FALSE); V <- list(); for (v in ms$variants) V[[paste(v$variant, v$mode)]] <- v
T <- read.csv(file.path(MET, "study_effects__primary__adjusted.csv"), stringsAsFactors = FALSE); PP <- read.csv(file.path(MET, "pair_level_pooled_rho.csv"), stringsAsFactors = FALSE)
CR <- read.csv(file.path(WS, "concordance_results.csv"), stringsAsFactors = FALSE); rownames(CR) <- CR$accession
pylist <- function(s) { m <- regmatches(s, gregexpr("'[^']*'", s))[[1]]; substring(m, 2, nchar(m) - 1) }                         # ast.literal_eval of a list of strings
H <- 134
draw <- function() {
  set_ls(1.1); tag("(a)"); tag("(b)", y = 56.5); tag("(c)", x = 88.0, y = 56.5)
  ## ------------------------------------------------------------------ (a) forest plot
  ax <- axes_new(57.0, 6.5, 70.0, 37.0, c(-0.22, 0.42), c(-0.5, 9.5))
  rect_mm(ax$X(-0.10) - 0.5 * PT, ax$top, ax$X(0.10) + 0.5 * PT, ax$bottom, "#ececea", NA); seg(ax$X(0), ax$top, ax$X(0), ax$bottom, INK, 0.5)           # axvspan: its 1 pt edge in the same color widens the band by half a point on each side
  ylab <- character(0); yy <- numeric(0); tx <- character(0)
  for (i in 1:5) { y <- c(8, 7, 6, 5, 4)[i]; r <- T[i, ]; lo <- r$E - 1.96 * r$SE; hi <- r$E + 1.96 * r$SE
    seg(ax$X(lo), ax$Y(y), ax$X(hi), ax$Y(y), BLUE, 1.1); marker(ax$X(r$E), ax$Y(y), "s", sqrt(19), fill = BLUE)
    ylab <- c(ylab, sprintf("%s (%d β-cells, %d pairs)", sub(" [(].*$", "", r$label), r$beta_cells, r$pairs)); yy <- c(yy, y); tx <- c(tx, sprintf("%s [%s, %s]", sg(r$E), sg(lo), sg(hi))) }
  pooled <- function(y, rec, name, col) {
    r <- Filter(function(e) startsWith(e$estimator, "REML+HK (primary)"), rec$estimators)[[1]]
    poly(c(ax$X(r$lo), ax$X(r$mu), ax$X(r$hi), ax$X(r$mu)), c(ax$Y(y), ax$Y(y + 0.30), ax$Y(y), ax$Y(y - 0.30)), col)
    ylab <<- c(ylab, name); yy <<- c(yy, y); tx <<- c(tx, sprintf("%s [%s, %s]", sg(r$mu), sg(r$lo), sg(r$hi)))
  }
  pooled(2.4, V[["primary adjusted"]], "Pooled, depth-adjusted (primary)", INK); pooled(1.4, V[["nondiabetic adjusted"]], "Pooled, non-diabetic donors only", GREYLINE); pooled(0.4, V[["primary unadjusted"]], "Pooled, not depth-adjusted", GREYLINE)
  seg(ax$left, ax$Y(3.2), ax$right, ax$Y(3.2), "#bdbdb8", 0.5)
  for (i in seq_along(yy)) mtext(ax$left + 1.03 * ax$w, ax$Y(yy[i]), tx[i], 6.5, 1, INK, "left", "center")
  mtext(ax$left + 1.03 * ax$w, ax$Y(9.0), list(list(list("E", 4), list(" [95% CI]", 2))), 6.5, 1, INK, "left", "center", plus = list(MATH))
  ax_spines(ax); tw <- ax_yticks(ax, yy, ylab, 6.5); th <- ax_xticks(ax, c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4), c(paste0(MINUS, "0.2"), paste0(MINUS, "0.1"), "0.0", "0.1", "0.2", "0.3", "0.4"), 6.5)
  ax_xlabel_lines(ax, list("Excess co-expression of the consensus pairs over matched", list(list("random pairs, "), list("E", 3), list(" (Fisher "), list("z", 3), list(")"))), th, 7, 1.15, plus = list(NULL, MATH))
  ## ------------------------------------------------------------------ (b) pair-level pooled rho
  ax <- axes_new(38.0, 60.5, 44.0, 60.0, c(-0.36, 0.52), c(-0.7, nrow(PP) - 0.3)); P <- PP[order(PP$pooled_rho), ]; yp <- seq_len(nrow(P)) - 1
  for (y in yp) seg(ax$left, ax$Y(y), ax$right, ax$Y(y), "#efefec", 0.5); seg(ax$X(0), ax$top, ax$X(0), ax$bottom, INK, 0.5)
  for (i in seq_len(nrow(P))) { y <- yp[i]; r <- P[i, ]; seg(ax$X(r$lo), ax$Y(y), ax$X(r$hi), ax$Y(y), "#8c8b86", 1.1)
    for (v in as.numeric(strsplit(r$rho_by_dataset, ";", fixed = TRUE)[[1]])) marker(ax$X(v), ax$Y(y), "o", sqrt(11), fill = BLUE, border = "white", bw_pt = 0.3)
    poly(c(ax$X(r$pooled_rho - 0.014), ax$X(r$pooled_rho), ax$X(r$pooled_rho + 0.014), ax$X(r$pooled_rho)), c(ax$Y(y), ax$Y(y + 0.34), ax$Y(y), ax$Y(y - 0.34)), INK) }
  ax_spines(ax); tw <- ax_yticks(ax, yp, as.list(sprintf("%s–%s (%d)", P$gene1, P$gene2, P$datasets)), 6.5, 3)                                   # plain italic text (digits of Tcf7l2, Cacna1c ... stay italic)
  th <- ax_xticks(ax, c(-0.2, 0, 0.2, 0.4), c(paste0(MINUS, "0.2"), "0.0", "0.2", "0.4"), 6.5); ax_xlabel(ax, NULL, th, runs = list(list("Pooled Spearman "), list("ρ", 3), list(" (95% CI)")))
  ## ------------------------------------------------------------------ (c) the three analyzable perturbation series against the INS-1 DKK-1 arm (small panels on shared axes)
  SER <- list(list("GSE60803", list(c("β-catenin shRNA", 1))), list("GSE167073", list(c("Fzd9", 3), c(" knockout", 1))), list("GSE99503", list(c("GSK3 inhibitor", 1))))
  XL <- c(-1.03, 0.22); YL <- c(-0.36, 0.58); PX0 <- 99.0; PW <- 22.0; PGAP <- 2.0; PTOP <- 71.0; PH <- 49.5
  for (k in 1:3) {
    acc <- SER[[k]][[1]]; ax <- axes_new(PX0 + (k - 1) * (PW + PGAP), PTOP, PW, PH, XL, YL); d <- read.csv(file.path(WS, sprintf("fold_changes_%s.csv", acc)), stringsAsFactors = FALSE, check.names = FALSE); r <- CR[acc, ]
    same <- d$gene %in% pylist(r$same_sign_genes)
    rect_mm(ax$X(XL[1]), ax$Y(0), ax$X(0), ax$Y(YL[1]), "#efeeea", NA); rect_mm(ax$X(0), ax$Y(YL[2]), ax$X(XL[2]), ax$Y(0), "#efeeea", NA)                  # lower left quadrant: both down; upper right: both up
    seg(ax$left, ax$Y(0), ax$right, ax$Y(0), INK, 0.5); seg(ax$X(0), ax$top, ax$X(0), ax$bottom, INK, 0.5)
    x <- d[["ins1_log2fc_DKK-1"]]; y <- d[["log2fc_signed"]]
    marker(ax$X(x[same]), ax$Y(y[same]), "o", sqrt(12), fill = BLUE, border = BLUE, bw_pt = 0.4); marker(ax$X(x[!same]), ax$Y(y[!same]), "o", sqrt(12), fill = "white", border = BLUE, bw_pt = 0.7)
    cxm <- ax$left + PW / 2; mtext(cxm, PTOP - 4.6, acc, 6.8, 2, INK, "center", "bottom"); mtext(cxm, PTOP - 1.6, list(lapply(SER[[k]][[2]], function(p) list(p[1], as.integer(p[2])))), 6.0, 1, INK2, "center", "bottom")
    ax_spines(ax); ax_xticks(ax, c(-1.0, -0.5, 0.0), c(paste0(MINUS, "1"), paste0(MINUS, "0.5"), "0"), 6.5); ax_yticks(ax, c(-0.2, 0.0, 0.2, 0.4), if (k == 1) c(paste0(MINUS, "0.2"), "0.0", "0.2", "0.4") else NULL, 6.5)
  }
  SUB2 <- list("2", 1, 0.86, 1.5); LOGSUB <- c(2, 4.7)                                       # mathtext subscript (86% size) and the extra ascent / descent (px) that a line with a subscript has in matplotlib
  mtext(PX0 + (3 * PW + 2 * PGAP) / 2, 130.2, list(list(list("INS-1, DKK-1 vs control (log"), SUB2, list(" fold change)"))), 7, 1, INK, "center", "bottom", plus = list(LOGSUB))
  mtext(89.2, PTOP + PH / 2, list("Perturbation vs control, sign-adjusted", list(list("(log"), SUB2, list(" fold change)"))), 7, 1, INK, "center", "center", rot = 90, ls = 1.15, plus = list(NULL, LOGSUB), mode = "default")
}
fig_save("Figure_8_islet_meta_final", H, draw, OUT)
