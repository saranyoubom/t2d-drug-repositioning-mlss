## make_fig7.R - R port of wp6a_assets/make_fig7_final.py: FINAL Fig. 7 (public islet cross-check, version 2) for the JTM manuscript. Plots saved results of wp4_islets/results only; no new analysis.
## (a) detection of the 17 genes in the author-labeled beta-cells of six datasets (numbers = % of beta-cells; blue shade grows with the percentage; gray = below the 10% rule, not evaluable; white = gene not in the
##     data file). (b) the 21 consensus pairs in the same six datasets: depth-adjusted partial Spearman rho (dot) and expression-matched random pairs (gray bar = null mean to 95th percentile); NE = not evaluable.
## No key strip; 170 mm wide, 300 dpi PNG and PDF. Usage: Rscript make_fig7.R   Environment: WP4_RES (folder with <tag>__primary__<mode>__<part>.csv and <tag>__primary__summary.json; default the archived
## wp4_islets/results), FIG_OUT (default results_R/figures). Output: Figure_7_islet_crosscheck_final.png / .pdf
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R")); suppressMessages(library(jsonlite))
RES <- Sys.getenv("WP4_RES", normalizePath(file.path(R_PORT, "..", "wp4_islets", "results"))); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures"))
GENES <- c("Hes1", "Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
DS <- list(c("gse81608", "GSE81608", "human"), c("gse84133_human", "GSE84133", "human"), c("gse86469", "GSE86469", "human"), c("gse81547", "GSE81547", "human"), c("emtab5061", "E-MTAB-5061", "human"), c("gse84133_mouse", "GSE84133", "mouse"))
rd <- function(tag, mode, part, variant = "primary") read.csv(file.path(RES, sprintf("%s__%s__%s__%s.csv", tag, variant, mode, part)), stringsAsFactors = FALSE, na.strings = c("", "NA"))
summ <- lapply(DS, function(d) fromJSON(file.path(RES, sprintf("%s__primary__summary.json", d[1])), simplifyVector = FALSE)$adjusted); names(summ) <- vapply(DS, `[`, "", 1)
cons <- rd("gse81608", "unadjusted", "pairs"); cons <- cons[cons$consensus_invitro %in% c("True", "TRUE", TRUE), c("gene1", "gene2")]; stopifnot(nrow(cons) == 21)
H <- 158; LEFT <- 24.0; PITCH <- 24.4; PW <- 21.9; ROW_A <- 3.05; ROW_B <- 3.35
draw <- function() {
  set_ls(1.1); tag("(a)"); tag("(b)", y = 66.5)
  ## ------------------------------------------------------------------ (a) detection heat map
  ng <- length(GENES); nd <- length(DS); ax <- axes_new(LEFT - 1.25, 11.8, nd * PITCH, ng * ROW_A, c(-0.5, nd - 0.5), c(ng - 0.5, -0.5))
  M <- matrix(NA_real_, ng, nd)
  for (j in seq_len(nd)) { d <- rd(DS[[j]][1], "unadjusted", "detection"); d$in_data <- d$in_data %in% c("True", "TRUE", TRUE); v <- ifelse(d$in_data, d$detected_pct, NA_real_); M[, j] <- v[match(GENES, d$gene_invitro)] }
  for (i in seq_len(ng)) for (j in seq_len(nd)) {
    v <- M[i, j]; x1 <- ax$X(j - 1 - 0.5); x2 <- ax$X(j - 1 + 0.5); y1 <- ax$Y(i - 1 - 0.5); y2 <- ax$Y(i - 1 + 0.5); xc <- ax$X(j - 1); yc <- ax$Y(i - 1)
    if (is.na(v)) { rect_mm(x1, y1, x2, y2, "white", "#d0d0cc", 0.4); txt(xc, yc, "absent", 6, 1, INK2, va = "center"); next }
    if (v >= 10) { col <- HEAT((v - 10) / 90); tcol <- if (luminance(col) < 0.30) "white" else INK } else { col <- "#e6e5e1"; tcol <- INK2 }
    rect_mm(x1, y1, x2, y2, col, "white", 0.7); txt(xc, yc, if (v >= 1) py_fmt_f(v, 0) else py_fmt_f(v, 1), 6.5, 1, tcol, va = "center")
  }
  ax_yticks(ax, 0:(ng - 1), as.list(GENES), 6.8, 3, len = 0, pad = 2.0)                                                                   # gene symbols in italics
  for (j in seq_len(nd)) txt(ax$X(j - 1), ax$top - 2.0 * PT, sprintf("%s\n%s\n(%d β-cells)", DS[[j]][2], DS[[j]][3], summ[[DS[[j]][1]]]$beta_cells), 6.8, 2, INK, va = "bottom", ls = 1.1)      # top tick labels: bottom of the box 2 pt above the axes
  ## ------------------------------------------------------------------ (b) the 21 consensus pairs
  nc <- nrow(cons); yp <- (nc - 1):0
  for (k in seq_along(DS)) {
    ax <- axes_new(LEFT + (k - 1) * PITCH, 78.6, PW, nc * ROW_B, c(-0.25, 0.72), c(-0.8, nc - 0.2)); cc <- rd(DS[[k]][1], "adjusted", "cons"); key <- paste(cc$gene1, cc$gene2)
    for (y in yp) seg(ax$left, ax$Y(y), ax$right, ax$Y(y), "#efefec", 0.5)                                                                # horizontal grid behind the data
    seg(ax$X(0), ax$top, ax$X(0), ax$bottom, INK, 0.5)
    for (i in seq_len(nc)) { y <- yp[i]; w <- match(paste(cons$gene1[i], cons$gene2[i]), key)
      if (!is.na(w)) { r <- cc[w, ]; seg(ax$X(r$null_mean), ax$Y(y), ax$X(r$null_95pct), ax$Y(y), NULLBAR, 2.8); marker(ax$X(r$rho), ax$Y(y), "o", sqrt(17), fill = BLUE, border = "white", bw_pt = 0.4) }
      else txt(ax$X(0.30), ax$Y(y), "NE", 6, 1, GREYLINE, va = "center") }
    ax_spines(ax); th <- ax_xticks(ax, c(0, 0.3, 0.6), c("0.0", "0.3", "0.6"), 6.5)
    ax_yticks(ax, yp, if (k == 1) as.list(paste0(cons$gene1, "–", cons$gene2)) else NULL, 6.8, 3)
    txt(ax$left + ax$w / 2, ax$top - 2.5 * PT, sprintf("%s\n%s\n%d evaluable", DS[[k]][2], DS[[k]][3], nrow(cc)), 6.8, 2, INK, va = "baseline", ls = 1.1)      # axes title: baseline of its last line 2.5 pt above the axes
  }
  mtext(LEFT + 3 * PITCH - 1.25, H - 0.8, list(list(list("Depth-adjusted partial Spearman "), list("ρ", 3))), 7, 1, INK, "center", "bottom", plus = list(MATH))
}
fig_save("Figure_7_islet_crosscheck_final", H, draw, OUT)
