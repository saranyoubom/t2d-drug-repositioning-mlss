## make_figA51.R - R port of wp6a_assets/make_figA51_final.py: FINAL Fig. A5.1 (positive control of the islet cross-check) for Additional file 5 of the JTM manuscript. Plots saved results of wp4_islets/positive_control only; no new analysis.
## (a) share of 200 simulations in which the pooled verdict of the five human datasets was "no meaningful excess co-expression" (false reassurance) after coupling of nominal within-donor correlation rho was added to the genes of
##     the evaluable consensus pairs: solid line = added only to cells in which a gene is detected (conservative, the primary model), dashed line = added to all cells (upper bound); dotted line = the 20% limit set in advance.
## (b) excess E of ribosomal protein pairs (dots, mean of 100 random sets of 11 pairs, whiskers = SD of the sets) and of the consensus pairs (squares, +/-1.96 SE) in each dataset; the shaded band is the +/-0.10 margin.
## No key strip and no in-image labels of the lines, dots or squares (the line and marker meanings are in the legend text); only the numbers 2.5% and 20% and the "+/-0.10 margin" label stand in the plot. Statistical parameters
## (rho, E, z) are italic. 170 mm wide, 300 dpi PNG and PDF. Usage: Rscript make_figA51.R   Environment: WP4_PC (folder with positive_control_pooled.csv, rp_control.csv, positive_control_by_dataset.csv; default the archived
## wp4_islets/positive_control), FIG_OUT (default results_R/figures). Output: Figure_A5-1_positive_control_final.png / .pdf
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R"))
PC <- Sys.getenv("WP4_PC", normalizePath(file.path(R_PORT, "..", "wp4_islets", "positive_control"))); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures"))
POOL <- read.csv(file.path(PC, "positive_control_pooled.csv"), stringsAsFactors = FALSE); RP <- read.csv(file.path(PC, "rp_control.csv"), stringsAsFactors = FALSE); rownames(RP) <- RP$tag
BY <- read.csv(file.path(PC, "positive_control_by_dataset.csv"), stringsAsFactors = FALSE)
H <- 66
draw <- function() {
  tag("(a)"); tag("(b)", x = 84.0)
  ## ------------------------------------------------------------------ (a) false reassurance against the nominal coupling
  ax <- axes_new(17.0, 6.5, 60.0, 47.0, c(-0.35, 4.75), c(-0.04, 1.06)); RHOS <- c(0.05, 0.10, 0.20, 0.30, 0.50); xs <- seq_along(RHOS) - 1
  det <- POOL[POOL$mode == "detected", ]; det <- det[match(RHOS, det$rho), ]; alln <- POOL[POOL$mode == "all", ]; alln <- alln[match(RHOS, alln$rho), ]
  dotted(ax$X(-0.35), ax$Y(0.20), ax$X(4.28), ax$Y(0.20), INK, 0.5)
  pline_dashed(ax$X(xs), ax$Y(alln$share_no_meaningful_excess), BLUE, 1.1, 3, 2); marker(ax$X(xs), ax$Y(alln$share_no_meaningful_excess), "s", 2.8, fill = "white", border = BLUE, bw_pt = 0.8)
  pline(ax$X(xs), ax$Y(det$share_no_meaningful_excess), BLUE, 1.2); marker(ax$X(xs), ax$Y(det$share_no_meaningful_excess), "o", 3.0, fill = BLUE, border = BLUE, bw_pt = 1.0)
  s30 <- det$share_no_meaningful_excess[RHOS == 0.30]
  txt(ax$X(3.0), ax$Y(s30 + 0.075), paste0(py_fmt_f(100 * s30, 1), "%"), 6.5, adj = c(0.5, 0), va = "bottom"); txt(ax$X(4.36), ax$Y(0.20), "20%", 6.5, adj = c(0, 0.5), va = "center")
  ax_spines(ax); th <- ax_xticks(ax, xs, c("0.05", "0.1", "0.2", "0.3", "0.5"), 6.5); tw <- ax_yticks(ax, c(0, 0.2, 0.4, 0.6, 0.8, 1.0), c("0", "20", "40", "60", "80", "100"), 6.5)
  ax_xlabel(ax, NULL, th, runs = list(list("Nominal coupling "), list("ρ", 3), list(" added to the consensus genes"))); ax_ylabel(ax, "Simulations reading \"no meaningful excess\" (%)", tw)
  ## ------------------------------------------------------------------ (b) ribosomal protein pairs against the consensus pairs
  ax <- axes_new(108.0, 6.5, 58.0, 47.0, c(-0.12, 0.56), c(-0.6, 5.6))
  ORDER <- list(c("gse81608", "GSE81608"), c("gse84133_human", "GSE84133 human"), c("gse86469", "GSE86469"), c("gse81547", "GSE81547"), c("emtab5061", "E-MTAB-5061"), c("gse84133_mouse", "GSE84133 mouse"))
  ys <- setNames(length(ORDER) - 1 - seq_along(ORDER) + 1, vapply(ORDER, `[`, "", 1))
  rect_mm(ax$X(-0.10) - 0.5 * PT, ax$top, ax$X(0.10) + 0.5 * PT, ax$bottom, "#ececea", NA); seg(ax$X(0), ax$top, ax$X(0), ax$bottom, INK, 0.5); seg(ax$left, ax$Y(0.5), ax$right, ax$Y(0.5), "#bdbdb8", 0.5)
  for (o in ORDER) {
    t <- o[1]; y <- ys[[t]]; b <- BY[BY$tag == t & BY$mode == "detected" & BY$rho == 0.3, ][1, ]; r <- RP[t, ]
    seg(ax$X(r$E_RP_mean - r$E_RP_sd), ax$Y(y + 0.17), ax$X(r$E_RP_mean + r$E_RP_sd), ax$Y(y + 0.17), BLUE, 1.0); marker(ax$X(r$E_RP_mean), ax$Y(y + 0.17), "o", sqrt(15), fill = BLUE)
    seg(ax$X(b$E_primary - 1.96 * b$SE_E), ax$Y(y - 0.17), ax$X(b$E_primary + 1.96 * b$SE_E), ax$Y(y - 0.17), BLUE, 0.9); marker(ax$X(b$E_primary), ax$Y(y - 0.17), "s", sqrt(15), fill = "white", border = BLUE, bw_pt = 0.8)
  }
  txt(ax$X(0.0), ax$top - 0.015 * ax$h, "±0.10 margin", 6.3, 1, INK2, adj = c(0.5, 0), va = "bottom")
  ax_spines(ax); tw <- ax_yticks(ax, unname(ys), vapply(ORDER, `[`, "", 2), 6.5); th <- ax_xticks(ax, c(0, 0.1, 0.2, 0.3, 0.4, 0.5), c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5"), 6.5)
  ax_xlabel(ax, NULL, th, runs = list(list("Excess co-expression "), list("E", 3), list(" (Fisher "), list("z", 3), list(")")))
}
fig_save("Figure_A5-1_positive_control_final", H, draw, OUT)
