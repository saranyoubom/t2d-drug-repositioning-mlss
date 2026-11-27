## make_figS1.R - R port of reanalysis_17genes/fig_s1.py: Fig. S1 of the 17-gene reanalysis (observed counts against the group-preserving permutation null, and significant pairs before and after adjustment for treatment).
## (a) left: pairs with pooled p < 0.05 in the four schemes, observed (diamond) against the 2000 residual permutations (violin = density, thin bar = 2.5th to 97.5th percentile, thick bar = quartiles; a constant null is a line);
##     right: consensus pairs (all four schemes; plus the same sign), observed 21. (b) significant pairs per scheme: pooled p, pooled BH q, adjusted p, adjusted BH q (df = 8; only |r| >= 0.90 can pass BH at q < 0.05).
## Style of the original (matplotlib defaults, Arial, 6.8 x 5.9 in at 300 dpi = 2040 x 1770 px), not that of Figs 7-10. Usage: Rscript make_figS1.R   Environment: RA_RES (folder with perm_arrays_17.csv and table_counts.json, default
## results_R/reanalysis_17genes), FIG_OUT (default results_R/figures). Output: fig_s1_17genes.png / .pdf
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R")); suppressMessages(library(jsonlite))
RA <- Sys.getenv("RA_RES", file.path(R_PORT, "results_R", "reanalysis_17genes")); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures"))
A <- read.csv(file.path(RA, "perm_arrays_17.csv"), check.names = FALSE); TC <- fromJSON(file.path(RA, "table_counts.json"), simplifyVector = FALSE)[["17"]]
NPG <- list(dark = "#3C5488", light = "#4DBBD5", salmon = "#F39B7F", red = "#E64B35", grey = "#C9CEE0"); K <- "black"
obs <- c(Raw = 66, `Hes1-normalized` = 120, `Kcnj11-normalized` = 84, `Ins1-normalized` = 80)
for (k in names(obs)) stopifnot("observed counts of the figure differ from table_counts.json" = obs[[k]] == TC[[k]]$pooled_p)
labels <- list(Raw = list(NULL, "Raw", "expression"), `Hes1-normalized` = list("Hes1", "-", "normalized"), `Kcnj11-normalized` = list("Kcnj11", "-", "normalized"), `Ins1-normalized` = list("Ins1", "-", "normalized"))
W <- 6.8 * 25.4 + 0.01; H <- 5.9 * 25.4 + 0.01                                                    # 2040 x 1770 px (the 0.01 mm guards against a floating-point truncation to 2039 or 1769 px)
G <- gridspec_mm(W, H, 2, 2, c(2.2, 1), c(1, 1), 0.55, 0.28, 0.1, 0.985, 0.94, 0.1)
segp <- function(x1, y1, x2, y2, col = K, lwpt = 0.8) segments(x1, y1, x2, y2, col = col, lwd = lw(lwpt), lend = "square")        # matplotlib Line2D / spine default: projecting caps
kde <- function(arr, width = 0.8, n = 100) {                                                       # violinplot: Gaussian KDE (Scott's factor n^(-1/5) x sd), 100 points over [min, max], widest half-width 0.5 x width
  h <- sd(arr) * length(arr)^(-1 / 5); coords <- seq(min(arr), max(arr), length.out = n); vals <- vapply(coords, function(x) mean(dnorm(x, arr, h)), 1)
  list(coords = coords, half = 0.5 * width * vals / max(vals))
}
auto_lim <- function(lo, hi) { d <- hi - lo; c(lo - 0.05 * d, hi + 0.05 * d) }                       # matplotlib default margins (5%)
draw <- function() {
  tick_labels <- function(ax, ticks, labs, size = 7.5) { w <- 0; for (i in seq_along(ticks)) { mtext(ax$left - 3.5 * PT, ax$Y(ticks[i]), labs[i], size, 1, K, "right", "center_baseline"); w <- max(w, wmm(labs[i], size)) }; w }
  frame <- function(ax, yt, ylabs, xs) { for (t in yt) segp(ax$left, ax$Y(t), ax$right, ax$Y(t), "#E6E6E6", 0.8); segp(ax$left, ax$bottom, ax$right, ax$bottom, K, 0.8); for (x in xs) seg(ax$X(x), ax$bottom, ax$X(x), ax$bottom + 3.5 * PT, K, 0.8) }
  violin <- function(ax, x, arr, observed, label_dx = 0.18) {
    if (sd(arr) * sqrt((length(arr) - 1) / length(arr)) > 0) {
      kv <- kde(arr); polygon(ax$X(c(x - kv$half, rev(x + kv$half))), ax$Y(c(kv$coords, rev(kv$coords))), col = adjustcolor(NPG$grey, alpha.f = 0.85), border = NA)
      q1 <- np_quantile(arr, 0.25); q3 <- np_quantile(arr, 0.75); lo <- np_quantile(arr, 0.025); hi <- np_quantile(arr, 0.975)
      segp(ax$X(x), ax$Y(lo), ax$X(x), ax$Y(hi), NPG$dark, 1.0); segp(ax$X(x), ax$Y(q1), ax$X(x), ax$Y(q3), NPG$dark, 3.2)
    } else segp(ax$X(x - 0.25), ax$Y(arr[1]), ax$X(x + 0.25), ax$Y(arr[1]), NPG$dark, 1.2)
    marker(ax$X(x), ax$Y(observed), "D", sqrt(42), fill = "white", border = "black", bw_pt = 1.3)
    mtext(ax$X(x + label_dx), ax$Y(observed), as.character(observed), 7.5, 1, K, "left", "center")
  }
  cat_label <- function(ax, x, gene, l1, l2) {                                                     # two-line category label: box 7 pt below the axis, boxes of the lines 2 pt apart (AnnotationBbox of a VPacker)
    ytop <- ax$bottom + 7 * PT; mtext(ax$X(x), ytop, list(if (is.null(gene)) list(list(l1, 1)) else list(list(gene, 3), list(l1, 1))), 7.5, 1, K, "center", "top")
    mtext(ax$X(x), ytop + sum(lpm(7.5)) + 2 * PT, list(l2), 7.5, 1, K, "center", "top")
  }
  ylab <- function(ax, tw, runs, size = 8) mtext(ax$left - 3.5 * PT - tw - 4 * PT, ax$top + ax$h / 2, list(runs), size, 2, K, "center", "bottom", rot = 90, plus = list(MATH))
  title <- function(ax, runs, size, pad) mtext(ax$left + ax$w / 2, ax$top - pad * PT, list(runs), size, 2, K, "center", "baseline")
  mtext(0.005 * W, (1 - 0.965) * H, "(a)", 10, 2, K, "left", "baseline"); mtext(0.005 * W, (1 - 0.485) * H, "(b)", 10, 2, K, "left", "baseline")
  ## ------------------------------------------------------------------ (a) left
  cl <- G(1, 1); ks <- names(obs); sdz <- vapply(ks, function(k) sd(A[[k]]) * sqrt((length(A[[k]]) - 1) / length(A[[k]])) == 0, TRUE)
  xl <- auto_lim(min(ifelse(sdz, 0.25, 0.4) * -1 + 0:3), max(ifelse(sdz, 0.25, 0.4) + 0:3)); ax <- axes_new(cl[["left"]], cl[["top"]], cl[["width"]], cl[["height"]], xl, c(55, 128))
  frame(ax, seq(60, 120, 10), NULL, 0:3); for (i in seq_along(ks)) violin(ax, i - 1, A[[ks[i]]], obs[[ks[i]]], if (sdz[i]) 0.32 else 0.18)
  for (i in seq_along(ks)) { lab <- labels[[ks[i]]]; cat_label(ax, i - 1, lab[[1]], lab[[2]], lab[[3]]) }
  tw <- tick_labels(ax, seq(60, 120, 10), as.character(seq(60, 120, 10))); ylab(ax, tw, list(list("Pairs with pooled "), list("p", 3), list(" < 0.05"))); title(ax, list(list("Observed counts vs group-preserving null")), 8.5, 6)
  ## ------------------------------------------------------------------ (a) right
  cl <- G(1, 2); ax <- axes_new(cl[["left"]], cl[["top"]], cl[["width"]], cl[["height"]], auto_lim(-0.4, 1.4), c(12, 34)); yt <- seq(12.5, 32.5, 2.5)
  frame(ax, yt, NULL, 0:1); violin(ax, 0, A$cons, 21); violin(ax, 1, A$cons_sign, 21)
  for (x in 0:1) mtext(ax$X(x), ax$bottom + 7 * PT, if (x == 0) list("All four", "schemes") else list("+ same", "sign"), 7.5, 1, K, "center", "top", ls = 1.2)
  tw <- tick_labels(ax, yt, py_fmt_f(yt, 1)); ylab(ax, tw, list(list("Consensus pairs"))); title(ax, list(list("Consensus counts")), 8.5, 6)
  ## ------------------------------------------------------------------ (b)
  cl <- G(2, 1, 2); w <- 0.19; hw <- w * 0.94 / 2; xl <- auto_lim(-1.5 * w - hw, 4 + 1.5 * w + hw); ax <- axes_new(cl[["left"]], cl[["top"]], cl[["width"]], cl[["height"]], xl, c(0, 175))
  groups <- c("Raw", "Hes1-normalized", "Kcnj11-normalized", "Ins1-normalized", "all"); glab <- list(labels$Raw, labels$`Hes1-normalized`, labels$`Kcnj11-normalized`, labels$`Ins1-normalized`, list(NULL, "All four", "schemes"))
  series <- list(list("pooled_p", list(list("pooled "), list("p", 3), list(" < 0.05")), NPG$dark), list("pooled_q", list(list("pooled BH "), list("q", 3), list(" < 0.05")), NPG$light),
                 list("adj_p", list(list("adjusted "), list("p", 3), list(" < 0.05")), NPG$salmon), list("adj_q", list(list("adjusted BH "), list("q", 3), list(" < 0.05")), NPG$red))
  yt <- seq(0, 160, 20); frame(ax, yt, NULL, 0:4)
  lp <- lpm(6.8)
  for (si in seq_along(series)) for (gi in seq_along(groups)) {
    xc <- (gi - 1) + (si - 1 - 1.5) * w; v <- TC[[groups[gi]]][[series[[si]][[1]]]]; rect_mm(ax$X(xc - hw), ax$Y(v), ax$X(xc + hw), ax$Y(0), series[[si]][[3]], NA)
    run_line(ax$X(xc) + (lp[1] - lp[2]) / 2, ax$Y(v + 2), norm_runs(as.character(v)), 6.8, K, 90, 0)                       # value above the bar: rotated text, bottom of its box 2 units above the bar, centred on it
  }
  for (gi in seq_along(groups)) cat_label(ax, gi - 1, glab[[gi]][[1]], glab[[gi]][[2]], glab[[gi]][[3]])
  tw <- tick_labels(ax, yt, as.character(yt)); ylab(ax, tw, list(list("Significant pairs")))
  title(ax, list(list("Significant pairs before and after adjustment for treatment ("), list("df", 3), list(" = 8); only |"), list("r", 3), list("| ≥ 0.90 can pass BH at "), list("q", 3), list(" < 0.05")), 8, 18)
  fs <- 7; hl <- 1.2 * fs * PT; hp <- 0.8 * fs * PT; cs <- 1.4 * fs * PT                                                        # legend: 4 columns, no frame, lower centre at the top of the axes (borderaxespad 0.5, borderpad 0.4)
  tws <- vapply(series, function(s) sum(vapply(norm_runs(s[[2]]), function(r) wmm(r$s, fs, r$face), 1)), 1); x <- ax$left + ax$w / 2 - (sum(hl + hp + tws) + 3 * cs) / 2
  yb <- ax$top - 0.5 * fs * PT - 0.4 * fs * PT - 8 * 25.4 / 300                                                                   # baseline: the descent of a mathtext line (8 px) above the padded bottom
  for (si in seq_along(series)) { rect_mm(x, yb - 0.7 * fs * PT, x + hl, yb, series[[si]][[3]], NA); run_line(x + hl + hp, yb, norm_runs(series[[si]][[2]]), fs, K, 0, 0); x <- x + hl + hp + tws[si] + cs }
}
fig_save("fig_s1_17genes", H, draw, OUT, W)
