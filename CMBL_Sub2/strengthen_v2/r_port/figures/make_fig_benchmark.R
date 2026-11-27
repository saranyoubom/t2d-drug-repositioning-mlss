## make_fig_benchmark.R - R port of the figure of benchmark_route3/summarize_benchmark.py (Benchmark_figure.png / .pdf): the route-3 benchmark of 37 public data sets, gene set A (17 differentially expressed genes per data set).
## (a) consensus pairs of the pooled analysis: dot = observed (red when the permutation P < 0.05), grey dot and bar = permutation null mean +/- 1.96 SD, data sets sorted by z; (b) permutation P of the observed counts (blue) and with injected
## within-group coupling (red); (c) share of the consensus pairs with q < 0.05 after treatment adjustment (at least one scheme). Style of the original (matplotlib defaults with Arial 7.5 pt and spines of 0.6 pt, 6.69 x 7.0 in at 300 dpi).
## Usage: Rscript make_fig_benchmark.R   Environment: ROUTE3_RES (folder with per_dataset_results.csv, default results_R/benchmark_route3), FIG_OUT (default results_R/figures). Output: Benchmark_figure.png / .pdf
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R"))
RES <- Sys.getenv("ROUTE3_RES", file.path(R_PORT, "results_R", "benchmark_route3")); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures")); K <- "black"
T <- read.csv(file.path(RES, "per_dataset_results.csv"), stringsAsFactors = FALSE); A <- T[T$set == "A", ]; A <- A[order(A$z), ]; n <- nrow(A)
W <- 6.69 * 25.4 + 0.01; H <- 7.0 * 25.4 + 0.01                                                  # 2007 x 2100 px (the 0.01 mm guards against a floating-point truncation)
G <- gridspec_mm(W, H, 2, 2, c(1, 1), c(1.25, 1), 0.38, 0.32, 0.09, 0.985, 0.965, 0.075)
segp <- function(x1, y1, x2, y2, col = K, lwpt = 0.6) segments(x1, y1, x2, y2, col = col, lwd = lw(lwpt), lend = "square")      # spines: projecting caps
draw <- function() {
  FS <- 7.5; lp <- lpm(FS)
  axis_frame <- function(ax, xt, yt) {                                                            # left and bottom spine (0.6 pt), tick marks 3.5 pt x 0.8 pt, tick labels 7.5 pt; returns the width of the widest y tick label and the height of the x tick label boxes
    segp(ax$left, ax$top, ax$left, ax$bottom, K, 0.6); segp(ax$left, ax$bottom, ax$right, ax$bottom, K, 0.6)
    for (x in xt) seg(ax$X(x), ax$bottom, ax$X(x), ax$bottom + 3.5 * PT, K, 0.8)
    yl <- mpl_tick_labels(yt); w <- 0; for (i in seq_along(yt)) { seg(ax$left, ax$Y(yt[i]), ax$left - 3.5 * PT, ax$Y(yt[i]), K, 0.8); mtext(ax$left - 7 * PT, ax$Y(yt[i]), yl[i], FS, 1, K, "right", "center_baseline"); w <- max(w, wmm(yl[i], FS)) }
    w
  }
  ylab <- function(ax, tw, s) mtext(ax$left - 7 * PT - tw - 4 * PT, ax$top + ax$h / 2, s, FS, 1, K, "center", "bottom", rot = 90)
  ttl <- function(ax, s) mtext(ax$left, ax$top - 6 * PT, s, FS, 1, K, "left", "baseline")
  mtext(0.005 * W, (1 - 0.99) * H, "(a)", 9, 2, K, "left", "top"); mtext(0.005 * W, (1 - 0.40) * H, "(b)", 9, 2, K, "left", "top"); mtext(0.52 * W, (1 - 0.40) * H, "(c)", 9, 2, K, "left", "top")
  ## ------------------------------------------------------------------ (a) observed against the permutation null, per data set
  cl <- G(1, 1, 2); x <- 0:(n - 1); lo <- A$null_mean - 1.96 * A$null_sd; hi <- A$null_mean + 1.96 * A$null_sd; ylo <- min(lo, A$null_mean, A$observed); yhi <- max(hi, A$null_mean, A$observed); d <- yhi - ylo
  ylim <- c(ylo - 0.05 * d, yhi + 0.05 * d); ax <- axes_new(cl[["left"]], cl[["top"]], cl[["width"]], cl[["height"]], c(-0.8, n - 0.2), ylim)
  yt <- mpl_auto_ticks(ylim[1], ylim[2], ax$h / PT, FS, "y"); tw <- axis_frame(ax, x, yt)
  for (i in x + 1) seg(ax$X(i - 1), ax$Y(lo[i]), ax$X(i - 1), ax$Y(hi[i]), "#C9C9C9", 3)                    # LineCollection: butt caps
  marker(ax$X(x), ax$Y(A$null_mean), "o", sqrt(8), fill = "#8A8A8A", border = "#8A8A8A", bw_pt = 1.0)      # scatter(color=...) without linewidths: edge of 1 pt (lines.markeredgewidth) in the same color
  for (i in x + 1) marker(ax$X(i - 1), ax$Y(A$observed[i]), "o", 4, fill = if (A$perm_p[i] < 0.05) "#E64B35" else "#3C5488")
  m <- lpm(5.5); for (i in x + 1) run_line(ax$X(i - 1) + (m[1] - m[2]) / 2, ax$bottom + 7 * PT, norm_runs(A$gse[i]), 5.5, K, 90, 1)      # rotated 90 degrees, the end of the label at the tick label position
  ylab(ax, tw, "Consensus pairs (pooled analysis)"); ttl(ax, "Gene set A: 17 differentially expressed genes per dataset; dot = observed, grey bar = permutation null mean ± 1.96 SD")
  ## ------------------------------------------------------------------ (b) permutation P of the observed and the spiked data sets
  cl <- G(2, 1); bins <- seq(0, 1, length.out = 11); pa <- T$perm_p[T$set == "A"]; ps <- T$perm_p[T$set == "A_spike"]
  hcount <- function(v) vapply(1:10, function(i) sum(v >= bins[i] & (if (i < 10) v < bins[i + 1] else v <= bins[i + 1])), 1)       # np.histogram: [left, right) bins, the last one closed
  ha <- hcount(pa); hs <- hcount(ps); xc_a <- bins[-11] + 0.05 - 0.036; xc_s <- bins[-11] + 0.05 + 0.036; bw <- 0.036
  xr <- c(min(xc_a - bw / 2), max(xc_s + bw / 2)); xlim <- c(xr[1] - 0.05 * diff(xr), xr[2] + 0.05 * diff(xr)); ylim <- c(0, 1.05 * max(ha, hs))
  ax <- axes_new(cl[["left"]], cl[["top"]], cl[["width"]], cl[["height"]], xlim, ylim); xt <- mpl_auto_ticks(xlim[1], xlim[2], ax$w / PT, FS, "x"); yt <- mpl_auto_ticks(ylim[1], ylim[2], ax$h / PT, FS, "y"); tw <- axis_frame(ax, xt, yt)
  for (i in 1:10) if (ha[i] > 0) rect_mm(ax$X(xc_a[i] - bw / 2), ax$Y(ha[i]), ax$X(xc_a[i] + bw / 2), ax$Y(0), "#3C5488", NA)
  for (i in 1:10) if (hs[i] > 0) rect_mm(ax$X(xc_s[i] - bw / 2), ax$Y(hs[i]), ax$X(xc_s[i] + bw / 2), ax$Y(0), "#E64B35", NA)
  xl <- mpl_tick_labels(xt); for (i in seq_along(xt)) mtext(ax$X(xt[i]), ax$bottom + 7 * PT, xl[i], FS, 1, K, "center", "top")
  mtext(ax$left + ax$w / 2, ax$bottom + 7 * PT + sum(lp) + 4 * PT, "Permutation P (null ≥ observed)", FS, 1, K, "center", "top"); ylab(ax, tw, "Datasets"); ttl(ax, "Observed (blue) and with injected within-group coupling (red)")
  ## ------------------------------------------------------------------ (c) share of consensus pairs retained after adjustment for treatment
  cl <- G(2, 2); v <- A$retained_any[!is.na(A$retained_any)]; hc <- hcount(v); xlim <- c(-0.05, 1.05); ylim <- c(0, 1.05 * max(hc))
  ax <- axes_new(cl[["left"]], cl[["top"]], cl[["width"]], cl[["height"]], xlim, ylim); xt <- mpl_auto_ticks(xlim[1], xlim[2], ax$w / PT, FS, "x"); yt <- mpl_auto_ticks(ylim[1], ylim[2], ax$h / PT, FS, "y"); tw <- axis_frame(ax, xt, yt)
  for (i in 1:10) if (hc[i] > 0) rect_mm(ax$X(bins[i]), ax$Y(hc[i]), ax$X(bins[i + 1]), ax$Y(0), "#4DBBD5", NA)
  xl <- mpl_tick_labels(xt); for (i in seq_along(xt)) mtext(ax$X(xt[i]), ax$bottom + 7 * PT, xl[i], FS, 1, K, "center", "top")
  mtext(ax$left + ax$w / 2, ax$bottom + 7 * PT + sum(lp) + 4 * PT, list("Consensus pairs with q < 0.05 after", "treatment adjustment (share, ≥ 1 scheme)"), FS, 1, K, "center", "top", ls = 1.2); ylab(ax, tw, "Datasets")
}
fig_save("Benchmark_figure", H, draw, OUT, W)
