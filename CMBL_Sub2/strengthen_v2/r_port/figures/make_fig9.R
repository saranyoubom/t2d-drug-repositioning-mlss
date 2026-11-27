## make_fig9.R - R port of wp6a_assets/make_fig9_final.py: FINAL Fig. 9 (benchmark part 2: calibration, prevalence and design guidance) for the JTM manuscript. Plots saved results of wp2_benchmark2 only; no new analysis.
## Six panels in three rows (calibration, application, design), 170 mm wide.
## (a) power of the three coupling tests (T_r2, T_p, T_eig) and of the consensus-count test when one latent factor is injected into 6 genes of decoupled public data (3 groups x 4 replicates), by within-group r;
##     leftmost points = type I error. (b) power of T_r2 by the number of coupled genes. (c) share of public data sets in which T_r2 detects within-group coupling (P < 0.05): gene set A (17 differentially
##     expressed genes, filled circles) and set B (17 random expressed genes, open squares), Wilson 95% intervals; dashed line = 5% expected without coupling. (d) observed adjusted q < 0.05 count divided
##     by its own null mean, by normalization scheme (gene set A, all 96 analyzed data sets). (e, f) design guidance by simulation (G = 17, one module of 6 genes, 3 groups): power of T_r2 and share of true
##     module pairs detected, by independent experiments per group.
## No key strip and no curve labels: markers, colors and shades are explained in the legend text. Statistical parameters are italic.
## Usage: Rscript make_fig9.R      Environment: WP2_RES = folder with calibration_tables.csv, validation_endpoints_long.csv, validation_per_series.csv, design_table.csv (default: the archived wp2_benchmark2;
## the R run is results_R/wp2), FIG_OUT = output folder (default results_R/figures). Output: Figure_9_benchmark_final.png / .pdf
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R"))
W2 <- Sys.getenv("WP2_RES", normalizePath(file.path(R_PORT, "..", "wp2_benchmark2"))); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures"))
rd <- function(name) read.csv(file.path(W2, name), stringsAsFactors = FALSE)
COL <- c(T_r2 = BLUE, T_p = AQUA, T_eig = GREYLINE, count = ORANGE); MK <- c(T_r2 = "o", T_p = "^", T_eig = "D", count = "s")      # blue circles T_r2, green triangles T_p, gray diamonds T_eig, orange squares count (explained in the legend text)
ORD3 <- c("#86b6ef", BLUE, "#0d366b"); XP <- c(`0` = 0, `0.3` = 1, `0.5` = 2, `0.7` = 3, `0.9` = 4); XT <- c("none", "0.3", "0.5", "0.7", "0.9"); YT <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
TBI <- list(list("T", 4), list("r", 4, 0.86, 1.5), list("2", 2, 0.86, 1.5))                  # T_r2 in bold italics for the bold panel titles (T and r italic, the index 2 upright)
## errorbar() of matplotlib: error lines (elinewidth 0.6 pt) with caps (2 x capsize wide, 1 pt thick), then the data line with markers (marker edge 1 pt in the line color), lowest first
err <- function(ax, x, y, lo, hi, col, shape, lwpt = 1.1, ms = 3.0, ew = 0.6, capsize = 1.3) {
  for (i in seq_along(x)) { X <- ax$X(x[i]); seg(X, ax$Y(lo[i]), X, ax$Y(hi[i]), col, ew); for (v in c(lo[i], hi[i])) seg(X - capsize * PT, ax$Y(v), X + capsize * PT, ax$Y(v), col, 1.0) }
  for (i in seq_len(length(x) - 1)) seg(ax$X(x[i]), ax$Y(y[i]), ax$X(x[i + 1]), ax$Y(y[i + 1]), col, lwpt)
  marker(ax$X(x), ax$Y(y), shape, ms, fill = col, border = col, bw_pt = 1.0)
}
H <- 150; LX <- 27; RX <- 106; AW <- 54; R1 <- 8; R2 <- 58.5; R3 <- 109; AH <- 31
draw <- function() {
  set_ls(1.1)
  for (p in list(list("(a)", 0.8, 0.6), list("(b)", 83, 0.6), list("(c)", 0.8, 51.1), list("(d)", 83, 51.1), list("(e)", 0.8, 101.6), list("(f)", 83, 101.6))) tag(p[[1]], p[[2]], p[[3]])
  P_share <- list(list("Share with ", 1), list("P", 3), list(" < 0.05", 1))
  lab_x <- function(ax, runs, th) ax_xlabel(ax, NULL, th, 7, runs = runs)
  frame <- function(ax, xt, xl, yt = YT, ylab = NULL, ytl = NULL, xtick_size = 6.5) {
    ax_spines(ax); th <- ax_xticks(ax, xt, xl, xtick_size); tw <- ax_yticks(ax, yt, if (is.null(ytl)) num(yt, 1) else ytl, 6.5); list(th = th, tw = tw)
  }
  ## ------------------------------------------------------------------ (a) power and type I error of the four tests
  T <- rd("calibration_tables.csv"); ax <- axes_new(LX, R1, AW, AH, c(-0.3, 4.3), c(-0.03, 1.08))
  dashed(ax$left, ax$Y(0.8), ax$right, ax$Y(0.8), INK, 0.5); dotted(ax$left, ax$Y(0.05), ax$right, ax$Y(0.05), INK, 0.5)
  for (s in c("T_eig", "T_p", "T_r2", "count")) {                                       # T_r2 drawn on top of the two statistics that almost coincide with it
    d <- rbind(T[T$statistic == s & T$m == 0, ], T[T$statistic == s & T$m == 6, ]); d <- d[order(d$r), ]; x <- XP[as.character(d$r)]
    err(ax, x, d$rate, d$ci_lo, d$ci_hi, COL[[s]], MK[[s]], lwpt = if (s %in% c("T_r2", "count")) 1.2 else 0.9, ms = if (s != "T_r2") 3.2 else 3.0)
  }
  f <- frame(ax, unname(XP), XT, ytl = num(YT, 1)); ax_xlabel(ax, NULL, f$th, runs = list(list("Injected within-group "), list("r", 3), list(" (6 genes)"))); ax_ylabel(ax, NULL, f$tw, runs = P_share); ax_title(ax, "Power and type I error")
  ## ------------------------------------------------------------------ (b) T_r2 by number of coupled genes
  ax <- axes_new(RX, R1, AW, AH, c(-0.3, 4.3), c(-0.03, 1.08)); null <- T[T$statistic == "T_r2" & T$m == 0, ][1, ]
  dashed(ax$left, ax$Y(0.8), ax$right, ax$Y(0.8), INK, 0.5); dotted(ax$left, ax$Y(0.05), ax$right, ax$Y(0.05), INK, 0.5)
  for (k in 1:3) {
    m <- c(3, 6, 10)[k]; d <- T[T$statistic == "T_r2" & T$m == m, ]; d <- d[order(d$r), ]; x <- c(0, XP[as.character(d$r)]); y <- c(null$rate, d$rate); lo <- c(null$ci_lo, d$ci_lo); hi <- c(null$ci_hi, d$ci_hi)
    err(ax, x, y, lo, hi, ORD3[k], "o")
  }
  f <- frame(ax, unname(XP), XT); ax_xlabel(ax, NULL, f$th, runs = list(list("Injected within-group "), list("r", 3))); ax_ylabel(ax, NULL, f$tw, runs = P_share)
  ax_title(ax, runs = c(TBI, list(list(" by number of coupled genes", 2))))
  ## ------------------------------------------------------------------ (c) prevalence of within-group coupling in public data
  L <- rd("validation_endpoints_long.csv"); L <- L[L$endpoint == "coupling_T_r2", ]
  ARMS <- list(c("Development arrays (37)", "Development\narrays"), c("Held-out GPL570", "Held-out\nhuman arrays"), c("Held-out GPL1261", "Held-out\nmouse arrays"), c("RNA-seq human", "RNA-seq\nhuman"), c("All validation sets pooled", "All validation\nsets"))
  ax <- axes_new(LX, R2, AW, AH, c(0, 1.0), c(-0.6, length(ARMS) - 0.4)); ypos <- setNames(length(ARMS) - 1 - seq_along(ARMS) + 1, vapply(ARMS, `[`, "", 1))
  dashed(ax$X(0.05), ax$top, ax$X(0.05), ax$bottom, INK, 0.5); seg(ax$left, ax$Y(0.5), ax$right, ax$Y(0.5), "#bdbdb8", 0.5)
  for (a in ARMS) for (g in list(list("A", 0.17, "o", BLUE), list("B", -0.17, "s", "white"))) {
    r <- L[L$arm == a[1] & L$gene_set == g[[1]], ][1, ]; y <- ypos[[a[1]]] + g[[2]]
    seg(ax$X(r$lo), ax$Y(y), ax$X(r$hi), ax$Y(y), BLUE, 1.0); marker(ax$X(r$estimate), ax$Y(y), g[[3]], sqrt(14), fill = g[[4]], border = BLUE, bw_pt = 0.8)
    txt(ax$X(1.03), ax$Y(y), sprintf("%d/%d", r$k, r$n), 6.3, adj = c(0, 0.5), va = "center")
  }
  ax_spines(ax); th <- ax_xticks(ax, c(0, 0.25, 0.5, 0.75, 1.0), c("0", "25%", "50%", "75%", "100%")); tw <- ax_yticks(ax, unname(ypos), vapply(ARMS, `[`, "", 2), 6.5, ls = 1.1)
  ax_xlabel(ax, "Data sets with within-group coupling", th); ax_title(ax, "Prevalence in public data")
  ## ------------------------------------------------------------------ (d) adjusted counts against their own null
  Pd <- rd("validation_per_series.csv"); Pd <- Pd[Pd$gene_set == "A", ]; ax <- axes_new(RX, R2, AW, AH, c(-0.55, 3.75), c(0.02, 250), ylog = TRUE); rng <- np_rng(1); FLOOR <- 0.03
  dashed(ax$left, ax$Y(1), ax$right, ax$Y(1), INK, 0.5); seg(ax$left, ax$Y(0.045), ax$right, ax$Y(0.045), "#bdbdb8", 0.5)
  for (i in 0:3) {
    s <- c("raw", "den1", "den2", "den3")[i + 1]; o <- Pd[[paste0("adjq_", s, "_obs")]]; nm <- Pd[[paste0("adjq_", s, "_null_mean")]]; r <- (o / nm)[!is.na(nm) & nm > 0]; r <- r[!is.na(r)]
    jit <- -0.2 + 0.4 * np_random(rng, length(r)); marker(ax$X(i + jit), ax$Y(pmax(r, FLOOR)), "o", sqrt(5), fill = BLUE, alpha = 0.5)
    q1 <- np_quantile(r, 0.25); md <- np_quantile(r, 0.5); q3 <- np_quantile(r, 0.75)
    seg(ax$X(i), ax$Y(max(q1, FLOOR)), ax$X(i), ax$Y(q3), INK, 0.9); seg(ax$X(i - 0.24), ax$Y(max(md, FLOOR)), ax$X(i + 0.24), ax$Y(max(md, FLOOR)), INK, 1.4)
    lab <- py_fmt_f(md, 2); xx <- ax$X(i + 0.30); yy <- ax$Y(max(md, FLOOR)); pad <- 0.1 * 6.3 * PT; hb <- sum(lpm(6.3)) / 2                  # white box behind the median value: text box (ascent + descent of "lp") plus pad 0.1 of the font size
    rect_mm(xx - pad, yy - hb - pad, xx + wmm(lab, 6.3) + pad, yy + hb + pad, adjustcolor("white", alpha.f = 0.85), NA); txt(xx, yy, lab, 6.3, adj = c(0, 0.5), va = "center")
  }
  ax_spines(ax); th <- ax_xticks(ax, 0:3, c("raw", "ratio to\ngene 1", "ratio to\ngene 2", "ratio to\ngene 3"), 6.5, ls = 1.1)
  tw <- ax_yticks(ax, c(FLOOR, 0.1, 1, 10, 100), c("0", "0.1", "1", "10", "100"), 6.5); ax_ylabel(ax, "Observed / null mean", tw)
  ax_title(ax, runs = list(list("Adjusted ", 2), list("q", 4), list(" < 0.05 counts against their null", 2)))
  ## ------------------------------------------------------------------ (e, f) design guidance
  D <- rd("design_table.csv"); NG <- c(4, 6, 8, 12, 20, 30)
  for (k in 1:2) {
    col <- c("power_T_r2", "pair_power")[k]; left <- c(LX, RX)[k]; ax <- axes_new(left, R3, AW, AH, c(log2(4) - 0.3, log2(30) + 0.35), c(-0.03, 1.08))
    dotted(ax$X(log2(4)), ax$top, ax$X(log2(4)), ax$bottom, INK, 0.5); if (k == 1) dashed(ax$left, ax$Y(0.8), ax$right, ax$Y(0.8), INK, 0.5)
    for (j in 1:4) { r <- c(0.3, 0.5, 0.7, 0.9)[j]; d <- D[D$r == r, ]; d <- d[order(d$n_per_group), ]; x <- log2(d$n_per_group); y <- d[[col]]
      for (i in seq_len(length(x) - 1)) seg(ax$X(x[i]), ax$Y(y[i]), ax$X(x[i + 1]), ax$Y(y[i + 1]), ORD4[j], 1.1); marker(ax$X(x), ax$Y(y), "o", 3.0, fill = ORD4[j], border = ORD4[j], bw_pt = 1.0) }
    ax_spines(ax); th <- ax_xticks(ax, log2(NG), as.character(NG), 6.5); tw <- ax_yticks(ax, YT, num(YT, 1), 6.5); ax_xlabel(ax, "Independent experiments per group", th)
    if (k == 1) { ax_ylabel(ax, NULL, tw, runs = P_share); ax_title(ax, runs = c(list(list("Power of ", 2)), TBI)) }
    else { ax_ylabel(ax, "Share of module pairs", tw); ax_title(ax, runs = list(list("True pairs detected (BH ", 2), list("q", 4), list(" < 0.05)", 2))) }
  }
}
fig_save("Figure_9_benchmark_final", H, draw, OUT)
