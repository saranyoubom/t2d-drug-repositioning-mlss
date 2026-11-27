## figstyle.R - R port of strengthen_v2/wp6a_assets/figstyle.py: shared style and helpers for the final figures of the JTM manuscript (Figs 7 to 10, A5.1, S1, benchmark figure).
## The figures are drawn with base R graphics on a canvas whose user coordinates are MILLIMETERS from the top-left corner of a 170 mm wide figure (x to the right, y downward), the same way the Python
## scripts place their axes ("ax_mm"): an axes box is given by left, top, width and height in mm and a data range; ticks, labels and titles follow the rules of the matplotlib style of the originals
## (tick length 2.5 pt, tick pad 3.5 pt, label pad 4 pt, title pad 3 pt, spines 0.6 pt, Arial, sizes in points at print size). Fonts: Arial through ragg/systemfonts (PNG, 300 dpi) and cairo_pdf (PDF,
## fonts embedded). No color key is ever drawn inside a figure (user rule of 4 Nov 2026, applied again on 7 Nov): no key strip, no marker-plus-text key and no text that only names a curve, line, symbol or
## shade; they are explained in the legend text. Statistical parameters (P, q, r, E, z, T, rho) are italic (rich() runs with face 3), so are gene symbols (font 3).
## Usage: source(this file); define draw <- function() {...} using the helpers; fig_save("stem", height_mm, draw, outdir). Run in R >= 4.3 with ragg installed.
suppressMessages(library(ragg))
FIGW <- 170                                                     # mm, JTM full page width
FONT <- "Arial"; PS <- 10                                       # device point size; text sizes are given in points and converted to cex = size / PS
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; NULLBAR <- "#c9c9c9"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"; GREYLINE <- "#6f6e69"; BLUE_TXT <- "#1c5cab"
RAMP <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#1c5cab", "#0d366b")                      # blue, light to dark (reference sequential hue)
ORD4 <- c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")                                             # ordinal ramp for r = 0.3, 0.5, 0.7, 0.9
MINUS <- "−"; PT <- 25.4 / 72                                                               # mm per point
.FIG <- new.env(); .FIG$ls <- 1.2; .FIG$tdy <- 0                                              # current matplotlib linespacing (default 1.2); vertical text correction (mm, see canvas)

## matplotlib LinearSegmentedColormap.from_list: nodes evenly spaced, 256-entry table, value x in [0, 1] maps to table entry floor(x * 256) (capped at 255)
HEAT <- function(x) {
  k <- pmin(floor(x * 256), 255); pos <- k / 255; nodes <- seq(0, 1, length.out = length(RAMP)); rgbm <- t(col2rgb(RAMP)) / 255
  vapply(seq_along(pos), function(i) { r <- vapply(1:3, function(j) approx(nodes, rgbm[, j], xout = pos[i])$y, 1); rgb(r[1], r[2], r[3]) }, "")
}
luminance <- function(col) { v <- as.numeric(col2rgb(col)) / 255; v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4); 0.2126 * v[1] + 0.7152 * v[2] + 0.0722 * v[3] }
lw <- function(pt) pt * 96 / 72                                                                   # R line width (1/96 inch) from points
cx <- function(size) size / PS
## numbers with a true minus sign; Python's '%+.3f' / '%.3f' (needs py_fmt_f of pycompat.R, source it first)
sg <- function(v, nd = 3) gsub("-", MINUS, paste0(ifelse(v < 0 | (v == 0 & 1 / v < 0), "-", "+"), py_fmt_f(abs(v), nd)), fixed = TRUE)
num <- function(v, nd = 3) gsub("-", MINUS, py_fmt_f(v, nd), fixed = TRUE)

## ---------------------------------------------------------------------------------------------------------------- canvas and devices
canvas <- function(H, grid = FALSE, W = FIGW) {      # grid = TRUE (PNG): map millimeters to pixels the way matplotlib's Agg renderer does: 300/25.4 px per mm from the left edge, and from the top edge of the float-size canvas (the image has
  .FIG$H <- H; .FIG$tdy <- 0                #   only the integer part of its size, so shapes sit delta = H x 300/25.4 - rows px higher). Text is not shifted by delta in matplotlib and its glyph bitmap is placed 0.6 px lower than the
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", lend = "butt", ljoin = "mitre", xpd = NA, family = FONT, bg = "white", ps = PS, cex = 1)      # nominal baseline: tdy (mm) is added to the y of every text
  plot.new()
  if (grid) { s <- 300 / 25.4; px <- dev.size("px"); delta <- H * s - px[2]; plot.window(xlim = c(0, px[1] / s), ylim = c(H, delta / s)); .FIG$tdy <- (delta + 0.6) / s }
  else plot.window(xlim = c(0, W), ylim = c(H, 0))
}
fig_save <- function(stem, H, draw, outdir, W = FIGW) {      # PNG at 300 dpi and PDF with embedded fonts, the same drawing (W: figure width in mm, default the 170 mm of the journal)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE); png <- file.path(outdir, paste0(stem, ".png")); pdf <- file.path(outdir, paste0(stem, ".pdf"))
  agg_png(png, width = W, height = H, units = "mm", res = 300, pointsize = PS, background = "white"); canvas(H, grid = TRUE, W = W); draw(); invisible(dev.off())
  cairo_pdf(pdf, width = W / 25.4, height = H / 25.4, family = FONT, pointsize = PS); canvas(H, W = W); draw(); invisible(dev.off())
  cat("saved", png, "|", pdf, "\n"); invisible(c(png, pdf))
}

## ---------------------------------------------------------------------------------------------------------------- primitives (mm coordinates, sizes in points)
txt <- function(x, y, s, size = 7, face = 1, col = INK, adj = c(0.5, 0.5), rot = 0, va = NULL, ls = .FIG$ls, plus = NULL) {      # with va (baseline, bottom, center, top, center_baseline) the text is placed by the layout rules of matplotlib (see mtext)
  if (!is.null(va)) return(invisible(mtext(x, y, as.list(strsplit(s, "\n", fixed = TRUE)[[1]]), size, face, col, c("left", "center", "right")[round(adj[1] * 2) + 1], va, rot, ls, plus)))
  text(x, y + .FIG$tdy, labels = s, cex = cx(size), font = face, col = col, adj = adj, srt = rot)
}
wmm <- function(s, size = 7, face = 1) strwidth(s, units = "user", cex = cx(size), font = face)              # width of a text in mm at print size
hmm <- function(s, size = 7, face = 1) abs(strheight(s, units = "user", cex = cx(size), font = face))      # height (cap height of the first line plus the line pitch of the others) in mm (the y axis runs downward, so strheight is negative)
set_ls <- function(ls) { .FIG$ls <- ls; par(lheight = (0.21 + 0.73 * ls) / 1.2) }                           # matplotlib linespacing ls: baseline pitch = descent of "lp" (0.21 em) + ls x ascent of "l" (0.73 em); R: 1.2 x lheight em
## a line of text made of runs in different styles (e.g. italic statistical symbols between upright words): runs = list(list(s, face, size_factor, dy_pt), ...); dy_pt > 0 lowers the run (subscript)
rich <- function(x, y, runs, size = 7, adj = c(0.5, 0.5), col = INK, rot = 0) {
  runs <- lapply(runs, function(r) list(s = r[[1]], face = if (length(r) >= 2) r[[2]] else 1, f = if (length(r) >= 3) r[[3]] else 1, dy = if (length(r) >= 4) r[[4]] else 0))
  w <- vapply(runs, function(r) wmm(r$s, size * r$f, r$face), 1); W <- sum(w); a <- (rot * pi / 180); ux <- cos(a); uy <- -sin(a); px <- sin(a); py <- cos(a)      # along the baseline (y downward) and perpendicular (down in the text frame)
  cap <- hmm("M", size); off <- -adj[1] * W; base <- adj[2] * cap                                                                                        # all runs share one baseline, shifted so that the requested vertical anchor sits at y
  for (i in seq_along(runs)) { r <- runs[[i]]; sh <- base + r$dy * PT; xi <- x + off * ux + sh * px; yi <- y + off * uy + sh * py + .FIG$tdy; text(xi, yi, labels = r$s, cex = cx(size * r$f), font = r$face, col = col, adj = c(0, 0), srt = rot); off <- off + w[i] }
  invisible(W)
}
## text with the layout rules of matplotlib (Text._get_layout), so that every text of the R figures stands where it stands in the matplotlib originals. matplotlib aligns text (va) and spaces the lines by the metrics of "lp":
## ascent of "l" and descent of "p" in whole pixels at 300 dpi (table LPT, measured with the renderer of the originals); a line with mathtext is 2 px higher and 2 px deeper (plus = list(MATH)), one with a subscript deeper still.
## Baseline of line i below the top of the text box: B1 = ascent; Bi = B(i-1) + d(i-1) + max(ls x ascent of "l", ls x ascent of the line); box height = Bn + dn; va = top / center / bottom / baseline / center_baseline (single line).
LP_SIZE <- c(5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11)
LP_ASC <- c(15, 15, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 30, 30, 30, 30, 30, 31, 31, 31, 31, 32, 32, 32, 34)
LP_D <- c(4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9)
lpm <- function(size) { k <- match(round(size, 1), LP_SIZE); v <- if (is.na(k)) c(round(0.745 * size * 300 / 72), round(0.205 * size * 300 / 72)) else c(LP_ASC[k], LP_D[k]); v * 25.4 / 300 }      # (ascent, descent) in mm
MATH <- c(2, 2)                                                                                                       # extra ascent and descent (px) of a line that contains mathtext (italic statistical symbols)
norm_runs <- function(l, face = 1) if (is.character(l)) list(list(s = l, face = face, f = 1, dy = 0)) else lapply(l, function(r) list(s = r[[1]], face = if (length(r) >= 2) r[[2]] else face, f = if (length(r) >= 3) r[[3]] else 1, dy = if (length(r) >= 4) r[[4]] else 0))
run_line <- function(xb, yb, runs, size, col = INK, rot = 0, hf = 0.5) {      # one line of runs on the baseline through (xb, yb); hf: 0 = the line starts at the point, 0.5 = centred on it, 1 = ends at it; rot degrees counter-clockwise
  w <- vapply(runs, function(r) wmm(r$s, size * r$f, r$face), 1); W <- sum(w); a <- rot * pi / 180; ex <- cos(a); ey <- -sin(a); nx <- -sin(a); ny <- -cos(a); off <- -hf * W      # (ex, ey) along the baseline, (nx, ny) up; y runs downward
  for (i in seq_along(runs)) { r <- runs[[i]]; text(xb + off * ex - r$dy * PT * nx, yb + off * ey - r$dy * PT * ny + .FIG$tdy, labels = r$s, cex = cx(size * r$f), font = r$face, col = col, adj = c(0, 0), srt = rot); off <- off + w[i] }
  invisible(W)
}
mtext <- function(x, y, lines, size = 7, face = 1, col = INK, ha = "center", va = "baseline", rot = 0, ls = .FIG$ls, plus = NULL, mode = "anchor") {   # lines: list of strings or run lists; plus: list of c(ascent, descent) px added per line
  lines <- lapply(as.list(lines), norm_runs, face = face); n <- length(lines); lp <- lpm(size); hf <- c(left = 0, center = 0.5, right = 1)[[ha]]; asc <- rep(lp[1], n); d <- rep(lp[2], n)
  if (!is.null(plus)) for (i in seq_len(n)) if (!is.null(plus[[i]])) { asc[i] <- asc[i] + plus[[i]][1] * 25.4 / 300; d[i] <- d[i] + plus[[i]][2] * 25.4 / 300 }
  B <- numeric(n); B[1] <- asc[1]; if (n > 1) for (i in 2:n) B[i] <- B[i - 1] + d[i - 1] + max(lp[1] * ls, asc[i] * ls)
  Ht <- B[n] + d[n]; bv <- asc[n] + (if (n > 1) B[n - 1] + d[n - 1] else 0)                                                   # bv: matplotlib's "baseline" quantity used by center_baseline
  topv <- switch(va, top = 0, center = Ht / 2, bottom = Ht, baseline = B[n], center_baseline = bv / 2, stop("va"))        # top of the text box above the anchor (baseline: the anchor is on the baseline of the last line)
  if (mode == "default" && rot != 0 && !(ha == "center" && va == "center")) stop("rotation_mode 'default' is implemented for centred text only")
  a <- rot * pi / 180; nx <- -sin(a); ny <- -cos(a)
  for (i in seq_len(n)) { v <- topv - B[i]; run_line(x + v * nx, y + v * ny, lines[[i]], size, col, rot, hf) }
  invisible(Ht)
}
seg <- function(x1, y1, x2, y2, col = INK, lwpt = 0.6, ...) segments(x1, y1, x2, y2, col = col, lwd = lw(lwpt), ...)
seg_dashed <- function(x1, y1, x2, y2, col = INK, lwpt = 0.5, on = 3.7, off = 1.6) {      # matplotlib dash pattern: on and off lengths in units of the line width, starting with a dash
  L <- sqrt((x2 - x1)^2 + (y2 - y1)^2); u <- c(x2 - x1, y2 - y1) / L; p <- (on + off) * lwpt * PT; d <- on * lwpt * PT; t <- seq(0, L, by = p)
  for (s in t) segments(x1 + u[1] * s, y1 + u[2] * s, x1 + u[1] * min(s + d, L), y1 + u[2] * min(s + d, L), col = col, lwd = lw(lwpt))
}
pline <- function(x, y, col = INK, lwpt = 1) lines(x, y, col = col, lwd = lw(lwpt), lend = "butt", ljoin = "round")                  # polyline like Line2D (butt caps, round joins)
pline_dashed <- function(x, y, col = INK, lwpt = 1, on = 3, off = 2) {      # matplotlib custom dash pattern (on, off in units of the line width, ls = (0, (on, off))) that runs on continuously along the polyline
  pon <- on * lwpt * PT; poff <- off * lwpt * PT; drawing <- TRUE; remain <- pon; sx <- x[1]; sy <- y[1]
  flush <- function() if (length(sx) > 1) lines(sx, sy, col = col, lwd = lw(lwpt), lend = "butt", ljoin = "round")
  for (i in 2:length(x)) {
    dx <- x[i] - x[i - 1]; dy <- y[i] - y[i - 1]; L <- sqrt(dx^2 + dy^2); done <- 0
    while (L - done > 1e-9) {
      step <- min(remain, L - done); done <- done + step; remain <- remain - step; px <- x[i - 1] + dx * done / L; py <- y[i - 1] + dy * done / L
      if (drawing) { sx <- c(sx, px); sy <- c(sy, py) }
      if (remain <= 1e-9) { if (drawing) flush(); drawing <- !drawing; remain <- if (drawing) pon else poff; sx <- px; sy <- py }
    }
  }
  if (drawing) flush()
}
dashed <- function(x1, y1, x2, y2, col = INK, lwpt = 0.5) seg_dashed(x1, y1, x2, y2, col, lwpt, 3.7, 1.6)
dotted <- function(x1, y1, x2, y2, col = INK, lwpt = 0.5) seg_dashed(x1, y1, x2, y2, col, lwpt, 1, 1.65)
rect_mm <- function(x1, y1, x2, y2, fill = NA, border = NA, lwpt = 0.5) rect(x1, y1, x2, y2, col = fill, border = border, lwd = lw(lwpt))
poly <- function(x, y, fill = NA, border = NA, lwpt = 0.5) polygon(x, y, col = fill, border = border, lwd = lw(lwpt))
## markers: d_pt = marker size in points (matplotlib: ms for plot(), sqrt(s) for scatter()); the edge of width bw_pt is centred on the outline, so the outer size is d + bw
marker <- function(x, y, shape = "o", d_pt = 3, fill = BLUE, border = NA, bw_pt = 0, alpha = NULL) {
  d <- d_pt * PT; n <- length(x); if (!is.null(alpha) && !is.na(fill)) fill <- adjustcolor(fill, alpha.f = alpha)
  for (i in seq_len(n)) {
    if (shape == "o") { th <- seq(0, 2 * pi, length.out = 49); px <- x[i] + d / 2 * cos(th); py <- y[i] + d / 2 * sin(th) }
    else if (shape == "s") { px <- x[i] + d / 2 * c(-1, 1, 1, -1); py <- y[i] + d / 2 * c(-1, -1, 1, 1) }
    else if (shape == "D") { r <- d * sqrt(2) / 2; px <- x[i] + r * c(0, 1, 0, -1); py <- y[i] + r * c(-1, 0, 1, 0) }                # unit square rotated by 45 degrees: diagonal = sqrt(2) * size
    else if (shape == "^") { px <- x[i] + d / 2 * c(0, -1, 1); py <- y[i] + d / 2 * c(-1, 1, 1) }                                  # vertices (0, 0.5), (-0.5, -0.5), (0.5, -0.5) in marker units, y up
    else stop("marker shape")
    polygon(px, py, col = fill, border = if (bw_pt > 0) border else NA, lwd = lw(max(bw_pt, 0.01)))
  }
}

## ---------------------------------------------------------------------------------------------------------------- axes in mm
## matplotlib GridSpec (get_grid_positions): returns cell(i, j, j2 = j) = c(left, top, width, height) in mm from the top-left corner of a W x H mm figure; wr, hr = width and height ratios; edges and spacings as fractions of the figure
gridspec_mm <- function(W, H, nrow, ncol, wr = rep(1, ncol), hr = rep(1, nrow), hspace = 0.2, wspace = 0.2, left = 0.125, right = 0.9, top = 0.88, bottom = 0.11) {
  ch <- hr * (top - bottom) / (nrow + hspace * (nrow - 1)) * nrow / sum(hr); sep_h <- hspace * (top - bottom) / (nrow + hspace * (nrow - 1))
  cw <- wr * (right - left) / (ncol + wspace * (ncol - 1)) * ncol / sum(wr); sep_w <- wspace * (right - left) / (ncol + wspace * (ncol - 1))
  tops <- numeric(nrow); p <- 0; for (i in seq_len(nrow)) { if (i > 1) p <- p + sep_h; tops[i] <- p; p <- p + ch[i] }
  lefts <- numeric(ncol); p <- 0; for (j in seq_len(ncol)) { if (j > 1) p <- p + sep_w; lefts[j] <- p; p <- p + cw[j] }
  function(i, j, j2 = j) c(left = (left + lefts[j]) * W, top = ((1 - top) + tops[i]) * H, width = (lefts[j2] + cw[j2] - lefts[j]) * W, height = ch[i] * H)
}
## matplotlib's default tick locator (AutoLocator = MaxNLocator(nbins = "auto", steps = c(1, 2, 2.5, 5, 10))) and ScalarFormatter, for figures whose ticks are not set by hand (axis = "x" or "y"; length_pt = axes length in points)
mpl_auto_ticks <- function(vmin, vmax, length_pt, size = 7.5, axis = "y") {
  nbins <- min(max(floor(length_pt / (size * if (axis == "y") 2 else 3)), 1), 9); dv <- abs(vmax - vmin); scale <- 10^floor(log10(dv / nbins)); steps <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20) * scale
  raw <- dv / nbins; istep <- if (any(steps >= raw)) which(steps >= raw)[1] else length(steps); ticks <- numeric(0)
  for (step in rev(steps[seq_len(istep)])) {
    best <- floor(vmin / step) * step; x <- vmin - best; d <- floor(x / step); m <- x - d * step; lo <- if (abs(m / step - 1) < 1e-10) d + 1 else d
    x <- vmax - best; d <- floor(x / step); m <- x - d * step; hi <- if (abs(m / step) < 1e-10) d else d + 1; ticks <- (lo:hi) * step + best
    if (sum(ticks <= vmax & ticks >= vmin) >= 2) break
  }
  ticks[ticks >= vmin - 1e-9 * dv & ticks <= vmax + 1e-9 * dv]
}
mpl_tick_labels <- function(ticks) {      # ScalarFormatter: the smallest number of decimals that reproduces all ticks to 1e-3 of the power of ten of their range; true minus sign
  rng <- diff(range(ticks)); if (rng == 0) rng <- max(abs(ticks)); if (rng == 0) rng <- 1; oom <- floor(log10(rng)); sig <- max(0, 3 - oom); thr <- 1e-3 * 10^oom
  while (sig >= 0 && max(abs(ticks - round(ticks, sig))) < thr) sig <- sig - 1
  gsub("-", MINUS, py_fmt_f(ticks, sig + 1), fixed = TRUE)
}
axes_new <- function(left, top, w, h, xlim, ylim, xlog = FALSE, ylog = FALSE) {
  tx <- function(v) if (xlog) log10(v) else v; ty <- function(v) if (ylog) log10(v) else v
  ax <- list(left = left, top = top, w = w, h = h, bottom = top + h, right = left + w, xlim = xlim, ylim = ylim)
  ax$X <- function(v) left + (tx(v) - tx(xlim[1])) / (tx(xlim[2]) - tx(xlim[1])) * w
  ax$Y <- function(v) top + h - (ty(v) - ty(ylim[1])) / (ty(ylim[2]) - ty(ylim[1])) * h
  ax
}
ax_spines <- function(ax, sides = c("left", "bottom"), lwpt = 0.6, col = INK) {
  if ("left" %in% sides) seg(ax$left, ax$top - lwpt * PT / 2, ax$left, ax$bottom + lwpt * PT / 2, col, lwpt); if ("bottom" %in% sides) seg(ax$left - lwpt * PT / 2, ax$bottom, ax$right + lwpt * PT / 2, ax$bottom, col, lwpt)
  if ("top" %in% sides) seg(ax$left, ax$top, ax$right, ax$top, col, lwpt); if ("right" %in% sides) seg(ax$right, ax$top, ax$right, ax$bottom, col, lwpt)
}
TICKLEN <- 2.5; TICKPAD <- 3.5; LABELPAD <- 4.0
ax_xticks <- function(ax, at, labels, size = 6.5, face = 1, col = INK, len = TICKLEN, pad = TICKPAD, y = NULL, ls = 1.2) {        # labels: character vector; \n gives several lines (matplotlib line spacing 1.2); returns the height of the tallest label box (mm), as matplotlib measures it
  y0 <- if (is.null(y)) ax$bottom else y; hmax <- 0
  for (i in seq_along(at)) { x <- ax$X(at[i]); if (len > 0) seg(x, y0, x, y0 + len * PT, col, 0.6); if (!is.null(labels)) { lab <- labels[[i]]; hmax <- max(hmax, mtext(x, y0 + (len + pad) * PT, as.list(strsplit(lab, "\n", fixed = TRUE)[[1]]), size, face, col, "center", "top", ls = ls)) } }
  invisible(hmax)
}
ax_yticks <- function(ax, at, labels, size = 6.5, face = 1, col = INK, len = TICKLEN, pad = TICKPAD, x = NULL, ls = 1.2) {       # tick labels: right-aligned, center_baseline; \n gives several lines; returns the width of the widest line (mm)
  x0 <- if (is.null(x)) ax$left else x; wmax <- 0
  for (i in seq_along(at)) { y <- ax$Y(at[i]); if (len > 0) seg(x0, y, x0 - len * PT, y, col, 0.6)
    if (!is.null(labels)) { lines <- strsplit(labels[[i]], "\n", fixed = TRUE)[[1]]; mtext(x0 - (len + pad) * PT, y, as.list(lines), size, face, col, "right", "center_baseline", ls = ls); wmax <- max(wmax, vapply(lines, wmm, 1, size = size, face = face)) } }
  invisible(wmax)
}
## axis labels: x label below the tick labels (top of the text at labelpad below their bottom), y label left of them (right end of the rotated text at labelpad left of their left edge)
ax_xlabel <- function(ax, s, tick_h, size = 7, adjx = 0.5, runs = NULL, math = !is.null(runs)) {      # x label: top of its text box at labelpad below the tick label boxes (tick_h from ax_xticks); math = the label has mathtext (italic symbols)
  y <- ax$bottom + (TICKLEN + TICKPAD) * PT + tick_h + LABELPAD * PT; xc <- ax$left + adjx * ax$w
  mtext(xc, y, list(if (is.null(runs)) s else runs), size, 1, INK, "center", "top", plus = if (math) list(MATH) else NULL)
}
ax_xlabel_lines <- function(ax, lines, tick_h, size = 7, ls = 1.15, plus = NULL) {      # a multi-line x label: lines = list of run lists (or strings), each centred; plus = list of c(ascent, descent) px added to the lp metrics of a line (mathtext lines)
  y <- ax$bottom + (TICKLEN + TICKPAD) * PT + tick_h + LABELPAD * PT
  mtext(ax$left + ax$w / 2, y, lines, size, 1, INK, "center", "top", ls = ls, plus = plus)
}
ax_ylabel <- function(ax, s, tick_w, size = 7, runs = NULL, math = !is.null(runs)) {      # y label: the right edge of its rotated text box (bottom, rotation_mode anchor) at labelpad left of the tick labels
  x <- ax$left - (TICKLEN + TICKPAD) * PT - tick_w - LABELPAD * PT; yc <- ax$top + ax$h / 2
  mtext(x, yc, list(if (is.null(runs)) s else runs), size, 1, INK, "center", "bottom", rot = 90, plus = if (math) list(MATH) else NULL)
}
ax_title <- function(ax, s = NULL, size = 7, face = 2, pad = 3.0, runs = NULL, adjx = 0.5) {         # baseline at pad points above the axes box
  y <- ax$top - pad * PT; xc <- ax$left + adjx * ax$w
  if (is.null(runs)) txt(xc, y, s, size, face, adj = c(0.5, 0)) else rich(xc, y, runs, size, adj = c(0.5, 0))
}
tag <- function(s, x = 0.8, y = 0.6) txt(x, y, s, 9, 2, INK, adj = c(0, 1))                          # bold panel tag, (x, y) mm from the top-left corner
text_mm <- function(x, y, s, ...) txt(x, y, s, ...)
