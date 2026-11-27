## make_fig10.R - R port of fig10_201126/make_fig10.py: Fig. 10 of the JTM manuscript, the candidate checks for pooled co-expression networks from small experiments, with the evidence behind each step.
## A flow of six steps; for each: the check, what the benchmark and the INS-1 data show (every number carries its status: simulation, benchmark, post hoc, INS-1) and how to read it. Plots no new analysis. Every number is
## read from a saved result file of the study (calibration_tables.csv, calibration_sensitivity_no_degenerate.csv, design_table.csv, validation_endpoints_long.csv, post_hoc/results_*/summary.json) or is an INS-1 number of the
## manuscript text; fig10_numbers.csv lists each number with its source, and the script checks that the INS-1 numbers occur in the manuscript text. Style of Figs 7-9 (figstyle.R): 170 mm wide, Arial >= 6 pt, PNG 300 dpi and PDF,
## no key strip, statistical parameters (P, q, r, D, T) in italics.
## Usage: Rscript make_fig10.R validation   (final figure; needs post_hoc/results_validation)      Rscript make_fig10.R dev   (layout test with development-series numbers; the file name carries _TEST and must not be used)
## Environment: WP2_RES (benchmark result folder, default the archived wp2_benchmark2), POSTHOC_RES (folder that holds results_validation, results_dev and results_ins1, default the archived ars_review_191126/post_hoc),
## FIG10_MANUSCRIPT (docx whose text must contain the INS-1 numbers), FIG_OUT (default results_R/figures).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "figstyle.R")); suppressMessages(library(jsonlite))
SOURCE <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "validation"; stopifnot(SOURCE %in% c("validation", "dev"))
V2 <- normalizePath(file.path(R_PORT, "..")); W2 <- Sys.getenv("WP2_RES", file.path(V2, "wp2_benchmark2")); PH <- Sys.getenv("POSTHOC_RES", file.path(V2, "ars_review_191126", "post_hoc")); OUT <- Sys.getenv("FIG_OUT", file.path(R_PORT, "results_R", "figures"))
MS_PATH <- Sys.getenv("FIG10_MANUSCRIPT", normalizePath(file.path(V2, "..", "..", "Manuscript", "Manuscript_JTM_v2_SUBMISSION_191126c.docx"), mustWork = FALSE))
rd <- function(p) read.csv(file.path(W2, p), stringsAsFactors = FALSE)

## ------------------------------------------------------------------ numbers and their sources
trace <- list(); put <- function(key, value, source) { trace[[length(trace) + 1]] <<- c(key, value, source); value }
cal <- rd("calibration_tables.csv"); cal_row <- function(r, m, stat) cal[cal$r == r & cal$m == m & cal$statistic == stat, ][1, ]
t1 <- cal_row(0.0, 0, "T_r2"); p_tr2 <- cal_row(0.7, 6, "T_r2"); c07 <- cal_row(0.7, 6, "count"); c09 <- cal_row(0.9, 6, "count"); c05 <- cal_row(0.5, 6, "count")
pct <- function(x) sprintf("%d%%", as.integer(round(100 * x))); f1 <- function(x) py_fmt_f(x, 1); q2 <- function(v) py_fmt_f(v, 2); kk <- function(a, b) sprintf("%d of %d", as.integer(a), as.integer(b))
N <- list()
N$typeI <- put("T_r2 type I error", paste0(f1(100 * t1$rate), "%"), paste0("calibration_tables.csv T_r2 r=0 m=0: ", kk(t1$n_reject, t1$n_tests)))
nd <- rd("calibration_sensitivity_no_degenerate.csv"); nd <- nd[nd$statistic == "T_r2", ][1, ]
N$typeI_nd <- put("T_r2 type I error without 2 series with constant genes (post hoc)", paste0(f1(100 * nd$type1), "%"), paste0("calibration_sensitivity_no_degenerate.csv: ", kk(nd$type1_k, nd$type1_n)))
N$power_tr2 <- put("T_r2 power r 0.7, 6 genes", pct(p_tr2$rate), paste0("calibration_tables.csv: ", kk(p_tr2$n_reject, p_tr2$n_tests)))
N$count07 <- put("count test power r 0.7", pct(c07$rate), paste0("calibration_tables.csv: ", kk(c07$n_reject, c07$n_tests)))
N$count09 <- put("count test power r 0.9", pct(c09$rate), paste0("calibration_tables.csv: ", kk(c09$n_reject, c09$n_tests)))
N$count05 <- put("count test power r 0.5", pct(c05$rate), paste0("calibration_tables.csv: ", kk(c05$n_reject, c05$n_tests)))
dt <- rd("design_table.csv"); d4 <- dt[dt$n_per_group == 4 & dt$r == 0.7, ][1, ]
N$d4_power <- put("design 4 per group: T_r2 power r 0.7", pct(d4$power_T_r2), "design_table.csv n=4 r=0.7")
N$d4_pairs <- put("design 4 per group: share of true module pairs detected", pct(d4$pair_power), "design_table.csv n=4 r=0.7")
N$d4_rbest <- put("design 4 per group: |r| the best single pair needs", q2(d4$r_best_pair), "design_table.csv n=4")
d6 <- dt[dt$n_per_group == 6 & dt$r == 0.5, ][1, ]; d8 <- dt[dt$n_per_group == 8 & dt$r == 0.5, ][1, ]
N$d6_8 <- put("design 6 and 8 per group: T_r2 power r 0.5", paste0(q2(d6$power_T_r2), " to ", q2(d8$power_T_r2)), "design_table.csv n=6, n=8 r=0.5")
vl <- rd("validation_endpoints_long.csv"); va <- vl[startsWith(vl$arm, "All validation sets pooled") & vl$endpoint == "coupling_T_r2" & vl$gene_set == "A", ][1, ]
N$flag_k <- c(as.integer(va$k), as.integer(va$n)); N$flag_pct <- put("T_r2 flagged, validation sets, gene set A", pct(va$estimate), paste0("validation_endpoints_long.csv: ", kk(va$k, va$n)))
ph_path <- file.path(PH, if (SOURCE == "validation") "results_validation" else "results_dev", "summary.json"); ph_path <- if (.Platform$OS.type == "windows") gsub("/", "\\", normalizePath(ph_path, mustWork = FALSE), fixed = TRUE) else ph_path
S <- fromJSON(ph_path, simplifyVector = FALSE); SA <- S$sets$A; SB <- S$sets$B; SAB <- S$sets[["A+B"]]
N$A_obs <- put("share carried by group means, gene set A (median of series)", q2(SA$obs_share_quartiles[[2]]), paste(ph_path, "sets.A.obs_share_quartiles"))
N$A_exp <- put("expected share without coupling, gene set A", q2(SA$null_share_quartiles[[2]]), paste(ph_path, "sets.A.null_share_quartiles"))
N$B_obs <- put("share carried by group means, gene set B", q2(SB$obs_share_quartiles[[2]]), paste(ph_path, "sets.B.obs_share_quartiles"))
N$B_exp <- put("expected share without coupling, gene set B", q2(SB$null_share_quartiles[[2]]), paste(ph_path, "sets.B.null_share_quartiles"))
N$above <- put("matrices above the 97.5th null percentile", kk(SAB$n_above_null_q975, SAB$n_share_defined), paste(ph_path, "sets.A+B"))
N$pc1_k <- put("T_r2 flagged after removal of PC1, gene set A", pct(SA$flagged_pc1 / SA$n_matrices), paste0(ph_path, " sets.A: ", kk(SA$flagged_pc1, SA$n_matrices)))
N$orig_k <- put("T_r2 flagged, original, same sets (gene set A)", pct(SA$flagged_original / SA$n_matrices), paste0(ph_path, " sets.A: ", kk(SA$flagged_original, SA$n_matrices)))
N$n_sets <- SA$n_matrices
if (SOURCE == "validation") stopifnot("post hoc validation results do not reproduce the benchmark prevalence" = (N$orig_k == N$flag_pct && SA$flagged_original == N$flag_k[1] && SA$n_matrices == N$flag_k[2] && N$flag_k[2] == 59))      # gate: the original flagged fraction recomputed by the post hoc script equals the benchmark's
I <- fromJSON(file.path(PH, "results_ins1", "summary.json"), simplifyVector = FALSE)
N$I_obs <- put("INS-1 share carried by group means", q2(I$share$obs_median), "results_ins1/summary.json"); N$I_exp <- put("INS-1 expected share", q2(I$share$null_mean), "results_ins1/summary.json")
## INS-1 numbers of the manuscript text (checked against the text below)
MS_NUM <- c(`null mean 22.1` = "22.1", `P = 0.68` = "0.68", `Tr2 P = 0.061` = "0.061", `same-run r` = "+0.149", `between-run r` = "−0.019", D = "+0.168", `D P` = "0.020", `row-layout simulation` = "+0.520", `observed same-run, 26 August` = "+0.536")
for (k in names(MS_NUM)) trace[[length(trace) + 1]] <- c(paste("INS-1", k), MS_NUM[[k]], "manuscript Results, paragraphs on the permutation null, the coupling test and the run structure")
html_unescape <- function(x) { if (length(x) == 0) return(x); x <- gsub("&lt;", "<", x, fixed = TRUE); x <- gsub("&gt;", ">", x, fixed = TRUE); x <- gsub("&quot;", "\"", x, fixed = TRUE); x <- gsub("&apos;", "'", x, fixed = TRUE)
  m <- gregexpr("&#x[0-9A-Fa-f]+;|&#[0-9]+;", x); for (e in unique(unlist(regmatches(x, m)))) { v <- if (startsWith(e, "&#x")) strtoi(substr(e, 4, nchar(e) - 1), 16L) else strtoi(substr(e, 3, nchar(e) - 1), 10L); x <- gsub(e, intToUtf8(v), x, fixed = TRUE) }
  gsub("&amp;", "&", x, fixed = TRUE) }
manuscript_text <- function(path) {
  tmp <- tempfile(); dir.create(tmp); unzip(path, files = "word/document.xml", exdir = tmp); f <- file.path(tmp, "word", "document.xml"); x <- readChar(f, file.info(f)$size, useBytes = TRUE); Encoding(x) <- "UTF-8"
  body <- substr(x, regexpr("<w:body>", x, fixed = TRUE)[1], nchar(x)); paras <- strsplit(body, "</w:p>", fixed = TRUE)[[1]]                                   # paragraphs; the text of a paragraph is the concatenation of its <w:t> runs
  paste(vapply(paras, function(p) { t <- regmatches(p, gregexpr("<w:t(?: [^>]*)?>[^<]*</w:t>", p, perl = TRUE))[[1]]; paste(html_unescape(sub("</w:t>$", "", sub("^<w:t(?: [^>]*)?>", "", t, perl = TRUE))), collapse = "") }, ""), collapse = "\n")
}
MST <- manuscript_text(MS_PATH); missing <- MS_NUM[!vapply(MS_NUM, function(v) grepl(v, MST, fixed = TRUE), TRUE)]
if (length(missing) > 0) stop("numbers not found in the manuscript text: ", paste(missing, collapse = ", "))

## ------------------------------------------------------------------ text of the figure (markup: <i> italic, <sub> subscript; the pieces are set as matplotlib's mathtext does: italic and upright in the math font, not bold)
TR2 <- "<i>T</i><sub><i>r</i>2</sub>"; rr <- "<i>r</i>"; PP <- "<i>P</i>"; QQ <- "<i>q</i>"; DD <- "<i>D</i>"
ROWS <- list(
  list("Design and genes",
       sprintf("With 4 replicates per group, %s detects one module of 6 genes at %s = 0.7 with power %s, but only %s of the module’s pairs reach %s < 0.05, and the best single pair would need |%s| ≥ %s. Six to eight independent experiments per group give power %s at %s = 0.5 (simulation).", TR2, rr, N$d4_power, N$d4_pairs, QQ, rr, N$d4_rbest, N$d6_8, rr),
       "Wells from one preparation are not independent experiments. The design limits what any pair-level result can show."),
  list("Permutation null of the consensus count",
       sprintf("The count detects injected within-group coupling in %s of tests (%s = 0.7, 6 genes; %s at %s = 0.9, %s at %s = 0.5) (benchmark, simulation). INS-1: 21 pairs, null mean 22.1, %s = 0.68.", N$count07, rr, N$count09, rr, N$count05, rr, PP),
       "A count inside its null does not show that coupling is absent."),
  list("Share carried by group means, against its expectation",
       sprintf("Expected without within-group coupling (post hoc, %d public data sets): genes selected for a group effect %s observed, %s expected; random expressed genes %s observed, %s expected. INS-1: %s observed, %s expected. The observed share was above its null in %s matrices.", as.integer(N$n_sets), N$A_obs, N$A_exp, N$B_obs, N$B_exp, N$I_obs, N$I_exp, N$above),
       "A share near 1 reflects the selected genes and the three groups of four, not co-regulation."),
  list(sprintf("Within-group coupling test, %s", TR2),
       sprintf("Power %s (%s = 0.7, 6 genes). Type I error %s, above the limit of 6%% set in the protocol (%s after a post hoc exclusion of 2 series with constant genes). Flagged in %s of 59 validation sets (gene set A); %s after removal of one shared axis (post hoc). INS-1: %s = 0.061.", N$power_tr2, rr, N$typeI, N$typeI_nd, N$flag_pct, N$pc1_k, PP),
       "A flag shows dependence within groups, not its source: shared sample-level variation or shared biology."),
  list("Treatment adjustment and ratio schemes",
       "After adjustment for treatment no pair was significant in all four schemes, in INS-1 or (median) in the public sets. Counts of the ratio schemes sat at their own null (median observed to expected 1.05 to 1.19).",
       "Ratio schemes share a denominator: their counts are not findings."),
  list("Run and reference structure",
       sprintf("INS-1: within-group %s was +0.149 among same-run pairs and −0.019 between runs (%s = +0.168, %s = 0.020); the row layout of one run alone gave +0.520 (observed +0.536).", rr, DD, PP),
       "Compatible with a shared reference, plate or layout, not with co-regulation."))
NOTE <- "Candidate checks. Whether they improve the conclusions drawn from small experiments, and whether they apply to toxicogenomic data, was not tested."

## ------------------------------------------------------------------ layout in millimeters
FS <- 6.5; LH <- 2.95; PAD <- 1.6; GAP <- 2.0; CX <- 9.5; CW <- 38.0; EX <- 49.5; EW <- 72.0; RX <- 123.5; RW <- 46.0; NBSP <- " "
mk_runs <- function(s, face = 1, size = FS) {      # markup to runs (s, face, size factor, drop in pt); subscripts at 86% size (SHRINK_FACTOR of the originals)
  toks <- regmatches(s, gregexpr("</?(i|sub)>|[^<]+|<", s, perl = TRUE))[[1]]; it <- 0; sb <- 0; runs <- list()
  for (t in toks) {
    if (t == "<i>") it <- it + 1 else if (t == "</i>") it <- it - 1 else if (t == "<sub>") sb <- sb + 1 else if (t == "</sub>") sb <- sb - 1
    else { f <- if (it > 0) 3 else if (sb > 0 || it > 0) 1 else face; r <- if (sb > 0) list(t, f, 0.86, 0.214 * size) else list(t, f); n <- length(runs)
      if (n > 0 && identical(runs[[n]][-1], r[-1])) runs[[n]][[1]] <- paste0(runs[[n]][[1]], t) else runs[[n + 1]] <- r }
  }
  runs
}
line_plus <- function(s) if (grepl("<sub>", s, fixed = TRUE)) c(2, 5.1) else if (grepl("<i>", s, fixed = TRUE)) c(2, 2) else NULL       # mathtext lines are 2 px higher and deeper than plain ones, with a subscript 5.1 px deeper (6.5 pt)
nb <- function(s) gsub("(\\S) ([=<≥]) (?=[−+-]?\\d)", paste0("\\1", NBSP, "\\2", NBSP), s, perl = TRUE)      # no break inside 'x = value' / 'x < value'
.FIG$mw <- function(s, face = 1) sum(vapply(mk_runs(s, face), function(r) wmm(r[[1]], FS * (if (length(r) >= 3) r[[3]] else 1), r[[2]]), 1))
wrap <- function(s, width, face = 1) {
  lines <- character(0); cur <- ""
  for (w in strsplit(nb(s), " ", fixed = TRUE)[[1]]) {
    trial <- if (nzchar(cur)) paste(cur, w) else w
    if (.FIG$mw(trial, face) <= width - 2 * 1.4 || !nzchar(cur)) cur <- trial else { lines <- c(lines, cur); cur <- w }
  }
  c(lines, cur)
}
agg_png(tempfile(fileext = ".png"), width = FIGW, height = 50, units = "mm", res = 300, pointsize = PS, background = "white"); canvas(50)          # probe device: text widths in mm
rows <- list(); y <- 9.0
for (i in seq_along(ROWS)) { tl <- wrap(ROWS[[i]][[1]], CW, 2); el <- wrap(ROWS[[i]][[2]], EW); rl <- wrap(ROWS[[i]][[3]], RW); n <- max(length(tl), length(el), length(rl)); h <- n * LH + 2 * PAD
  rows[[i]] <- list(i = i, y = y, h = h, tl = tl, el = el, rl = rl); y <- y + h + GAP }
nl <- wrap(NOTE, FIGW - 2 * 1.0); nh <- length(nl) * LH + 2 * PAD; H <- y + 1.0 + nh + 1.5; invisible(dev.off())

roundrect <- function(x, y, w, h, r, fill, border, lwpt) {      # matplotlib FancyBboxPatch 'round' (pad 0): straight sides, corners are quadratic Bezier curves with the control point at the corner
  x0 <- x; y0 <- y; x1 <- x + w; y1 <- y + h; t <- seq(0, 1, length.out = 13)
  q <- function(p0, c1, p2) cbind((1 - t)^2 * p0[1] + 2 * t * (1 - t) * c1[1] + t^2 * p2[1], (1 - t)^2 * p0[2] + 2 * t * (1 - t) * c1[2] + t^2 * p2[2])
  P <- rbind(q(c(x1 - r, y0), c(x1, y0), c(x1, y0 + r)), q(c(x1, y1 - r), c(x1, y1), c(x1 - r, y1)), q(c(x0 + r, y1), c(x0, y1), c(x0, y1 - r)), q(c(x0, y0 + r), c(x0, y0), c(x0 + r, y0)))
  polygon(P[, 1], P[, 2], col = fill, border = border, lwd = lw(lwpt))
}
put_line <- function(x, y, s, size = FS, face = 1, col = INK) mtext(x, y, list(mk_runs(s, face, size)), size, face, col, "left", "center", plus = if (is.null(line_plus(s))) NULL else list(line_plus(s)))
draw <- function() {
  for (cc in list(list(CX, "Check"), list(EX, "What the benchmark and the INS-1 data show"), list(RX, "How to read it"))) txt(cc[[1]] + 1.4, 4.8, cc[[2]], 7, 2, INK, adj = c(0, 0.5), va = "center")
  prev <- NA; th <- seq(0, 2 * pi, length.out = 73)
  for (r in rows) {
    for (cc in list(list(CX, CW, r$tl, "#e8f1fc", 2), list(EX, EW, r$el, "white", 1), list(RX, RW, r$rl, "#f4f3f0", 1))) {
      roundrect(cc[[1]], r$y, cc[[2]], r$h, 1.2, cc[[4]], GREYLINE, 0.5)
      for (k in seq_along(cc[[3]])) put_line(cc[[1]] + 1.4, r$y + PAD + (k - 0.5) * LH, cc[[3]][k], FS, cc[[5]])
    }
    cy <- r$y + r$h / 2; polygon(4.6 + 2.6 * cos(th), cy + 2.6 * sin(th), col = BLUE, border = NA); mtext(4.6, cy + 0.1, as.character(r$i), 7.5, 2, "white", "center", "center")
    if (!is.na(prev)) seg(4.6, prev + 2.6, 4.6, cy - 2.6, GREYLINE, 0.8)
    prev <- cy
  }
  ny <- H - 1.0 - nh; roundrect(1.0, ny, FIGW - 2.0, nh, 1.2, "#f4f3f0", GREYLINE, 0.5)
  for (k in seq_along(nl)) put_line(2.4, ny + PAD + (k - 0.5) * LH, nl[k], FS, 1)
}
stem <- paste0("Figure_10_candidate_checks", if (SOURCE == "validation") "" else "_TEST")
fig_save(stem, H, draw, OUT)
csv_field <- function(s) if (grepl("[,\"\r\n]", s)) paste0("\"", gsub("\"", "\"\"", s, fixed = TRUE), "\"") else s
txt_out <- paste0(vapply(c(list(c("number", "value", "source")), trace), function(r) paste(vapply(r, csv_field, ""), collapse = ","), ""), "\r\n", collapse = "")
nf <- file.path(OUT, sprintf("fig10_numbers%s.csv", if (SOURCE == "validation") "" else "_TEST")); con <- file(nf, "wb"); writeBin(charToRaw(enc2utf8(txt_out)), con); close(con)
cat(sprintf("figure height %.1f mm, %d rows; numbers traced: %d\n", H, length(rows), length(trace)))
