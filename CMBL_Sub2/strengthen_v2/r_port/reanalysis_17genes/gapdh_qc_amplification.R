## gapdh_qc_amplification.R - R port of reanalysis_17genes/gapdh_qc_amplification.py: QC of the amplification curves of the 26 August run (all target and Gapdh wells): logistic fit
## base + top / (1 + exp(-(x - mid) / k)) to every curve (scipy curve_fit, Levenberg-Marquardt of MINPACK, start [max, 25, 1.5, 0]; here minpack.lm::nls.lm, the same MINPACK routine and tolerances),
## plateau RFU (mean of cycles 36-40), Cq and sample of the well. Writes qc26_wells.csv and prints the Gapdh wells, the plateau summary and the Gapdh melt peaks.
suppressMessages(library(minpack.lm))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "cfx.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run <- "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Kcnj11-Cacna1d-Cacna1c-Ptbp1"
amp <- cfx_sheet_rows(run, "Quantification Amplification Results"); cq <- cfx_sheet_rows(run, "Quantification Cq Results")[["0"]]; mp <- cfx_sheet_rows(run, "Melt Curve Peak Results")[["0"]]
txt <- function(x) if (is.null(x)) NA_character_ else as.character(x)
cqd <- list(); for (r in cq[-1]) if (!is.null(r[[2]]) && nzchar(txt(r[[2]]))) cqd[[paste(txt(r[[4]]), txt(r[[2]]), sep = "\t")]] <- list(sample = r[[6]], cq = r[[8]])   # (Target, Well) -> (Sample, Cq)
sig <- function(x, top, mid, k, base) base + top / (1 + exp(-(x - mid) / k))
fit_curve <- function(cyc, y) {                         # returns list(p = 4 parameters, r2); NaN when the fit does not converge (scipy: RuntimeError)
  ctl <- nls.lm.control(maxfev = 20000, maxiter = 20000)
  fit <- tryCatch(suppressWarnings(nls.lm(par = c(max(y), 25, 1.5, 0), fn = function(p) y - sig(cyc, p[1], p[2], p[3], p[4]), control = ctl)), error = function(e) NULL)
  if (is.null(fit) || !(fit$info %in% 1:4) || !all(is.finite(fit$par))) return(list(p = rep(NA_real_, 4), r2 = NA_real_))
  p <- fit$par; pred <- sig(cyc, p[1], p[2], p[3], p[4]); list(p = p, r2 = 1 - sum((y - pred)^2) / sum((y - mean(y))^2))
}
rows <- list()
for (gene in names(amp)) {
  if (gene == "Run Information") next
  data <- amp[[gene]]; hdr <- data[[1]][-1]; wells <- vapply(hdr[-1], function(x) if (is.null(x)) NA_character_ else as.character(x), ""); cyc <- vapply(data[-1], function(r) as.numeric(r[[2]]), 1)
  for (jj in seq_along(wells)) {
    w <- wells[jj]; if (is.na(w)) next; col <- jj + 2L
    y <- vapply(data[-1], function(r) if (is.numeric(r[[col]])) r[[col]] else NA_real_, 1); ok <- !is.na(y)
    ft <- fit_curve(cyc[ok], y[ok]); wl <- if (nchar(w) == 3) w else paste0(substr(w, 1, 1), "0", substr(w, 2, nchar(w))); cc <- cqd[[paste(gene, wl, sep = "\t")]]
    rows[[length(rows) + 1]] <- data.frame(gene = gene, well = w, sample = if (is.null(cc) || is.null(cc$sample)) NA_character_ else as.character(cc$sample), cq = if (is.null(cc) || is.null(cc$cq)) NA_real_ else as.numeric(cc$cq),
                                           plateau = mean(tail(y, 5), na.rm = TRUE), top = ft$p[1], mid = ft$p[2], width = ft$p[3], r2 = ft$r2, stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)
g <- d[d$gene == "Gapdh", ]; cat("Gapdh wells (26 Aug): sample, Cq, plateau RFU (cycles 36-40 mean), logistic mid / slope-width, fit R2\n"); print(format(g, digits = 4), row.names = FALSE)
tgt <- d[d$gene != "Gapdh", ]
cat(sprintf("\nplateau RFU, all target wells: median %.0f; Gapdh by sample:", median(tgt$plateau)), "\n")
for (s in py_sort(unique(g$sample))) { v <- g$plateau[g$sample == s]; cat(sprintf("  '%s': mean %.0f min %.0f max %.0f\n", s, mean(v), min(v), max(v))) }
cat("\nGapdh melt peaks:\n"); for (r in mp[-1]) if (identical(txt(r[[4]]), "Gapdh")) cat("  ", txt(r[[2]]), txt(r[[6]]), "Tm", txt(r[[7]]), "height", txt(r[[8]]), "\n")
py_write_csv_repr(d, file.path(OUT, "qc26_wells.csv"))
