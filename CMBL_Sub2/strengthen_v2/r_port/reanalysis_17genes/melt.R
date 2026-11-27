## melt.R - R port of reanalysis_17genes/melt.py: all wells of the five CFX runs (Cq, sample, melt peaks) as cfx_long.json, and per target gene the number of wells, wells with Cq >= 35 or no Cq,
## wells without or with several melt peaks and the mean and SD of the first-peak melt temperature as cfx_summary.json. Needs cfx.R. Prints the same table.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "cfx.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
truthy <- function(x) !is.null(x) && !(is.character(x) && !nzchar(x)) && !(is.numeric(x) && x == 0)
header_ix <- function(h) { k <- vapply(h, function(x) if (is.null(x)) "" else as.character(x), ""); setNames(seq_along(k)[nzchar(k)], k[nzchar(k)]) }
GENE_MAP <- c(GAPDH = "Gapdh", `GLUT-2` = "Glut2", `Ins-1` = "Ins1", ITPR1 = "Itpr1", RYR1 = "Ryr1", RYR2 = "Ryr2", RYR3 = "Ryr3")
rows <- list()
for (run in cfx_runs()) {
  cq <- list(); mt <- list()
  sh <- cfx_sheet_rows(run, "Quantification Cq Results")
  for (t in names(sh)) {
    if (startsWith(t, "Run")) next
    r <- sh[[t]]; ix <- header_ix(r[[1]])
    for (x in r[-1]) { w <- x[[ix[["Well"]]]]; if (truthy(w)) cq[[w]] <- list(target = x[[ix[["Target"]]]], sample = x[[ix[["Sample"]]]], cq = x[[ix[["Cq"]]]], content = x[[ix[["Content"]]]]) }
  }
  sh <- cfx_sheet_rows(run, "Melt Curve Peak Results")
  for (t in names(sh)) {
    if (startsWith(t, "Run")) next
    r <- sh[[t]]; ix <- header_ix(r[[1]])
    for (x in r[-1]) {
      w <- x[[ix[["Well"]]]]; if (!truthy(w)) next
      tm <- x[[ix[["Melt Temperature"]]]]
      if (is.numeric(tm)) mt[[w]] <- c(mt[[w]], list(list(tm, x[[ix[["Peak Height"]]]]))) else if (is.null(mt[[w]])) mt[[w]] <- list()
    }
  }
  for (w in names(cq)) {
    v <- cq[[w]]; g <- trimws(sub("^[0-9]+\\.\\s*", "", if (is.null(v$target)) "None" else as.character(v$target)))
    suf <- regmatches(g, regexec("_(Ctrl|DAPT|DKK1)$", g))[[1]]; g <- sub("_(Ctrl|DAPT|DKK1)$", "", g)
    sm <- if (truthy(v$sample)) v$sample else if (length(suf)) suf[2] else v$sample
    if (g %in% names(GENE_MAP)) g <- GENE_MAP[[g]]
    rows[[length(rows) + 1]] <- list(run = run, well = w, gene = g, sample = sm, cq = v$cq, content = v$content, melt = if (is.null(mt[[w]])) list() else mt[[w]])
  }
}
write_text_bytes(py_json_typed(rows, indent = NULL), file.path(OUT, "cfx_long.json"))
genes <- vapply(rows, function(r) r$gene, "")
cat(length(rows), "wells;", length(unique(genes)), "targets:", paste0("[", paste(sprintf("'%s'", py_sort(unique(genes))), collapse = ", "), "]"), "\n")
cat(sprintf("%-10s %5s %9s %6s %5s %8s %6s  %s\n", "gene", "wells", "Cq>=35/NA", "nomelt", "multi", "Tm mean", "Tm SD", "runs"))
nanf <- function(x, w, d) if (is.null(x) || is.na(x) || x == 0) formatC("nan", width = w) else formatC(x, format = "f", digits = d, width = w)      # Python: x or float('nan')
summ <- list()
for (g in py_sort(unique(genes))) {
  rs <- rows[genes == g]; n <- length(rs)
  late <- sum(vapply(rs, function(r) !is.numeric(r$cq) || r$cq >= 35, TRUE)); nomelt <- sum(vapply(rs, function(r) length(r$melt) == 0, TRUE)); multi <- sum(vapply(rs, function(r) length(r$melt) > 1, TRUE))
  tms <- vapply(rs[vapply(rs, function(r) length(r$melt) > 0, TRUE)], function(r) r$melt[[1]][[1]], 1)
  summ[[g]] <- list(n = n, late = late, nomelt = nomelt, multi = multi, tm = if (length(tms)) mean(tms) else NULL, sd = if (length(tms) > 1) np_std(tms) else NULL,
                    runs = as.list(py_sort(unique(substr(vapply(rs, function(r) r$run, ""), 1, 16)))))
  cat(sprintf("%-10s %5d %9d %6d %5d %s %s  %s\n", g, n, late, nomelt, multi, nanf(summ[[g]]$tm, 8, 2), nanf(summ[[g]]$sd, 6, 2), paste(unlist(summ[[g]]$runs), collapse = ",")))
}
write_text_bytes(py_json_typed(summ, indent = NULL), file.path(OUT, "cfx_summary.json"))
