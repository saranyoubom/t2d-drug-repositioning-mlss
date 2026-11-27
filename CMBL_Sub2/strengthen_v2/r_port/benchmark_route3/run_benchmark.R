## run_benchmark.R - R port of benchmark_route3/run_benchmark.py. Run the pre-specified benchmark (PROTOCOL sections 2-5) on the frozen dataset list. One JSON line per processed series in results_per_dataset.jsonl;
## every series, analysable or not, is logged with its status. Series are processed in the frozen order; a series that is not analysable (fewer than 34 DE genes, unreadable file, fewer than 12 usable samples,
## platform not annotated) is logged and the next eligible series is used, until 20 series per platform are analysable or the list is exhausted. Downloaded files are untrusted data (kept in <data dir>, read only as text/numbers).
## Usage: Rscript run_benchmark.R <data dir> [max_per_platform=20] [dry]    Environment: RB_OUT (output folder, default results_R/benchmark_route3), RB_FROZEN (frozen_datasets.csv, default the archived list),
## RB_ONLY (comma-separated series to run, for checks on a few series; the duplicate rule and the per-platform counters then see only these series).
## 'dry' only decides analysability (DE gene count and file integrity, no correlation is computed) and writes dry_run_per_dataset.jsonl.
suppressMessages({ library(data.table); library(curl); library(jsonlite) })
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "bench_pipeline.R"))
args <- commandArgs(trailingOnly = TRUE)
OUT <- Sys.getenv("RB_OUT", file.path(R_PORT, "results_R", "benchmark_route3")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DATA <- if (length(args) >= 1) args[1] else file.path(OUT, "data"); MAXP <- if (length(args) >= 2) as.integer(args[2]) else 20L; DRY <- length(args) >= 3 && args[3] == "dry"
ONLY <- setdiff(strsplit(Sys.getenv("RB_ONLY", ""), ",", fixed = TRUE)[[1]], ""); FROZEN <- Sys.getenv("RB_FROZEN", normalizePath(file.path(R_PORT, "..", "..", "benchmark_route3", "frozen_datasets.csv")))
for (d in c("matrices", "platforms")) dir.create(file.path(DATA, d), showWarnings = FALSE, recursive = TRUE)
FTP <- "https://ftp.ncbi.nlm.nih.gov/geo/"; NPERM <- 2000L; MIN_DE <- 34L; NGENE <- 17L
## RB_ARCHIVED (optional, added 25 Nov 2026, restricted 26 Nov 2026 after the independent audit, audit_261126 C5a): path of an archived results jsonl, used only for the series listed in RB_ARCHIVED_ONLY
## (comma-separated GSE ids) and only when the download of that series' own matrix fails (GSE31851 answered HTTP 403 on 25 Nov 2026). Such a series then takes its archived line verbatim, so the frozen loop
## and its per-platform counters stay as in the archive. Every such use is listed in <OUT>/archived_records_used.txt, flagged in the console line, and printed by compare_route3.R as a MISSING line (not
## computed by R), which the agreement report never explains. A failed platform-annotation download never uses the archive. Without both variables a failed series is logged not analysable.
ARCHIVED <- Sys.getenv("RB_ARCHIVED", ""); ARCHIVED_ONLY <- setdiff(trimws(strsplit(Sys.getenv("RB_ARCHIVED_ONLY", ""), ",", fixed = TRUE)[[1]]), "")
append_line <- function(txt, path) { con <- file(path, "ab"); writeBin(charToRaw(enc2utf8(paste0(txt, .NL))), con); close(con) }   # one line ended as Python's text-mode file ends it (os.linesep: CRLF on Windows, as in the archived jsonl). The 25 Nov version wrote LF on the wrong premise that the archive is LF (audit C5a)
py_json_ascii <- function(s) {   # Python's json.dumps (ensure_ascii, the default in run_benchmark.py) writes every non-ASCII character as a \uXXXX escape (lowercase hex; surrogate pairs above U+FFFF); added 26 Nov 2026 (audit C3: the beta of GSE29881)
  cp <- utf8ToInt(enc2utf8(s)); if (!length(cp) || all(cp < 128L)) return(s)
  paste(vapply(cp, function(c) if (c < 128L) intToUtf8(c) else if (c < 65536L) sprintf("\\u%04x", c) else { c <- c - 65536L; sprintf("\\u%04x\\u%04x", 55296L + c %/% 1024L, 56320L + c %% 1024L) }, ""), collapse = "")
}
archived_record <- function(gse, e) {
  if (!nzchar(ARCHIVED) || !(gse %in% ARCHIVED_ONLY) || !grepl(sprintf("download failed .*/%s/matrix/", gse), conditionMessage(e))) return(NULL)   # only a listed series, only its own matrix download
  for (line in readLines(ARCHIVED, encoding = "UTF-8", warn = FALSE)) {
    if (!nzchar(trimws(line))) next; j <- fromJSON(line, simplifyVector = FALSE)
    if (identical(j$gse, gse)) {
      append_line(sprintf("%s: %s; archived record used (status %s)", gse, conditionMessage(e), j$status), file.path(OUT, "archived_records_used.txt"))
      attr(j, "raw") <- line; return(j)
    }
  }
  NULL
}
PANDAS_NA <- c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null")   # pandas default na_values

http_status <- function(url, path = NULL) { h <- new_handle(useragent = "geo-benchmark/1.0", timeout = 180); tryCatch(if (is.null(path)) curl_fetch_memory(url, h) else curl_fetch_disk(url, path, h), error = function(e) NULL) }
fetch <- function(url, path) {
  if (file.exists(path) && file.size(path) > 0) return(path)
  for (k in 0:3) {
    r <- http_status(url, paste0(path, ".part"))
    if (!is.null(r) && r$status_code == 200) { file.rename(paste0(path, ".part"), path); return(path) }
    if (!is.null(r) && r$status_code == 404) { unlink(paste0(path, ".part")); stop(structure(class = c("http404", "error", "condition"), list(message = paste("HTTP 404", url), call = NULL))) }
    Sys.sleep(3 + 3 * k)
  }
  stop("download failed ", url)
}
series_path <- function(gse, gpl) {
  n <- as.integer(substring(gse, 4)); folder <- if (n >= 1000) sprintf("GSE%dnnn", n %/% 1000) else "GSEnnn"; base <- sprintf("%sseries/%s/%s/matrix/", FTP, folder, gse); path <- file.path(DATA, "matrices", sprintf("%s_%s.txt.gz", gse, gpl))
  tryCatch(fetch(paste0(base, gse, "_series_matrix.txt.gz"), path), http404 = function(e) {
    idx <- rawToChar(curl_fetch_memory(base, new_handle(timeout = 60))$content); hit <- regmatches(idx, gregexpr("href=\"[^\"]*series_matrix\\.txt\\.gz\"", idx))[[1]]; hit <- sub("^href=\"", "", sub("\"$", "", hit)); hit <- hit[grepl(gpl, hit, fixed = TRUE)]
    fetch(paste0(base, hit[1]), path)
  })
}
read_table_block <- function(path, begin, end) {                    # the text table between two marker lines of a GEO file
  lines <- readLines(gzfile(path, "rt"), encoding = "UTF-8", warn = FALSE); b <- which(startsWith(lines, begin))[1]; e <- which(startsWith(lines, end))[1]
  fread(text = lines[(b + 1):(e - 1)], sep = "\t", quote = "\"", header = TRUE, colClasses = "character", na.strings = PANDAS_NA, data.table = FALSE, fill = TRUE)
}
read_matrix <- function(path) read_table_block(path, "!series_matrix_table_begin", "!series_matrix_table_end")
annotation <- function(gpl) {
  n <- as.integer(substring(gpl, 4)); folder <- if (n >= 1000) sprintf("GPL%dnnn", n %/% 1000) else "GPLnnn"
  path <- fetch(sprintf("%splatforms/%s/%s/annot/%s.annot.gz", FTP, folder, gpl, gpl), file.path(DATA, "platforms", paste0(gpl, ".annot.gz")))
  t <- read_table_block(path, "!platform_table_begin", "!platform_table_end"); col <- names(t)[grepl("symbol", tolower(names(t)), fixed = TRUE)][1]
  out <- data.frame(probe = t$ID, symbol = t[[col]], stringsAsFactors = FALSE); out <- out[!is.na(out$probe) & !is.na(out$symbol), ]
  out$symbol <- trimws(sub(" /// .*$", "", out$symbol)); out <- out[out$symbol != "" & out$symbol != "---", ]; setNames(out$symbol, out$probe)
}
ANNOT <- list()
gene_matrix <- function(gse, gpl, gsms) {                            # 12 samples x genes on the log2 scale (probe with the highest mean linear expression per gene) and the scale decision
  M <- read_matrix(series_path(gse, gpl)); miss <- setdiff(gsms, names(M)); if (length(miss)) stop(structure(class = c("ValueError", "error", "condition"), list(message = sprintf("samples missing: [%s]", paste(sprintf("'%s'", head(miss, 3)), collapse = ", ")), call = NULL)))
  X <- vapply(gsms, function(g) suppressWarnings(as.numeric(M[[g]])), numeric(nrow(M))); rownames(X) <- M[[1]]
  q99 <- unname(quantile(X[!is.na(X)], 0.99, type = 7)); is_log <- q99 < 30; lin <- if (is_log) 2^X else X
  ok <- rowSums(!is.na(lin) & lin > 0) == ncol(lin); lin <- lin[ok, , drop = FALSE]
  if (is.null(ANNOT[[gpl]])) ANNOT[[gpl]] <<- annotation(gpl)
  sym <- unname(ANNOT[[gpl]][rownames(lin)]); keep <- !is.na(sym); lin <- lin[keep, , drop = FALSE]; sym <- sym[keep]
  mn <- rowMeans(lin); o <- order(sym, -mn, seq_along(sym), method = "radix"); best <- o[!duplicated(sym[o])]       # one probe per gene symbol: the highest mean, first in file order on ties
  list(L2 = log2(t(lin[best, , drop = FALSE])), genes = sym[best], is_log = is_log, q99 = q99)
}
## One-way ANOVA across the 3 groups for every gene (columns of L), BH FDR. The sums follow numpy's order on purpose: a probe that is constant in all 12 samples (the detection floor, thousands of probes on
## GPL1261) has ssw = 0 exactly, and numpy's grand mean (the 12 rows added one after the other: the log2 matrix is C-contiguous, so the sum over axis 0 is sequential, not pairwise; checked bit for bit on
## GSE20677 and GSE49079, 25 Nov 2026; the first version of this port assumed the 8-way pairwise order and was right on GSE307491 and GSE49079 only because their floor values behave alike under both orders)
## differs from the value by about one ulp, so ssb > 0 and F = Inf, p = 0 (R's long-double colMeans would give 0/0 = NaN, p = 1). That artefact of the Python pipeline feeds the BH adjustment and so the
## number of DE genes of the archived results; the port reproduces it.
## ANOVA_EXACT=1 (added 26 Nov 2026 for the impact check O6 of the independent audit; default off, archived behaviour unchanged): a gene that is constant in all 12 samples gets p = 1, the value of exact
## arithmetic (ssb = ssw = 0, F = 0/0, which this function already maps to p = 1); every other gene keeps numpy's order.
ANOVA_EXACT <- identical(Sys.getenv("ANOVA_EXACT", ""), "1")
seq_sum_rows <- function(M) { s <- M[1, ]; for (i in seq_len(nrow(M))[-1]) s <- s + M[i, ]; s }
anova_fdr <- function(L) {
  gm <- t(vapply(1:B_K, function(g) { ix <- which(B_GROUP == g); s <- L[ix[1], ]; for (i in ix[-1]) s <- s + L[i, ]; s / length(ix) }, numeric(ncol(L)))); grand <- seq_sum_rows(L) / B_N
  d2 <- (gm - matrix(grand, B_K, ncol(L), byrow = TRUE))^2; ssb <- 4 * (d2[1, ] + d2[2, ] + d2[3, ]); ssw <- seq_sum_rows((L - gm[B_GROUP, , drop = FALSE])^2)
  F <- (ssb / (B_K - 1)) / (ssw / (B_N - B_K)); p <- pf(F, B_K - 1, B_N - B_K, lower.tail = FALSE); p[is.nan(p)] <- 1
  if (ANOVA_EXACT) p[colSums(L != matrix(L[1, ], nrow(L), ncol(L), byrow = TRUE)) == 0] <- 1
  list(p = p, fdr = p.adjust(p, "BH"))
}
spike <- function(Llog, den, rng) {                                  # positive control: shared latent factor added to 6 of the 14 non-denominator genes (log scale)
  cand <- setdiff(seq_len(ncol(Llog)), den); pick <- cand[np_choice_noreplace_idx(rng, length(cand), 6) + 1]
  gm <- t(vapply(1:B_K, function(g) colMeans(Llog[B_GROUP == g, , drop = FALSE]), numeric(ncol(Llog)))); sdv <- sqrt(colSums((Llog - gm[B_GROUP, , drop = FALSE])^2) / (B_N - B_K)); z <- np_standard_normal(rng, B_N)
  out <- Llog; out[, pick] <- out[, pick] + 1.5 * (z %o% sdv[pick]); list(out = out, pick = pick)
}
run_set <- function(X, den, seed) B_analyse(X, den, nperm = NPERM, seed = seed)
process <- function(row) {
  gse <- row$gse; gpl <- row$platform; n <- as.integer(substring(gse, 4)); gsms <- strsplit(row$gsms, ";", fixed = TRUE)[[1]]
  rec <- list(gse = gse, platform = gpl, order = as.integer(row$order), how = row$how, groups = row$groups)
  gmx <- gene_matrix(gse, gpl, gsms); L2 <- gmx$L2; genes <- gmx$genes; rec$n_genes_mapped <- length(genes); rec$scale <- if (gmx$is_log) "log2" else "linear"; rec$q99 <- gmx$q99
  mu <- colMeans(L2); expressed <- mu > median(mu); af <- anova_fdr(L2); fdr <- af$fdr; de <- which(expressed & fdr < 0.05); rec$n_expressed <- sum(expressed); rec$n_de <- length(de)
  if (length(de) < MIN_DE) { rec$status <- sprintf("not analysable: %d DE genes (< %d)", length(de), MIN_DE); return(rec) }
  if (DRY) { rec$status <- "analysable (dry run)"; return(rec) }
  de_sorted <- de[order(genes[de], method = "radix")]
  drawA <- de_sorted[np_choice_noreplace_idx(np_rng(2026 + n), length(de_sorted), NGENE) + 1]
  ex <- which(expressed); ex_sorted <- ex[order(genes[ex], method = "radix")]; drawB <- ex_sorted[np_choice_noreplace_idx(np_rng(c(2026 + n, 1)), length(ex_sorted), NGENE) + 1]
  for (tag in c("A", "B")) {
    draw <- if (tag == "A") drawA else drawB; res <- run_set(2^L2[, draw, drop = FALSE], 1:3, 2026 + n); res$genes <- as.list(genes[draw]); res$n_de_in_set <- sum(fdr[draw] < 0.05); rec[[tag]] <- res
  }
  LA <- log(2^L2[, drawA, drop = FALSE]); sp <- spike(LA, 1:3, np_rng(c(2026 + n, 2))); res <- run_set(exp(sp$out), 1:3, 2026 + n)
  res$spiked_genes <- as.list(genes[drawA[sp$pick]]); rec$A_spike <- res; rec$status <- "analysed"; rec
}
main <- function() {
  rows <- read.csv(FROZEN, stringsAsFactors = FALSE, check.names = FALSE); if (length(ONLY)) rows <- rows[rows$gse %in% ONLY, , drop = FALSE]
  done <- c(GPL570 = 0L, GPL1261 = 0L); fn <- file.path(OUT, if (DRY) "dry_run_per_dataset.jsonl" else "results_per_dataset.jsonl"); if (file.exists(fn)) file.remove(fn); used <- list()
  for (i in seq_len(nrow(rows))) {
    row <- as.list(rows[i, ]); gpl <- row$platform; if (done[[gpl]] >= MAXP) next
    t0 <- Sys.time(); gs <- strsplit(row$gsms, ";", fixed = TRUE)[[1]]; clash <- py_sort(unique(unlist(used[gs[gs %in% names(used)]])))
    rec <- if (length(clash)) list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = paste0("not analysable: duplicate of ", clash[1])) else
      tryCatch(process(row), error = function(e) { a <- archived_record(row$gse, e); if (!is.null(a)) a else
        list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = sprintf("not analysable: %s %s", class(e)[1], substr(conditionMessage(e), 1, 120))) })
    if (rec$status %in% c("analysed", "analysable (dry run)")) { done[[gpl]] <- done[[gpl]] + 1L; for (g in gs) used[[g]] <- row$gse }
    append_line(if (!is.null(attr(rec, "raw"))) attr(rec, "raw") else py_json_ascii(py_json_typed(rec, indent = NULL)), fn)
    cat(sprintf("%s %s (#%s): %s | DE %s | done %s | %.0fs%s\n", gpl, row$gse, row$order, rec$status, if (is.null(rec$n_de)) "-" else rec$n_de, paste0("{GPL570: ", done[["GPL570"]], ", GPL1261: ", done[["GPL1261"]], "}"), as.numeric(difftime(Sys.time(), t0, units = "secs")),
                if (!is.null(attr(rec, "raw"))) " | ARCHIVED RECORD COPIED, NOT COMPUTED" else ""))
  }
  cat("finished", paste0("{GPL570: ", done[["GPL570"]], ", GPL1261: ", done[["GPL1261"]], "}"), "\n")
}
if (sys.nframe() == 0L) main()
