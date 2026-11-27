## bench_data_io.R - R port of wp2_benchmark2/bench_data_io.py (+ ncbi_http.py): importable loaders of benchmark part 2 (WP2, protocol v2 section 3), no side effects at source time.
## The loaders of benchmark_route3/run_benchmark.R (fetch, series_path, read_matrix, annotation, gene_matrix, anova_fdr) plus the code path that turns one frozen series into the two 12 x 17 analysis
## matrices of the protocol (gene set A = 17 random differentially expressed genes, B = 17 random expressed genes). Development series, held-out arrays and RNA-seq series use the SAME functions
## (prepare_from_log -> prepare_array / prepare_rnaseq). Downloaded files are untrusted data: read as text or numbers only, never executed. NCBI politeness: throttle() keeps requests >= 0.34 s apart.
## Needs rng/numpy_rng.R, pycompat.R and bench_pipeline.R (B_N, B_K, B_GROUP) sourced first. Usage: io_configure(matrices, platforms).
suppressMessages({ library(data.table); library(curl) })
FTP <- "https://ftp.ncbi.nlm.nih.gov/geo/"; MIN_DE <- 34L; NGENE <- 17L; NPERM_IO <- 2000L
PANDAS_NA <- c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null")   # pandas default na_values
IO <- new.env(); IO$matrices <- NULL; IO$platforms <- NULL; IO$annot <- list(); IO$last <- 0
MIN_GAP <- as.numeric(Sys.getenv("WP2_MIN_GAP", "0.34")); RETRY_SCALE <- as.numeric(Sys.getenv("WP2_RETRY_SCALE", "1"))     # at most 3 requests per second by default; NCBI answers 403 to bursts: raise both for a slower, gentler run
io_configure <- function(matrices, platforms) { IO$matrices <- matrices; IO$platforms <- platforms; for (d in c(matrices, platforms)) dir.create(d, showWarnings = FALSE, recursive = TRUE) }
throttle <- function() { wait <- IO$last + MIN_GAP - as.numeric(Sys.time()); if (wait > 0) Sys.sleep(wait); IO$last <- as.numeric(Sys.time()) }
USER_AGENT <- "geo-benchmark/1.0"
.http404 <- function(url) structure(class = c("http404", "error", "condition"), list(message = paste("HTTP 404", url), call = NULL))
http_fetch <- function(url, path) { h <- new_handle(useragent = USER_AGENT, timeout = 180); tryCatch(curl_fetch_disk(url, path, h), error = function(e) NULL) }
fetch <- function(url, path) {         # download to path (atomic .part rename) unless it exists; retries on transient errors, 404 is raised at once
  if (file.exists(path) && file.size(path) > 0) return(path)
  last <- "no answer"
  for (k in 0:3) {
    throttle(); r <- http_fetch(url, paste0(path, ".part")); last <- if (is.null(r)) "no answer" else paste("HTTP", r$status_code)
    if (!is.null(r) && r$status_code == 200) { file.rename(paste0(path, ".part"), path); return(path) }
    if (!is.null(r) && r$status_code == 404) { unlink(paste0(path, ".part")); stop(.http404(url)) }
    unlink(paste0(path, ".part")); Sys.sleep((3 + 3 * k) * RETRY_SCALE)
  }
  stop("download failed (", last, ") ", url)
}
http_text <- function(url, tries = 4) {  # GET as text, 404 raised, other errors retried
  for (k in 0:(tries - 1)) {
    throttle(); r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = USER_AGENT, timeout = 90)), error = function(e) NULL)
    if (!is.null(r) && r$status_code == 200) return(rawToChar(r$content)); if (!is.null(r) && r$status_code == 404) stop(.http404(url)); Sys.sleep(3 + 3 * k)
  }
  stop("download failed ", url)
}
series_path <- function(gse, gpl) {
  n <- as.integer(substring(gse, 4)); folder <- if (n >= 1000) sprintf("GSE%dnnn", n %/% 1000) else "GSEnnn"; base <- sprintf("%sseries/%s/%s/matrix/", FTP, folder, gse); path <- file.path(IO$matrices, sprintf("%s_%s.txt.gz", gse, gpl))
  tryCatch(fetch(paste0(base, gse, "_series_matrix.txt.gz"), path), http404 = function(e) {
    idx <- http_text(base); hit <- regmatches(idx, gregexpr("href=\"[^\"]*series_matrix\\.txt\\.gz\"", idx))[[1]]; hit <- sub("^href=\"", "", sub("\"$", "", hit)); hit <- hit[grepl(gpl, hit, fixed = TRUE)]
    fetch(paste0(base, hit[1]), path)
  })
}
## The text table between two marker lines of a GEO file, streamed to a temporary plain file in chunks (memory: one chunk, not the whole file; this laptop has about 190 MB free), then only the
## wanted columns are read: select(header names) returns the names to keep (the first column always). Returns a data frame of character columns; the temporary file is removed.
table_block_cols <- function(path, begin, end, select) {
  con <- gzfile(path, "rt"); tmp <- tempfile(fileext = ".tsv"); out <- file(tmp, "w"); state <- 0L
  tryCatch({
    repeat {
      ch <- readLines(con, n = 20000L, warn = FALSE, encoding = "UTF-8"); if (!length(ch)) break
      if (state == 0L) { ib <- which(startsWith(ch, begin))[1]; if (is.na(ib)) next; ch <- ch[-seq_len(ib)]; state <- 1L }
      if (state == 1L) { ie <- which(startsWith(ch, end))[1]; if (!is.na(ie)) { ch <- ch[seq_len(ie - 1L)]; state <- 2L }; if (length(ch)) writeLines(ch, out) }
      if (state == 2L) break
    }
  }, finally = { close(con); close(out) })
  on.exit(unlink(tmp), add = TRUE); if (state != 2L) stop("table markers not found in ", basename(path))
  hdr <- gsub("^\"+|\"+$", "", strsplit(readLines(tmp, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]); want <- unique(c(1L, which(hdr %in% select(hdr))))
  fread(tmp, sep = "\t", quote = "\"", header = TRUE, select = want, colClasses = "character", na.strings = PANDAS_NA, data.table = FALSE, fill = TRUE)
}
read_matrix <- function(path, gsms) table_block_cols(path, "!series_matrix_table_begin", "!series_matrix_table_end", function(h) gsms)     # probe ids and the selected samples (missing samples are not returned)
annotation <- function(gpl) {
  n <- as.integer(substring(gpl, 4)); folder <- if (n >= 1000) sprintf("GPL%dnnn", n %/% 1000) else "GPLnnn"
  path <- fetch(sprintf("%splatforms/%s/%s/annot/%s.annot.gz", FTP, folder, gpl, gpl), file.path(IO$platforms, paste0(gpl, ".annot.gz")))
  t <- table_block_cols(path, "!platform_table_begin", "!platform_table_end", function(h) c("ID", h[grepl("symbol", tolower(h), fixed = TRUE)][1])); col <- names(t)[grepl("symbol", tolower(names(t)), fixed = TRUE)][1]
  out <- data.frame(probe = t$ID, symbol = t[[col]], stringsAsFactors = FALSE); out <- out[!is.na(out$probe) & !is.na(out$symbol), ]
  out$symbol <- trimws(sub(" /// .*$", "", out$symbol)); out <- out[out$symbol != "" & out$symbol != "---", ]; setNames(out$symbol, out$probe)
}
## sequential (outer-axis) numpy mean of the rows of M: s = M[1] + M[2] + ... then / n (R's colMeans accumulates in long double)
np_colmeans <- function(M) { s <- M[1, ]; for (i in seq_len(nrow(M))[-1]) s <- s + M[i, ]; s / nrow(M) }
gene_matrix <- function(gse, gpl, gsms) {   # 12 samples x genes on the log2 scale (probe with the highest mean linear expression per gene) and the scale decision
  M <- read_matrix(series_path(gse, gpl), gsms); miss <- setdiff(gsms, names(M))
  if (length(miss)) stop(structure(class = c("ValueError", "error", "condition"), list(message = sprintf("samples missing: [%s]", paste(sprintf("'%s'", head(miss, 3)), collapse = ", ")), call = NULL)))
  .zde <- function() stop(structure(class = c("ZeroDivisionError", "error", "condition"), list(message = "division by zero", call = NULL)))   # Python: a matrix without usable rows ends in a ZeroDivisionError (logged as such in the archived runs)
  if (!nrow(M)) .zde()
  X <- vapply(gsms, function(g) suppressWarnings(as.numeric(M[[g]])), numeric(nrow(M))); rownames(X) <- M[[1]]
  q99 <- np_quantile(X[!is.na(X)], 0.99); is_log <- q99 < 30; lin <- if (is_log) 2^X else X
  ok <- rowSums(!is.na(lin) & lin > 0) == ncol(lin); lin <- lin[ok, , drop = FALSE]
  if (is.null(IO$annot[[gpl]])) IO$annot[[gpl]] <- annotation(gpl)
  sym <- unname(IO$annot[[gpl]][rownames(lin)]); keep <- !is.na(sym); lin <- lin[keep, , drop = FALSE]; sym <- sym[keep]
  if (!nrow(lin)) .zde()
  mn <- rowMeans(lin); o <- order(sym, -mn, seq_along(sym), method = "radix"); best <- o[!duplicated(sym[o])]       # one probe per gene symbol: the highest mean, first in file order on ties
  list(L2 = log2(t(lin[best, , drop = FALSE])), genes = sym[best], is_log = is_log, q99 = q99)
}
## One-way ANOVA across the 3 groups for every gene (columns of L), BH FDR, with numpy's summation order. The log2 matrix of gene_matrix() is C-contiguous (12 x genes), and numpy's sum over axis 0 of such an
## array adds the rows one after the other (checked bit for bit on GSE20677 and GSE49079: 22,189 of 22,189 probes equal to the sequential sum, the 8-way pairwise sum differs in 5,811). A probe that is constant
## in all 12 samples (the detection floor) has ssw = 0 exactly, while the sequentially summed grand mean can differ from the value by one ulp (the sum of 12 equal values is rounded), so ssb > 0, F = Inf, p = 0
## (artefact of the Python pipeline, reproduced on purpose); with an exact mean ssb = 0 and F = 0/0 = NaN, p = 1.
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
## shared code path of part 1 (run_benchmark.process) after the log-scale gene matrix exists: expressed genes (mean log expression above the median), differential expression (ANOVA, BH FDR < 0.05),
## analysability (>= MIN_DE DE genes) and the two random draws of 17 genes (seed 2026 + GSE number; set B seed [2026 + n, 1]). xfun(draw) gives the 12 x 17 linear matrix (draw: 1-based column indices).
## Returns list(rec, sets): sets = list(A = list(X, genes, n_de_in_set, idx), B = ...) or NULL when not analysable / dry run.
prepare_from_log <- function(rec, L2, genes, n, xfun, dry = FALSE) {
  mu <- np_colmeans(L2); expressed <- mu > median(mu); af <- anova_fdr(L2); fdr <- af$fdr; de <- which(expressed & fdr < 0.05); rec$n_expressed <- as.integer(sum(expressed)); rec$n_de <- as.integer(length(de))
  if (length(de) < MIN_DE) { rec$status <- sprintf("not analysable: %d DE genes (< %d)", length(de), MIN_DE); return(list(rec = rec, sets = NULL)) }
  if (dry) { rec$status <- "analysable (dry run)"; return(list(rec = rec, sets = NULL)) }
  de_sorted <- de[order(genes[de], method = "radix")]; drawA <- de_sorted[np_choice_noreplace_idx(np_rng(2026 + n), length(de_sorted), NGENE) + 1]
  ex <- which(expressed); ex_sorted <- ex[order(genes[ex], method = "radix")]; drawB <- ex_sorted[np_choice_noreplace_idx(np_rng(c(2026 + n, 1)), length(ex_sorted), NGENE) + 1]
  sets <- list()
  for (tag in c("A", "B")) { draw <- if (tag == "A") drawA else drawB; sets[[tag]] <- list(X = as.matrix(xfun(draw)), genes = genes[draw], n_de_in_set = as.integer(sum(fdr[draw] < 0.05)), idx = draw) }
  rec$status <- "prepared"; list(rec = rec, sets = sets)
}
prepare_array <- function(row, dry = FALSE) {   # one frozen microarray series (row of a frozen csv, a list) -> list(rec, sets); same record fields as part 1
  gse <- row$gse; gpl <- row$platform; n <- as.integer(substring(gse, 4)); gsms <- strsplit(row$gsms, ";", fixed = TRUE)[[1]]
  rec <- list(gse = gse, platform = gpl, order = as.integer(row$order), how = row$how, groups = row$groups)
  gm <- gene_matrix(gse, gpl, gsms); rec$n_genes_mapped <- as.integer(length(gm$genes)); rec$scale <- if (gm$is_log) "log2" else "linear"; rec$q99 <- gm$q99
  L2 <- gm$L2; prepare_from_log(rec, L2, gm$genes, n, function(d) 2^L2[, d, drop = FALSE], dry)
}
## ---------------------------------------------------------------------------------------------------------------- RNA-seq (protocol 3.5)
RNASEQ_BUILD <- c(human = "GRCh38.p13", mouse = "GRCm39")
COUNT_URL <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_raw_counts_%s_NCBI.tsv.gz"
counts_path <- function(root, gse, build) file.path(root, "v2_data", "rnaseq", "counts", sprintf("%s_%s.tsv.gz", gse, build))
fetch_counts <- function(root, gse, org) {   # list(path, why): the NCBI-generated raw count file (cached) or the reason why there is none (run_rnaseq.fetch_counts; http_get: 404 and 403 are never retried)
  build <- RNASEQ_BUILD[[org]]; path <- counts_path(root, gse, build); dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(path) && file.size(path) > 0) return(list(path = path, why = NULL))
  url <- sprintf(COUNT_URL, gse, gse, build)
  for (k in 0:3) {
    throttle(); r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = USER_AGENT, timeout = 180)), error = function(e) NULL)
    if (!is.null(r) && r$status_code %in% c(404L, 403L)) return(list(path = NULL, why = sprintf("no NCBI-generated counts (HTTP %d)", r$status_code)))
    if (!is.null(r) && r$status_code == 200) {
      if (!identical(r$content[1:2], as.raw(c(0x1f, 0x8b)))) return(list(path = NULL, why = "blocked: the reply was a web page, not a count file"))
      writeBin(r$content, paste0(path, ".part")); file.rename(paste0(path, ".part"), path); return(list(path = path, why = NULL))
    }
    Sys.sleep(3 + 3 * k)
  }
  stop("download failed ", url)
}
read_counts_file <- function(path, gsms) {     # genes x 12 matrix of the selected samples (by accession) and the gene ids (run_rnaseq.read_counts)
  T <- fread(path, sep = "\t", header = TRUE, colClasses = list(character = 1L), data.table = FALSE, na.strings = PANDAS_NA)
  miss <- setdiff(gsms, names(T)); if (length(miss)) stop(structure(class = c("ValueError", "error", "condition"), list(message = sprintf("%d of 12 selected samples are missing in the count file", length(miss)), call = NULL)))
  C <- vapply(gsms, function(g) suppressWarnings(as.numeric(T[[g]])), numeric(nrow(T)))
  if (anyNA(C)) stop(structure(class = c("ValueError", "error", "condition"), list(message = "non-numeric entries in the count file", call = NULL)))
  list(C = C, ids = as.character(T[[1]]))
}
rnaseq_cpm <- function(counts) {   # counts: genes x 12 matrix of raw counts. CPM on all genes of the file; universe = genes with count >= 1 in all 12 samples
  C <- as.matrix(counts); storage.mode(C) <- "double"; lib <- colSums(C); cpm <- C / matrix(lib, nrow(C), ncol(C), byrow = TRUE) * 1e6; keep <- which(rowSums(C >= 1) == ncol(C))
  list(cpm = t(cpm[keep, , drop = FALSE]), L2 = log2(t(cpm[keep, , drop = FALSE]) + 1.0), keep = keep)
}
prepare_rnaseq <- function(rec, counts, ids, symbols, n, dry = FALSE) {   # genes need count >= 1 in all 12 samples; expressed and DE as in part 1 on log2(CPM + 1); X for the analysis = CPM of the drawn genes
  r <- rnaseq_cpm(counts); keep <- r$keep
  genes <- vapply(keep, function(i) if (!is.null(symbols[[i]]) && !is.na(symbols[[i]]) && nzchar(as.character(symbols[[i]]))) as.character(symbols[[i]]) else as.character(ids[[i]]), "")
  C <- as.matrix(counts); storage.mode(C) <- "double"
  rec$n_genes_mapped <- as.integer(length(genes)); rec$scale <- "counts->CPM"; rec$n_genes_total <- as.integer(nrow(C)); rec$min_lib <- min(colSums(C))
  cpm <- r$cpm; prepare_from_log(rec, r$L2, genes, n, function(d) cpm[, d, drop = FALSE], dry)
}
