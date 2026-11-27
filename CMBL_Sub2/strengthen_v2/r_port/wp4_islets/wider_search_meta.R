## wider_search_meta.R - R port of wp4_islets/wider_search_meta.py: metadata of the candidates that passed the reading of titles and summaries (B2 of PRESPEC_WP4b...): sample titles and characteristics from the GEO series matrix headers (only
## lines starting with '!' are read; the values table is never parsed), the supplementary file lists, and for ArrayExpress the SDRF and the study record. Used to write the contrast of each kept series BEFORE the freeze (B3).
## Downloads go to <WP4_DATA>/wider/<accession>/ (untrusted data, only read), each logged in attempts_R.csv. Output (default results_R/wp4/wider_search_meta): sample_records_<accession>.tsv (one row per sample: title, characteristics).
## Usage: Rscript wider_search_meta.R ACC [ACC ...]
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "wp4_http.R")); suppressMessages(library(jsonlite))
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); ROOT <- file.path(DATA, "wider"); OUT <- Sys.getenv("WP4_WSM_OUT", file.path(R_PORT, "results_R", "wp4", "wider_search_meta"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); ATT <- file.path(dirname(OUT), "attempts_R.csv"); G <- "https://ftp.ncbi.nlm.nih.gov/geo/series/"
geo_dir <- function(acc) { n <- as.integer(substring(acc, 4)); sprintf("%sGSE%dnnn/%s/", G, n %/% 1000L, acc) }
KEYS <- c("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1", "!Sample_characteristics_ch1", "!Sample_description", "!Sample_platform_id", "!Series_overall_design")
tsv_field <- function(s) if (is.na(s)) "" else if (grepl("[\t\"\r\n]", s)) paste0("\"", gsub("\"", "\"\"", s, fixed = TRUE), "\"") else s      # pandas to_csv(sep='\t'): quotes only fields with the separator, a quote or a line break
write_tsv_py <- function(df, path) { lines <- c(paste(vapply(names(df), tsv_field, ""), collapse = "\t"), if (nrow(df)) do.call(paste, c(lapply(df, function(x) vapply(x, tsv_field, "")), sep = "\t"))); write_text_bytes(paste0(paste(lines, collapse = "\r\n"), "\r\n"), path) }
fetch <- function(acc, item, url, path, max_mb, src, outcome, note) {
  tryCatch({ n <- download(url, path, max_mb = max_mb); attempts_record(ATT, acc, item, url, n, src, outcome, note); TRUE },
           error = function(e) { h <- head_info(url); attempts_record(ATT, acc, item, url, h$size, "", "failed or refused", substr(conditionMessage(e), 1, 200)); cat(acc, item, "FAILED", substr(conditionMessage(e), 1, 120), "\n"); FALSE })
}
split_tab <- function(l) { p <- strsplit(paste0(l, "\tEND"), "\t", fixed = TRUE)[[1]]; p[-length(p)] }
geo_meta <- function(acc) {
  d <- file.path(ROOT, acc); dir.create(d, showWarnings = FALSE, recursive = TRUE); base <- geo_dir(acc)
  fetch(acc, "series matrix", paste0(base, sprintf("matrix/%s_series_matrix.txt.gz", acc)), file.path(d, sprintf("%s_series_matrix.txt.gz", acc)), 300, "GEO sample records (authors)", "downloaded (metadata read only)", "size checked by HEAD; only lines starting with ! are read")
  fetch(acc, "suppl filelist", paste0(base, "suppl/filelist.txt"), file.path(d, "filelist.txt"), 300, "GEO sample records (authors)", "downloaded (metadata read only)", "size checked by HEAD; only lines starting with ! are read")
  rows <- list(); p <- file.path(d, sprintf("%s_series_matrix.txt.gz", acc))
  if (file.exists(p)) { con <- gzfile(p, "rt", encoding = "UTF-8"); repeat {
    l <- readLines(con, n = 1, warn = FALSE); if (!length(l) || startsWith(l, "!series_matrix_table_begin")) break
    if (any(startsWith(l, KEYS))) { pr <- split_tab(l); k <- pr[1]; rows[[k]] <- c(rows[[k]], list(gsub("^\"+|\"+$", "", pr[-1]))) } }; close(con) }
  if (is.null(rows[["!Sample_title"]])) return(list(T = NULL, extra = rows))
  n <- length(rows[["!Sample_title"]][[1]]); first <- function(k, dflt) if (is.null(rows[[k]])) dflt else rows[[k]][[1]]
  T <- list(gsm = first("!Sample_geo_accession", rep("", n)), title = first("!Sample_title", NULL), source = first("!Sample_source_name_ch1", rep("", n)))
  for (i in seq_along(rows[["!Sample_characteristics_ch1"]])) T[[sprintf("ch%d", i - 1L)]] <- rows[["!Sample_characteristics_ch1"]][[i]]
  if (!is.null(rows[["!Sample_platform_id"]])) T[["platform"]] <- rows[["!Sample_platform_id"]][[1]]
  fl <- file.path(d, "filelist.txt"); files <- if (file.exists(fl)) head(readLines(fl, warn = FALSE, encoding = "UTF-8"), 40) else character(0)
  list(T = as.data.frame(T, check.names = FALSE, stringsAsFactors = FALSE), extra = list(design = if (is.null(rows[["!Series_overall_design"]])) "" else substr(rows[["!Series_overall_design"]][[1]][1], 1, 600), files = files))
}
NA_STR <- c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null")      # the strings pandas reads as missing
ae_meta <- function(acc) {
  d <- file.path(ROOT, acc); dir.create(d, showWarnings = FALSE, recursive = TRUE); nn <- sub("^.*-", "", acc); base <- sprintf("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/%s/%s/Files/", substring(nn, nchar(nn) - 2L), acc); out <- list()
  for (fn in c(sprintf("%s.sdrf.txt", acc), sprintf("%s.idf.txt", acc))) if (fetch(acc, fn, paste0(base, fn), file.path(d, fn), 50, "ArrayExpress SDRF/IDF (authors)", "downloaded (metadata)", "")) out[[fn]] <- file.path(d, fn)
  s <- fromJSON(get_text(paste0("https://www.ebi.ac.uk/biostudies/api/v1/studies/", acc)), simplifyVector = FALSE); files <- list()
  walk <- function(x) { if (is.list(x) && !is.null(names(x))) { if (is.list(x$files)) for (f in x$files) { f1 <- if (is.null(names(f)) && length(f)) f[[1]] else f; files[[length(files) + 1L]] <<- list(f1$path, f1$size) }; for (v in x) walk(v) } else if (is.list(x)) for (v in x) walk(v) }
  walk(s); list(out = out, files = files)
}
sdrf_frame <- function(p) {      # pandas read_csv(sep='\t', dtype=str): missing-value strings become empty, repeated column names get the suffixes .1, .2
  X <- read.delim(p, sep = "\t", header = FALSE, quote = "\"", colClasses = "character", na.strings = character(0), comment.char = "", fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
  nm <- as.character(unlist(X[1, ])); X <- X[-1, , drop = FALSE]; cnt <- list(); for (i in seq_along(nm)) { k <- nm[i]; c0 <- if (is.null(cnt[[k]])) 0L else cnt[[k]]; if (c0 > 0L) nm[i] <- sprintf("%s.%d", k, c0); cnt[[k]] <- c0 + 1L }
  names(X) <- nm; X[] <- lapply(X, function(v) { v[v %in% NA_STR] <- ""; v }); rownames(X) <- NULL; X
}
if (sys.nframe() == 0L) for (acc in commandArgs(trailingOnly = TRUE)) {
  cat(strrep("=", 100), "\n", acc, "\n", sep = "")
  if (startsWith(acc, "GSE")) {
    r <- geo_meta(acc); if (is.null(r$T)) { cat("no sample records\n"); next }
    write_tsv_py(r$T, file.path(OUT, sprintf("sample_records_%s.tsv", acc))); cat("samples", nrow(r$T), "| platform", if (!is.null(r$T$platform)) r$T$platform[1] else "?", "| design:", substr(r$extra$design, 1, 400), "\n"); cat("suppl files:", head(r$extra$files, 12), sep = "\n  ")
  } else {
    r <- ae_meta(acc)
    if (length(r$out) && any(grepl("sdrf", names(r$out)))) { S <- sdrf_frame(r$out[[sprintf("%s.sdrf.txt", acc)]]); write_tsv_py(S, file.path(OUT, sprintf("sample_records_%s.tsv", acc))); cat("samples", nrow(S), "\n") }
    cat("files in the study record:", length(r$files), "\n")
  }
}
