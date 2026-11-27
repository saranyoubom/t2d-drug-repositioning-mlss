## wp4_http.R - R port of wp4_islets/wp4_http.py (+ the download list of the fetch_*.py scripts): polite plain-HTTPS helpers and the fetch of the raw data of WP4. Rules of the work package: plain public HTTPS from NCBI / EBI with a plain
## tool User-Agent, at most 3 requests per second, nothing sent about the user, 403/404 and challenge pages are reported and NEVER worked around, a file above 300 MB is not fetched without an explicit decision (the size is read
## with a HEAD request first), downloads are untrusted data written to the folder given and only read, never executed.
## Usage: Rscript wp4_http.R [dataset ...]   datasets: GSE81608 GSE84133 GSE86469 GSE81547 E-MTAB-5061 (default: all); data folder WP4_DATA (default results_R/wp4_data); log in <data>/fetch.log
suppressMessages(library(curl))
USER_AGENT <- "geo-benchmark/1.0"; MIN_GAP <- 0.34; .H <- new.env(); .H$last <- 0
throttle <- function() { w <- .H$last + MIN_GAP - as.numeric(Sys.time()); if (w > 0) Sys.sleep(w); .H$last <- as.numeric(Sys.time()) }
is_challenge <- function(path) { b <- readBin(path, "raw", 2000); s <- rawToChar(b[b != as.raw(0)]); Encoding(s) <- "latin1"; s <- tolower(s); grepl("<html", s, fixed = TRUE) || grepl("<!doctype html", s, fixed = TRUE) || grepl("recaptcha", s, fixed = TRUE) }
head_info <- function(url, tries = 3) {      # status and size (bytes) without downloading
  last <- "no answer"
  for (k in 0:(tries - 1)) {
    throttle(); r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = USER_AGENT, nobody = TRUE, timeout = 60, followlocation = TRUE)), error = function(e) NULL)
    if (!is.null(r)) { if (r$status_code %in% c(403, 404, 405, 400)) return(list(status = r$status_code, size = NA_real_)); if (r$status_code == 200) { h <- parse_headers_list(r$headers); cl <- h[["content-length"]]; return(list(status = 200L, size = if (is.null(cl)) NA_real_ else as.numeric(cl))) }; last <- paste("HTTP", r$status_code) }
    Sys.sleep(3 + 3 * k)
  }
  list(status = last, size = NA_real_)
}
download <- function(url, path, max_mb = 300, tries = 3) {      # streams url into path (atomic .part rename); refuses above max_mb; 403 / 404 / challenge pages stop with an error
  if (file.exists(path) && file.size(path) > 0) return(file.size(path))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); h <- head_info(url)
  if (is.numeric(h$status) && h$status >= 400) stop(sprintf("HTTP %s (HEAD) %s", h$status, url))
  if (!is.na(h$size) && h$size > max_mb * 1e6) stop(sprintf("size %.1f MB above the limit of %d MB: %s", h$size / 1e6, max_mb, url))
  last <- "no answer"
  for (k in 0:(tries - 1)) {
    throttle(); part <- paste0(path, ".part"); r <- tryCatch(curl_fetch_disk(url, part, new_handle(useragent = USER_AGENT, timeout = 3600, followlocation = TRUE)), error = function(e) NULL)
    if (!is.null(r) && r$status_code == 200) {
      if (is_challenge(part) && !grepl("[.]html?$", path)) { file.rename(part, paste0(path, ".html_challenge")); stop("HTML challenge page instead of data: ", url) }
      file.rename(part, path); return(file.size(path))
    }
    last <- if (is.null(r)) "no answer" else paste("HTTP", r$status_code); unlink(part); if (!is.null(r) && r$status_code %in% c(403, 404)) stop(sprintf("%s %s", last, url)); Sys.sleep(3 + 3 * k)
  }
  stop(sprintf("download failed (%s) %s", last, url))
}
get_text <- function(url, tries = 3, timeout = 120, maxbytes = 20e6) {      # body of a plain GET as UTF-8 text (wp4_http.py get_text): 403 / 404 stop at once, other failures are retried after 3, 6, 9 s
  for (k in 0:(tries - 1)) {
    throttle(); r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = USER_AGENT, timeout = timeout, followlocation = TRUE)), error = function(e) NULL)
    if (!is.null(r)) {
      if (r$status_code %in% c(403, 404)) stop(sprintf("HTTP %d %s", r$status_code, url))
      if (r$status_code == 200) { if (length(r$content) > maxbytes) stop(sprintf("response larger than %.0f bytes: %s", maxbytes, url)); s <- iconv(rawToChar(r$content), "UTF-8", "UTF-8", sub = "?"); Encoding(s) <- "UTF-8"; return(s) }
    }
    Sys.sleep(3 + 3 * k)
  }
  stop("request failed ", url)
}
url_query <- function(l) paste(vapply(names(l), function(k) paste0(utils::URLencode(k, reserved = TRUE), "=", utils::URLencode(as.character(l[[k]]), reserved = TRUE)), ""), collapse = "&")      # urllib.parse.urlencode (spaces as %20, not +)
ATT_COLS <- c("when", "dataset", "item", "url", "size_bytes", "label_source", "outcome", "reason")
attempts_record <- function(path, dataset, item, url, size_bytes, label_source, outcome, reason = "") {      # attempts_log.py record(): append-only record of an access attempt (csv.writer, CRLF), one row per call
  new <- !file.exists(path); dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); q <- function(s) if (grepl("[\",\r\n]", s)) paste0("\"", gsub("\"", "\"\"", s, fixed = TRUE), "\"") else s
  row <- c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), dataset, item, url, if (is.null(size_bytes) || is.na(size_bytes)) "" else format(size_bytes, scientific = FALSE), label_source, outcome, reason)
  con <- file(path, "ab"); writeBin(charToRaw(enc2utf8(paste0(if (new) paste0(paste(ATT_COLS, collapse = ","), "\r\n") else "", paste(vapply(row, q, ""), collapse = ","), "\r\n"))), con); close(con)
}
## the raw data of WP4 (names as in parse_*.py: <data>/<dataset>/<file>)
WP4_FILES <- list(
  GSE81608 = c(series_matrix.txt.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81608/matrix/GSE81608_series_matrix.txt.gz", human_islets_rpkm.txt.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81608/suppl/GSE81608_human_islets_rpkm.txt.gz"),
  GSE84133 = c(GSE84133_RAW.tar = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84133/suppl/GSE84133_RAW.tar"),
  GSE86469 = c(GSE86469_series_matrix.txt.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86469/matrix/GSE86469_series_matrix.txt.gz", GSE86469_RSEM_expected_counts.csv.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86469/suppl/GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz"),
  GSE81547 = c(GSE81547_series_matrix.txt.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/matrix/GSE81547_series_matrix.txt.gz", GSE81547_RAW.tar = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/suppl/GSE81547_RAW.tar"),
  `E-MTAB-5061` = c(E.MTAB.5061.sdrf.txt = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/061/E-MTAB-5061/Files/E-MTAB-5061.sdrf.txt", E.MTAB.5061.processed.1.zip = "https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-5061/E-MTAB-5061.processed.1.zip"),
  GSE60803 = c(GSE60803_series_matrix.txt.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60803/matrix/GSE60803_series_matrix.txt.gz", GPL11532.annot.gz = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11532/annot/GPL11532.annot.gz"),
  GSE99503 = c(GSE99503_counts_GRCh38.p13_NCBI.tsv.gz = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE99503&format=file&file=GSE99503_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
  GSE167073 = c(GSE167073_series_matrix.txt.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167073/matrix/GSE167073_series_matrix.txt.gz", GPL17400_data.txt = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17400&targ=self&form=text&view=data"))
names(WP4_FILES[["E-MTAB-5061"]]) <- c("E-MTAB-5061.sdrf.txt", "E-MTAB-5061.processed.1.zip")
wp4_fetch <- function(datasets = names(WP4_FILES), root = NULL, log = NULL) {
  for (d in datasets) for (f in names(WP4_FILES[[d]])) {
    path <- file.path(root, d, f); t0 <- Sys.time()
    msg <- tryCatch(sprintf("ok %s %s %d bytes (%.0f s)", d, f, download(WP4_FILES[[d]][[f]], path), as.numeric(difftime(Sys.time(), t0, units = "secs"))), error = function(e) sprintf("FAILED %s %s: %s", d, f, conditionMessage(e)))
    cat(msg, "\n"); if (!is.null(log)) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = log, append = TRUE)
  }
}
if (sys.nframe() == 0L) {
  R_PORT <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), ".."))
  root <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); dir.create(root, showWarnings = FALSE, recursive = TRUE); a <- commandArgs(trailingOnly = TRUE)
  wp4_fetch(if (length(a)) a else names(WP4_FILES), root, file.path(root, "fetch.log")); cat("fetch finished\n", file = file.path(root, "fetch.log"), append = TRUE)
}
