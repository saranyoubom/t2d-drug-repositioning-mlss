## wp2_common.R - shared set-up of the R port of strengthen_v2/wp2_benchmark2 (protocol v2 section 3). Source it after defining R_PORT (folder of the port).
## ARC = archived Python results (compared against, never written); OUT = results of the R run (env WP2_OUT); DATA = downloaded GEO/NCBI files, untrusted data read as text/numbers only (env WP2_DATA_ROOT).
suppressMessages(library(jsonlite))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "agree.R"))
ARC <- Sys.getenv("WP2_ARCHIVE", normalizePath(file.path(R_PORT, "..", "wp2_benchmark2"))); P1 <- normalizePath(file.path(R_PORT, "..", "..", "benchmark_route3"))
OUT <- Sys.getenv("WP2_OUT", file.path(R_PORT, "results_R", "wp2")); dir.create(file.path(OUT, "compact"), showWarnings = FALSE, recursive = TRUE)
DATA <- Sys.getenv("WP2_DATA_ROOT", file.path(R_PORT, "results_R", "wp2_data")); dir.create(DATA, showWarnings = FALSE, recursive = TRUE)
options(py_json_nan = TRUE)                                               # json.dump writes NaN as NaN
NL <- .NL                                                                   # Python text-mode files end lines with os.linesep (CRLF on Windows)
## JSON lines written by Python (NaN and Infinity literals): the literals are turned into null (NA in R) before parsing
read_jsonl <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE); lines <- lines[nzchar(lines)]
  lapply(lines, function(l) fromJSON(gsub("(?<=[:,\\[ ])(-?(NaN|Infinity)|NA)(?=[,\\]}])", "null", l, perl = TRUE), simplifyVector = FALSE))      # (a bare NA in a line of an earlier R run is read as null too)
}
nullNA <- function(l) lapply(l, function(v) if (is.null(v)) NA_real_ else v)
## rows (list of named lists) -> data.frame, columns in first-appearance order, integer + double columns promote to double like pandas
rows_df <- function(rows) { d <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE))); rownames(d) <- NULL; d }
write_csv_crlf <- function(df, path) py_write_csv_repr(df, path, eol = "\r\n")        # pandas to_csv(index=False) on Windows
write_md <- function(lines, path, final_newline = FALSE) write_text_bytes(paste0(paste(lines, collapse = NL), if (final_newline) NL else ""), path)
sha256_file <- function(path) { if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(path, algo = "sha256", file = TRUE)); h <- system2("certutil", c("-hashfile", shQuote(normalizePath(path)), "SHA256"), stdout = TRUE); tolower(gsub("\\s", "", h[2])) }
