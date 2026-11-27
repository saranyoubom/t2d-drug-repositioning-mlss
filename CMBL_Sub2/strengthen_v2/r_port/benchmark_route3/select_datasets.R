## select_datasets.R - R port of benchmark_route3/select_datasets.py. Blind selection of the benchmark datasets (PROTOCOL sections 3 and 9). Only sample METADATA (series-matrix header) are read here;
## no expression value is touched. Writes eligibility.csv (every screened series with its status or reason) and frozen_datasets.csv (the ordered eligible list with the chosen groups and the 12 sample accessions),
## then FROZEN_SHA256.txt. Downloaded headers are cached in <data dir>/headers (untrusted data, never executed). Usage: Rscript select_datasets.R <data dir>
## Output goes to RB_OUT (default results_R/benchmark_route3), never over the archived, fingerprinted list. The candidate list comes from NCBI GEO (esearch/esummary) and is shuffled with Python's
## random.Random(2026) (rng/pyrandom.R); it changes when GEO receives new series, so a rerun today need not give the list of 4 Nov 2026. To test the screening alone on the archived series set
## RB_SELECT_FROM=<archived eligibility.csv>: the (platform, order, gse) rows of that file are then screened in that order instead of querying GEO.
suppressMessages({ library(curl); library(jsonlite) })
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "rng", "pyrandom.R"))
args <- commandArgs(trailingOnly = TRUE); OUT <- Sys.getenv("RB_OUT", file.path(R_PORT, "results_R", "benchmark_route3")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DATA <- if (length(args) >= 1) args[1] else file.path(OUT, "data"); dir.create(file.path(DATA, "headers"), showWarnings = FALSE, recursive = TRUE); FROM <- Sys.getenv("RB_SELECT_FROM", "")
EU <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; FTP <- "https://ftp.ncbi.nlm.nih.gov/geo/series/"; PLATFORMS <- c("GPL570", "GPL1261"); N_SCREEN <- 150L
QUERY <- "%s[ACCN] AND gse[ETYP] AND 12:36[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine)"
KEYWORDS <- c("treatment", "agent", "compound", "drug", "stimulus", "stimulation", "condition", "genotype", "knockdown", "knockout", "sirna", "shrna", "transfection", "infection", "cytokine", "ligand", "inhibitor", "vehicle", "dose", "protocol", "variation")
CONTROL <- "(?i)\\b(control|ctrl|vehicle|dmso|untreated|mock|nt|scrambled|wild|wt|naive|none|sham|pbs|saline|nontarget)\\b"
TITLE_TAIL <- "(?i)[\\s_\\-\\.,#]*(rep(licate)?|sample|replic)?[\\s_\\-\\.#]*\\d{1,2}\\s*$"
get_url <- function(url, tries = 4) {
  for (k in 0:(tries - 1)) {
    r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = "geo-benchmark/1.0", timeout = 90)), error = function(e) NULL)
    if (!is.null(r) && r$status_code == 200) return(rawToChar(r$content)); if (!is.null(r) && r$status_code == 404) stop(structure(class = c("http404", "error", "condition"), list(message = paste("HTTP 404", url), call = NULL)))
    Sys.sleep(2 + 2 * k)
  }
  stop("failed ", url)
}
candidates <- function(gpl) {
  term <- sprintf(QUERY, gpl); j <- fromJSON(get_url(paste0(EU, "esearch.fcgi?db=gds&term=", URLencode(term, reserved = TRUE), "&retmax=5000&retmode=json")), simplifyVector = FALSE); Sys.sleep(0.4)
  ids <- unlist(j$esearchresult$idlist); accs <- character(0)
  for (i in seq(1, length(ids), by = 150)) {
    chunk <- ids[i:min(i + 149, length(ids))]; s <- fromJSON(get_url(paste0(EU, "esummary.fcgi?db=gds&id=", paste(chunk, collapse = ","), "&retmode=json")), simplifyVector = FALSE)$result; Sys.sleep(0.4)
    accs <- c(accs, unlist(lapply(chunk, function(x) if (!is.null(s[[x]]) && startsWith(as.character(s[[x]]$accession), "GSE")) s[[x]]$accession)))
  }
  accs <- unique(accs); accs <- accs[order(as.integer(substring(accs, 4)))]; accs <- py_shuffle(py_random(2026), accs)       # sorted by accession number, then random.Random(2026).shuffle
  list(count = as.integer(j$esearchresult$count), accs = accs)
}
matrix_url <- function(gse, gpl) {
  n <- as.integer(substring(gse, 4)); folder <- if (n >= 1000) sprintf("GSE%dnnn", n %/% 1000) else "GSEnnn"; base <- sprintf("%s%s/%s/matrix/", FTP, folder, gse); u <- paste0(base, gse, "_series_matrix.txt.gz")
  st <- tryCatch(curl_fetch_memory(u, new_handle(useragent = "geo-benchmark/1.0", timeout = 90, nobody = TRUE))$status_code, error = function(e) NA)
  if (isTRUE(st == 200)) return(u)
  idx <- get_url(base); names_ <- regmatches(idx, gregexpr("href=\"[^\"]*series_matrix\\.txt\\.gz\"", idx))[[1]]; names_ <- sub("^href=\"", "", sub("\"$", "", names_)); hit <- names_[grepl(gpl, names_, fixed = TRUE)]
  if (length(hit)) paste0(base, hit[1]) else NULL
}
header <- function(gse, gpl) {                                   # the !Sample_*, !Series_title and !Series_geo_accession lines before the data table (streamed, the download stops there)
  path <- file.path(DATA, "headers", sprintf("%s_%s.txt", gse, gpl)); if (file.exists(path)) return(strsplit(readChar(path, file.size(path), useBytes = TRUE), "\n", fixed = TRUE)[[1]])
  url <- matrix_url(gse, gpl); if (is.null(url)) return(NULL)
  con <- gzcon(curl(url, open = "rb", handle = new_handle(useragent = "geo-benchmark/1.0", timeout = 90))); on.exit(try(close(con), silent = TRUE)); lines <- character(0)
  repeat { l <- readLines(con, n = 1L, warn = FALSE, encoding = "UTF-8"); if (!length(l) || startsWith(l, "!series_matrix_table_begin")) break; l <- sub("[\r\n]+$", "", l); if (any(startsWith(l, c("!Sample_", "!Series_title", "!Series_geo_accession")))) lines <- c(lines, l) }
  write_text_bytes(paste(lines, collapse = "\n"), path); lines
}
parse_header <- function(lines) {
  rows <- list(); for (l in lines) { parts <- strsplit(l, "\t", fixed = TRUE)[[1]]; vals <- gsub("^\"+|\"+$", "", parts[-1]); rows[[parts[1]]] <- c(rows[[parts[1]]], list(vals)) }
  first <- function(k) if (is.null(rows[[k]])) character(0) else rows[[k]][[1]]
  gsm <- first("!Sample_geo_accession"); title <- first("!Sample_title"); plat <- unique(first("!Sample_platform_id")); keys <- character(0); per <- replicate(length(gsm), list(), simplify = FALSE)
  for (line in rows[["!Sample_characteristics_ch1"]]) for (i in seq_len(min(length(line), length(gsm)))) {
    v <- line[i]; pos <- regexpr(": ", v, fixed = TRUE); if (pos < 0) { k <- "unnamed"; val <- v } else { k <- substr(v, 1, pos - 1); val <- substring(v, pos + 2) }
    k <- tolower(trimws(k)); if (!(k %in% keys)) keys <- c(keys, k); per[[i]][[k]] <- trimws(val)
  }
  list(gsm = gsm, title = title, plat = plat, keys = keys, per = per)
}
groups_from <- function(labels) { d <- list(); for (i in seq_along(labels)) { lab <- labels[i]; if (!is.na(lab) && nzchar(lab)) d[[lab]] <- c(d[[lab]], i) }; d }
screen <- function(lines, gpl) {
  p <- parse_header(lines); n <- length(p$gsm)
  if (!(n >= 12 && n <= 36)) return(list(status = sprintf("ineligible: %d samples", n), info = NULL))
  if (!(length(p$plat) == 1 && p$plat == gpl)) return(list(status = sprintf("ineligible: platforms [%s]", paste(sprintf("'%s'", py_sort(p$plat)), collapse = ", ")), info = NULL))
  chosen <- NULL; how <- NULL
  for (k in p$keys) if (any(vapply(KEYWORDS, function(w) grepl(w, k, fixed = TRUE), TRUE))) {
    d <- groups_from(vapply(p$per, function(x) if (is.null(x[[k]])) NA_character_ else x[[k]], "")); q <- d[vapply(d, length, 1L) >= 4]
    if (length(q) >= 3) { chosen <- q; how <- paste0("characteristic: ", k); break }
  }
  if (is.null(chosen)) { d <- groups_from(sub(TITLE_TAIL, "", p$title, perl = TRUE)); q <- d[vapply(d, length, 1L) >= 4]; if (length(q) >= 3) { chosen <- q; how <- "title" } }
  if (is.null(chosen)) return(list(status = "ineligible: no grouping with >= 3 groups of >= 4 samples", info = NULL))
  labels <- names(chosen); ctrl <- labels[grepl(CONTROL, gsub("_", " ", labels, fixed = TRUE), perl = TRUE)]; first <- if (length(ctrl)) ctrl[1] else py_sort(labels)[1]
  sel <- c(first, head(labels[labels != first], 2))
  list(status = "eligible", info = list(how = how, n_groups = length(chosen), groups = sel, gsms = lapply(sel, function(l) p$gsm[head(chosen[[l]], 4)])))
}
sha256_file <- function(path) { if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(path, algo = "sha256", file = TRUE)); h <- system2("certutil", c("-hashfile", shQuote(normalizePath(path)), "SHA256"), stdout = TRUE); tolower(gsub("\\s", "", h[2])) }
main <- function() {
  elig <- list(); frozen <- list(); src <- if (nzchar(FROM)) read.csv(FROM, stringsAsFactors = FALSE) else NULL
  for (gpl in PLATFORMS) {
    if (is.null(src)) { cand <- candidates(gpl); cat(sprintf("%s: %d series returned; screening the first %d of the shuffled list\n", gpl, cand$count, N_SCREEN)); accs <- head(cand$accs, N_SCREEN) }
    else { s <- src[src$platform == gpl, ]; accs <- s$gse[order(s$order)]; nsel <- suppressWarnings(as.integer(Sys.getenv("RB_SELECT_N", ""))); if (!is.na(nsel)) accs <- head(accs, nsel)      # RB_SELECT_N: only the first n per platform
           cat(sprintf("%s: screening the %d series of %s in their archived order\n", gpl, length(accs), basename(FROM))) }
    for (pos in seq_along(accs)) {
      gse <- accs[pos]
      r <- tryCatch({ lines <- header(gse, gpl); if (is.null(lines)) list(status = "ineligible: no series matrix for the platform", info = NULL) else screen(lines, gpl) }, error = function(e) list(status = sprintf("ineligible: error %s", class(e)[1]), info = NULL))
      info <- r$info
      elig[[length(elig) + 1]] <- list(platform = gpl, order = pos, gse = gse, status = r$status, how = if (is.null(info)) "" else info$how, n_groups = if (is.null(info)) "" else info$n_groups, groups = if (is.null(info)) "" else paste(info$groups, collapse = "|"))
      if (!is.null(info)) frozen[[length(frozen) + 1]] <- list(platform = gpl, order = pos, gse = gse, how = info$how, n_groups = info$n_groups, groups = paste(info$groups, collapse = "|"), gsms = paste(unlist(info$gsms), collapse = ";"))
      if (pos %% 25 == 0) cat(sprintf("  %s: %d screened, %d eligible\n", gpl, pos, sum(vapply(frozen, function(f) f$platform == gpl, TRUE))))
      Sys.sleep(0.3)
    }
  }
  dictwriter_csv(elig, names(elig[[1]]), file.path(OUT, "eligibility.csv")); dictwriter_csv(frozen, names(frozen[[1]]), file.path(OUT, "frozen_datasets.csv"))
  h <- sha256_file(file.path(OUT, "frozen_datasets.csv")); write_text_bytes(sprintf("%s *frozen_datasets.csv\n%s\n", h, format(Sys.time(), "frozen %Y-%m-%d %H:%M:%S %z")), file.path(OUT, "FROZEN_SHA256.txt"))
  cat("eligible per platform:", paste0("{", paste(sprintf("'%s': %d", PLATFORMS, vapply(PLATFORMS, function(g) sum(vapply(frozen, function(f) f$platform == g, TRUE)), 1L)), collapse = ", "), "}"), "| sha256", h, "\n")
}
main()
