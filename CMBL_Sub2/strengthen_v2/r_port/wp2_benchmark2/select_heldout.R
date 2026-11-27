## select_heldout.R - R port of wp2_benchmark2/select_heldout.py (+ the metadata helpers that select_rnaseq.py imports from it): protocol v2 section 3.5, held-out microarrays (WP2 deliverable D).
## Blind selection: only sample METADATA (series-matrix headers) are read, no expression value. Copy of the selection logic of benchmark_route3/select_datasets (same query, numeric sort and shuffle with
## random.Random(2026), same parsing, perturbation keywords, title fallback, control-like first group, 3 groups of 4 samples) with ONE change: the screening window is positions 151-450 of the same shuffled
## lists (GPL570 and GPL1261) instead of 1-150. Added rules: (1) the re-run query must reproduce the part 1 lists (positions 1-150 equal to benchmark_route3/eligibility.csv), otherwise the script stops;
## (2) a series that shares at least one sample (any sample of the series) with one of the 37 analysed development series is excluded and listed with that reason; (3) NCBI requests at most 3 per second.
## Writes eligibility_heldout.csv, frozen_heldout_arrays.csv (and the SHA-256 of the frozen list to <OUT>/PROTOCOL_v2_SHA256_R.txt: the archived fingerprint file is never touched).
## The candidate list (esearch) grows with GEO, so a query today does not reproduce the list of 5 Nov 2026 and rule (1) stops a live run; WP2_SELECT_FROM=<compact/heldout_arrays_screened_order.json> screens the
## archived window instead (the same series in the same order) and checks the result against the fingerprint of the archived frozen list.
## Usage: Rscript select_heldout.R [--check-only | --download-only]      Headers are cached in <data>/v2_data/arrays_heldout/headers (untrusted data, never executed).
if (!exists("HERE") || !file.exists(file.path(HERE, "wp2_common.R"))) { HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, "..")) }
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "rng", "pyrandom.R")); source(file.path(HERE, "bench_data_io.R")); suppressMessages(library(curl))
EU <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; FTP_SERIES <- "https://ftp.ncbi.nlm.nih.gov/geo/series/"; PLATFORMS <- c("GPL570", "GPL1261"); WINDOW <- c(150L, 450L)       # python slice of the shuffled list: positions 151-450
QUERY <- "%s[ACCN] AND gse[ETYP] AND 12:36[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine)"
KEYWORDS <- c("treatment", "agent", "compound", "drug", "stimulus", "stimulation", "condition", "genotype", "knockdown", "knockout", "sirna", "shrna", "transfection", "infection", "cytokine", "ligand", "inhibitor", "vehicle", "dose", "protocol", "variation")
CONTROL <- "(?i)\\b(control|ctrl|vehicle|dmso|untreated|mock|nt|scrambled|wild|wt|naive|none|sham|pbs|saline|nontarget)\\b"
TITLE_TAIL <- "(?i)[\\s_\\-\\.,#]*(rep(licate)?|sample|replic)?[\\s_\\-\\.#]*\\d{1,2}\\s*$"
candidates <- function(gpl) {
  term <- sprintf(QUERY, gpl); j <- fromJSON(http_text(paste0(EU, "esearch.fcgi?", "db=gds&term=", URLencode(term, reserved = TRUE), "&retmax=5000&retmode=json")), simplifyVector = FALSE)
  ids <- unlist(j$esearchresult$idlist); accs <- character(0)
  for (i in seq(1, length(ids), by = 150)) {
    chunk <- ids[i:min(i + 149, length(ids))]; s <- fromJSON(http_text(paste0(EU, "esummary.fcgi?db=gds&id=", paste(chunk, collapse = ","), "&retmode=json")), simplifyVector = FALSE)$result
    accs <- c(accs, unlist(lapply(chunk, function(x) if (!is.null(s[[x]]) && startsWith(as.character(s[[x]]$accession), "GSE")) s[[x]]$accession)))
  }
  accs <- unique(accs); accs <- accs[order(as.integer(substring(accs, 4)))]; accs <- py_shuffle(py_random(2026), accs)
  list(count = as.integer(j$esearchresult$count), accs = accs)
}
matrix_url_sh <- function(gse, gpl) {
  n <- as.integer(substring(gse, 4)); folder <- if (n >= 1000) sprintf("GSE%dnnn", n %/% 1000) else "GSEnnn"; base <- sprintf("%s%s/%s/matrix/", FTP_SERIES, folder, gse); u <- paste0(base, gse, "_series_matrix.txt.gz")
  throttle(); st <- tryCatch(curl_fetch_memory(u, new_handle(useragent = USER_AGENT, timeout = 90, nobody = TRUE))$status_code, error = function(e) NA)
  if (isTRUE(st == 200)) return(u)
  idx <- tryCatch(http_text(base), http404 = function(e) NULL); if (is.null(idx)) return(NULL)
  nm <- regmatches(idx, gregexpr("href=\"[^\"]*series_matrix\\.txt\\.gz\"", idx))[[1]]; nm <- sub("^href=\"", "", sub("\"$", "", nm)); hit <- nm[grepl(gpl, nm, fixed = TRUE)]
  if (length(hit)) paste0(base, hit[1]) else NULL
}
stream_header <- function(url) {       # the !Sample_*, !Series_title and !Series_geo_accession lines before the data table (the download stops there)
  throttle(); con <- gzcon(curl(url, open = "rb", handle = new_handle(useragent = USER_AGENT, timeout = 90))); on.exit(try(close(con), silent = TRUE)); lines <- character(0)
  repeat { l <- readLines(con, n = 1L, warn = FALSE, encoding = "UTF-8"); if (!length(l) || startsWith(l, "!series_matrix_table_begin")) break; l <- sub("[\r\n]+$", "", l); if (any(startsWith(l, c("!Sample_", "!Series_title", "!Series_geo_accession")))) lines <- c(lines, l) }
  lines
}
header <- function(gse, gpl, hdir) {   # metadata lines of the series matrix (cached in hdir)
  path <- file.path(hdir, sprintf("%s_%s.txt", gse, gpl)); if (file.exists(path)) return(strsplit(readChar(path, file.size(path), useBytes = TRUE), "\n", fixed = TRUE)[[1]])
  url <- matrix_url_sh(gse, gpl); if (is.null(url)) return(NULL)
  lines <- stream_header(url); write_text_bytes(paste(lines, collapse = "\n"), path); lines
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
pick_groups <- function(p) {           # perturbation key (first qualifying in file order) or title fallback; 3 groups of >= 4 samples; control-like group first
  chosen <- NULL; how <- NULL
  for (k in p$keys) if (any(vapply(KEYWORDS, function(w) grepl(w, k, fixed = TRUE), TRUE))) {
    d <- groups_from(vapply(p$per, function(x) if (is.null(x[[k]])) NA_character_ else x[[k]], "")); q <- d[vapply(d, length, 1L) >= 4]
    if (length(q) >= 3) { chosen <- q; how <- paste0("characteristic: ", k); break }
  }
  if (is.null(chosen)) { d <- groups_from(sub(TITLE_TAIL, "", p$title, perl = TRUE)); q <- d[vapply(d, length, 1L) >= 4]; if (length(q) >= 3) { chosen <- q; how <- "title" } }
  if (is.null(chosen)) return(NULL)
  labels <- names(chosen); ctrl <- labels[grepl(CONTROL, gsub("_", " ", labels, fixed = TRUE), perl = TRUE)]; first <- if (length(ctrl)) ctrl[1] else py_sort(labels)[1]; sel <- c(first, head(labels[labels != first], 2))
  list(how = how, n_groups = length(chosen), groups = sel, gsms = lapply(sel, function(l) p$gsm[head(chosen[[l]], 4)]))
}
screen <- function(lines, gpl) {
  p <- parse_header(lines); n <- length(p$gsm)
  if (!(n >= 12 && n <= 36)) return(list(status = sprintf("ineligible: %d samples", n), info = NULL))
  if (!(length(p$plat) == 1 && p$plat == gpl)) return(list(status = sprintf("ineligible: platforms [%s]", paste(sprintf("'%s'", py_sort(p$plat)), collapse = ", ")), info = NULL))
  info <- pick_groups(p); if (is.null(info)) return(list(status = "ineligible: no grouping with >= 3 groups of >= 4 samples", info = NULL)); list(status = "eligible", info = info)
}
all_gsm <- function(lines) { for (l in lines) if (startsWith(l, "!Sample_geo_accession")) return(gsub("^\"+|\"+$", "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])); character(0) }
dev_sample_sets <- function(root) {    # GSM accessions (all samples of the series header) of the 37 analysed development series -> named character vector gsm -> gse; the headers come from the cached series matrices
  ref <- read_jsonl(file.path(P1, "results_per_dataset.jsonl")); out <- character(0)
  for (r in ref) if (identical(r$status, "analysed")) {
    hp <- file.path(root, "bench_data", "headers", sprintf("%s_%s.txt", r$gse, r$platform))
    if (!file.exists(hp) && !file.exists(file.path(root, "bench_data", "matrices", sprintf("%s_%s.txt.gz", r$gse, r$platform)))) { cat("note: no header and no matrix of", r$gse, "(not downloadable); its samples cannot be checked for overlap\n"); next }
    lines <- if (file.exists(hp)) strsplit(readChar(hp, file.size(hp), useBytes = TRUE), "\n", fixed = TRUE)[[1]] else {
      mp <- file.path(root, "bench_data", "matrices", sprintf("%s_%s.txt.gz", r$gse, r$platform)); con <- gzfile(mp, "rt"); on.exit(try(close(con), silent = TRUE), add = TRUE); l <- character(0)
      repeat { x <- readLines(con, n = 1L, warn = FALSE, encoding = "UTF-8"); if (!length(x) || startsWith(x, "!series_matrix_table_begin")) break; if (any(startsWith(x, c("!Sample_", "!Series_title", "!Series_geo_accession")))) l <- c(l, x) }; close(con); l }
    for (g in all_gsm(lines)) if (!(g %in% names(out))) out[g] <- r$gse
  }
  out
}
main <- function() {
  a <- commandArgs(trailingOnly = TRUE); root <- DATA; hdir <- file.path(root, "v2_data", "arrays_heldout", "headers"); dir.create(hdir, showWarnings = FALSE, recursive = TRUE)
  check_only <- "--check-only" %in% a; dl_only <- "--download-only" %in% a; from <- Sys.getenv("WP2_SELECT_FROM", ""); dev <- if (!dl_only && !check_only) dev_sample_sets(root) else character(0)
  p1 <- list(); er <- read.csv(file.path(P1, "eligibility.csv"), stringsAsFactors = FALSE); for (g in PLATFORMS) p1[[g]] <- er$gse[er$platform == g]
  lists <- list()
  if (nzchar(from)) { lists <- fromJSON(from, simplifyVector = FALSE); lists <- lapply(lists, unlist); cat("screening the archived window of", basename(from), "\n") }
  else for (gpl in PLATFORMS) {
    cand <- candidates(gpl); lists[[gpl]] <- cand$accs; same <- identical(cand$accs[1:150], p1[[gpl]])
    cat(sprintf("%s: %d series returned (%d accessions); positions 1-150 identical to part 1: %s\n", gpl, cand$count, length(cand$accs), if (same) "True" else "False"))
    if (!same) stop(sprintf("the re-run query does not reproduce the part 1 list for %s: stop (document and decide before continuing)", gpl))
  }
  if (check_only) return(invisible())
  elig <- list(); frozen <- list()
  for (gpl in PLATFORMS) {
    accs <- if (nzchar(from)) lists[[gpl]] else lists[[gpl]][(WINDOW[1] + 1):min(WINDOW[2], length(lists[[gpl]]))]
    for (k in seq_along(accs)) {
      gse <- accs[k]; pos <- WINDOW[1] + k
      if (dl_only) { tryCatch(header(gse, gpl, hdir), error = function(e) cat("header failed", gse, class(e)[1], "\n")); if (k %% 25 == 0) cat(sprintf("  %s: headers cached up to position %d\n", gpl, pos)); next }
      r <- tryCatch({ lines <- header(gse, gpl, hdir)
        if (is.null(lines)) list(status = "ineligible: no series matrix for the platform", info = NULL)
        else { s <- screen(lines, gpl); if (!is.null(s$info)) { g_all <- all_gsm(lines); hit <- g_all[g_all %in% names(dev)]; shared <- py_sort(unique(dev[hit])); if (length(shared)) s <- list(status = sprintf("excluded: shares %d samples with development series %s", length(hit), paste(shared, collapse = ",")), info = NULL) }; s }
      }, error = function(e) list(status = sprintf("ineligible: error %s", class(e)[1]), info = NULL))
      info <- r$info
      elig[[length(elig) + 1]] <- list(platform = gpl, order = as.integer(pos), gse = gse, status = r$status, how = if (is.null(info)) "" else info$how, n_groups = if (is.null(info)) "" else as.integer(info$n_groups), groups = if (is.null(info)) "" else paste(info$groups, collapse = "|"))
      if (!is.null(info)) frozen[[length(frozen) + 1]] <- list(platform = gpl, order = as.integer(pos), gse = gse, how = info$how, n_groups = as.integer(info$n_groups), groups = paste(info$groups, collapse = "|"), gsms = paste(unlist(info$gsms), collapse = ";"))
      if (k %% 25 == 0) cat(sprintf("  %s: positions %d-%d screened, %d eligible\n", gpl, WINDOW[1] + 1, pos, sum(vapply(frozen, function(f) f$platform == gpl, TRUE))))
    }
  }
  if (dl_only) { cat("headers cached\n"); return(invisible()) }
  dictwriter_csv(elig, names(elig[[1]]), file.path(OUT, "eligibility_heldout.csv")); dictwriter_csv(frozen, names(frozen[[1]]), file.path(OUT, "frozen_heldout_arrays.csv"))
  h <- sha256_file(file.path(OUT, "frozen_heldout_arrays.csv")); write_text_bytes(sprintf("%s *wp2_benchmark2/frozen_heldout_arrays.csv\nWP2 held-out microarray list (R port) frozen: %s\n", h, format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), file.path(OUT, "PROTOCOL_v2_SHA256_R.txt"))
  cat("eligible per platform:", paste(sprintf("%s=%d", PLATFORMS, vapply(PLATFORMS, function(g) sum(vapply(frozen, function(f) f$platform == g, TRUE)), 1L)), collapse = ", "), "| sha256", h, "\n")
  st <- table(vapply(elig, function(e) e$status, "")); cat("status counts:", paste(sprintf("'%s': %d", names(st), st), collapse = ", "), "\n")
}
if (sys.nframe() == 0L) main()
