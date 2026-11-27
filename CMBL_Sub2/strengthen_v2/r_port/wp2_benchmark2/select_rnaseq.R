## select_rnaseq.R - R port of wp2_benchmark2/select_rnaseq.py: protocol v2 section 3.5, held-out RNA-seq (WP2 deliverable E): candidate lists, pattern check of the NCBI-generated count download,
## metadata screening and freezing. Query (verbatim from the protocol): ("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[ETYP]
## AND 12:36[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine). The set has more than 30000 records, above the 10000 that
## one esearch call can return, so it is retrieved with the same query split by n_samples (12, 13, ... 36; the union is the same set). A GEO series record has UID 200000000 + GSE number, so accessions come
## from the UIDs; the organism comes from running the query once per organism; series found under both organisms are excluded as mixed. Per organism the accessions are sorted numerically and shuffled with
## random.Random(2026); the first 300 are screened on metadata by the rules of select_heldout.R (12-36 samples, ONE platform for the whole series, same perturbation keys / title fallback, 3 groups of 4,
## control-like group first).
## Steps (Rscript select_rnaseq.R <step>):
##   lists          build and save the shuffled lists (<OUT>/compact/rnaseq_lists.json); a live query: GEO has grown since 5 Nov 2026, so the lists differ from the frozen ones (compact/rnaseq_lists.json of the archive)
##   pattern-check  status and content type (never the content) of the NCBI-generated count URL for 20 series per organism from positions 301 and later of the shuffled lists, plus the annotation file;
##                  written to rnaseq_pattern_check.csv (a statement about the NCBI service on the day it is run)
##   screen [human,mouse]   screen the first 300 per organism of the lists in WP2_RNASEQ_LISTS (default: the archived frozen lists) on metadata, write eligibility_heldout_rnaseq.csv and frozen_heldout_rnaseq.csv
##                  and the SHA-256 to <OUT>/PROTOCOL_v2_SHA256_R.txt BEFORE any count is read
if (!exists("HERE") || !file.exists(file.path(HERE, "wp2_common.R"))) { HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, "..")) }
source(file.path(HERE, "select_heldout.R"))
TERMS <- "\"Expression profiling by high throughput sequencing\"[DataSet Type] AND gse[ETYP] AND %d:%d[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine)"
ORGS <- c(human = "Homo sapiens", mouse = "Mus musculus"); N_SCREEN <- 300L
ANNOT_URL <- c(human = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz", mouse = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Mouse.GRCm39.annot.tsv.gz")
F_LISTS <- Sys.getenv("WP2_RNASEQ_LISTS", file.path(ARC, "compact", "rnaseq_lists.json"))
uids <- function(org_name, lo, hi) {
  term <- paste0("\"", org_name, "\"[Organism] AND ", sprintf(TERMS, lo, hi)); j <- fromJSON(http_text(paste0(EU, "esearch.fcgi?db=gds&term=", URLencode(term, reserved = TRUE), "&retmax=10000&retmode=json")), simplifyVector = FALSE)$esearchresult
  n <- as.integer(j$count); ids <- unlist(j$idlist); if (n > length(ids)) stop(sprintf("n_samples %d-%d for %s has %d records, more than one call returns: split further", lo, hi, org_name, n)); ids
}
build_lists <- function() {
  out <- list(); raw <- list()
  for (key in names(ORGS)) { s <- numeric(0); for (n in 12:36) s <- union(s, as.numeric(uids(ORGS[[key]], n, n)) - 200000000); raw[[key]] <- s; cat(key, "records", length(s), "\n") }
  mixed <- intersect(raw$human, raw$mouse)
  for (key in names(ORGS)) { accs <- sort(setdiff(raw[[key]], mixed)); lst <- sprintf("GSE%d", accs); out[[key]] <- py_shuffle(py_random(2026), lst) }
  out$n_records <- lapply(raw, function(v) as.integer(length(v))); out$n_mixed <- as.integer(length(mixed)); out$queried <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  dir.create(file.path(OUT, "compact"), showWarnings = FALSE, recursive = TRUE); write_json_typed(out, file.path(OUT, "compact", "rnaseq_lists.json"), indent = NULL); cat("lists saved; mixed-organism series excluded:", length(mixed), "\n"); out
}
probe <- function(url, nbytes = 2L) {      # status, content type and (if gzip) the first two bytes only: the content is never read
  h <- new_handle(useragent = USER_AGENT, timeout = 60); throttle()
  tryCatch({
    con <- curl(url, open = "rb", handle = h); on.exit(try(close(con), silent = TRUE)); d <- handle_data(h); hd <- parse_headers_list(d$headers); ct <- if (is.null(hd[["content-type"]])) "" else hd[["content-type"]]
    first <- readBin(con, "raw", if (grepl("html", ct)) 60L else nbytes)
    kind <- if (length(first) >= 2 && identical(first[1:2], as.raw(c(0x1f, 0x8b)))) "gzip file" else if (any(first[1:min(60, length(first))] == as.raw(0x3c))) "HTML (challenge or error page)" else "other"
    list(status = as.integer(d$status_code), content_type = ct, kind = kind)
  }, error = function(e) { m <- regmatches(conditionMessage(e), regexpr("[0-9]{3}", conditionMessage(e))); list(status = if (length(m)) as.integer(m) else -1L, content_type = "", kind = if (length(m)) "HTTP error" else class(e)[1]) })
}
pattern_check <- function(lists) {
  rows <- list()
  for (key in names(ORGS)) {
    for (i in seq_along(lists[[key]][(N_SCREEN + 1):(N_SCREEN + 20)])) { acc <- lists[[key]][N_SCREEN + i]; p <- probe(sprintf(COUNT_URL, acc, acc, RNASEQ_BUILD[[key]])); rows[[length(rows) + 1]] <- c(list(organism = key, position = as.integer(N_SCREEN + i), what = paste("counts", acc)), p) }
    p <- probe(ANNOT_URL[[key]]); rows[[length(rows) + 1]] <- c(list(organism = key, position = 0L, what = "annotation"), p)
  }
  dictwriter_csv(rows, names(rows[[1]]), file.path(OUT, "rnaseq_pattern_check.csv"))
  st <- vapply(rows, function(r) r$status, 1L); kd <- vapply(rows, function(r) r$kind, ""); cat(sprintf("pattern check: %d probes; plain gzip files %d, HTML pages %d, HTTP 404 %d\n", length(rows), sum(kd == "gzip file"), sum(grepl("HTML", kd)), sum(st == 404L))); rows
}
header_rna <- function(gse, hdir) {          # metadata lines; NULL if the series has no single series-matrix file (several platforms)
  path <- file.path(hdir, paste0(gse, ".txt"))
  if (file.exists(path)) { t <- readChar(path, file.size(path), useBytes = TRUE); return(if (identical(t, "MULTI")) NULL else strsplit(t, "\n", fixed = TRUE)[[1]]) }
  n <- as.integer(substring(gse, 4)); folder <- if (n >= 1000) sprintf("GSE%dnnn", n %/% 1000) else "GSEnnn"; url <- sprintf("%s%s/%s/matrix/%s_series_matrix.txt.gz", FTP_SERIES, folder, gse, gse)
  throttle(); st <- tryCatch(curl_fetch_memory(url, new_handle(useragent = USER_AGENT, timeout = 90, nobody = TRUE))$status_code, error = function(e) NA)
  if (isTRUE(st == 404L)) { write_text_bytes("MULTI", path); return(NULL) }
  if (!isTRUE(st == 200L)) stop(structure(class = c("HTTPError", "error", "condition"), list(message = paste("HTTP", st), call = NULL)))
  lines <- stream_header(url); write_text_bytes(paste(lines, collapse = "\n"), path); lines
}
screen_rna <- function(lines) {
  p <- parse_header(lines); n <- length(p$gsm)
  if (!(n >= 12 && n <= 36)) return(list(status = sprintf("ineligible: %d samples", n), info = NULL))
  if (length(p$plat) != 1) return(list(status = sprintf("ineligible: platforms [%s]", paste(sprintf("'%s'", py_sort(p$plat)), collapse = ", ")), info = NULL))
  info <- pick_groups(p); if (is.null(info)) return(list(status = "ineligible: no grouping with >= 3 groups of >= 4 samples", info = NULL)); list(status = "eligible", info = info)
}
screen_lists <- function(lists, orgs = names(ORGS)) {
  root <- DATA; hdir <- file.path(root, "v2_data", "rnaseq", "headers"); dir.create(hdir, showWarnings = FALSE, recursive = TRUE); dev <- dev_sample_sets(root); elig <- list(); frozen <- list()
  arr <- Sys.getenv("WP2_FROZEN_ARRAYS", file.path(ARC, "frozen_heldout_arrays.csv")); held <- character(0)
  if (file.exists(arr)) for (r in split(read.csv(arr, stringsAsFactors = FALSE), seq_len(nrow(read.csv(arr))))) {
    p <- file.path(root, "v2_data", "arrays_heldout", "headers", sprintf("%s_%s.txt", r$gse, r$platform))
    if (file.exists(p)) for (g in all_gsm(strsplit(readChar(p, file.size(p), useBytes = TRUE), "\n", fixed = TRUE)[[1]])) if (!(g %in% names(held))) held[g] <- r$gse
  }
  for (key in orgs) for (pos in seq_len(min(N_SCREEN, length(lists[[key]])))) {
    gse <- lists[[key]][pos]
    r <- tryCatch({ lines <- header_rna(gse, hdir)
      if (is.null(lines)) list(status = "ineligible: several platforms (no single series matrix)", info = NULL)
      else { s <- screen_rna(lines); if (!is.null(s$info)) { g_all <- all_gsm(lines); shared <- py_sort(unique(c(dev[g_all[g_all %in% names(dev)]], held[g_all[g_all %in% names(held)]]))); if (length(shared)) s <- list(status = sprintf("excluded: shares samples with series %s", paste(shared, collapse = ",")), info = NULL) }; s }
    }, error = function(e) list(status = sprintf("ineligible: error %s", class(e)[1]), info = NULL))
    info <- r$info
    elig[[length(elig) + 1]] <- list(organism = key, order = as.integer(pos), gse = gse, status = r$status, how = if (is.null(info)) "" else info$how, n_groups = if (is.null(info)) "" else as.integer(info$n_groups), groups = if (is.null(info)) "" else paste(info$groups, collapse = "|"))
    if (!is.null(info)) frozen[[length(frozen) + 1]] <- list(organism = key, order = as.integer(pos), gse = gse, how = info$how, n_groups = as.integer(info$n_groups), groups = paste(info$groups, collapse = "|"), gsms = paste(unlist(info$gsms), collapse = ";"))
    if (pos %% 25 == 0) cat(sprintf("  %s: %d screened, %d eligible\n", key, pos, sum(vapply(frozen, function(f) f$organism == key, TRUE))))
  }
  dictwriter_csv(elig, names(elig[[1]]), file.path(OUT, "eligibility_heldout_rnaseq.csv")); dictwriter_csv(frozen, names(frozen[[1]]), file.path(OUT, "frozen_heldout_rnaseq.csv"))
  h <- sha256_file(file.path(OUT, "frozen_heldout_rnaseq.csv")); write_text_bytes(sprintf("%s *wp2_benchmark2/frozen_heldout_rnaseq.csv\nWP2 held-out RNA-seq list (R port) frozen BEFORE any count was read: %s\n", h, format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), file.path(OUT, "PROTOCOL_v2_SHA256_R_rnaseq.txt"))
  cat("eligible per organism:", paste(sprintf("%s=%d", orgs, vapply(orgs, function(k) sum(vapply(frozen, function(f) f$organism == k, TRUE)), 1L)), collapse = ", "), "| sha256", h, "\n")
}
if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE); step <- if (length(a)) a[1] else "lists"
  if (step == "lists") build_lists() else if (step == "pattern-check") pattern_check(lapply(fromJSON(F_LISTS, simplifyVector = FALSE), function(v) if (is.list(v) && length(v) > 5) unlist(v) else v))
  else if (step == "screen") { L <- lapply(fromJSON(F_LISTS, simplifyVector = FALSE)[names(ORGS)], unlist); screen_lists(L, if (length(a) >= 2) strsplit(a[2], ",", fixed = TRUE)[[1]] else names(ORGS)) }
  else stop("lists | pattern-check | screen [human,mouse]")
}
