## wider_search.R - R port of wp4_islets/wider_search.py: B1 of PRESPEC_WP4b_positive_control_and_wider_search.md (queries fixed before this ran): three sources, each run once, every returned record listed with the source(s) that returned it.
##  S1  NCBI GEO DataSets (db=gds):  gse[ETYP] AND (NOTCH OR WNT) AND (CELL)
##  S2  EBI BioStudies, ArrayExpress collection (api/v1/arrayexpress/search), the same three term groups as one query (E-GEOD-n records are the GEO series GSEn and are merged with S1)
##  S3  PubMed: (CELL) AND (NOTCH OR WNT) AND (transcriptom* OR "RNA-seq" OR "RNA sequencing" OR microarray OR "expression profiling" OR "gene expression"), PMIDs capped at 1500, linked to GEO by elink (pubmed_gds)
## Plain HTTPS, at most 3 requests per second, plain tool User-Agent, nothing about the user is sent; a refusal or challenge page is reported, not worked around. No expression value is read.
## Output (default results_R/wp4/wider_search_scout): queries.json, candidates_all.csv (accession, sources, pmids, entry type, title, organism, assay, n_samples, platform, summary). LIVE QUERY: the databases grow, so a run
## today returns the series of the archived run of 6 Nov 2026 and possibly new ones; the analyses read the ARCHIVED, frozen list (wp4_islets/wider_search/candidates_all.csv, frozen_candidates.csv). compare_scout.R lists the differences.
## (The Python script fetched the esummary records of the S3 links twice, once without using them; the R port fetches them once.)   Usage: Rscript wider_search.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "wp4_http.R")); suppressMessages(library(jsonlite))
OUT <- Sys.getenv("WP4_WS_SCOUT_OUT", file.path(R_PORT, "results_R", "wp4", "wider_search_scout")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE); ATT <- file.path(dirname(OUT), "attempts_R.csv")
CELL <- "(islet OR islets OR \"beta cell\" OR \"beta cells\" OR \"beta-cell\" OR \"INS-1\" OR INS1 OR \"INS-1E\" OR MIN6 OR EndoC OR betaTC OR \"BRIN-BD11\" OR RINm5F OR insulinoma OR \"pancreatic beta\")"
NOTCH <- "(Notch OR \"gamma-secretase\" OR \"gamma secretase\" OR DAPT OR LY411575 OR RO4929097 OR \"compound E\" OR \"MK-0752\" OR RBPJ OR Jagged OR Jag1 OR Dll1 OR Dll4 OR Notch1 OR Notch2)"
WNT <- "(Wnt OR DKK1 OR Dickkopf OR \"beta-catenin\" OR CTNNB1 OR TCF7L2 OR TCF4 OR IWP OR \"IWP-2\" OR IWR OR XAV939 OR LGK974 OR \"ICG-001\" OR CHIR99021 OR GSK3 OR Wnt3a OR Axin OR SFRP OR Frizzled)"
PERT <- sprintf("(%s OR %s)", NOTCH, WNT); Q1 <- sprintf("gse[ETYP] AND %s AND %s", PERT, CELL); Q2 <- sprintf("%s AND %s", CELL, PERT)
Q3 <- sprintf("%s AND %s AND (transcriptom* OR \"RNA-seq\" OR \"RNA sequencing\" OR microarray OR \"expression profiling\" OR \"gene expression\")", CELL, PERT)
EU <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; BS <- "https://www.ebi.ac.uk/biostudies/api/v1/"
qlog <- list(date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), terms = list(CELL = CELL, NOTCH = NOTCH, WNT = WNT), S1 = list(query = Q1), S2 = list(query = Q2), S3 = list(query = Q3, cap = 1500L))
cand <- new.env(); jget <- function(url) fromJSON(get_text(url), simplifyVector = FALSE); blank <- function(v) is.null(v) || length(v) == 0 || is.na(v[1]) || identical(v, "")
add <- function(acc, src, ...) {
  d <- if (exists(acc, envir = cand, inherits = FALSE)) get(acc, envir = cand) else list(accession = acc, sources = character(0), pmids = character(0), entrytype = NA, title = NA, organism = NA, assay = NA, n_samples = NA, platform = NA, summary = NA)
  d$sources <- union(d$sources, src); kw <- list(...)
  for (k in names(kw)) { v <- kw[[k]]; if (k == "pmid") d$pmids <- union(d$pmids, as.character(v)) else if (!blank(v) && blank(d[[k]])) d[[k]] <- v }
  assign(acc, d, envir = cand)
}
gds_records <- function(uids) {      # esummary of GEO DataSets uids -> list(uid, acc, kw) for series (GSE); curated datasets (GDS) are mapped to their series
  out <- list()
  for (k in if (length(uids)) seq(1, length(uids), by = 100) else integer(0)) {
    r <- jget(paste0(EU, "esummary.fcgi?", url_query(list(db = "gds", id = paste(uids[k:min(k + 99, length(uids))], collapse = ","), retmode = "json"))))$result
    for (uid in unlist(r$uids)) {
      x <- r[[uid]]; acc <- if (is.null(x$accession)) "" else x$accession; if (startsWith(acc, "GDS") && !blank(x$gse)) acc <- paste0("GSE", x$gse)
      if (!startsWith(acc, "GSE")) next
      out[[length(out) + 1]] <- list(uid = uid, acc = acc, kw = list(entrytype = x$entrytype, title = x$title, organism = x$taxon, assay = x$gdstype, n_samples = x$n_samples, platform = x$gpl, summary = substr(gsub("\n", " ", if (is.null(x$summary)) "" else x$summary, fixed = TRUE), 1, 700)))
    }
  }
  out
}
## ---- S1 GEO
u <- paste0(EU, "esearch.fcgi?", url_query(list(db = "gds", term = Q1, retmax = 5000, retmode = "json"))); j <- jget(u)$esearchresult; ids1 <- as.character(unlist(j$idlist)); qlog$S1$count <- as.integer(j$count); qlog$S1$ids <- length(ids1)
for (rec in gds_records(ids1)) do.call(add, c(list(rec$acc, "S1"), rec$kw))
attempts_record(ATT, "NCBI GEO DataSets", "wider search S1", "esearch.fcgi (db=gds)", NA, "E-utilities", "queried", paste("count", j$count)); cat("S1 GEO: count", j$count, "| series", sum(vapply(ls(cand), function(a) "S1" %in% get(a, envir = cand)$sources, TRUE)), "\n")
## ---- S3 PubMed -> GEO
u <- paste0(EU, "esearch.fcgi?", url_query(list(db = "pubmed", term = Q3, retmax = 5000, retmode = "json"))); j <- jget(u)$esearchresult; allp <- as.character(unlist(j$idlist)); pm <- allp[order(as.numeric(allp))]; pm <- pm[seq_len(min(1500L, length(pm)))]
qlog$S3$count <- as.integer(j$count); qlog$S3$pmids <- length(pm); link <- new.env()
for (k in if (length(pm)) seq(1, length(pm), by = 50) else integer(0)) {
  batch <- pm[k:min(k + 49, length(pm))]; r <- jget(paste0(EU, "elink.fcgi?dbfrom=pubmed&db=gds&linkname=pubmed_gds&retmode=json&", paste0("id=", batch, collapse = "&")))
  for (lset in r$linksets) for (db in lset$linksetdbs) if (identical(db$linkname, "pubmed_gds")) for (lu in unlist(db$links)) { key <- as.character(lu); assign(key, union(if (exists(key, envir = link, inherits = FALSE)) get(key, envir = link) else character(0), as.character(unlist(lset$ids))), envir = link) }
}
uids <- sort(ls(link), method = "radix"); qlog$S3$linked_gds_uids <- length(uids)
for (rec in gds_records(uids)) { do.call(add, c(list(rec$acc, "S3"), rec$kw)); for (p in get(rec$uid, envir = link)) add(rec$acc, "S3", pmid = p) }
attempts_record(ATT, "PubMed -> GEO", "wider search S3", "esearch.fcgi (db=pubmed) + elink pubmed_gds", NA, "E-utilities", "queried", sprintf("PMIDs %d, linked GEO uids %d", length(pm), length(uids)))
cat("S3 PubMed: count", j$count, "| PMIDs", length(pm), "| linked GEO records", length(uids), "\n")
## ---- S2 BioStudies / ArrayExpress
hits <- list(); page <- 1L; total <- NULL
repeat {
  r <- tryCatch(jget(paste0(BS, "arrayexpress/search?", url_query(list(query = Q2, pageSize = 100, page = page)))), error = function(e) e)
  if (inherits(r, "error")) { qlog$S2$error <- substr(conditionMessage(r), 1, 300); cat("S2 FAILED (reported, not worked around):", substr(conditionMessage(r), 1, 200), "\n"); break }
  total <- r$totalHits; hits <- c(hits, r$hits); if (length(hits) >= (if (is.null(total)) 0 else total) || !length(r$hits)) break
  page <- page + 1L
}
qlog$S2["count"] <- list(if (is.null(total)) NULL else as.integer(total)); qlog$S2$hits <- length(hits); ae <- character(0)
for (h in hits) { a <- h$accession; acc <- if (startsWith(a, "E-GEOD-")) paste0("GSE", substring(a, nchar("E-GEOD-") + 1L)) else a; add(acc, "S2", title = h$title); if (!startsWith(acc, "GSE")) ae <- c(ae, acc) }
for (a in ae) tryCatch({      # non-GEO ArrayExpress studies: description from the study record
  s <- jget(paste0(BS, "studies/", a)); attrs <- list(); for (x in s$section$attributes) if (!is.null(x$name)) attrs[[x$name]] <- x$value
  title <- attrs[["Title"]]; if (blank(title)) { title <- NULL; for (x in s$attributes) if (identical(x$name, "Title")) { title <- x$value; break } }
  desc <- attrs[["Description"]]; if (blank(desc)) desc <- attrs[["Abstract"]]; if (blank(desc)) desc <- ""; assay <- attrs[["Study type"]]; if (blank(assay)) assay <- attrs[["Assay"]]
  add(a, "S2", entrytype = "ArrayExpress", title = title, summary = substr(gsub("\n", " ", desc, fixed = TRUE), 1, 700), organism = attrs[["Organism"]], assay = assay)
}, error = function(e) add(a, "S2", summary = sprintf("(study record not retrievable: %s)", substr(conditionMessage(e), 1, 80))))
need <- Filter(function(a) startsWith(a, "GSE") && (blank(get(a, envir = cand)$title) || blank(get(a, envir = cand)$organism)), ls(cand))      # GEO series reached only through ArrayExpress (E-GEOD) need their GEO metadata
for (a in need) { ids <- as.character(unlist(jget(paste0(EU, "esearch.fcgi?", url_query(list(db = "gds", term = paste0(a, "[ACCN] AND gse[ETYP]"), retmode = "json"))))$esearchresult$idlist)); for (rec in gds_records(ids)) do.call(add, c(list(rec$acc, "S2"), rec$kw)) }
attempts_record(ATT, "EBI BioStudies", "wider search S2 (arrayexpress collection)", paste0(BS, "arrayexpress/search"), NA, "BioStudies API", "queried", paste("hits", if (is.null(total)) "None" else total)); cat("S2 BioStudies: hits", if (is.null(total)) "None" else total, "| non-GEO studies", length(ae), "\n")
accs <- sort(ls(cand), method = "radix"); g <- function(d, k) { v <- d[[k]]; if (blank(v)) NA else v }
D <- do.call(rbind, lapply(accs, function(a) { d <- get(a, envir = cand); data.frame(accession = a, sources = paste(sort(d$sources, method = "radix"), collapse = "+"), pmids = paste(d$pmids[order(as.numeric(d$pmids))], collapse = ";"), entrytype = as.character(g(d, "entrytype")), title = as.character(g(d, "title")),
                                                                                  organism = as.character(g(d, "organism")), assay = as.character(g(d, "assay")), n_samples = suppressWarnings(as.numeric(g(d, "n_samples"))), platform = as.character(g(d, "platform")), summary = as.character(g(d, "summary")), stringsAsFactors = FALSE) }))
if (!anyNA(D$n_samples)) D$n_samples <- as.integer(D$n_samples)
py_write_csv_repr(D, file.path(OUT, "candidates_all.csv"), eol = "\r\n"); tb <- table(D$sources); tb <- tb[order(-as.integer(tb), match(names(tb), unique(D$sources)))]
qlog$unique_records <- nrow(D); qlog$by_source_combination <- setNames(as.list(as.integer(tb)), names(tb)); write_json_typed(qlog, file.path(OUT, "queries.json")); cat(nrow(D), "unique records;", paste(names(tb), as.integer(tb), collapse = ", "), "\n")
