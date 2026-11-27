## perturbation_scout.R - R port of wp4_islets/perturbation_scout.py: protocol 4.6, step 1: run the two fixed GEO DataSets queries (PRESPEC_WP4_choices.md section D) once through the NCBI E-utilities (db=gds) and list EVERY candidate
## series with title, organism, assay type, number of samples, platform and the first 400 characters of the summary, plus an automatic stem-cell-design flag (regex on title + summary). Nothing is dropped here; the keep / drop
## decision with its reason is written afterwards from title, summary and design text only (no expression value is read at this step; perturbation_decisions.R).
## Plain HTTPS, at most 3 requests per second, plain tool User-Agent, nothing about the user is sent. Output (default results_R/wp4/perturbation_scout): candidates.csv, queries.json. LIVE QUERY: GEO grows, so a run today lists the
## series of the archived run of 5-6 Nov 2026 and possibly new ones; the analyses read the ARCHIVED list (wp4_islets/perturbation/candidates.csv), which is frozen. compare_scout.R lists the differences.
## Usage: Rscript perturbation_scout.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "wp4_http.R")); suppressMessages(library(jsonlite))
OUT <- Sys.getenv("WP4_SCOUT_OUT", file.path(R_PORT, "results_R", "wp4", "perturbation_scout")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE); ATT <- file.path(dirname(OUT), "attempts_R.csv")
CELL <- "(islet OR islets OR \"beta cell\" OR \"beta cells\" OR \"beta-cell\" OR INS-1 OR INS1 OR MIN6 OR \"pancreatic beta\")"
QUERIES <- list(notch = paste0("gse[ETYP] AND (Notch OR \"gamma-secretase\" OR \"gamma secretase\" OR DAPT) AND ", CELL),
                wnt = paste0("gse[ETYP] AND (DKK1 OR Dickkopf OR \"Wnt inhibitor\" OR \"beta-catenin\" OR TCF7L2 OR TCF4 OR \"Wnt signaling\" OR \"Wnt signalling\") AND ", CELL))
EU <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; STEM <- "differentiat|iPSC|\\bESC\\b|embryonic stem|pluripoten|stem cell|stem-cell|progenitor|reprogram"
na1 <- function(x) if (is.null(x) || length(x) == 0) NA else x
rows <- list(); qlog <- list()
for (name in names(QUERIES)) {
  term <- QUERIES[[name]]; url <- paste0(EU, "esearch.fcgi?", url_query(list(db = "gds", term = term, retmax = 500, retmode = "json")))
  j <- fromJSON(get_text(url), simplifyVector = FALSE)$esearchresult; ids <- as.character(unlist(j$idlist)); qlog[[name]] <- list(term = term, count = as.integer(j$count), ids = as.list(ids), translation = na1(j$querytranslation))
  attempts_record(ATT, "NCBI GEO DataSets", paste("query", name), substr(url, 1, 250), NA, "E-utilities esearch db=gds", "queried", paste("count", j$count)); cat(name, "count", j$count, "\n")
  for (k in if (length(ids)) seq(1, length(ids), by = 100) else integer(0)) {
    u <- paste0(EU, "esummary.fcgi?", url_query(list(db = "gds", id = paste(ids[k:min(k + 99, length(ids))], collapse = ","), retmode = "json"))); res <- fromJSON(get_text(u), simplifyVector = FALSE)$result
    for (uid in unlist(res$uids)) {
      r <- res[[uid]]; acc <- na1(r$accession); if (is.na(acc)) acc <- ""
      d <- if (is.null(rows[[acc]])) list(accession = acc, entrytype = na1(r$entrytype), title = na1(r$title), organism = na1(r$taxon), assay = na1(r$gdstype), n_samples = na1(r$n_samples), platform = na1(r$gpl), queries = character(0),
                                          summary = gsub("\n", " ", if (is.null(r$summary)) "" else r$summary, fixed = TRUE)) else rows[[acc]]
      d$queries <- union(d$queries, name); rows[[acc]] <- d
    }
  }
}
D <- do.call(rbind, lapply(rows, function(d) data.frame(accession = d$accession, entrytype = as.character(d$entrytype), title = as.character(d$title), organism = as.character(d$organism), assay = as.character(d$assay), n_samples = suppressWarnings(as.numeric(d$n_samples)),
                                                       platform = as.character(d$platform), queries = paste(sort(d$queries, method = "radix"), collapse = "+"), summary = d$summary, stringsAsFactors = FALSE)))
D$stem_cell_flag <- grepl(STEM, paste(ifelse(is.na(D$title), "", D$title), ifelse(is.na(D$summary), "", D$summary)), perl = TRUE, ignore.case = TRUE); D$summary <- substr(D$summary, 1, 400)
D <- D[order(D$accession, method = "radix"), ]; rownames(D) <- NULL; if (!anyNA(D$n_samples)) D$n_samples <- as.integer(D$n_samples)
py_write_csv_repr(D, file.path(OUT, "candidates.csv"), eol = "\r\n"); write_json_typed(qlog, file.path(OUT, "queries.json"))
cat(nrow(D), "unique records;", paste(names(table(D$entrytype)), as.integer(table(D$entrytype)), collapse = ", "), "\n")
