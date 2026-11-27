## wp12_annotations.R - R port of wp12_database_annotations/wp12_annotations.py: WP12 (9 Nov 2026): annotations from the updated databases, as pre-specified in PRESPEC_WP12_annotations.md.
##   part 1  hub-gene target annotation from the DrugCentral file of 25 Sep 2026 (record counts, drugs, TDL)
##   part 2  score-quartile classes of the MLSS combinations (153 reported, 171 of the 25 Sep 2026 release)
##   part 3  regulatory / market status of the 19 MLSS drugs: DrugCentral approval-list membership and openFDA Drugs@FDA (retrieved on the day of the run)
## Usage: Rscript wp12_annotations.R retrieve   (network: openFDA drugsfda, cached in <OUT>/raw/ under the file name sha1(url).json, safe to repeat; the cache of the archived run is read first)
##        Rscript wp12_annotations.R analyze    (offline: reads saved files only; writes <OUT>/*.csv and summary.json)
## Environment: WP12_OUT (default results_R/wp12), WP12_W10 (folder with the WP10 results: approval_tiers_old_vs_new.csv, drugcentral_summary.json, mlss_comparison.csv, mlss_new_all_combinations.csv; default the
## archived wp10_versions/results), WP_PROJECT (project root holding Database/ and Analysis_outputs/), WP12_RAW_ARCHIVE (archived cache of openFDA responses). Read only with respect to every other folder.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- HERE
source(file.path(R_PORT, "pycompat.R")); suppressMessages({ library(jsonlite); library(data.table); library(digest); library(curl) })
options(py_json_nan = TRUE)
ARCH12 <- normalizePath(file.path(R_PORT, "..", "wp12_database_annotations")); P <- Sys.getenv("WP_PROJECT", normalizePath(file.path(R_PORT, "..", "..", "..")))
W10 <- Sys.getenv("WP12_W10", normalizePath(file.path(R_PORT, "..", "wp10_versions", "results"))); RES <- Sys.getenv("WP12_OUT", file.path(R_PORT, "results_R", "wp12")); RAW <- file.path(RES, "raw"); dir.create(RAW, showWarnings = FALSE, recursive = TRUE)
RAW_ARCH <- Sys.getenv("WP12_RAW_ARCHIVE", file.path(ARCH12, "results", "raw"))
BASE <- "https://api.fda.gov/drug/drugsfda.json"; HUBS <- c("GLP1R", "TCF7L2", "CACNA1D", "PTBP1"); CONTEXT <- c("WNT5A", "CACNA1C"); MIN_GAP <- 0.5; .last <- 0
quote_plus <- function(s) {          # urllib.parse.quote_plus: letters, digits and _.-~ stay, space -> +, everything else percent-encoded (UTF-8 bytes, upper-case hex)
  vapply(s, function(x) { b <- charToRaw(enc2utf8(x)); paste(vapply(b, function(by) { ch <- rawToChar(by); if (grepl("[A-Za-z0-9_.~-]", ch)) ch else if (ch == " ") "+" else sprintf("%%%02X", as.integer(by)) }, ""), collapse = "") }, "", USE.NAMES = FALSE)
}
url_for <- function(name, skip = 0, limit = 100) paste0(BASE, "?search=", quote_plus(sprintf("products.active_ingredients.name:\"%s\"", name)), sprintf("&limit=%d&skip=%d", limit, skip))
cache_name <- function(url) paste0(digest(url, algo = "sha1", serialize = FALSE), ".json")
fetch <- function(url, tries = 8) {
  for (d in c(RAW, RAW_ARCH)) { p <- file.path(d, cache_name(url)); if (file.exists(p)) return(fromJSON(p, simplifyVector = FALSE)$response) }
  path <- file.path(RAW, cache_name(url)); delay <- 3
  for (k in seq_len(tries)) {
    gap <- as.numeric(Sys.time()) - .last; if (gap < MIN_GAP) Sys.sleep(MIN_GAP - gap); .last <<- as.numeric(Sys.time()); err <- NULL; resp <- NULL
    r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = "research-screen/1.0", timeout = 120)), error = function(e) { err <<- paste0(class(e)[1], ": ", substr(conditionMessage(e), 1, 100)); NULL })
    if (!is.null(r)) {
      if (r$status_code == 200) resp <- fromJSON(rawToChar(r$content), simplifyVector = FALSE)
      else if (r$status_code == 404) resp <- list(meta = list(results = list(total = 0L, skip = 0L, limit = 100L)), results = list(), note = "HTTP 404: no match")
      else { err <- paste("HTTP", r$status_code); if (!(r$status_code %in% c(429, 500, 502, 503, 504))) stop(err) }
    }
    if (is.null(err)) { write_json_typed(list(url = url, retrieved = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), response = resp), path, indent = NULL); return(resp) }
    cat(sprintf("  retry %d/%d after %s\n", k, tries, err)); Sys.sleep(delay)
  }
  stop("failed: ", url)
}
read_csv_dict <- function(path) read.csv(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = "UTF-8-BOM")
drug_list <- function() read_csv_dict(file.path(W10, "approval_tiers_old_vs_new.csv"))
query_name <- function(drug) { d <- tolower(drug); ifelse(startsWith(d, "(s)-"), substring(d, 5), d) }
listing <- function(name) {            # all applications for an active-ingredient name (pages of 100, at most 6 pages)
  out <- list(); skip <- 0
  repeat { r <- fetch(url_for(name, skip)); total <- if (is.null(r$meta$results$total)) 0L else r$meta$results$total; out <- c(out, r$results); skip <- skip + 100; if (skip >= total || skip >= 600) break }
  list(total = as.integer(total), res = out)
}
retrieve <- function() {
  names_ <- py_sort(unique(query_name(drug_list()$drug))); cat(sprintf("queries for %d names\n", length(names_)))
  for (n in names_) { l <- listing(n); cat(sprintf("  %-14s total applications %d (retrieved %d)\n", n, l$total, length(l$res))) }
}
# --------------------------------------------------------------------------- analysis
read_dc <- function(kind) {            # kind 'reported': Analysis_outputs/R_17genes/drug.target.interaction.csv; kind 'new': Database/drug.target.interaction.tsv.gz (dump of 25 Sep 2026)
  if (kind == "reported") read.csv(file.path(P, "CMBL_Sub2", "R_17genes", "drug.target.interaction.csv"), stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = "UTF-8-BOM", quote = "\"")
  else as.data.frame(fread(cmd = paste("gzip -dc", shQuote(normalizePath(file.path(P, "Database", "drug.target.interaction.tsv.gz")))), sep = "\t", quote = "\"", colClasses = "character", na.strings = NULL, data.table = FALSE, header = TRUE))
}
annot <- function(d, g) {
  has <- vapply(strsplit(gsub(" ", "", d$GENE, fixed = TRUE), "|", fixed = TRUE), function(x) g %in% x, TRUE); sel <- d[has, , drop = FALSE]; hum <- sel[sel$ORGANISM == "Homo sapiens", , drop = FALSE]
  tdl <- py_sort(unique(unlist(strsplit(sel$TDL, "|", fixed = TRUE)))); tdl <- tdl[nzchar(tdl)]
  list(records_all_species = nrow(sel), records_human = nrow(hum), distinct_drugs_human = length(unique(tolower(hum$DRUG_NAME))), moa_records_human = as.integer(sum(hum$MOA == "1")),
       TDL_values = if (length(tdl)) paste(tdl, collapse = "|") else if (!nrow(sel)) "no record in the file" else "none stated")
}
part1 <- function() {
  summ <- fromJSON(file.path(W10, "drugcentral_summary.json"), simplifyVector = FALSE); old <- read_dc("reported"); new <- read_dc("new"); out <- list()
  for (g in c(HUBS, CONTEXT)) {
    a <- annot(old, g); b <- annot(new, g); r <- list(gene = g, role = if (g %in% HUBS) "hub (degree >= 5)" else "context (degree 4)")
    for (k in names(a)) r[[paste0("reported_file_", k)]] <- a[[k]]; for (k in names(b)) r[[paste0("file_25Sep2026_", k)]] <- b[[k]]
    r$MLSS_drugs_reported_file <- if (is.null(summ$drugs_per_hub_old[[g]])) "" else summ$drugs_per_hub_old[[g]]; r$MLSS_drugs_25Sep2026_file <- if (is.null(summ$drugs_per_hub_new[[g]])) "" else summ$drugs_per_hub_new[[g]]
    out[[length(out) + 1]] <- r
  }
  dictwriter_csv(out, names(out[[1]]), file.path(RES, "hub_annotation.csv"), fmt = "repr"); out
}
classes <- function(scores) {
  t75 <- np_quantile(scores, 0.75); t50 <- np_quantile(scores, 0.50); t25 <- np_quantile(scores, 0.25)
  list(th = list(t75 = t75, t50 = t50, t25 = t25), cl = function(s) if (s >= t75) "Q1" else if (s >= t50) "Q2" else if (s >= t25) "Q3" else "Q4")
}
part2 <- function() {
  comp <- read_csv_dict(file.path(W10, "mlss_comparison.csv")); allc <- read_csv_dict(file.path(W10, "mlss_new_all_combinations.csv")); rownames(allc) <- allc$Combination
  has_o <- comp$MLSS_old != ""; has_n <- comp$MLSS_new != ""
  old <- setNames(as.numeric(comp$MLSS_old[has_o]), comp$combination[has_o]); new <- setNames(as.numeric(comp$MLSS_new[has_n]), comp$combination[has_n])
  ro <- setNames(as.integer(comp$rank_old[comp$rank_old != ""]), comp$combination[comp$rank_old != ""]); rn <- setNames(as.integer(comp$rank_new[comp$rank_new != ""]), comp$combination[comp$rank_new != ""])
  co <- classes(old); cn_ <- classes(new); rows <- list(); size_o <- list(); size_n <- list(); changes <- list()
  for (c in names(new)[order(rn[names(new)], method = "radix")]) {
    a <- allc[c, ]; cross <- a$Targets_A != a$Targets_B; cl_o <- if (c %in% names(old)) co$cl(old[[c]]) else ""; cl_n <- cn_$cl(new[[c]])
    if (nzchar(cl_o)) size_o[[cl_o]] <- (if (is.null(size_o[[cl_o]])) 0L else size_o[[cl_o]]) + 1L; size_n[[cl_n]] <- (if (is.null(size_n[[cl_n]])) 0L else size_n[[cl_n]]) + 1L
    if (nzchar(cl_o) && cl_o != cl_n) changes[[length(changes) + 1]] <- list(c, cl_o, cl_n)
    rows[[length(rows) + 1]] <- list(combination = c, in_reported_153 = as.integer(c %in% names(old)), MLSS = round(new[[c]], 6), rank_reported = if (c %in% names(ro)) ro[[c]] else "", rank_25Sep2026 = rn[[c]], class_reported_153 = cl_o, class_171 = cl_n,
                                     cross_target_pair = as.integer(cross), targets = sprintf("%s + %s", a$Targets_A, a$Targets_B), reported_top15 = as.integer(c %in% names(ro) && ro[[c]] <= 15))
  }
  dictwriter_csv(rows, names(rows[[1]]), file.path(RES, "mlss_quartile_classes.csv"), fmt = "repr")
  top15 <- lapply(Filter(function(r) r$reported_top15 == 1L, rows), function(r) list(r$combination, r$rank_reported, r$class_reported_153)); srt <- function(l) l[order(names(l), method = "radix")]
  list(n_reported = length(old), n_171 = length(new), thresholds_reported_153 = co$th, thresholds_171 = cn_$th, class_sizes_reported_153 = srt(size_o), class_sizes_171 = srt(size_n), reported_combinations_changing_class_when_171_used = changes,
       top15_reported_classes = top15, distinct_scores_reported_153 = length(unique(round(old, 9))), max_abs_score_difference_shared = max(abs(old[names(old)] - new[names(old)])))
}
brands_of <- function(a) vapply(a$products, function(p) if (is.null(p$brand_name)) "" else p$brand_name, "")
part3 <- function() {
  dl <- drug_list(); apps_rows <- list(); drug_rows <- list(); cache <- list()
  for (i in seq_len(nrow(dl))) {
    r <- as.list(dl[i, ]); d <- r$drug; q <- query_name(d); if (is.null(cache[[q]])) cache[[q]] <- listing(q); total <- cache[[q]]$total; res <- cache[[q]]$res
    nb <- Filter(function(a) startsWith(a$application_number, "NDA") || startsWith(a$application_number, "BLA"), res)
    first_ap <- function(a) { ds <- vapply(Filter(function(s) identical(s$submission_type, "ORIG") && identical(s$submission_status, "AP"), a$submissions), function(s) s$submission_status_date, ""); ds <- py_sort(ds); if (length(ds)) ds[1] else "" }
    stat <- list()
    for (a in res) for (p in a$products) if (any(vapply(p$active_ingredients, function(x) grepl(q, tolower(if (is.null(x$name)) "" else x$name), fixed = TRUE), TRUE))) {
      k <- if (is.null(p$marketing_status)) "not stated" else p$marketing_status; stat[[k]] <- (if (is.null(stat[[k]])) 0L else stat[[k]]) + 1L
    }
    fa <- py_sort(Filter(nzchar, vapply(nb, first_ap, "")))
    for (a in nb) {
      sts <- py_sort(unique(vapply(a$products, function(p) if (is.null(p$marketing_status)) "not stated" else p$marketing_status, "")))
      apps_rows[[length(apps_rows) + 1]] <- list(drug = d, application = a$application_number, sponsor = if (is.null(a$sponsor_name)) "" else a$sponsor_name,
                                                 brands = paste(py_sort(unique(brands_of(a))), collapse = "|"), marketing_status_of_products = paste(sts, collapse = "|"), first_ORIG_AP_date = first_ap(a))
    }
    marketed <- as.integer(sum(unlist(stat[names(stat) %in% c("Prescription", "Over-the-counter")]))); tf <- strsplit(gsub(" ", "", gsub("^[()]+|[()]+$", "", r$FDA_EMA_PMDA_new)), ",", fixed = TRUE)[[1]]
    stat_sorted <- if (length(stat)) stat[order(names(stat), method = "radix")] else setNames(list(), character(0))
    drug_rows[[length(drug_rows) + 1]] <- list(drug = d, queried_name = q, applications_total_in_Drugs_at_FDA = total, NDA_or_BLA_applications = as.integer(length(nb)), products_by_marketing_status = py_json_typed(stat_sorted, indent = NULL),
                                               products_not_discontinued = marketed, US_status = if (total == 0) "no FDA application" else if (marketed) "listed as marketed (not discontinued)" else "application(s) present, all products discontinued or tentative",
                                               first_US_NDA_BLA_approval = if (length(fa)) fa[1] else "", NDA_BLA_brands = paste(py_sort(unique(unlist(lapply(nb, brands_of)))), collapse = "|"),
                                               DrugCentral_FDA = tf[1], DrugCentral_EMA = tf[2], DrugCentral_PMDA = tf[3])
  }
  dictwriter_csv(drug_rows, names(drug_rows[[1]]), file.path(RES, "market_status_by_drug.csv"), fmt = "repr"); dictwriter_csv(apps_rows, names(apps_rows[[1]]), file.path(RES, "market_status_nda_bla_applications.csv"), fmt = "repr")
  rd <- if (length(list.files(RAW, pattern = "[.]json$"))) RAW else RAW_ARCH
  dates <- vapply(list.files(rd, pattern = "[.]json$", full.names = TRUE), function(f) fromJSON(f, simplifyVector = FALSE)$retrieved, "")
  list(drug_rows = drug_rows, apps_rows = apps_rows, meta = list(responses_saved = as.integer(length(dates)), first_retrieved = min(dates), last_retrieved = max(dates)))
}
analyze <- function() {
  h <- part1(); q <- part2(); p3 <- part3(); out <- list(hubs = h, quartiles = q, market = p3$drug_rows, market_retrieval = p3$meta)
  write_json_typed(out, file.path(RES, "summary.json")); cat("analysis written to", RES, "\n")
  for (r in p3$drug_rows) cat(sprintf("  %-18s apps %3s NDA/BLA %2s  %-60s first NDA/BLA approval %s  brands %s\n", r$drug, r$applications_total_in_Drugs_at_FDA, r$NDA_or_BLA_applications, r$US_status, r$first_US_NDA_BLA_approval, substr(r$NDA_BLA_brands, 1, 60)))
}
if (sys.nframe() == 0L) { cmd <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else ""; if (cmd == "retrieve") retrieve() else if (cmd == "analyze") analyze() else cat("usage: Rscript wp12_annotations.R retrieve | analyze\n") }
