## wp9_safety.R - R port of wp9_safety/fda_fetch.py, wp9_screen.py, wp9_fallback.py and freeze_raw.py: WP9: pharmacovigilance safety screen (an annotation, NOT a prediction) of the MLSS drugs and
## cross-target drug pairs. Pre-specification: PRESPEC_WP9_safety_screen.md (fingerprinted before any report count of an MLSS drug was retrieved); every rule below is from it.
## Usage: Rscript wp9_safety.R <command>
##   retrieve           needs the network; cached and resumable; at most 850 new requests per run (openFDA allows 1000 per day without a key)
##   analyze            reads only the saved responses (offline); writes drug_counts.csv, signals_single.csv, pairs_all77.csv, pairs_long.csv, pairs_top15.csv and summary.json
##   fallback-retrieve  SECONDARY, POST HOC (spec section 9, 7 Nov 15:22): drugs with fewer than 100 reports through substance_name are searched again through the free-text field (at most 120 requests per run)
##   fallback-analyze   offline; writes fallback_signals_single.csv, fallback_pairs.csv, fallback_summary.json
##   freeze             fingerprint of all saved responses (RAW_SHA256.txt)
## Every response is stored as JSON in <raw>/ (file name = sha1 of the URL) with its URL and retrieval time; the cache of the archived run (wp9_safety/results/raw, WP9_RAW_ARCHIVE) is read first, so the
## analysis runs offline and on exactly the responses of the archived run. Environment: WP9_OUT (default results_R/wp9). Inputs (read only): wp3_mlss/direction_table.csv, mlss_scenario_scores.csv (archived).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- HERE
source(file.path(R_PORT, "pycompat.R")); suppressMessages({ library(jsonlite); library(digest); library(curl) })
V2 <- normalizePath(file.path(R_PORT, "..")); RES <- Sys.getenv("WP9_OUT", file.path(R_PORT, "results_R", "wp9")); RAW <- file.path(RES, "raw"); dir.create(RAW, showWarnings = FALSE, recursive = TRUE)
RAW_ARCH <- Sys.getenv("WP9_RAW_ARCHIVE", file.path(V2, "wp9_safety", "results", "raw")); options(py_json_nan = TRUE)
## ---------------------------------------------------------------------------------------------------------------- fda_fetch
BASE <- "https://api.fda.gov/drug/event.json"; MIN_GAP <- 0.4; ST <- new.env(); ST$last <- 0; ST$requests <- 0L; ST$cached <- 0L; ST$offline <- FALSE
quote_plus <- function(s) vapply(s, function(x) { b <- charToRaw(enc2utf8(x)); paste(vapply(b, function(by) { ch <- rawToChar(by); if (grepl("[A-Za-z0-9_.~-]", ch)) ch else if (ch == " ") "+" else sprintf("%%%02X", as.integer(by)) }, ""), collapse = "") }, "", USE.NAMES = FALSE)
url_for <- function(search = NULL, count = NULL, limit = 1) { parts <- character(0); if (!is.null(search) && nzchar(search)) parts <- c(parts, paste0("search=", quote_plus(search))); if (!is.null(count) && nzchar(count)) parts <- c(parts, paste0("count=", count)); parts <- c(parts, sprintf("limit=%d", limit)); paste0(BASE, "?", paste(parts, collapse = "&")) }
cache_file <- function(url) paste0(digest(url, algo = "sha1", serialize = FALSE), ".json")
fetch <- function(url, tries = 9) {
  for (d in c(RAW, RAW_ARCH)) { p <- file.path(d, cache_file(url)); if (file.exists(p)) { ST$cached <- ST$cached + 1L; return(fromJSON(p, simplifyVector = FALSE)$response) } }
  if (ST$offline) stop("offline and not saved: ", url)
  delay <- 3; last_err <- NULL; path <- file.path(RAW, cache_file(url))
  for (k in seq_len(tries)) {
    gap <- as.numeric(Sys.time()) - ST$last; if (gap < MIN_GAP) Sys.sleep(MIN_GAP - gap); ST$last <- as.numeric(Sys.time())
    r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = "research-screen/1.0", timeout = 120)), error = function(e) { last_err <<- paste0(class(e)[1], ": ", substr(conditionMessage(e), 1, 120)); NULL })
    if (!is.null(r)) {
      if (r$status_code == 200 || r$status_code == 404) {
        resp <- if (r$status_code == 200) fromJSON(rawToChar(r$content), simplifyVector = FALSE) else list(meta = list(results = list(total = 0L)), results = list(), note = "HTTP 404: no match")
        ST$requests <- ST$requests + 1L; write_json_typed(list(url = url, retrieved = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), response = resp), path, indent = NULL); return(resp)
      }
      last_err <- paste("HTTP", r$status_code); if (!(r$status_code %in% c(429, 500, 502, 503, 504))) stop(last_err)
    }
    cat(sprintf("  retry %d/%d after %s (wait %.0f s)\n", k, tries, last_err, delay)); Sys.sleep(delay); delay <- min(delay * 2, 120)
  }
  stop("giving up on ", url, ": ", last_err)
}
total <- function(search = NULL) fetch(url_for(search, NULL, 1))$meta$results$total
last_updated <- function() fetch(url_for(NULL, NULL, 1))$meta$last_updated
## ---------------------------------------------------------------------------------------------------------------- wp9_screen
D <- "patient.drug.openfda.substance_name"; Rr <- "patient.reaction.reactionmeddrapt.exact"; MIN_DRUG <- 100L; MIN_PAIR <- 20L; MIN_A <- 3L; MAX_NEW <- 850L; Z <- 1.959963984540054
GROUPS <- list(
  list("HYPO", "hypoglycemia", c("HYPOGLYCAEMIA", "HYPOGLYCAEMIC COMA", "HYPOGLYCAEMIC UNCONSCIOUSNESS", "BLOOD GLUCOSE DECREASED")),
  list("HYPER", "hyperglycemia, loss of glycemic control", c("HYPERGLYCAEMIA", "BLOOD GLUCOSE INCREASED", "DIABETES MELLITUS INADEQUATE CONTROL", "DIABETIC KETOACIDOSIS", "GLYCOSYLATED HAEMOGLOBIN INCREASED")),
  list("PANC", "pancreatic events", c("PANCREATITIS", "PANCREATITIS ACUTE", "PANCREATITIS CHRONIC", "PANCREATITIS NECROTISING", "LIPASE INCREASED", "AMYLASE INCREASED")),
  list("GI", "gastrointestinal intolerance", c("NAUSEA", "VOMITING", "DIARRHOEA", "ABDOMINAL PAIN", "ABDOMINAL PAIN UPPER", "DECREASED APPETITE", "CONSTIPATION")),
  list("BILI", "biliary events", c("CHOLELITHIASIS", "CHOLECYSTITIS", "CHOLECYSTITIS ACUTE", "BILIARY COLIC")),
  list("HYPOT", "hypotension, syncope", c("HYPOTENSION", "ORTHOSTATIC HYPOTENSION", "SYNCOPE", "PRESYNCOPE", "BLOOD PRESSURE DECREASED")),
  list("OEDEMA", "edema", c("OEDEMA PERIPHERAL", "OEDEMA", "PERIPHERAL SWELLING")),
  list("CARD", "conduction and rhythm", c("BRADYCARDIA", "ATRIOVENTRICULAR BLOCK", "ATRIOVENTRICULAR BLOCK COMPLETE", "ELECTROCARDIOGRAM QT PROLONGED", "TORSADE DE POINTES", "VENTRICULAR TACHYCARDIA")))
CODES <- vapply(GROUPS, function(g) g[[1]], ""); PT <- setNames(lapply(GROUPS, function(g) g[[3]]), CODES)
term <- function(name) ifelse(name == "(S)-nitrendipine", "nitrendipine", tolower(name)); dq <- function(t) sprintf("%s:\"%s\"", D, t); gq <- function(c) sprintf("%s:(%s)", Rr, paste(sprintf("\"%s\"", PT[[c]]), collapse = " OR ")); pq <- function(a, b) paste0(dq(a), " AND ", dq(b))
rd_csv <- function(p) read.csv(p, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = "UTF-8-BOM")
load_mlss <- function() {
  dt <- rd_csv(file.path(V2, "wp3_mlss", "direction_table.csv")); dt <- dt[dt$used_in_mlss == "True", ]; cls <- setNames(dt$target_gene, dt$drug)
  sc <- rd_csv(file.path(V2, "wp3_mlss", "mlss_scenario_scores.csv")); ab <- strsplit(sc$Combination, " + ", fixed = TRUE)
  rows <- lapply(seq_len(nrow(sc)), function(i) list(rank = as.integer(as.numeric(sc$Rank[i])), pair = sc$Combination[i], mlss = as.numeric(sc$reported_MLSS[i]), A = ab[[i]][1], B = ab[[i]][2], cA = cls[[ab[[i]][1]]], cB = cls[[ab[[i]][2]]]))
  cross <- Filter(function(r) r$cA != r$cB, rows); cross <- cross[order(vapply(cross, function(r) r$rank, 1L))]
  stopifnot(length(cls) == 18, length(rows) == 153, length(cross) == 77, identical(vapply(cross, function(r) r$rank, 1L), 1:77)); list(cls = cls, rows = rows, cross = cross)
}
pair_keys <- function(cross) { seen <- character(0); keys <- list(); for (r in cross) { k <- py_sort(c(term(r$A), term(r$B))); kk <- paste(k, collapse = "\t"); if (!(kk %in% seen)) { seen <- c(seen, kk); keys[[length(keys) + 1]] <- k } }; keys }
budget <- function(max_new = MAX_NEW) if (ST$requests >= max_new) { cat(sprintf("request budget of this run reached (%d new requests); run again later (the daily limit without a key is 1000)\n", ST$requests)); quit(status = 2) }
retrieve <- function() {
  m <- load_mlss(); cls <- m$cls; cross <- m$cross; terms <- py_sort(unique(term(names(cls)))); keys <- pair_keys(cross); stopifnot(length(terms) == 17, length(keys) == 70)
  cat("database: all reports", total(NULL), ", last updated", last_updated(), "\n")
  for (c in CODES) for (t in PT[[c]]) { n <- total(sprintf("%s:\"%s\"", Rr, t)); cat(sprintf("  PT %-7s %-38s %9d\n", c, t, n)); if (n == 0) { cat("PT without any report: fix the spelling BEFORE any drug count is retrieved and log it in section 9 of the pre-specification\n"); quit(status = 3) }; budget() }
  for (c in CODES) { cat("group", c, total(gq(c)), "\n"); budget() }
  for (t in terms) { n <- total(dq(t)); cat(sprintf("drug %-14s %9d reports\n", t, n)); if (n >= MIN_DRUG) for (c in CODES) { total(paste0(dq(t), " AND ", gq(c))); budget() }; budget() }
  for (i in seq_along(keys)) { a <- keys[[i]][1]; b <- keys[[i]][2]; n <- total(pq(a, b)); cat(sprintf("pair %2d/70 %-14s + %-14s %8d reports%s\n", i, a, b, n, if (n >= MIN_PAIR) "" else "  (not evaluable)")); if (n >= MIN_PAIR) for (c in CODES) { total(paste0(pq(a, b), " AND ", gq(c))); budget() }; budget() }
  cat(sprintf("retrieval complete: %d new requests in this run, %d answered from the saved files\n", ST$requests, ST$cached))
}
ror <- function(a, n_d, E, N) {
  a <- as.numeric(a); n_d <- as.numeric(n_d); E <- as.numeric(E); N <- as.numeric(N)            # doubles: Python integers do not overflow, R integers (2^31) would in the products below
  cells <- c(a, n_d - a, E - a, N - n_d - E + a); if (min(cells) <= 0) cells <- cells + 0.5
  lor <- log(cells[1] * cells[4] / (cells[2] * cells[3])); se <- sqrt(1 / cells[1] + 1 / cells[2] + 1 / cells[3] + 1 / cells[4]); c(exp(lor), exp(lor - Z * se), exp(lor + Z * se))
}
cp <- function(a, n) c(if (a == 0) 0.0 else qbeta(0.025, a, n - a + 1), if (a == n) 1.0 else qbeta(0.975, a + 1, n - a))
write_csv <- function(name, header, data, dir = RES) {       # csv.writer: default dialect (CRLF, minimal quoting), str() of every value (strings are formatted by the caller)
  lines <- c(paste(vapply(header, .dw_quote, ""), collapse = ","), vapply(data, function(r) paste(vapply(r, function(v) .dw_quote(.dw_val(v, "repr")), ""), collapse = ","), ""))
  write_text_bytes(paste0(paste(lines, collapse = "\r\n"), "\r\n"), file.path(dir, name))
}
f <- function(x, nd) py_fmt_f(x, nd)
analyze <- function() {
  ST$offline <- TRUE; m <- load_mlss(); cls <- m$cls; cross <- m$cross; N <- total(NULL); LU <- last_updated(); E <- setNames(vapply(CODES, function(c) as.integer(total(gq(c))), 1L), CODES)
  tm <- term(names(cls)); entries <- split(names(cls), tm); terms <- py_sort(names(entries)); entries <- entries[terms]; tclass <- setNames(vapply(terms, function(t) cls[[entries[[t]][1]]], ""), terms)
  nd <- setNames(vapply(terms, function(t) as.integer(total(dq(t))), 1L), terms); ev <- nd >= MIN_DRUG; a_dg <- list(); for (t in terms[ev]) for (c in CODES) a_dg[[paste(t, c)]] <- as.integer(total(paste0(dq(t), " AND ", gq(c))))
  warn <- character(0); sig <- list(); single_rows <- list()
  for (t in terms[ev]) for (c in CODES) {
    a <- a_dg[[paste(t, c)]]; if (a > min(nd[[t]], E[[c]])) warn <- c(warn, sprintf("a_dg > min(n_d, E_g): %s %s", t, c)); o <- ror(a, nd[[t]], E[[c]], N); s <- a >= MIN_A && o[2] > 1; sig[[paste(t, c)]] <- s
    single_rows[[length(single_rows) + 1]] <- list(t, tclass[[t]], c, a, nd[[t]], f(a / nd[[t]], 5), E[[c]], f(o[1], 3), f(o[2], 3), f(o[3], 3), as.integer(s))
  }
  write_csv("signals_single.csv", c("term", "target_gene", "group", "reports_in_group", "reports_of_drug", "proportion", "reports_in_group_all_drugs", "ROR", "ROR_lo95", "ROR_hi95", "signal"), single_rows)
  write_csv("drug_counts.csv", c("term", "mlss_entries", "target_gene", "reports", "evaluable"), lapply(terms, function(t) list(t, paste(entries[[t]], collapse = ";"), tclass[[t]], nd[[t]], as.integer(ev[[t]]))))
  issig <- function(t, c) isTRUE(sig[[paste(t, c)]])
  nab <- list(); aab <- list(); for (k in pair_keys(cross)) { kk <- paste(k, collapse = "\t"); nab[[kk]] <- as.integer(total(pq(k[1], k[2]))); if (nab[[kk]] >= MIN_PAIR) for (c in CODES) aab[[paste(kk, c)]] <- as.integer(total(paste0(pq(k[1], k[2]), " AND ", gq(c)))) }
  pair_rows <- list(); long_rows <- list()
  for (r in cross) {
    tA <- term(r$A); tB <- term(r$B); k <- py_sort(c(tA, tB)); kk <- paste(k, collapse = "\t"); n_ab <- nab[[kk]]; pev <- n_ab >= MIN_PAIR; cev <- ev[[tA]] && ev[[tB]]
    sA <- CODES[vapply(CODES, function(c) issig(tA, c), TRUE)]; sB <- CODES[vapply(CODES, function(c) issig(tB, c), TRUE)]; flags <- character(0)
    if (n_ab > min(nd[[tA]], nd[[tB]])) warn <- c(warn, sprintf("n_AB > min(n_A, n_B): (%s)", paste(sprintf("'%s'", k), collapse = ", ")))
    if (pev && cev) for (c in CODES) {
      a <- aab[[paste(kk, c)]]; if (a > min(a_dg[[paste(tA, c)]], a_dg[[paste(tB, c)]], n_ab)) warn <- c(warn, sprintf("a_AB > min(a_A, a_B, n_AB): (%s) %s", paste(sprintf("'%s'", k), collapse = ", "), c))
      p <- a / n_ab; ci <- cp(a, n_ab); pA <- if (nd[[tA]] > n_ab) (a_dg[[paste(tA, c)]] - a) / (nd[[tA]] - n_ab) else NaN; pB <- if (nd[[tB]] > n_ab) (a_dg[[paste(tB, c)]] - a) / (nd[[tB]] - n_ab) else NaN
      fl <- a >= MIN_A && isTRUE(ci[1] > max(pA, pB)); if (fl) flags <- c(flags, c)
      long_rows[[length(long_rows) + 1]] <- list(r$rank, r$pair, c, a, n_ab, f(p, 5), f(ci[1], 5), f(ci[2], 5), f(pA, 5), f(pB, 5), as.integer(fl))
    }
    status <- if (pev && cev) "evaluable" else if (pev) sprintf("component drug below %d reports", MIN_DRUG) else sprintf("fewer than %d co-reports", MIN_PAIR)
    pair_rows[[length(pair_rows) + 1]] <- list(r$rank, r$pair, f(r$mlss, 4), tA, tB, n_ab, status, if (ev[[tA]]) (if (length(sA)) paste(sA, collapse = " ") else "-") else "n/e", if (ev[[tB]]) (if (length(sB)) paste(sB, collapse = " ") else "-") else "n/e", if (pev && cev) (if (length(flags)) paste(flags, collapse = " ") else "-") else "n/a")
  }
  hdr <- c("rank", "pair", "MLSS", "term_A", "term_B", "reports_listing_both", "status", "signals_drug_A", "signals_drug_B", "co_report_excess")
  write_csv("pairs_all77.csv", hdr, pair_rows); write_csv("pairs_top15.csv", hdr, Filter(function(p) p[[1]] <= 15, pair_rows))
  write_csv("pairs_long.csv", c("rank", "pair", "group", "reports_both_in_group", "reports_both", "p_both", "p_both_lo95", "p_both_hi95", "p_A_without_B", "p_B_without_A", "co_report_excess"), long_rows)
  n_sig <- sum(unlist(sig)); per_group <- list()
  for (c in CODES) {
    drugs_c <- py_sort(terms[vapply(terms, function(t) issig(t, c), TRUE)]); has <- function(r) issig(term(r$A), c) || issig(term(r$B), c)
    per_group[[c]] <- list(E_g = E[[c]], signalled_drugs = as.list(drugs_c), pairs_top15_with_signalled_drug = as.integer(sum(vapply(cross, function(r) r$rank <= 15 && has(r), TRUE))), pairs_all77_with_signalled_drug = as.integer(sum(vapply(cross, has, TRUE))))
  }
  flagged <- lapply(Filter(function(p) !(p[[10]] %in% c("-", "n/a")), pair_rows), function(p) list(p[[1]], p[[2]], p[[10]]))
  summary <- list(last_updated = LU, N_reports = N, terms = length(terms), evaluable_drugs = as.integer(sum(ev)), not_evaluable_drugs = as.list(terms[!ev]), drug_reports = as.list(nd), single_signals_total = as.integer(n_sig), single_tests_total = length(single_rows), per_group = per_group,
                  pairs_total = length(cross), pairs_with_20_or_more_co_reports = as.integer(sum(vapply(pair_rows, function(p) !startsWith(p[[7]], "fewer"), TRUE))), pairs_evaluable = as.integer(sum(vapply(pair_rows, function(p) p[[7]] == "evaluable", TRUE))),
                  pairs_evaluable_top15 = as.integer(sum(vapply(pair_rows, function(p) p[[1]] <= 15 && p[[7]] == "evaluable", TRUE))), co_report_flags = as.integer(sum(vapply(flagged, function(x) length(strsplit(x[[3]], " ", fixed = TRUE)[[1]]), 1L))),
                  co_report_flagged_pair_count = length(flagged), co_report_flagged_pairs = flagged, pair_group_tests = length(long_rows), warnings = as.list(warn))
  write_json_typed(summary, file.path(RES, "summary.json"))
  cat(sprintf("N reports %d (last updated %s); drugs evaluable %d of %d; single-drug signals %d of %d tests; pairs evaluable %d of %d; co-report flags %d (in %d pairs) of %d pair-group tests; warnings %d\n", N, LU, summary$evaluable_drugs, length(terms), n_sig, length(single_rows),
              summary$pairs_evaluable, length(cross), summary$co_report_flags, length(flagged), length(long_rows), length(warn))); for (w in warn) cat("WARNING", w, "\n")
}
## ---------------------------------------------------------------------------------------------------------------- wp9_fallback (SECONDARY, POST HOC)
TEXT <- "patient.drug.medicinalproduct"; NAMES <- list(albiglutide = c("TANZEUM", "ALBIGLUTIDE"), bepridil = c("VASCOR", "BEPRIDIL"), nilvadipine = "NILVADIPINE", nitrendipine = "NITRENDIPINE"); MAX_NEW_FB <- 120L
tq <- function(t) sprintf("%s:(%s)", TEXT, paste(sprintf("\"%s\"", NAMES[[t]]), collapse = " OR "))
partners <- function(cross, t) { out <- list(); for (r in cross) { a <- term(r$A); b <- term(r$B); if (t %in% c(a, b)) out[[length(out) + 1]] <- list(r = r, p = if (a == t) b else a) }; out }
low_terms <- function(terms) terms[vapply(terms, function(t) total(dq(t)) < MIN_DRUG, TRUE)]
fallback_retrieve <- function() {
  m <- load_mlss(); cross <- m$cross; terms <- py_sort(unique(term(names(m$cls)))); low <- low_terms(terms); cat(sprintf("drugs below %d reports through substance_name: %s\n", MIN_DRUG, paste(low, collapse = ", ")))
  for (t in low) {
    if (is.null(NAMES[[t]])) { cat("no fallback name was fixed in advance for", t, "(not searched)\n"); next }
    n <- total(tq(t)); cat(sprintf("fallback %-12s %8d reports (free text %s)\n", t, n, paste(NAMES[[t]], collapse = ", "))); budget(MAX_NEW_FB); if (n < MIN_DRUG) next
    for (c in CODES) { total(paste0(tq(t), " AND ", gq(c))); budget(MAX_NEW_FB) }
    ps <- partners(cross, t); ps <- ps[order(vapply(ps, function(x) x$r$rank, 1L))]; seen <- character(0)
    for (x in ps) { key <- paste(x$r$rank, x$p); if (key %in% seen) next; seen <- c(seen, key); nab <- total(paste0(tq(t), " AND ", dq(x$p))); cat(sprintf("  pair %-12s + %-14s %7d co-reports\n", t, x$p, nab)); budget(MAX_NEW_FB)
      if (nab >= MIN_PAIR) for (c in CODES) { total(paste0(tq(t), " AND ", dq(x$p), " AND ", gq(c))); budget(MAX_NEW_FB) } }
  }
  cat(sprintf("fallback retrieval complete: %d new requests, %d from the saved files\n", ST$requests, ST$cached))
}
fallback_analyze <- function() {
  ST$offline <- TRUE; m <- load_mlss(); cls <- m$cls; cross <- m$cross; terms <- py_sort(unique(term(names(cls)))); N <- total(NULL); E <- setNames(vapply(CODES, function(c) as.integer(total(gq(c))), 1L), CODES)
  nd <- setNames(vapply(terms, function(t) as.integer(total(dq(t))), 1L), terms); ev <- nd >= MIN_DRUG; a_dg <- list(); for (t in terms[ev]) for (c in CODES) a_dg[[paste(t, c)]] <- as.integer(total(paste0(dq(t), " AND ", gq(c))))
  sig_primary <- list(); for (k in names(a_dg)) { t <- sub(" [A-Z]+$", "", k); c <- sub("^.* ", "", k); o <- ror(a_dg[[k]], nd[[t]], E[[c]], N); sig_primary[[k]] <- a_dg[[k]] >= MIN_A && o[2] > 1 }
  single <- list(); pairs_out <- list(); summary <- list(note = "SECONDARY, POST HOC, free-text names; not merged into the primary counts", drugs = list())
  for (t in low_terms(terms)) {
    if (is.null(NAMES[[t]])) { summary$drugs[[t]] <- list(searched = FALSE); next }
    n_t <- total(tq(t)); info <- list(searched = TRUE, names = as.list(NAMES[[t]]), reports = n_t, evaluable = n_t >= MIN_DRUG); summary$drugs[[t]] <- info; if (n_t < MIN_DRUG) next
    a_t <- setNames(vapply(CODES, function(c) as.integer(total(paste0(tq(t), " AND ", gq(c)))), 1L), CODES); sig_t <- list()
    tg <- cls[[names(cls)[term(names(cls)) == t][1]]]
    for (c in CODES) { o <- ror(a_t[[c]], n_t, E[[c]], N); s <- a_t[[c]] >= MIN_A && o[2] > 1; sig_t[[c]] <- s; single[[length(single) + 1]] <- list(t, tg, c, a_t[[c]], n_t, f(a_t[[c]] / n_t, 5), E[[c]], f(o[1], 3), f(o[2], 3), f(o[3], 3), as.integer(s)) }
    summary$drugs[[t]]$signals <- as.list(CODES[vapply(CODES, function(c) isTRUE(sig_t[[c]]), TRUE)])
    for (x in partners(cross, t)) {
      r <- x$r; p <- x$p; nab <- total(paste0(tq(t), " AND ", dq(p))); pev <- nab >= MIN_PAIR; cev <- ev[[p]]; sP <- CODES[vapply(CODES, function(c) isTRUE(sig_primary[[paste(p, c)]]), TRUE)]; flags <- character(0)
      if (pev && cev) for (c in CODES) {
        xx <- total(paste0(tq(t), " AND ", dq(p), " AND ", gq(c))); ci <- cp(xx, nab); pA <- if (n_t > nab) (a_t[[c]] - xx) / (n_t - nab) else NaN; pB <- if (nd[[p]] > nab) (a_dg[[paste(p, c)]] - xx) / (nd[[p]] - nab) else NaN
        if (xx >= MIN_A && isTRUE(ci[1] > max(pA, pB))) flags <- c(flags, c)
      }
      status <- if (pev && cev) "evaluable" else if (pev) sprintf("partner drug below %d reports", MIN_DRUG) else sprintf("fewer than %d co-reports", MIN_PAIR)
      sigt <- CODES[vapply(CODES, function(c) isTRUE(sig_t[[c]]), TRUE)]
      pairs_out[[length(pairs_out) + 1]] <- list(r$rank, r$pair, f(r$mlss, 4), t, p, nab, status, if (length(sigt)) paste(sigt, collapse = " ") else "-", if (cev) (if (length(sP)) paste(sP, collapse = " ") else "-") else "n/e", if (pev && cev) (if (length(flags)) paste(flags, collapse = " ") else "-") else "n/a")
    }
  }
  pairs_out <- pairs_out[order(vapply(pairs_out, function(p) p[[1]], 1L))]
  write_csv("fallback_signals_single.csv", c("term", "target_gene", "group", "reports_in_group", "reports_of_drug", "proportion", "reports_in_group_all_drugs", "ROR", "ROR_lo95", "ROR_hi95", "signal"), single)
  write_csv("fallback_pairs.csv", c("rank", "pair", "MLSS", "fallback_term", "partner_term", "reports_listing_both", "status", "signals_fallback_drug", "signals_partner_drug", "co_report_excess"), pairs_out)
  summary$pairs <- length(pairs_out); summary$pairs_evaluable <- as.integer(sum(vapply(pairs_out, function(p) p[[7]] == "evaluable", TRUE))); summary$pairs_evaluable_top15 <- as.integer(sum(vapply(pairs_out, function(p) p[[1]] <= 15 && p[[7]] == "evaluable", TRUE)))
  summary$pairs_top15 <- as.integer(sum(vapply(pairs_out, function(p) p[[1]] <= 15, TRUE))); summary$co_report_flags <- as.integer(sum(vapply(pairs_out, function(p) if (p[[10]] %in% c("-", "n/a")) 0L else length(strsplit(p[[10]], " ", fixed = TRUE)[[1]]), 1L)))
  write_json_typed(summary, file.path(RES, "fallback_summary.json")); cat("fallback analysis written to", RES, "\n")
}
## ---------------------------------------------------------------------------------------------------------------- freeze_raw: sha256 over 'file name:sha256 of the file' lines of all saved responses, files sorted by name
freeze <- function() {
  rd <- if (length(list.files(RAW, pattern = "[.]json$"))) RAW else RAW_ARCH; files <- py_sort(list.files(rd, pattern = "[.]json$")); lines <- paste0(vapply(files, function(fn) sprintf("%s:%s\n", fn, digest(file.path(rd, fn), algo = "sha256", file = TRUE)), ""), collapse = "")
  h <- digest(lines, algo = "sha256", serialize = FALSE); recs <- lapply(file.path(rd, files), function(p) fromJSON(p, simplifyVector = FALSE)); times <- vapply(recs, function(r) r$retrieved, ""); lu <- unique(Filter(Negate(is.null), lapply(recs, function(r) r$response$meta$last_updated)))
  txt <- sprintf("%s *wp9_safety/results/raw (%d files)\nfirst retrieval %s, last retrieval %s, openFDA last_updated %s; written %s BEFORE the analysis\n", h, length(files), min(times), max(times), paste(py_sort(unlist(lu)), collapse = ","), format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  write_text_bytes(txt, file.path(RES, "RAW_SHA256.txt")); cat(txt)
}
if (sys.nframe() == 0L) {
  cmd <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else ""
  if (cmd == "retrieve") retrieve() else if (cmd == "analyze") analyze() else if (cmd == "fallback-retrieve") fallback_retrieve() else if (cmd == "fallback-analyze") fallback_analyze() else if (cmd == "freeze") freeze() else cat("usage: Rscript wp9_safety.R retrieve | analyze | fallback-retrieve | fallback-analyze | freeze\n")
}
