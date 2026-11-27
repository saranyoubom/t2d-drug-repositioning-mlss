## wp9_check_numbers.R - R port of wp9_safety/check_wp9_numbers.py: independent rebuild of every WP9 number from the saved openFDA responses (the raw cache of the archived run, wp9_safety/results/raw, which wp9_safety.R reads) and from the
## pre-specification text. It has its own parser of the saved responses, its own query builders and its own formulas; wp9_safety.R is sourced only for the group definitions, the fallback names, the URL quoting and the exact binomial
## interval, and its output files (WP9_OUT, default results_R/wp9) are compared with the rebuilt values. The pre-specification, README and wording are the archived texts of wp9_safety/.
## Prints OK or FAIL per check, 'N checks, M failures', and writes check_wp9_last_run.txt (into WP9_CHECK_OUT, default WP9_OUT) last.   Usage: Rscript wp9_check_numbers.R
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "wp9_safety.R"))      # PT, NAMES, quote_plus, cp, V2, RES, RAW, RAW_ARCH, BASE, D, Rr
ARC9 <- file.path(V2, "wp9_safety"); CHK <- Sys.getenv("WP9_CHECK_OUT", RES); dir.create(CHK, showWarnings = FALSE, recursive = TRUE); RAWD <- if (length(list.files(RAW_ARCH, pattern = "[.]json$"))) RAW_ARCH else RAW
count <- 0L; fails <- character(0); ok <- function(c, m) { count <<- count + 1L; cat(if (isTRUE(c)) "OK   " else "FAIL ", m, "\n", sep = ""); if (!isTRUE(c)) fails <<- c(fails, m) }
close_ <- function(a, b, tol = 6e-4) abs(as.numeric(a) - as.numeric(b)) <= tol * max(1, abs(as.numeric(b))); fmtc <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
rd_text <- function(p) { b <- readBin(p, "raw", file.size(p)); s <- rawToChar(b); Encoding(s) <- "UTF-8"; gsub("\r\n", "\n", s, fixed = TRUE) }; rcsv <- function(p) read.csv(p, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = "UTF-8-BOM")
inn <- function(s, x) grepl(enc2utf8(s), x, fixed = TRUE, useBytes = TRUE); jstr <- function(...) paste(..., sep = "")
## 1. groups from the pre-specification text
spec <- rd_text(file.path(ARC9, "PRESPEC_WP9_safety_screen.md")); mm <- regmatches(spec, gregexpr("(?m)^\\| ([A-Z]+) \\| [^|]+ \\| ([^|]+) \\|$", spec, perl = TRUE))[[1]]; groups <- list()
for (x in mm) { g <- regmatches(x, regexec("^\\| ([A-Z]+) \\| [^|]+ \\| ([^|]+) \\|$", x, perl = TRUE))[[1]]; groups[[g[2]]] <- trimws(strsplit(g[3], ";", fixed = TRUE)[[1]]) }
ok(length(groups) == 8, sprintf("pre-specification table has 8 event groups (%d)", length(groups)))
ok(identical(sort(names(groups)), sort(names(PT))) && all(vapply(names(groups), function(c) identical(groups[[c]], PT[[c]]), TRUE)), "group definitions of wp9_safety.R equal the pre-specification table term by term")
ok(sum(vapply(groups, length, 1L)) == 40, sprintf("the groups hold 40 preferred terms (%d)", sum(vapply(groups, length, 1L))))
## 2. saved responses
resp <- list(); for (p in list.files(RAWD, pattern = "[.]json$", full.names = TRUE)) { rec <- fromJSON(p, simplifyVector = FALSE); resp[[rec$url]] <- rec$response$meta$results$total; if (paste0(digest::digest(rec$url, algo = "sha1", serialize = FALSE), ".json") != basename(p)) ok(FALSE, paste0("file name is not the sha1 of its url: ", basename(p))) }
U <- function(search) paste0(BASE, "?", if (!is.null(search) && nzchar(search)) paste0("search=", quote_plus(search), "&") else "", "limit=1"); T_ <- function(search) { v <- resp[[U(search)]]; if (is.null(v)) stop("no saved response for ", search); as.numeric(v) }
drug <- function(t) paste0(D, ":\"", t, "\""); grp <- function(c) paste0(Rr, ":(", paste(sprintf("\"%s\"", groups[[c]]), collapse = " OR "), ")")
## 3. MLSS inputs
dt <- rcsv(file.path(V2, "wp3_mlss", "direction_table.csv")); dt <- dt[dt$used_in_mlss == "True", ]; cls <- setNames(dt$target_gene, dt$drug); sc <- rcsv(file.path(V2, "wp3_mlss", "mlss_scenario_scores.csv")); cross <- list()
for (i in seq_len(nrow(sc))) { ab <- strsplit(sc$Combination[i], " + ", fixed = TRUE)[[1]]; if (cls[[ab[1]]] != cls[[ab[2]]]) cross[[length(cross) + 1L]] <- list(rank = as.integer(as.numeric(sc$Rank[i])), name = sc$Combination[i], a = ab[1], b = ab[2]) }
cross <- cross[order(vapply(cross, function(z) z$rank, 1L))]; ok(length(cls) == 18 && length(cross) == 77 && identical(vapply(cross, function(z) z$rank, 1L), 1:77), "18 drugs, 77 cross-target pairs with ranks 1 to 77")
tm <- function(n) ifelse(n == "(S)-nitrendipine", "nitrendipine", tolower(n)); terms <- py_sort(unique(tm(names(cls)))); ok(length(terms) == 17, "17 distinct search terms")
## 4. rebuild and compare
S <- fromJSON(file.path(RES, "summary.json"), simplifyVector = FALSE); N <- T_(NULL); ok(N == S$N_reports, sprintf("N all reports %d equals summary", N))
E <- vapply(names(groups), function(c) T_(grp(c)), 1); ok(all(vapply(names(groups), function(c) E[[c]] == S$per_group[[c]]$E_g, TRUE)), "event-group report totals equal summary: " %+% paste(sprintf("%s %d", names(groups), E), collapse = ", "))
nd <- vapply(terms, function(t) T_(drug(t)), 1); ok(all(vapply(terms, function(t) nd[[t]] == S$drug_reports[[t]], TRUE)), "drug report totals equal summary")
ev <- nd >= 100; ok(sum(ev) == S$evaluable_drugs, sprintf("%d of 17 drugs evaluable (100 reports or more)", sum(ev)))
a_dg <- list(); for (t in terms[ev]) for (c in names(groups)) a_dg[[paste(t, c, sep = "\r")]] <- T_(paste0(drug(t), " AND ", grp(c))); adg <- function(t, c) a_dg[[paste(t, c, sep = "\r")]]
ok(all(vapply(names(a_dg), function(k) { s <- strsplit(k, "\r", fixed = TRUE)[[1]]; a_dg[[k]] <= min(nd[[s[1]]], E[[s[2]]]) }, TRUE)), "every drug-group count is at most the drug total and the group total")
lor_ci <- function(a, n, e) { cells <- c(a, n - a, e - a, N - n - e + a); if (min(cells) <= 0) cells <- cells + 0.5; lo_r <- log(cells[1] * cells[4] / (cells[2] * cells[3])); se <- sqrt(sum(1 / cells)); c(exp(lo_r), exp(lo_r - 1.959963984540054 * se), exp(lo_r + 1.959963984540054 * se)) }
rows <- rcsv(file.path(RES, "signals_single.csv")); ok(nrow(rows) == sum(ev) * 8, sprintf("signals_single.csv has %d rows (evaluable drugs x 8 groups)", nrow(rows))); bad <- 0L; sig <- list()
for (i in seq_len(nrow(rows))) { r <- rows[i, ]; t <- r$term; c <- r$group; a <- adg(t, c); o <- lor_ci(a, nd[[t]], E[[c]]); s <- a >= 3 && o[2] > 1; sig[[paste(t, c, sep = "\r")]] <- s
  if (!(as.integer(r$reports_in_group) == a && as.integer(r$reports_of_drug) == nd[[t]] && close_(r$ROR, o[1]) && close_(r$ROR_lo95, o[2]) && close_(r$ROR_hi95, o[3]) && as.integer(r$signal) == as.integer(s))) bad <- bad + 1L }
sg <- function(t, c) isTRUE(sig[[paste(t, c, sep = "\r")]])
ok(bad == 0, sprintf("every row of signals_single.csv (counts, ROR, interval, signal) equals the independent rebuild (%d rows differ)", bad)); ok(sum(unlist(sig)) == S$single_signals_total, sprintf("%d single-drug signals of %d tests equal summary", sum(unlist(sig)), nrow(rows)))
nab <- list(); aab <- list(); pk <- function(a, b) paste(py_sort(c(tm(a), tm(b))), collapse = "\t")
for (z in cross) { k <- py_sort(c(tm(z$a), tm(z$b))); kk <- paste(k, collapse = "\t"); if (is.null(nab[[kk]])) { nab[[kk]] <- T_(paste0(drug(k[1]), " AND ", drug(k[2]))); if (nab[[kk]] >= 20) for (c in names(groups)) aab[[paste(kk, c, sep = "\r")]] <- T_(paste0(drug(k[1]), " AND ", drug(k[2]), " AND ", grp(c))) } }
ok(length(nab) == 70, "70 distinct search-term pairs retrieved"); ok(all(vapply(cross, function(z) nab[[pk(z$a, z$b)]] <= min(nd[[tm(z$a)]], nd[[tm(z$b)]]), TRUE)), "co-report counts never exceed either drug total")
pr <- rcsv(file.path(RES, "pairs_all77.csv")); lg <- rcsv(file.path(RES, "pairs_long.csv")); bad <- 0L; nflag <- 0L; nlong <- 0L; nevaluable <- 0L
for (z in cross) {
  tA <- tm(z$a); tB <- tm(z$b); kk <- pk(z$a, z$b); n_ab <- nab[[kk]]; r <- pr[as.integer(pr$rank) == z$rank, ][1, ]
  status <- if (n_ab >= 20 && ev[[tA]] && ev[[tB]]) "evaluable" else if (n_ab >= 20) "component drug below 100 reports" else "fewer than 20 co-reports"
  sgn <- function(t) if (ev[[t]]) { s <- names(groups)[vapply(names(groups), function(c) sg(t, c), TRUE)]; if (length(s)) paste(s, collapse = " ") else "-" } else "n/e"; sA <- sgn(tA); sB <- sgn(tB); flags <- character(0)
  if (status == "evaluable") { nevaluable <- nevaluable + 1L
    for (c in names(groups)) { x <- aab[[paste(kk, c, sep = "\r")]]; ci <- cp(x, n_ab); pA <- (adg(tA, c) - x) / (nd[[tA]] - n_ab); pB <- (adg(tB, c) - x) / (nd[[tB]] - n_ab); fl <- x >= 3 && ci[1] > max(pA, pB)
      if (fl) { flags <- c(flags, c); nflag <- nflag + 1L }; nlong <- nlong + 1L; row <- lg[as.integer(lg$rank) == z$rank & lg$group == c, ]
      if (!(nrow(row) == 1 && as.integer(row$reports_both_in_group) == x && as.integer(row$reports_both) == n_ab && close_(row$p_both_lo95, ci[1], 1e-3) && close_(row$p_both_hi95, ci[2], 1e-3) && as.integer(row$co_report_excess) == as.integer(fl))) bad <- bad + 1L } }
  fl_s <- if (status == "evaluable") (if (length(flags)) paste(flags, collapse = " ") else "-") else "n/a"
  if (!(as.integer(r$reports_listing_both) == n_ab && r$status == status && r$signals_drug_A == sA && r$signals_drug_B == sB && r$co_report_excess == fl_s && r$pair == z$name)) bad <- bad + 1L
}
ok(bad == 0, sprintf("every row of pairs_all77.csv and pairs_long.csv equals the independent rebuild (%d rows differ)", bad))
ok(nevaluable == S$pairs_evaluable && nlong == S$pair_group_tests && nflag == S$co_report_flags && length(S$co_report_flagged_pairs) == S$co_report_flagged_pair_count, sprintf("%d evaluable pairs, %d pair-group tests and %d co-report flags (in %d pairs) equal summary", nevaluable, nlong, nflag, S$co_report_flagged_pair_count))
ok(length(S$warnings) == 0, sprintf("analysis raised no consistency warning (%d)", length(S$warnings))); top <- rcsv(file.path(RES, "pairs_top15.csv")); ok(identical(as.integer(top$rank), 1:15), "pairs_top15.csv holds ranks 1 to 15")
for (c in names(groups)) { t15 <- sum(vapply(cross, function(z) z$rank <= 15 && (sg(tm(z$a), c) || sg(tm(z$b), c)), TRUE)); a77 <- sum(vapply(cross, function(z) sg(tm(z$a), c) || sg(tm(z$b), c), TRUE))
  ok(t15 == S$per_group[[c]]$pairs_top15_with_signalled_drug && a77 == S$per_group[[c]]$pairs_all77_with_signalled_drug, sprintf("group %s: pairs with a signalled drug, top 15 %d, all 77 %d", c, t15, a77)) }
## 5. fingerprints
pl <- strsplit(rd_text(file.path(V2, "PROTOCOL_v2_SHA256.txt")), "\n", fixed = TRUE)[[1]]; hashes <- vapply(strsplit(pl[endsWith(pl, "*wp9_safety/PRESPEC_WP9_safety_screen.md")], " *", fixed = TRUE), function(p) p[1], "")
ok(length(hashes) >= 2 && hashes[length(hashes)] == sha256_file(file.path(ARC9, "PRESPEC_WP9_safety_screen.md")), sprintf("pre-specification file equals its latest fingerprint (%d fingerprints recorded)", length(hashes)))
files <- py_sort(list.files(RAWD)); h <- digest::digest(paste0(vapply(files, function(f) sprintf("%s:%s\n", f, sha256_file(file.path(RAWD, f))), ""), collapse = ""), algo = "sha256", serialize = FALSE)
## 27 Nov 2026: RAW_SHA256.txt is read from the archived folder (wp9_safety/results), as the Python check reads it; results_R/wp9 held only a copy left by an earlier run
rs <- strsplit(trimws(rd_text(if (file.exists(file.path(ARC9, "results", "RAW_SHA256.txt"))) file.path(ARC9, "results", "RAW_SHA256.txt") else file.path(RES, "RAW_SHA256.txt"))), "\\s+")[[1]][1]; ok(h == rs && any(startsWith(pl, rs)), sprintf("raw files (%d) equal the fingerprint written before the analysis (RAW_SHA256.txt and PROTOCOL_v2_SHA256.txt)", length(files))); ok(all(endsWith(files, ".json")), "raw folder holds only json files")
## 6. secondary, post hoc free-text fallback (pre-specification section 9)
FBN <- list(albiglutide = c("TANZEUM", "ALBIGLUTIDE"), bepridil = c("VASCOR", "BEPRIDIL"), nilvadipine = "NILVADIPINE", nitrendipine = "NITRENDIPINE")
for (t in names(FBN)) ok(inn(sprintf("%s: %s", t, paste(FBN[[t]], collapse = ", ")), spec), sprintf("fallback names of %s equal the pre-specification text (%s)", t, paste(FBN[[t]], collapse = ", ")))
ok(identical(lapply(NAMES, as.character), lapply(FBN, as.character)), "wp9_safety.R fallback names equal the checker list"); low <- terms[nd < 100]; ok(identical(py_sort(low), py_sort(names(FBN))), paste0("drugs below 100 reports through substance_name: ", paste(py_sort(low), collapse = ", ")))
tx <- function(t) paste0("patient.drug.medicinalproduct:(", paste(sprintf("\"%s\"", FBN[[t]]), collapse = " OR "), ")")
F <- fromJSON(file.path(RES, "fallback_summary.json"), simplifyVector = FALSE); fsig <- rcsv(file.path(RES, "fallback_signals_single.csv")); fpair <- rcsv(file.path(RES, "fallback_pairs.csv")); bad <- 0L; n_ev_pairs <- 0L; n_flags <- 0L; fsigs <- list(); n_fb_rows <- 0L
for (t in py_sort(names(FBN))) {
  n_t <- T_(tx(t)); info <- F$drugs[[t]]; if (!(info$reports == n_t && identical(info$evaluable, n_t >= 100))) bad <- bad + 1L; if (n_t < 100) next
  a_t <- vapply(names(groups), function(c) T_(paste0(tx(t), " AND ", grp(c))), 1)
  for (c in names(groups)) { o <- lor_ci(a_t[[c]], n_t, E[[c]]); s <- a_t[[c]] >= 3 && o[2] > 1; fsigs[[paste(t, c, sep = "\r")]] <- s; row <- fsig[fsig$term == t & fsig$group == c, ]
    if (!(nrow(row) == 1 && as.integer(row$reports_in_group) == a_t[[c]] && close_(row$ROR, o[1]) && close_(row$ROR_lo95, o[2]) && close_(row$ROR_hi95, o[3]) && as.integer(row$signal) == as.integer(s))) bad <- bad + 1L }
  for (z in cross) { ta <- tm(z$a); tb <- tm(z$b); if (!(t %in% c(ta, tb))) next; p <- if (ta == t) tb else ta; n_ab <- T_(paste0(tx(t), " AND ", drug(p))); pev <- n_ab >= 20; cev <- nd[[p]] >= 100
    status <- if (pev && cev) "evaluable" else if (pev) "partner drug below 100 reports" else "fewer than 20 co-reports"; flags <- character(0)
    if (pev && cev) { n_ev_pairs <- n_ev_pairs + 1L; for (c in names(groups)) { x <- T_(paste0(tx(t), " AND ", drug(p), " AND ", grp(c))); ci <- cp(x, n_ab); pA <- (a_t[[c]] - x) / (n_t - n_ab); pB <- (adg(p, c) - x) / (nd[[p]] - n_ab); if (x >= 3 && ci[1] > max(pA, pB)) { flags <- c(flags, c); n_flags <- n_flags + 1L } } }
    sA <- { s <- names(groups)[vapply(names(groups), function(c) isTRUE(fsigs[[paste(t, c, sep = "\r")]]), TRUE)]; if (length(s)) paste(s, collapse = " ") else "-" }
    sP <- if (nd[[p]] >= 100) { s <- names(groups)[vapply(names(groups), function(c) sg(p, c), TRUE)]; if (length(s)) paste(s, collapse = " ") else "-" } else "n/e"
    row <- fpair[as.integer(fpair$rank) == z$rank & fpair$fallback_term == t, ]; n_fb_rows <- n_fb_rows + 1L
    if (!(nrow(row) == 1 && as.integer(row$reports_listing_both) == n_ab && row$status == status && row$signals_fallback_drug == sA && row$signals_partner_drug == sP && row$co_report_excess == (if (status == "evaluable") (if (length(flags)) paste(flags, collapse = " ") else "-") else "n/a"))) bad <- bad + 1L }
}
ok(bad == 0, sprintf("every secondary (fallback) row equals the independent rebuild (%d differ; %d pair rows)", bad, n_fb_rows)); ok(F$pairs_evaluable == n_ev_pairs && F$co_report_flags == n_flags && F$pairs == nrow(fpair) && nrow(fpair) == n_fb_rows, sprintf("secondary: %d pair rows, %d evaluable pairs and %d co-report flags equal fallback_summary", n_fb_rows, n_ev_pairs, n_flags))
## 7. numbers stated in README.md, rebuilt from the saved responses
RD <- rd_text(file.path(ARC9, "README.md")); has <- function(txt, msg) ok(inn(txt, RD), paste0("README states ", msg, ": \"", txt, "\""))
has(fmtc(N), "the number of reports in the database"); has(S$last_updated, "the snapshot date"); has(sprintf("%d responses", length(files)), "the number of saved responses"); has(sprintf("**%d of %d terms are evaluable**", sum(ev), length(terms)), "evaluable terms")
has(sprintf("**%d of %d drug-group tests are signals**", sum(unlist(sig)), nrow(rows)), "signal count"); for (t in terms) if (ev[[t]]) has(fmtc(nd[[t]]), paste0("report total of ", t))
bad_terms <- terms[!ev]; ok(identical(py_sort(bad_terms), c("albiglutide", "bepridil", "nilvadipine", "nitrendipine")), "not evaluable terms are albiglutide, bepridil, nilvadipine, nitrendipine")
inv <- Filter(function(z) tm(z$a) %in% bad_terms || tm(z$b) %in% bad_terms, cross); has(sprintf("**%d of the 77 pairs contain such a drug, and %d of the top 15**", length(inv), sum(vapply(inv, function(z) z$rank <= 15, TRUE))), "pairs with an invisible drug")
nn <- vapply(cross, function(z) nab[[pk(z$a, z$b)]], 1); has(sprintf("**%d have no report listing both drugs, %d have 1 to 19, and %d have 20 or more**", sum(nn == 0), sum(nn >= 1 & nn <= 19), sum(nn >= 20)), "co-report distribution")
evr <- Filter(function(z) nab[[pk(z$a, z$b)]] >= 20, cross); evn <- vapply(evr, function(z) nab[[pk(z$a, z$b)]], 1); evk <- vapply(evr, function(z) z$rank, 1L)
has(paste0("ranks ", paste(evk[-length(evk)], collapse = ", "), " and ", evk[length(evk)]), "ranks of the evaluable pairs"); has(sprintf("%s to %s co-reports, median %d", format(min(evn), scientific = FALSE), fmtc(max(evn)), as.integer(np_median(evn))), "range and median of the co-reports")
top15n <- vapply(Filter(function(z) z$rank <= 15, cross), function(z) nab[[pk(z$a, z$b)]], 1); has(sprintf("**Of the top 15 pairs, %d have fewer than 20 co-reports (%d have none", sum(top15n < 20), sum(top15n == 0)), "top 15 co-reports")
has(sprintf("**%d flags in %d of the %d evaluable pairs, among %d pair-group tests; none in the top 15**", S$co_report_flags, S$co_report_flagged_pair_count, S$pairs_evaluable, S$pair_group_tests), "flag counts")
for (fp in S$co_report_flagged_pairs) has(sprintf("rank %d %s (%s)", fp[[1]], fp[[2]], gsub(" ", ", ", fp[[3]], fixed = TRUE)), "flagged pair")
ord <- c("PANC", "HYPO", "HYPER", "GI", "BILI", "HYPOT", "CARD", "OEDEMA"); gc <- vapply(names(groups), function(c) sum(vapply(terms, function(t) sg(t, c), TRUE)), 1L); gc <- gc[order(-gc, match(names(gc), ord))]
has(paste0("By group, drugs with a signal: ", paste(sprintf("%s %d", names(gc), gc), collapse = ", ")), "drugs with a signal per group")
has(paste(vapply(names(groups), function(c) sprintf("%s %d / %d", c, S$per_group[[c]]$pairs_top15_with_signalled_drug, S$per_group[[c]]$pairs_all77_with_signalled_drug), ""), collapse = ", "), "pairs with a signalled drug (top 15 / all 77)")
by_drug <- vapply(terms[ev], function(t) sum(vapply(names(groups), function(c) sg(t, c), TRUE)), 1L); d8 <- py_sort(names(by_drug)[by_drug == 8L]); has(paste0("8 of 8 groups for ", paste(d8[-length(d8)], collapse = ", "), " and ", d8[length(d8)]), "drugs with 8 of 8")
ok(identical(py_sort(names(by_drug)[by_drug == 5L]), c("dulaglutide", "exenatide", "isradipine", "liraglutide", "nisoldipine", "semaglutide")) && identical(py_sort(names(by_drug)[by_drug == 3L]), c("clevidipine", "lixisenatide")) && identical(py_sort(names(by_drug)[by_drug == 2L]), "nimodipine"), "drugs with 5, 3 and 2 signalled groups as stated")
fmt <- function(t, c, d) { o <- lor_ci(adg(t, c), nd[[t]], E[[c]]); sprintf("%s (%s to %s)", py_fmt_f(o[1], d), py_fmt_f(o[2], d), py_fmt_f(o[3], d)) }
fa <- fmt("liraglutide", "PANC", 1); has(paste0("ROR ", strsplit(fa, " (", fixed = TRUE)[[1]][1], " (", strsplit(fa, " (", fixed = TRUE)[[1]][2]), "liraglutide PANC ROR"); has(fa, "liraglutide PANC interval")
has(sprintf("glucagon HYPO %d of %s, ROR %s", adg("glucagon", "HYPO"), fmtc(nd[["glucagon"]]), fmt("glucagon", "HYPO", 1)), "glucagon HYPO"); has(paste0("amlodipine OEDEMA ROR ", fmt("amlodipine", "OEDEMA", 2)), "amlodipine OEDEMA"); has(paste0("HYPOT ", fmt("amlodipine", "HYPOT", 2)), "amlodipine HYPOT")
has(paste0("HYPOT ", fmt("isradipine", "HYPOT", 1)), "isradipine HYPOT"); has(paste0("OEDEMA ", fmt("isradipine", "OEDEMA", 1)), "isradipine OEDEMA"); has(paste0("HYPOT ", fmt("liraglutide", "HYPOT", 2)), "liraglutide HYPOT"); has(paste0("OEDEMA ", fmt("liraglutide", "OEDEMA", 2)), "liraglutide OEDEMA")
has(sprintf("%s of %s reports", fmtc(adg("liraglutide", "PANC")), fmtc(nd[["liraglutide"]])), "liraglutide PANC counts")
for (i in seq_len(nrow(top))) { r <- top[i, ]; ca <- if (length(strsplit(trimws(r$signals_drug_A), "\\s+")[[1]]) == 8) "all 8 groups" else r$signals_drug_A; cb <- if (length(strsplit(trimws(r$signals_drug_B), "\\s+")[[1]]) == 8) "all 8 groups" else r$signals_drug_B
  fl <- if (r$co_report_excess == "n/a") "n/a" else if (r$co_report_excess == "-") "none" else r$co_report_excess; has(sprintf("| %s | %s | %s | %s | %s | %s | %s |", r$rank, r$pair, r$MLSS, r$reports_listing_both, ca, cb, fl), sprintf("top-15 row %s", r$rank)) }
## secondary
fb_n <- vapply(names(FBN), function(t) T_(tx(t)), 1); for (t in names(FBN)) has(fmtc(fb_n[[t]]), paste0("free-text report total of ", t))
for (z in list(c("albiglutide", "HYPO"), c("albiglutide", "HYPER"), c("albiglutide", "PANC"), c("albiglutide", "GI"), c("bepridil", "HYPOT"), c("nilvadipine", "OEDEMA"), c("nilvadipine", "CARD"), c("nitrendipine", "HYPOT"), c("nitrendipine", "BILI"), c("nitrendipine", "OEDEMA"))) {
  rr <- fsig[fsig$term == z[1] & fsig$group == z[2], ]; has(sprintf("%s (%s)", z[2], py_fmt_f(as.numeric(rr$ROR), 2L)), sprintf("free-text ROR %s %s", z[1], z[2])) }
has(sprintf("CARD (%s)", py_fmt_f(as.numeric(fsig[fsig$term == "bepridil" & fsig$group == "CARD", ]$ROR), 1L)), "bepridil CARD"); has(sprintf("CARD (%s)", py_fmt_f(as.numeric(fsig[fsig$term == "nitrendipine" & fsig$group == "CARD", ]$ROR), 1L)), "nitrendipine CARD")
fs <- function(t, c) isTRUE(fsigs[[paste(t, c, sep = "\r")]]); ok(!fs("nitrendipine", "HYPO") && all(vapply(setdiff(names(groups), "HYPO"), function(c) fs("nitrendipine", c), TRUE)), "nitrendipine has signals in 7 of 8 groups, all but HYPO")
has(sprintf("Of the %d pairs that contain one of these four drugs, **one has 20 or more co-reports (albiglutide + amlodipine, rank 21, 86 co-reports; excess flag GI)**", nrow(fpair)), "secondary pair counts")
big <- fpair[as.integer(fpair$reports_listing_both) >= 20, ]; ok(nrow(big) == 1 && identical(c(big$rank, big$pair, big$reports_listing_both, big$co_report_excess), c("21", "albiglutide + amlodipine", "86", "GI")), "the only secondary pair with 20 or more co-reports is rank 21 albiglutide + amlodipine, 86, flag GI")
s15 <- fpair[as.integer(fpair$rank) <= 15, ]; key15 <- paste(s15$rank, s15$fallback_term, sep = "|"); v15 <- setNames(as.integer(s15$reports_listing_both), key15); excl <- c("4|nitrendipine", "8|nitrendipine", "11|nitrendipine", "12|nitrendipine", "15|albiglutide")
ok(v15[["4|nitrendipine"]] == 5 && v15[["11|nitrendipine"]] == 4 && v15[["15|albiglutide"]] == 3 && max(v15[!(names(v15) %in% excl)]) == 0, "secondary co-reports in the top 15: nitrendipine + liraglutide 5, nitrendipine + semaglutide 4, albiglutide + nifedipine 3, all others 0")
## 8. WORDING_safety_screen.md: every number of sections 1 to 7 is a number of the README, and the framing words stay inside negations
WD <- rd_text(file.path(ARC9, "WORDING_safety_screen.md")); body <- strsplit(WD, "## 8. What not to write", fixed = TRUE)[[1]][1]
tok <- function(x) { x <- gsub("(?<=\\d),(?=\\d{3})", "", x, perl = TRUE); unique(regmatches(x, gregexpr("(?<![\\w.])\\d+(?:\\.\\d+)?", x, perl = TRUE))[[1]]) }; rd_tok <- tok(RD); NONRESULT <- c("256", "350", "5.5")
tb <- tok(body); missing <- py_sort(tb[(nchar(tb) >= 2 | grepl(".", tb, fixed = TRUE)) & !(tb %in% rd_tok) & !(tb %in% NONRESULT)]); ok(!length(missing), sprintf("every number of the wording (sections 1 to 7) is a number of the README (not found: %s)", paste(missing, collapse = ", ")))
sents <- unlist(strsplit(body, "(?<=[.;:])\\s+|\n", perl = TRUE)); bad_sent <- substr(sents[grepl("\\b(predict\\w*|safe|unsafe|risk|validated|confirms)\\b", sents, ignore.case = TRUE, perl = TRUE) & !grepl("\\b(not|no|never|neither|nor|without)\\b", sents, ignore.case = TRUE, perl = TRUE)], 1, 90)
ok(!length(bad_sent), sprintf("framing words (predict, safe, unsafe, risk, validated, confirms) only inside negations in sections 1 to 7 (offending: %s)", paste(bad_sent, collapse = " | ")))
res_txt <- regmatches(body, regexec("(?s)## 2\\. Results.*?\n(.*?)\n\n## 3", body, perl = TRUE))[[1]][2]; ok(lengths(regmatches(res_txt, gregexpr("\\S+", res_txt, perl = TRUE))) <= 330, "the Results text is at most 330 words")
cat(sprintf("\n%d checks, %d failures\n", count, length(fails)))
write_text_bytes(gsub("\n", .NL, sprintf("%d checks, %d failures%s\n", count, length(fails), if (length(fails)) paste0(": ", paste(fails, collapse = "; ")) else ""), fixed = TRUE), file.path(CHK, "check_wp9_last_run.txt"))
