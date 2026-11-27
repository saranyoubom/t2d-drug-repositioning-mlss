## wp10_string_sensitivity.R - R port of wp10_versions/string_sensitivity.py. WP10 part A: does the STRING classification of Fig. 5 change with the current STRING release (12.5) instead of 12.0?
## Reproduces the classification of R_17genes/Fig_5abcd_Gene_Interaction_Network_Analysis.r from the reported inputs (read only): string_interactions.csv (STRING 12.0 web export, rat, minimum score 0.4),
## Fig_4d_consensus_pairs.csv, Fig_4e_Hes1_significant_correlations_only.csv. Categories (R code): Known Validated = in STRING and four-method consensus; Novel Discovery = not in STRING but consensus;
## STRING Only = in STRING, not consensus; Not Validated = neither, over the pairs with |r| >= 0.7 in the Hes1-normalized significant correlations (120 pairs). Pair names: the two gene names sorted and joined by a hyphen.
## NOTE: the correlation matrix names the gene Glut2, the STRING export Slc2a2, and the R code matches names literally, so no pair with Glut2 can be 'in STRING' in the reported analysis; both versions are reported:
## 'as reported' (literal names) and 'alias' (Slc2a2 read as Glut2).
## Usage: Rscript wp10_string_sensitivity.R gate | retrieve | analyze   (gate and analyze read only saved files; retrieve needs the network).
## Environment: R_OUT_DIR (output folder, default results_R/wp10), WP10_RETRIEVED_DIR (folder with string_v12_5_network.tsv and string_version.json for analyze, default R_OUT_DIR).
suppressMessages({ library(jsonlite); library(curl) })
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "pycompat.R"))
RDIR <- Sys.getenv("R17_DIR", normalizePath(file.path(HERE, "..", "..", "R_17genes"))); RES <- Sys.getenv("R_OUT_DIR", file.path(HERE, "results_R", "wp10")); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
RETR <- Sys.getenv("WP10_RETRIEVED_DIR", RES); STABLE <- "https://version-12-5.string-db.org"; ALIAS <- c(Slc2a2 = "Glut2")
al <- function(x) ifelse(x %in% names(ALIAS), ALIAS[x], x)
pname <- function(a, b) vapply(seq_along(a), function(i) paste(py_sort(c(a[i], b[i])), collapse = "-"), "")
read_csv_bom <- function(f, header = TRUE) read.csv(f, header = header, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8-BOM", na.strings = character(0))
load_inputs <- function() {
  old <- read_csv_bom(file.path(RDIR, "string_interactions.csv"))
  cp <- read_csv_bom(file.path(RDIR, "Fig_4d_consensus_pairs.csv"))$Consensus_Pairs; cons <- unique(vapply(strsplit(cp, "_", fixed = TRUE), function(s) pname(s[1], s[2]), ""))
  m <- read_csv_bom(file.path(RDIR, "Fig_4e_Hes1_significant_correlations_only.csv"), header = FALSE); hdr <- unlist(m[1, -1], use.names = FALSE); pairs <- list()
  for (i in 2:nrow(m)) for (j in seq_along(hdr)) { v <- m[i, j + 1]; if (!(v %in% c("", "NA")) && m[i, 1] != hdr[j] && abs(as.numeric(v)) >= 0.7) pairs[[pname(m[i, 1], hdr[j])]] <- as.numeric(v) }
  list(old = old, cons = cons, pairs = pairs)
}
string_set <- function(a, b, alias) { if (alias) { a <- al(a); b <- al(b) }; unique(pname(a, b)) }
classify <- function(pairs, cons, sset) vapply(names(pairs), function(p) { ins <- p %in% sset; val <- p %in% cons; if (ins && val) "Known Validated" else if (!ins && val) "Novel Discovery" else if (ins) "STRING Only" else "Not Validated" }, "")
counts <- function(cl) { cats <- c("Known Validated", "Novel Discovery", "STRING Only", "Not Validated"); setNames(lapply(cats, function(k) as.integer(sum(cl == k))), cats) }
with_cov <- function(d, npairs) { d$in_STRING <- d$`Known Validated` + d$`STRING Only`; d$coverage_percent <- round(100.0 * d$in_STRING / npairs, 1); d }
gate <- function() {
  x <- load_inputs(); pairs <- x$pairs; cons <- x$cons; a <- x$old$node1; b <- x$old$node2
  out <- list(`pairs_abs_r_ge_0.7` = length(pairs), consensus_pairs = length(cons), string_export_edges = length(a), string_export_distinct_pairs = length(string_set(a, b, FALSE)))
  for (nm in c("as_reported", "alias")) out[[nm]] <- with_cov(counts(classify(pairs, cons, string_set(a, b, nm == "alias"))), length(pairs))
  out$pairs_with_Glut2 <- sum(vapply(names(pairs), function(p) "Glut2" %in% strsplit(p, "-", fixed = TRUE)[[1]], TRUE))
  sl <- a == "Slc2a2" | b == "Slc2a2"; out$export_edges_with_Slc2a2_among_matrix_genes <- as.list(py_sort(pname(a[sl], b[sl])))
  write_json_typed(out, file.path(RES, "string_gate.json")); cat(py_json_typed(out), "\n"); rep <- out$as_reported
  ok <- length(pairs) == 120 && rep$`Known Validated` == 5 && rep$`Novel Discovery` == 16 && rep$`STRING Only` == 29 && rep$`Not Validated` == 70 && rep$in_STRING == 34
  cat(sprintf("GATE (reported numbers 120 pairs; 5 / 16 / 29 / 70; 34 in STRING = 28.3%%): %s\n", if (ok) "PASS" else "FAIL")); ok
}
GENES <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Slc2a2", "Ins1", "Glp1r", "Wnt2", "Wnt2b", "Hey1", "Ryr1")
get_url <- function(url, tries = 6) { delay <- 3
  for (k in seq_len(tries)) {
    r <- tryCatch(curl_fetch_memory(url, new_handle(useragent = "research-version-check/1.0", timeout = 90)), error = function(e) NULL)
    if (!is.null(r) && r$status_code == 200) return(rawToChar(r$content)); cat(sprintf("  retry %d\n", k)); Sys.sleep(delay); delay <- min(delay * 2, 60)
  }
  stop("giving up on ", url)
}
retrieve <- function() {
  v <- get_url(paste0(STABLE, "/api/json/version")); write_text_bytes(v, file.path(RES, "string_version.json")); cat("version endpoint:", trimws(v), "\n")
  ids <- paste(vapply(GENES, function(g) URLencode(g, reserved = TRUE), ""), collapse = "%0d")
  t <- get_url(paste0(STABLE, "/api/tsv/network?identifiers=", ids, "&species=10116&required_score=400&caller_identity=wp10_version_check")); write_text_bytes(t, file.path(RES, "string_v12_5_network.tsv"))
  cat("rows (with header):", length(strsplit(trimws(t), "\n", fixed = TRUE)[[1]]), "\n")
}
analyze <- function() {
  x <- load_inputs(); pairs <- x$pairs; cons <- x$cons; oa <- x$old$node1; ob <- x$old$node2
  lines <- strsplit(trimws(paste(readLines(file.path(RETR, "string_v12_5_network.tsv"), encoding = "UTF-8", warn = FALSE), collapse = "\n")), "\n", fixed = TRUE)[[1]]; hdr <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  cc <- lapply(lines[-1], function(l) strsplit(l, "\t", fixed = TRUE)[[1]]); na <- vapply(cc, function(c) c[match("preferredName_A", hdr)], ""); nb <- vapply(cc, function(c) c[match("preferredName_B", hdr)], ""); sc <- vapply(cc, function(c) as.numeric(c[match("score", hdr)]), 1)
  nscore <- setNames(sc, pname(na, nb)); nscore <- nscore[!duplicated(names(nscore), fromLast = TRUE)]; oscore <- setNames(as.numeric(x$old$combined_score), pname(oa, ob)); oscore <- oscore[!duplicated(names(oscore), fromLast = TRUE)]
  out <- list(version_response = fromJSON(file.path(RETR, "string_version.json"), simplifyVector = FALSE), edges_v12_0 = length(oscore), edges_v12_5 = length(nscore)); cls <- list()
  for (nm in c("as_reported", "alias")) {
    al_ <- nm == "alias"; so <- string_set(oa, ob, al_); sn <- string_set(na, nb, al_); c0 <- classify(pairs, cons, so); c1 <- classify(pairs, cons, sn); cls[[nm]] <- list(c0, c1)
    out[[nm]] <- list(v12_0 = with_cov(counts(c0), length(pairs)), v12_5 = with_cov(counts(c1), length(pairs))); pn <- names(pairs)
    out[[nm]]$pairs_gained <- as.list(py_sort(pn[pn %in% sn & !(pn %in% so)])); out[[nm]]$pairs_lost <- as.list(py_sort(pn[pn %in% so & !(pn %in% sn)]))
    ch <- pn[c0 != c1]; out[[nm]]$category_changes <- setNames(lapply(ch, function(p) list(c0[[p]], c1[[p]])), ch)
  }
  g17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
  e17 <- function(a, b) { a2 <- al(a); b2 <- al(b); unique(pname(a2[a2 %in% g17 & b2 %in% g17], b2[a2 %in% g17 & b2 %in% g17])) }
  a0 <- e17(oa, ob); a1 <- e17(na, nb); out$edges_among_17_genes <- list(v12_0 = length(a0), v12_5 = length(a1), gained = as.list(py_sort(setdiff(a1, a0))), lost = as.list(py_sort(setdiff(a0, a1)))); out$pairs <- length(pairs)
  pn <- py_sort(names(pairs)); sfmt <- function(s, p) if (p %in% names(s)) sprintf("%.3f", s[[p]]) else ""
  rows <- lapply(pn, function(p) list(p, .json_num(pairs[[p]]), as.character(as.integer(p %in% cons)), sfmt(oscore, p), sfmt(nscore, p), cls$as_reported[[1]][[p]], cls$as_reported[[2]][[p]], cls$alias[[1]][[p]], cls$alias[[2]][[p]]))
  hdr_out <- c("pair", "r_Hes1_normalized", "four_method_consensus", "score_v12_0_literal_names", "score_v12_5_literal_names", "category_as_reported_v12_0", "category_as_reported_v12_5", "category_alias_v12_0", "category_alias_v12_5")
  write_text_bytes(paste0(paste(c(paste(hdr_out, collapse = ","), vapply(rows, function(r) paste(vapply(r, .dw_quote, ""), collapse = ","), "")), collapse = "\r\n"), "\r\n"), file.path(RES, "string_classification.csv"))
  write_json_typed(out, file.path(RES, "string_summary.json")); cat(substr(py_json_typed(out[names(out) != "version_response"]), 1, 3500), "\n")
}
if (sys.nframe() == 0L) {                                          # as a script: gate | retrieve | analyze; sourced by another script (wp13_gapdh_versions.R) only the functions are defined
  cmd <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""
  if (cmd == "gate") quit(status = if (gate()) 0L else 1L) else if (cmd == "retrieve") retrieve() else if (cmd == "analyze") analyze() else cat("usage: gate | retrieve | analyze\n")
}
