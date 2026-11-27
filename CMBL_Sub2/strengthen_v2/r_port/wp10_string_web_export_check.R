## wp10_string_web_export_check.R - R port of wp10_versions/string_web_export_check.py: WP10 addendum (8 Nov 2026): the user's own STRING 12.5 web export of 8 Nov 2026 compared with the two STRING tables already used.
## Inputs (read only):
##   Database/STRING12.5/string_interactions_21_nodes.tsv   user's web export, 8 Nov 2026; network URL in STRING12.5.txt (version-12-5.string-db.org, networkId beSoaw2hyqcw)
##   R_17genes/string_interactions.csv                      STRING 12.0 web export that Fig. 5 used (original download date not recorded; copy date 2 Nov 2026)
##   results_R/wp10/string_v12_5_network.tsv                12.5 network retrieved through the API by wp10_string_sensitivity.R retrieve (rat 10116, 17 + 4 genes, minimum score 0.4), and string_summary.json of its analyze step
## Questions: (1) do the web export and the API network of 12.5 agree (edges, scores)? (2) what changed from 12.0 to 12.5 in the web exports? (3) does the Fig. 5 classification computed from the web export equal the API-based result?
## Usage: Rscript wp10_string_web_export_check.R   (writes string_web_export_check.json into R_OUT_DIR, default results_R/wp10; the other WP10 files are read from there, WP10_RETRIEVED_DIR overrides)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "wp10_string_sensitivity.R"))      # the functions (pname, string_set, classify, counts, with_cov, load_inputs, ALIAS), RES, RETR
P <- normalizePath(file.path(HERE, "..", "..", "..")); WEB <- file.path(P, "Database", "STRING12.5", "string_interactions_21_nodes.tsv"); URLF <- file.path(P, "Database", "STRING12.5", "STRING12.5.txt"); OUTF <- file.path(RES, "string_web_export_check.json")
G17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
read_web <- function(path) { x <- read.delim(path, sep = "\t", header = FALSE, quote = "\"", colClasses = "character", comment.char = "", na.strings = character(0), fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE); hdr <- sub("^#+", "", as.character(unlist(x[1, ]))); x <- x[-1, , drop = FALSE]; names(x) <- hdr; rownames(x) <- NULL; x }
web <- read_web(WEB); old <- load_inputs()$old
lines <- strsplit(trimws(paste(readLines(file.path(RETR, "string_v12_5_network.tsv"), encoding = "UTF-8", warn = FALSE), collapse = "\n")), "\n", fixed = TRUE)[[1]]; ah <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
api <- as.data.frame(do.call(rbind, lapply(lines[-1], function(l) { p <- strsplit(l, "\t", fixed = TRUE)[[1]]; length(p) <- length(ah); p })), stringsAsFactors = FALSE); names(api) <- ah
out <- list(network_url_file = trimws(paste(readLines(URLF, warn = FALSE, encoding = "UTF-8"), collapse = "\n")))
summarize <- function(rows, a = "node1", b = "node2", sc = "combined_score") {
  nodes <- py_sort(unique(c(rows[[a]], rows[[b]]))); keys <- pname(rows[[a]], rows[[b]]); v <- as.numeric(rows[[sc]]); by <- lapply(split(v, keys), function(x) sort(unique(x))); list(nodes = nodes, by = by)
}
w <- summarize(web); o <- summarize(old); a <- summarize(api, "preferredName_A", "preferredName_B", "score")
out$rows <- list(web_12_5 = nrow(web), old_12_0 = nrow(old), api_12_5 = nrow(api)); out$nodes <- list(web_12_5 = as.list(w$nodes), old_12_0 = as.list(o$nodes), api_12_5 = as.list(a$nodes))
out$distinct_pairs <- list(web_12_5 = length(w$by), old_12_0 = length(o$by), api_12_5 = length(a$by))
two <- function(s) as.list(py_sort(names(s$by)[vapply(s$by, length, 1L) > 1])); out$pairs_with_two_different_scores_in_one_file <- list(web_12_5 = two(w), old_12_0 = two(o), api_12_5 = two(a))
mn <- function(s) vapply(s$by, function(x) x[1], 1); W <- mn(w); O <- mn(o); A <- mn(a)
out$min_score_in_file <- list(web_12_5 = min(W), old_12_0 = min(O), api_12_5 = min(A))
cwa <- intersect(names(W), names(A))
out$web_vs_api_12_5 <- list(pairs_only_in_web = as.list(py_sort(setdiff(names(W), names(A)))), pairs_only_in_api = as.list(py_sort(setdiff(names(A), names(W)))), pairs_in_both = length(cwa), max_abs_score_difference_in_common_pairs = py_round(if (length(cwa)) max(abs(W[cwa] - A[cwa])) else 0, 4L))
common <- intersect(names(W), names(O)); dd <- W[common] - O[common]; od <- order(dd, common, method = "radix"); dd <- dd[od]; cm <- common[od]
out$web_12_5_vs_old_12_0 <- list(pairs_gained_in_12_5 = as.list(py_sort(setdiff(names(W), names(O)))), pairs_lost_in_12_5 = as.list(py_sort(setdiff(names(O), names(W)))), pairs_in_both = length(common),
                                 score_change_in_common_pairs = list(n_changed = as.integer(sum(abs(dd) > 0.0005)), min = py_round(dd[[1]], 3L), max = py_round(dd[[length(dd)]], 3L), min_pair = cm[1], max_pair = cm[length(cm)]))
x <- load_inputs(); pairs <- x$pairs; cons <- x$cons; res <- list()
for (nm in c("as_reported", "alias")) { cl <- classify(pairs, cons, string_set(web$node1, web$node2, nm == "alias")); res[[nm]] <- with_cov(counts(cl), length(pairs)) }
out$classification_from_web_export_12_5 <- res
## 27 Nov 2026: the classification is read from the output of this port (RES); RETR holds the retrieved network (the frozen retrieval when run_all_v2.R runs it)
sj <- fromJSON(file.path(RES, "string_summary.json"), simplifyVector = FALSE); out$classification_from_api_12_5_saved <- list(as_reported = sj$as_reported$v12_5, alias = sj$alias$v12_5)
same_cls <- function(p, q) identical(names(unlist(p)), names(unlist(q))) && identical(as.numeric(unlist(p)), as.numeric(unlist(q))); out$classification_web_equals_api <- same_cls(res$as_reported, sj$as_reported$v12_5) && same_cls(res$alias, sj$alias$v12_5)
e17 <- function(a, b) { a2 <- al(a); b2 <- al(b); k <- a2 %in% G17 & b2 %in% G17; unique(pname(a2[k], b2[k])) }
ew <- e17(web$node1, web$node2); ea <- e17(api$preferredName_A, api$preferredName_B); eo <- e17(old$node1, old$node2)
out$edges_among_17_genes <- list(web_12_5 = length(ew), api_12_5 = length(ea), old_12_0 = length(eo), web_equals_api = setequal(ew, ea), gained_vs_12_0 = as.list(py_sort(setdiff(ew, eo))), lost_vs_12_0 = as.list(py_sort(setdiff(eo, ew))))
write_json_typed(out, OUTF); cat("written", OUTF, "\n")
