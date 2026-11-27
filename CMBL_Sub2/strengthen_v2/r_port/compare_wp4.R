## compare_wp4.R - agreement of the R port of wp4_islets with the archived Python results. Usage: Rscript compare_wp4.R [meta | ...]   (WP4_OUT = folder of the R meta-analysis, default results_R/wp4/meta)
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(R_PORT, "agree.R"))
ARC <- normalizePath(file.path(R_PORT, "..", "wp4_islets")); secs <- commandArgs(trailingOnly = TRUE); if (!length(secs)) secs <- "meta"
same_bytes <- function(a, b, label) { x <- readBin(a, "raw", file.size(a)); y <- readBin(b, "raw", file.size(b)); cat(sprintf("%-6s %-46s %s\n", if (identical(x, y)) "SAME" else "BYTES", label, if (identical(x, y)) "byte-identical" else sprintf("bytes differ (%d vs %d)", length(x), length(y)))) }
if ("meta" %in% secs) {
  OUT <- Sys.getenv("WP4_OUT", file.path(R_PORT, "results_R", "wp4", "meta")); A <- file.path(ARC, "meta"); cat("== wp4_meta_analysis.R ==\n")
  for (f in list.files(A, pattern = "[.]csv$")) if (file.exists(file.path(OUT, f))) agree_csv(file.path(A, f), file.path(OUT, f), tol = 1e-9, label = f)
  agree_json(file.path(A, "meta_summary.json"), file.path(OUT, "meta_summary.json"), tol = 1e-9, label = "meta_summary.json")
  for (f in c("meta_report.md", "meta_summary.json", list.files(A, pattern = "[.]csv$"))) if (file.exists(file.path(OUT, f))) same_bytes(file.path(A, f), file.path(OUT, f), f)
}
if ("results" %in% secs) {      # run_datasets.R: every <tag>__<variant>__* file of results_R/wp4/results against the archived results/ (tags: WP4_TAGS, comma separated; default all tags that have a summary in the R run)
  OUT <- Sys.getenv("WP4_RES_OUT", file.path(R_PORT, "results_R", "wp4", "results")); A <- file.path(ARC, "results"); cat("== run_datasets.R ==\n"); tol <- as.numeric(Sys.getenv("WP4_TOL", "1e-9"))
  tags <- Sys.getenv("WP4_TAGS", ""); tags <- if (nzchar(tags)) strsplit(tags, ",", fixed = TRUE)[[1]] else unique(sub("__.*$", "", list.files(OUT, pattern = "__summary[.]json$")))
  for (tag in tags) for (f in sort(list.files(A, pattern = paste0("^", tag, "__")))) if (file.exists(file.path(OUT, f))) { if (grepl("[.]csv$", f)) agree_csv(file.path(A, f), file.path(OUT, f), tol = tol, label = f) else agree_json(file.path(A, f), file.path(OUT, f), tol = tol, label = f) } else cat("MISSING", f, "\n")
  if (Sys.getenv("WP4_INDEX", "") == "1") { agree_csv(file.path(A, "index.csv"), file.path(OUT, "index.csv"), tol = tol, label = "index.csv (six datasets)"); same_bytes(file.path(A, "index.csv"), file.path(OUT, "index.csv"), "index.csv") }
}
cmp_dir <- function(title, folder, files, env, default_sub, tol = 1e-9) {      # CSV / JSON files of one WP4 sub-folder against the archive
  OUT <- Sys.getenv(env, file.path(R_PORT, "results_R", "wp4", default_sub)); A <- file.path(ARC, folder); cat("==", title, "==\n")
  for (f in files) if (!file.exists(file.path(OUT, f))) cat("MISSING", f, "\n") else if (grepl("[.]csv$", f)) agree_csv(file.path(A, f), file.path(OUT, f), tol = tol, label = f) else if (grepl("[.]json$", f)) agree_json(file.path(A, f), file.path(OUT, f), tol = tol, label = f) else same_bytes(file.path(A, f), file.path(OUT, f), f)
}
if ("perturbation" %in% secs) cmp_dir("perturbation_gse60803.R", "perturbation", c("gse60803_fold_changes.csv", "gse60803_concordance.json"), "WP4_PERT_OUT", "perturbation")
if ("wider" %in% secs) cmp_dir("wider_search_analyze.R", "wider_search", c("fold_changes_GSE60803.csv", "fold_changes_GSE99503.csv", "fold_changes_GSE167073.csv", "concordance_results.csv", "concordance_meta.csv", "concordance_report.md"), "WP4_WS_OUT", "wider_search")
if ("meta_graph" %in% secs) cmp_dir("gene_graph_null.R", "meta", c("gene_graph_by_dataset.csv", "gene_graph_meta.csv", "gene_graph_report.md"), "WP4_OUT", "meta")
if ("metacell" %in% secs) cmp_dir("metacell.R", "metacell", c("metacell_runs.csv", "metacell_summary.csv", "metacell_gene_evaluability.csv", "metacell_pair_evaluability.csv"), "WP4_MC_OUT", "metacell")
if ("scout" %in% secs) {      # the live-query scripts against the archived runs of 5-6 Nov 2026: records, fields (a field that differs only in the ORDER of its ';'-separated items is a server-side re-ordering), query texts, counts; the sample records byte for byte
  suppressMessages(library(jsonlite)); W4 <- file.path(R_PORT, "results_R", "wp4"); cat("== perturbation_scout.R, wider_search.R, wider_search_meta.R ==\n")
  fieldcmp <- function(rp, ap, label) {
    if (!file.exists(rp)) { cat("MISSING", label, "\n"); return(invisible()) }
    rd <- function(p) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = character(0)); x <- rd(rp); y <- rd(ap); ox <- setdiff(x$accession, y$accession); oy <- setdiff(y$accession, x$accession); cm <- intersect(x$accession, y$accession)
    nd <- 0L; no <- 0L; msgs <- character(0); sset <- function(s) sort(trimws(strsplit(s, ";", fixed = TRUE)[[1]]), method = "radix")
    for (col in setdiff(intersect(names(x), names(y)), "accession")) for (a in cm) { u <- x[[col]][match(a, x$accession)]; v <- y[[col]][match(a, y$accession)]
      if (!identical(u, v)) { if (identical(sset(u), sset(v))) { no <- no + 1L; msgs <- c(msgs, sprintf("   order only  %s %s: R '%s' | archive '%s'", a, col, substr(u, 1, 50), substr(v, 1, 50))) } else { nd <- nd + 1L; msgs <- c(msgs, sprintf("   DIFFERENT   %s %s: R '%s' | archive '%s'", a, col, substr(u, 1, 50), substr(v, 1, 50))) } } }
    cat(sprintf("%-6s %-26s %d records (archive %d); only in R %d, only in archive %d; %d fields differ in content, %d only in the order of the items\n", if (!length(ox) && !length(oy) && nd == 0) "AGREE" else "DIFFER", label, nrow(x), nrow(y), length(ox), length(oy), nd, no)); for (m in msgs) cat(m, "\n")
  }
  fieldcmp(file.path(W4, "perturbation_scout", "candidates.csv"), file.path(ARC, "perturbation", "candidates.csv"), "perturbation candidates.csv"); fieldcmp(file.path(W4, "wider_search_scout", "candidates_all.csv"), file.path(ARC, "wider_search", "candidates_all.csv"), "wider candidates_all.csv")
  same_bytes(file.path(W4, "perturbation_scout", "queries.json"), file.path(ARC, "perturbation", "queries.json"), "perturbation queries.json")
  rj <- fromJSON(file.path(W4, "wider_search_scout", "queries.json"), simplifyVector = FALSE); aj <- fromJSON(file.path(ARC, "wider_search", "queries.json"), simplifyVector = FALSE)
  cat(sprintf("%-6s wider queries.json: query texts (terms, S1-S3) %s; S1 count %d vs %d, ids %d vs %d; S2 hits %d vs %d; S3 PMIDs %d vs %d (cap %d), linked GEO records %d vs %d; unique records %d vs %d\n", if (identical(rj$terms, aj$terms) && identical(rj$S1$query, aj$S1$query) && identical(rj$S2$query, aj$S2$query) && identical(rj$S3$query, aj$S3$query)) "AGREE" else "DIFFER",
              if (identical(rj$terms, aj$terms) && identical(rj$S1$query, aj$S1$query) && identical(rj$S2$query, aj$S2$query) && identical(rj$S3$query, aj$S3$query)) "identical" else "DIFFER", rj$S1$count, aj$S1$count, rj$S1$ids, aj$S1$ids, rj$S2$hits, aj$S2$hits, rj$S3$pmids, aj$S3$pmids, rj$S3$cap, rj$S3$linked_gds_uids, aj$S3$linked_gds_uids, rj$unique_records, aj$unique_records))
  for (f in list.files(file.path(ARC, "wider_search"), pattern = "^sample_records_.*[.]tsv$")) if (file.exists(file.path(W4, "wider_search_meta", f))) same_bytes(file.path(ARC, "wider_search", f), file.path(W4, "wider_search_meta", f), f) else cat("MISSING", f, "\n")
}
if ("positive" %in% secs) {      # positive_control.R: the CSV and Markdown files, and the replicate arrays inj_<tag>.rds against the archived inj_<tag>.npz
  cmp_dir("positive_control.R", "positive_control", c("positive_control_by_dataset.csv", "positive_control_pooled.csv", "rp_control.csv", "margin_sensitivity.csv", "positive_control_report.md", paste0("rp_", c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061", "gse84133_mouse"), ".json")), "WP4_PC_OUT", "positive_control")
  source(file.path(R_PORT, "oracle", "read_npz.R")); OUT <- Sys.getenv("WP4_PC_OUT", file.path(R_PORT, "results_R", "wp4", "positive_control"))
  for (t in c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061", "gse84133_mouse")) { f <- file.path(OUT, sprintf("inj_%s.rds", t)); if (!file.exists(f)) next; a <- read_npz(file.path(ARC, "positive_control", sprintf("inj_%s.npz", t))); b <- readRDS(f)
    d <- max(abs(a$Erep - b$Erep)); dp <- max(abs(a$Pemp - b$Pemp)); cat(sprintf("%-6s inj_%s: Erep %s values, largest difference %.2g; Pemp %.2g; E_primary %.2g; SE %.2g\n", if (d < 1e-9 && dp < 1e-9) "AGREE" else "DIFFER", t, length(a$Erep), d, dp, abs(a$E_primary - b$E_primary), abs(a$SE - b$SE))) }
}
