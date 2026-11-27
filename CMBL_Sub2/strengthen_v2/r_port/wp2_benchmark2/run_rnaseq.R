## run_rnaseq.R - R port of wp2_benchmark2/run_rnaseq.py: protocol v2 section 3.5, held-out RNA-seq arm (WP2 deliverable E): counts, dry run and full run on the frozen list (frozen_heldout_rnaseq.csv,
## fingerprinted before any count was read). Same structure and code path as run_heldout.R (bench_data_io.R prepare_rnaseq -> prepare_from_log -> bench2_core.R analyse_matrix).
## Usage: Rscript run_rnaseq.R dry | full
##   For every frozen series in order (per organism): download the NCBI-generated raw count matrix <GSE>_raw_counts_<build>_NCBI.tsv.gz (human GRCh38.p13, mouse GRCm39; cached under
##   <data>/v2_data/rnaseq/counts, untrusted data read as numbers only). A 404 means the series has no NCBI-generated counts and is recorded as such. A reply that is a web page instead of a count
##   file is recorded as blocked and is NOT worked around: after 5 such replies in a row the organism is stopped ('not done'). The 12 selected samples are taken by GSM accession; CPM on all genes of
##   the file, genes need count >= 1 in all 12 samples, expressed and DE as in part 1 on log2(CPM + 1); a series needs >= 34 DE genes. Gene ids (Entrez GeneID) are the gene names: the draws of the
##   17 genes are random, so no annotation is needed. Stop at 20 analysable series per organism or when the list is exhausted. Duplicate rule as in part 1: a series sharing a selected sample with an
##   earlier analysable series is excluded.
##   dry   decides analysability only (counts available, samples present, DE genes >= 34); writes heldout_rnaseq_dry.jsonl; no correlation is computed
##   full  same loop with the analysis of gene sets A and B; one JSON line per series in heldout_rnaseq_results.jsonl; every decision is checked against the dry run (a mismatch is logged)
## WP2_FROZEN_RNASEQ = the frozen list (default: the archived file); WP2_ONLY = comma-separated series for checks on a few series.
if (!exists("HERE") || !file.exists(file.path(HERE, "wp2_common.R"))) { HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, "..")) }
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "coexp_checks.R")); source(file.path(R_PORT, "cc_fast.R")); source(file.path(R_PORT, "bench_pipeline.R")); source(file.path(HERE, "bench_data_io.R")); source(file.path(HERE, "bench2_core.R"))
MAXP <- 20L; MAX_BLOCKED_IN_A_ROW <- 5L
read_counts <- read_counts_file      # genes x 12 matrix of the selected samples (by accession) and the gene ids (bench_data_io.R)
main <- function(dry) {
  rows <- read.csv(Sys.getenv("WP2_FROZEN_RNASEQ", file.path(ARC, "frozen_heldout_rnaseq.csv")), stringsAsFactors = FALSE, check.names = FALSE); only <- setdiff(strsplit(Sys.getenv("WP2_ONLY", ""), ",", fixed = TRUE)[[1]], ""); if (length(only)) rows <- rows[rows$gse %in% only, , drop = FALSE]
  orgs <- py_sort(unique(rows$organism)); dry_path <- file.path(OUT, "heldout_rnaseq_dry.jsonl"); ref <- list(); if (!dry && file.exists(dry_path)) for (j in read_jsonl(dry_path)) ref[[j$gse]] <- j
  fn <- if (dry) dry_path else file.path(OUT, "heldout_rnaseq_results.jsonl"); if (file.exists(fn)) file.remove(fn)
  done <- setNames(rep(0L, length(orgs)), orgs); used <- list(); blocked <- setNames(rep(0L, length(orgs)), orgs); stopped <- setNames(rep(FALSE, length(orgs)), orgs); mism <- 0L; tally <- setNames(lapply(orgs, function(o) list()), orgs)
  for (i in seq_len(nrow(rows))) {
    row <- as.list(rows[i, ]); org <- row$organism; if (done[[org]] >= MAXP || stopped[[org]]) next
    t0 <- Sys.time(); gs <- strsplit(row$gsms, ";", fixed = TRUE)[[1]]; clash <- py_sort(unique(unlist(used[gs[gs %in% names(used)]]))); sets <- NULL
    rec <- list(gse = row$gse, platform = org, organism = org, order = as.integer(row$order), how = row$how, groups = row$groups)
    if (length(clash)) rec$status <- paste0("not analysable: duplicate of ", clash[1])
    else {
      fc <- fetch_counts(DATA, row$gse, org)
      if (is.null(fc$path)) {
        rec$status <- paste0("not analysable: ", fc$why); blocked[[org]] <- if (startsWith(fc$why, "blocked")) blocked[[org]] + 1L else 0L
        if (blocked[[org]] >= MAX_BLOCKED_IN_A_ROW) { stopped[[org]] <- TRUE; rec$status <- paste0(rec$status, sprintf(" (organism stopped after %d in a row: not done, no work-around)", MAX_BLOCKED_IN_A_ROW)) }
      } else {
        blocked[[org]] <- 0L
        r <- tryCatch({ cc <- read_counts(fc$path, gs); prepare_rnaseq(rec, cc$C, cc$ids, NULL, as.integer(substring(row$gse, 4)), dry) },
                      error = function(e) { rec$status <- sprintf("not analysable: %s %s", class(e)[1], substr(conditionMessage(e), 1, 120)); list(rec = rec, sets = NULL) })
        rec <- r$rec; sets <- r$sets
      }
    }
    if (!is.null(sets) && !dry) {
      tryCatch({
        for (tag in c("A", "B")) { res <- analyse_matrix(sets[[tag]]$X, sets[[tag]]$genes, as.integer(substring(row$gse, 4))); rec[[tag]] <- list(flat = res$flat, detail = res$detail, genes = as.list(sets[[tag]]$genes), n_de_in_set = sets[[tag]]$n_de_in_set) }
        rec$status <- "analysed"
      }, error = function(e) { rec$status <<- sprintf("not analysable: analysis error %s %s", class(e)[1], substr(conditionMessage(e), 1, 120)); rec$A <<- NULL; rec$B <<- NULL })
    }
    if (rec$status %in% c("analysed", "analysable (dry run)")) { done[[org]] <- done[[org]] + 1L; for (g in gs) used[[g]] <- row$gse }
    key <- substr(rec$status, 1, 60); tally[[org]][[key]] <- (if (is.null(tally[[org]][[key]])) 0L else tally[[org]][[key]]) + 1L
    if (length(ref)) {
      r0 <- ref[[row$gse]]; a0 <- if (is.null(r0)) NA else identical(r0$status, "analysable (dry run)"); a1 <- identical(rec$status, "analysed")
      if (is.null(r0) || !identical(a0, a1) || !identical(r0$n_de, rec$n_de)) { mism <- mism + 1L; rec$dry_run_mismatch <- TRUE }
    }
    cat(py_json_typed(rec, indent = NULL), "\n", file = fn, sep = "", append = TRUE)
    cat(sprintf("%s %s (#%s): %s | DE %s | done {%s} | %.0fs\n", org, row$gse, row$order, substr(rec$status, 1, 110), if (is.null(rec$n_de)) "-" else rec$n_de, paste(sprintf("'%s': %d", names(done), done), collapse = ", "), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  cat("finished {", paste(sprintf("'%s': %d", names(done), done), collapse = ", "), "} dry-run mismatches:", mism, "\n")
  for (o in orgs) cat("status counts", o, py_json_typed(tally[[o]], indent = NULL), "\n")
}
if (sys.nframe() == 0L) { a <- commandArgs(trailingOnly = TRUE); mode <- if (length(a)) a[1] else ""; if (mode %in% c("dry", "full")) main(mode == "dry") else stop("usage: run_rnaseq.R dry|full") }
