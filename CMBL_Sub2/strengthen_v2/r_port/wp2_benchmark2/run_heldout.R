## run_heldout.R - R port of wp2_benchmark2/run_heldout.py: dry run and full run on the frozen held-out list (protocol v2 section 3.5, WP2 deliverable D), same code path as the development series
## (bench_data_io.R prepare_array -> bench2_core.R analyse_matrix; the loop, the duplicate rule and the stopping rule are those of run_benchmark).
## Usage: Rscript run_heldout.R arrays dry | arrays full
##   dry   decides analysability only (file integrity, number of DE genes >= 34), writes heldout_arrays_dry.jsonl; no correlation is computed
##   full  processes the frozen list in order until 20 series per platform are analysable or the list is exhausted; one JSON line per series (status or reason) in heldout_arrays_results.jsonl with the
##         flat endpoints and details of gene sets A and B; every decision is checked against the dry run (a mismatch is logged, never silently accepted)
## Duplicate rule (as in part 1, deviation D1): a series sharing at least one of its 12 selected samples with an earlier analysable held-out series is excluded. Downloads are throttled and kept under
## <data>/v2_data/arrays_heldout/matrices (untrusted data, read as numbers only). WP2_FROZEN_ARRAYS = the frozen list (default: the archived file); WP2_ONLY = comma-separated series for checks on a few series.
if (!exists("HERE") || !file.exists(file.path(HERE, "wp2_common.R"))) { HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, "..")) }
source(file.path(HERE, "wp2_common.R")); source(file.path(R_PORT, "coexp_checks.R")); source(file.path(R_PORT, "cc_fast.R")); source(file.path(R_PORT, "bench_pipeline.R")); source(file.path(HERE, "bench_data_io.R")); source(file.path(HERE, "bench2_core.R"))
MAXP <- 20L
arrays <- function(dry) {
  io_configure(file.path(DATA, "v2_data", "arrays_heldout", "matrices"), file.path(DATA, "bench_data", "platforms"))
  rows <- read.csv(Sys.getenv("WP2_FROZEN_ARRAYS", file.path(ARC, "frozen_heldout_arrays.csv")), stringsAsFactors = FALSE, check.names = FALSE); only <- setdiff(strsplit(Sys.getenv("WP2_ONLY", ""), ",", fixed = TRUE)[[1]], ""); if (length(only)) rows <- rows[rows$gse %in% only, , drop = FALSE]
  plats <- py_sort(unique(rows$platform)); dry_path <- file.path(OUT, "heldout_arrays_dry.jsonl"); ref <- list()
  if (!dry && file.exists(dry_path)) { for (j in read_jsonl(dry_path)) ref[[j$gse]] <- j }
  fn <- if (dry) dry_path else file.path(OUT, "heldout_arrays_results.jsonl"); if (file.exists(fn)) file.remove(fn)
  done <- setNames(rep(0L, length(plats)), plats); used <- list(); mism <- 0L
  for (i in seq_len(nrow(rows))) {
    row <- as.list(rows[i, ]); gpl <- row$platform; if (done[[gpl]] >= MAXP) next
    t0 <- Sys.time(); gs <- strsplit(row$gsms, ";", fixed = TRUE)[[1]]; clash <- py_sort(unique(unlist(used[gs[gs %in% names(used)]]))); sets <- NULL
    if (length(clash)) rec <- list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = paste0("not analysable: duplicate of ", clash[1]))
    else {
      r <- tryCatch(prepare_array(row, dry), error = function(e) list(rec = list(gse = row$gse, platform = gpl, order = as.integer(row$order), status = sprintf("not analysable: %s %s", class(e)[1], substr(conditionMessage(e), 1, 120))), sets = NULL))
      rec <- r$rec; sets <- r$sets
    }
    if (!is.null(sets) && !dry) {
      ok <- tryCatch({
        for (tag in c("A", "B")) { res <- analyse_matrix(sets[[tag]]$X, sets[[tag]]$genes, as.integer(substring(row$gse, 4))); rec[[tag]] <- list(flat = res$flat, detail = res$detail, genes = as.list(sets[[tag]]$genes), n_de_in_set = sets[[tag]]$n_de_in_set) }
        rec$status <- "analysed"; TRUE
      }, error = function(e) { rec$status <<- sprintf("not analysable: analysis error %s %s", class(e)[1], substr(conditionMessage(e), 1, 120)); rec$A <<- NULL; rec$B <<- NULL; FALSE })
    }
    if (rec$status %in% c("analysed", "analysable (dry run)")) { done[[gpl]] <- done[[gpl]] + 1L; for (g in gs) used[[g]] <- row$gse }
    if (length(ref)) {
      r0 <- ref[[row$gse]]; a0 <- if (is.null(r0)) NA else identical(r0$status, "analysable (dry run)"); a1 <- identical(rec$status, "analysed")
      if (is.null(r0) || !identical(a0, a1) || !identical(r0$n_de, rec$n_de)) { mism <- mism + 1L; rec$dry_run_mismatch <- TRUE }
    }
    cat(py_json_typed(rec, indent = NULL), "\n", file = fn, sep = "", append = TRUE)
    cat(sprintf("%s %s (#%s): %s | DE %s | done {%s} | %.0fs\n", gpl, row$gse, row$order, rec$status, if (is.null(rec$n_de)) "-" else rec$n_de, paste(sprintf("'%s': %d", names(done), done), collapse = ", "), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  cat("finished {", paste(sprintf("'%s': %d", names(done), done), collapse = ", "), "} dry-run mismatches:", mism, "\n")
}
if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE); arm <- if (length(a) >= 1) a[1] else ""; mode <- if (length(a) >= 2) a[2] else ""
  if (arm == "arrays" && mode %in% c("dry", "full")) arrays(mode == "dry") else stop("usage: run_heldout.R arrays dry|full")
}
