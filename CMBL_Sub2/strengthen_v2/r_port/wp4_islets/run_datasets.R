## run_datasets.R - R port of wp4_islets/run_datasets.py: protocol 4.2 / 4.3: cell-level analysis of every public islet dataset, unadjusted and depth-adjusted, with the expression-matched null (wp4_core.R analyse).
## Datasets (all frozen before analysis): gse81608, gse86469, gse81547, emtab5061 (parsed by wp4_parse.R, human) and GSE84133 human and mouse. Variants: 'primary' (all frozen beta-cells) and, for every dataset with donors that have type 2
## diabetes, 'nondiabetic' (those donors removed, same rules). Output in <WP4_OUT> (default results_R/wp4/results): <tag>__<variant>__<mode>__{pairs,cons,detection}.csv, <tag>__<variant>__summary.json (both modes), index_<tag>.csv
## and index.csv (one row per dataset, variant and mode). One dataset at a time is held in memory.
## Usage: Rscript run_datasets.R [tag ...]   tags: gse81608 gse84133_human gse86469 gse81547 emtab5061 gse84133_mouse (default: all); WP4_DATA = folder of the raw data (default results_R/wp4_data)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp4_core.R")); source(file.path(HERE, "wp4_parse.R")); options(py_json_nan = TRUE)
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); OUT <- Sys.getenv("WP4_OUT", file.path(R_PORT, "results_R", "wp4", "results")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOGS <- file.path(dirname(OUT), "logs"); dir.create(LOGS, showWarnings = FALSE, recursive = TRUE)
TAGS <- c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061", "gse84133_mouse")
KEEP <- c("tag", "variant", "mode", "kind", "depth_covariate", "beta_cells", "donors", "t2d_beta_cells", "df_for_p", "genes_evaluable", "consensus_pairs_evaluable", "other_pairs_evaluable", "median_rho_consensus", "median_rho_other", "null_median_rho_mean",
          "null_median_rho_95pct", "empirical_p_median_ge_observed", "empirical_p_mean_ge_observed", "mean_z_observed", "mean_z_null", "E", "SE_E", "empirical_p_E", "consensus_positive", "consensus_sig_positive_q05", "consensus_sig_negative_q05", "edges_q05_rho02")
is_t2d <- function(d) { s <- tolower(trimws(d)); startsWith(s, "t2d") | startsWith(s, "type 2") | startsWith(s, "type ii") }
run <- function(tag) {
  t0 <- Sys.time(); cells <- load_cells(tag, DATA, LOGS); t2d <- is_t2d(cells$disease); rows <- list(); variants <- list(list("primary", cells)); if (any(t2d) && sum(!t2d) > 0) variants[[2]] <- list("nondiabetic", cells_subset(cells, !t2d))
  for (v in variants) {
    vname <- v[[1]]; cv <- v[[2]]; res <- analyse(cv, 0.10, c("unadjusted", "adjusted")); n_t2d <- as.integer(sum(is_t2d(cv$disease)))
    for (mode in names(res)) {
      r <- res[[mode]]; for (part in c("pairs", "cons", "detection")) py_write_csv_repr(r[[part]], file.path(OUT, sprintf("%s__%s__%s__%s.csv", tag, vname, mode, part)), eol = "\r\n")
      s <- r$summary; s$tag <- tag; s$variant <- vname; s$t2d_beta_cells <- n_t2d; s$kind <- cv$kind; s$depth_covariate <- if (cv$kind %in% c("umi", "count")) "log10 library size" else "log10 detected genes"; rows[[length(rows) + 1]] <- s
      cat(sprintf("%-16s %-11s %-10s cells %5d donors %2d | evaluable genes %2d pairs %2d | median rho %s | null %s | E %s SE %s\n", tag, vname, mode, s$beta_cells, s$donors, s$genes_evaluable, s$consensus_pairs_evaluable,
                  if (!is.null(s$median_rho_consensus)) sprintf("%.3f", s$median_rho_consensus) else "-", if (!is.null(s$null_median_rho_mean)) sprintf("%.3f", s$null_median_rho_mean) else "-", if (!is.null(s$E)) sprintf("%.3f", s$E) else "-", if (!is.null(s$SE_E)) sprintf("%.3f", s$SE_E) else "-"))
    }
    write_json_typed(lapply(res, function(r) r$summary), file.path(OUT, sprintf("%s__%s__summary.json", tag, vname)))
  }
  I <- do.call(rbind, lapply(rows, function(s) { x <- lapply(KEEP, function(k) if (is.null(s[[k]])) NA else s[[k]]); names(x) <- KEEP; as.data.frame(x, stringsAsFactors = FALSE) }))
  py_write_csv_repr(I, file.path(OUT, sprintf("index_%s.csv", tag)), eol = "\r\n"); cat(sprintf("%s done in %.0f s\n", tag, as.numeric(difftime(Sys.time(), t0, units = "secs")))); invisible(I)
}
a <- commandArgs(trailingOnly = TRUE); tags <- if (length(a)) a else TAGS
if (sys.nframe() == 0L) {
  for (tag in tags) run(tag)
  parts <- lapply(TAGS[file.exists(file.path(OUT, sprintf("index_%s.csv", TAGS)))], function(t) read.csv(file.path(OUT, sprintf("index_%s.csv", t)), stringsAsFactors = FALSE, check.names = FALSE))
  py_write_csv_repr(do.call(rbind, parts), file.path(OUT, "index.csv"), eol = "\r\n"); cat("index.csv written with", length(parts), "datasets\n")
}
