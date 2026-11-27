## prep_exact_calib.R - 27 Nov 2026, impact check O6, step 2 (after `ANOVA_EXACT=1 WP2_OUT=results_R/o6_anova_exact/wp2_exact gate_part1.R extract`).
## Compares the 12 x 17 matrices of the exact-rule extract with those of the R run with numpy's order (results_R/wp2/compact/dev_sets.json, the inputs of the R calibration and of dev_core_results.jsonl).
## Where a matrix is bit-identical (identical() on the numbers and the gene list), the result of the R run is reused, because every later step depends only on the matrix and the GSE number
## (seeds 2026 + GSE number, SeedSequence([2026, GSE number, condition, replicate])):
##   calib_series/<gse>.rds  copied when the set A matrix is identical (calibrate_t2.R run then computes only the other series);
##   dev_core_results.jsonl  the line of a set copied when that matrix is identical (gate_part1.R consistency then computes only the missing sets).
## Step "order" (after the consistency run) rewrites dev_core_results.jsonl in the frozen order of dev_prepared.json (A then B), as the R run has it.
## Output: results_R/o6_anova_exact/wp2_exact/matrix_comparison.csv; prints a summary. Never writes into results_R/wp2.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
OLD <- file.path(R_PORT, "results_R", "wp2", "compact"); NEW <- file.path(R_PORT, "results_R", "o6_anova_exact", "wp2_exact", "compact")
step <- commandArgs(trailingOnly = TRUE); step <- if (length(step)) step[1] else "prep"
mat12 <- function(rows) do.call(rbind, lapply(rows, function(r) as.numeric(unlist(r))))
key_of <- function(l) { j <- fromJSON(gsub("(?<=[:\\[,] )(NaN|-?Infinity)(?=[,}\\]])", "null", l, perl = TRUE), simplifyVector = FALSE); paste(j$gse, j$tag) }   # NaN literals as written by the core step
if (step == "prep") {
  Zo <- fromJSON(file.path(OLD, "dev_sets.json"), simplifyVector = FALSE); Po <- fromJSON(file.path(OLD, "dev_prepared.json"), simplifyVector = FALSE)
  Zn <- fromJSON(file.path(NEW, "dev_sets.json"), simplifyVector = FALSE); Pn <- fromJSON(file.path(NEW, "dev_prepared.json"), simplifyVector = FALSE)
  ## the R calibration (calib_series/*.rds, 25 Nov 00:40 to 02:48) and dev_core_results.jsonl (25 Nov 00:30) were computed from the archived Python matrices (oracle/dev_sets_archive.json), before the R
  ## extract existed; the R extract equals them to 1.2e-15 (R_PORT_PLAN.md, WP2 validation chain); the largest difference of every reused matrix from the archived one is recorded here
  Za <- fromJSON(file.path(R_PORT, "oracle", "dev_sets_archive.json"), simplifyVector = FALSE)
  core_old <- readLines(file.path(OLD, "dev_core_results.jsonl"), warn = FALSE, encoding = "UTF-8"); core_old <- core_old[nzchar(core_old)]; names(core_old) <- vapply(core_old, key_of, "")
  if (file.exists(file.path(NEW, "dev_core_results.jsonl"))) stop("dev_core_results.jsonl already exists in the exact-rule folder: remove it first (prep must start from an empty file)")
  dir.create(file.path(NEW, "calib_series"), showWarnings = FALSE)
  rows <- list(); keep <- character(0)
  for (g in union(names(Pn), names(Po))) for (tag in c("A", "B")) {
    k <- sprintf("%s__%s", g, tag); a <- Zo[[k]]; b <- Zn[[k]]
    same_mat <- !is.null(a) && !is.null(b) && identical(mat12(a), mat12(b))
    ga <- if (!is.null(Po[[g]])) unlist(Po[[g]][[paste0("genes_", tag)]]) else NULL; gb <- if (!is.null(Pn[[g]])) unlist(Pn[[g]][[paste0("genes_", tag)]]) else NULL
    same_genes <- !is.null(ga) && !is.null(gb) && identical(ga, gb)
    rows[[length(rows) + 1]] <- data.frame(gse = g, set = tag, platform = if (!is.null(Pn[[g]])) Pn[[g]]$platform else Po[[g]]$platform, in_numpy_run = !is.null(a), in_exact_run = !is.null(b),
      identical_matrix = same_mat, identical_genes = same_genes, n_genes_shared = if (!is.null(ga) && !is.null(gb)) length(intersect(ga, gb)) else NA_integer_,
      max_abs_diff_exact_to_archive = if (!is.null(b) && !is.null(Za[[k]])) max(abs(mat12(b) - mat12(Za[[k]]))) else NA_real_,
      n_de_numpy = if (!is.null(Po[[g]])) Po[[g]]$n_de else NA_integer_, n_de_exact = if (!is.null(Pn[[g]])) Pn[[g]]$n_de else NA_integer_, stringsAsFactors = FALSE)
    if (same_mat && same_genes) {
      kk <- paste(g, tag); if (!kk %in% names(core_old)) stop("no core line for ", kk); keep <- c(keep, core_old[[kk]])
      if (tag == "A") { src <- file.path(OLD, "calib_series", paste0(g, ".rds")); if (!file.exists(src)) stop("missing ", src); file.copy(src, file.path(NEW, "calib_series", paste0(g, ".rds")), overwrite = FALSE) }
    }
  }
  M <- do.call(rbind, rows); write.csv(M, file.path(dirname(NEW), "matrix_comparison.csv"), row.names = FALSE)
  con <- file(file.path(NEW, "dev_core_results.jsonl"), "wb"); writeLines(enc2utf8(keep), con, sep = "\n", useBytes = TRUE); close(con)
  cat(sprintf("series: numpy run %d, exact run %d; matrices identical (numbers and genes): set A %d, set B %d\n", length(Po), length(Pn), sum(M$identical_matrix & M$identical_genes & M$set == "A"), sum(M$identical_matrix & M$identical_genes & M$set == "B")))
  same <- M$identical_matrix & M$identical_genes
  cat(sprintf("largest absolute difference of a reused matrix from the archived Python matrix: %s (largest value in those matrices %s)\n", format(max(M$max_abs_diff_exact_to_archive[same], na.rm = TRUE), digits = 3),
    format(max(vapply(which(same), function(i) max(abs(mat12(Zn[[sprintf("%s__%s", M$gse[i], M$set[i])]]))), 1)), digits = 6)))
  d <- M[!same, ]; cat("not identical:\n"); print(d, row.names = FALSE)
  cat(sprintf("calibration files copied: %d; core lines copied: %d\n", length(list.files(file.path(NEW, "calib_series"), "[.]rds$")), length(keep)))
} else if (step == "order") {
  Pn <- fromJSON(file.path(NEW, "dev_prepared.json"), simplifyVector = FALSE); f <- file.path(NEW, "dev_core_results.jsonl")
  L <- readLines(f, warn = FALSE, encoding = "UTF-8"); L <- L[nzchar(L)]; names(L) <- vapply(L, key_of, "")
  want <- unlist(lapply(names(Pn), function(g) paste(g, c("A", "B")))); if (!setequal(want, names(L)) || anyDuplicated(names(L))) stop("core lines do not match the 2 x ", length(Pn), " sets")
  con <- file(f, "wb"); writeLines(enc2utf8(unname(L[want])), con, sep = "\n", useBytes = TRUE); close(con); cat(sprintf("dev_core_results.jsonl rewritten in frozen order: %d lines\n", length(want)))
} else stop("prep | order")
