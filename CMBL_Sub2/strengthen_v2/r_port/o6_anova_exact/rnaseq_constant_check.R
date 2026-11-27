## rnaseq_constant_check.R - 26 Nov 2026, impact check O6: ANOVA_EXACT changes the p value of a gene only when the gene is constant in all 12 samples, so a series without such genes has the same
## DE count under both rules. This script counts, for every frozen held-out RNA-seq series whose NCBI count file is cached, the genes of the analysis universe (count >= 1 in all 12 samples,
## log2(CPM + 1)) that are constant in all 12 samples. Read-only on the cached count files (untrusted data, read as numbers). Output: results_R/o6_anova_exact/rnaseq_constant_genes.csv
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "bench_pipeline.R")); source(file.path(R_PORT, "wp2_benchmark2", "bench_data_io.R"))
rows <- read.csv(file.path(R_PORT, "..", "wp2_benchmark2", "frozen_heldout_rnaseq.csv"), stringsAsFactors = FALSE, check.names = FALSE)
root <- file.path(R_PORT, "results_R", "wp2_data"); out <- list()
for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]; path <- counts_path(root, r$gse, RNASEQ_BUILD[[r$organism]]); if (!file.exists(path)) next
  gsms <- strsplit(r$gsms, ";", fixed = TRUE)[[1]]
  res <- tryCatch({ cf <- read_counts_file(path, gsms); rc <- rnaseq_cpm(cf$C); L <- rc$L2
    n_const <- sum(colSums(L != matrix(L[1, ], nrow(L), ncol(L), byrow = TRUE)) == 0); list(n_genes = ncol(L), n_const = n_const, note = "") },
    error = function(e) list(n_genes = NA, n_const = NA, note = substr(conditionMessage(e), 1, 80)))
  out[[length(out) + 1]] <- data.frame(organism = r$organism, gse = r$gse, n_genes = res$n_genes, n_constant = res$n_const, note = res$note)
  cat(sprintf("%s %s: %s genes, %s constant %s\n", r$organism, r$gse, res$n_genes, res$n_const, res$note))
}
d <- do.call(rbind, out); dir.create(file.path(R_PORT, "results_R", "o6_anova_exact"), showWarnings = FALSE, recursive = TRUE)
write.csv(d, file.path(R_PORT, "results_R", "o6_anova_exact", "rnaseq_constant_genes.csv"), row.names = FALSE)
cat(sprintf("series read: %d; with at least one constant gene: %d; largest count: %s\n", sum(!is.na(d$n_const)), sum(d$n_constant > 0, na.rm = TRUE), max(d$n_constant, na.rm = TRUE)))
