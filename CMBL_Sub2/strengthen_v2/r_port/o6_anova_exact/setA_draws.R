## setA_draws.R - 26 Nov 2026, impact check O6: for the analysed series whose DE count changes under ANOVA_EXACT, draw gene sets A and B under both rules (prepare_array, the code path of the
## development and held-out series) and report how many of the 17 genes are shared. No correlation is computed. Read-only on the cached GEO files.
## Output: results_R/o6_anova_exact/setA_draws.txt
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "bench_pipeline.R")); source(file.path(R_PORT, "wp2_benchmark2", "bench_data_io.R"))
D <- file.path(R_PORT, "results_R", "wp2_data"); ARC <- normalizePath(file.path(R_PORT, ".."))
route3 <- read.csv(file.path(ARC, "..", "benchmark_route3", "frozen_datasets.csv"), stringsAsFactors = FALSE, check.names = FALSE)
held <- read.csv(file.path(ARC, "wp2_benchmark2", "frozen_heldout_arrays.csv"), stringsAsFactors = FALSE, check.names = FALSE)
todo <- list(list(route3, "bench_data/matrices", c("GSE15892", "GSE28477", "GSE53183")), list(held, "v2_data/arrays_heldout/matrices", "GSE68424"))
out <- character(0); say <- function(s) { out <<- c(out, s); cat(s, "\n") }
for (td in todo) {
  io_configure(file.path(D, td[[2]]), file.path(D, "bench_data", "platforms"))
  for (g in td[[3]]) {
    row <- as.list(td[[1]][td[[1]]$gse == g, ][1, ])
    ANOVA_EXACT <- FALSE; a <- prepare_array(row); ANOVA_EXACT <- TRUE; b <- prepare_array(row); ANOVA_EXACT <- FALSE
    for (tag in c("A", "B")) say(sprintf("%s %s set %s: DE %d -> %d; genes shared %d of 17; n_de_in_set %d -> %d; numpy [%s] | exact [%s]", row$platform, g, tag, a$rec$n_de, b$rec$n_de,
      length(intersect(a$sets[[tag]]$genes, b$sets[[tag]]$genes)), a$sets[[tag]]$n_de_in_set, b$sets[[tag]]$n_de_in_set, paste(a$sets[[tag]]$genes, collapse = ","), paste(b$sets[[tag]]$genes, collapse = ",")))
  }
}
writeLines(out, file.path(R_PORT, "results_R", "o6_anova_exact", "setA_draws.txt"), useBytes = TRUE)
