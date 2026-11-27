## wp4_stats.R - small statistics and reader helpers shared by the WP4 perturbation and wider-search scripts: scipy-equivalent Spearman correlation and exact binomial test, Python list repr for CSV cells, and a reader for the text
## table of a GEO series matrix. Source wp4_core.R (pycompat.R, data.table) first.
suppressMessages(library(data.table))
py_spearmanr <- function(x, y) {      # scipy.stats.spearmanr: rho = Pearson correlation of the average ranks; two-sided p from the t distribution with n - 2 degrees of freedom (no exact p)
  n <- length(x); rho <- cor(rank(x), rank(y)); t <- rho * sqrt(max((n - 2) / ((rho + 1) * (1 - rho)), 0)); list(rho = rho, p = 2 * pt(abs(t), n - 2, lower.tail = FALSE))
}
py_binom_p <- function(k, n, p = 0.5) binom.test(k, n, p)$p.value        # scipy.stats.binomtest(k, n, p).pvalue (two-sided, the same method of small p-values as R: relative error 1 + 1e-7)
py_list_repr <- function(v) if (!length(v)) "[]" else paste0("[", paste0("'", v, "'", collapse = ", "), "]")        # repr of a Python list of strings, as pandas writes it into a CSV cell
read_series_matrix <- function(path) {      # sample titles and accessions (header lines) and the table between the markers as a data frame (first column ID_REF, then one numeric column per sample)
  con <- gzfile(path, "rt"); L <- readLines(con, warn = FALSE, encoding = "UTF-8"); close(con); fields <- function(l) gsub("^\"+|\"+$", "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  b <- which(startsWith(L, "!series_matrix_table_begin"))[1]; e <- which(startsWith(L, "!series_matrix_table_end"))[1]; hdr <- L[seq_len(b - 1L)]
  tab <- fread(text = paste(L[(b + 1L):(e - 1L)], collapse = "\n"), sep = "\t", header = TRUE, quote = "\"", data.table = FALSE, colClasses = list(character = 1))
  names(tab) <- gsub("^\"+|\"+$", "", names(tab)); tab[[1]] <- gsub("^\"+|\"+$", "", tab[[1]])
  list(title = fields(hdr[startsWith(hdr, "!Sample_title")][1]), gsm = fields(hdr[startsWith(hdr, "!Sample_geo_accession")][1]), table = tab)
}
