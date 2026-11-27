## freeze_lists.R - R port of wp4_islets/freeze_lists.py and read_geo_meta.py: build the frozen dataset table and the frozen beta-cell lists of protocol 4.3 from METADATA ONLY (GEO series matrices, ArrayExpress SDRF), before any
## expression value of these datasets is read. Output: frozen_islet_datasets.csv (order, accession, decision, label source, beta cells, donors, reason) and frozen_beta_cells.csv (one row per author-labelled beta cell: dataset,
## cell_key = join key to the expression file, gsm, donor, disease, label, label_source). Label rules (fixed here): the label is the cell-type annotation written by the data authors in the public record; no label is derived from
## expression by us; cells without that label are not beta-cells. Input: <WP4_DATA>/GSE86469/GSE86469_series_matrix.txt.gz, GSE81547/GSE81547_series_matrix.txt.gz, E-MTAB-5061/E-MTAB-5061.sdrf.txt (wp4_http.R fetches them).
## Output folder: WP4_FREEZE_OUT (default results_R/wp4/frozen). The frozen files that the analyses read are the archived ones shipped with the repository (wp4_core.R: FROZEN()); this script rebuilds them for the check.
## Usage: Rscript freeze_lists.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R"))
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); OUT <- Sys.getenv("WP4_FREEZE_OUT", file.path(R_PORT, "results_R", "wp4", "frozen")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sample_table <- function(path) {      # the sample-level header block of a GEO series matrix (!Sample_ lines) as a data frame: columns title, geo_accession, ..., and 'ch:<key>' for a characteristics line whose values all share one 'key: value' key
  con <- gzfile(path, "rt"); L <- readLines(con, warn = FALSE, encoding = "UTF-8"); close(con); b <- which(startsWith(L, "!series_matrix_table_begin"))[1]; if (!is.na(b)) L <- L[seq_len(b - 1L)]; L <- L[startsWith(L, "!Sample_")]
  keys <- vapply(strsplit(L, "\t", fixed = TRUE), function(p) trimws(p[1]), ""); vals <- lapply(strsplit(L, "\t", fixed = TRUE), function(p) gsub("^\"+|\"+$", "", trimws(p[-1]))); cols <- list()
  for (k in unique(keys)) { ix <- which(keys == k); for (i in seq_along(ix)) {
    v <- vals[[ix[i]]]; name <- paste0(substring(k, nchar("!Sample_") + 1L), if (length(ix) == 1L) "" else sprintf("_%d", i))
    if (k == "!Sample_characteristics_ch1") { kk <- unique(sub(": .*$", "", v[grepl(": ", v, fixed = TRUE)])); if (length(kk) == 1L) { name <- paste0("ch:", kk); v <- ifelse(grepl(": ", v, fixed = TRUE), sub("^[^:]*: ", "", v), v) } }
    cols[[name]] <- v } }
  as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
}
rows <- list(); summ <- list()
## 1. GSE86469 (Lawlor 2017): 'cell type' = Beta in the GEO sample characteristics; join key = sample title (column name of the expression csv)
T <- sample_table(file.path(DATA, "GSE86469", "GSE86469_series_matrix.txt.gz")); B <- T[T[["ch:cell type"]] == "Beta", ]
rows[[1]] <- data.frame(dataset = "GSE86469", cell_key = B$title, gsm = B$geo_accession, donor = B[["ch:islet unos id"]], disease = ifelse(startsWith(B[["ch:disease"]], "Type 2"), "T2D", "non-diabetic"), label = "Beta", label_source = "GEO sample characteristic \"cell type\" (authors)", stringsAsFactors = FALSE)
summ[[1]] <- data.frame(order = 1L, dataset = "GSE86469", source = "GEO ftp (series matrix + RSEM expected counts csv.gz)", decision = "attempt", label_source = "GEO sample characteristic \"cell type\" (authors)", beta_cells = nrow(B), donors = length(unique(B[["ch:islet unos id"]])), reason = "", stringsAsFactors = FALSE)
## 2. GSE81547 (Enge 2017): 'inferred_cell_type' = beta in the GEO sample characteristics; join key = GSM accession (file name prefix inside GSE81547_RAW.tar)
T <- sample_table(file.path(DATA, "GSE81547", "GSE81547_series_matrix.txt.gz")); B <- T[T[["ch:inferred_cell_type"]] == "beta", ]; don <- sprintf("age%s_%s", B[["ch:donor_age"]], B[["ch:gender"]])
rows[[2]] <- data.frame(dataset = "GSE81547", cell_key = B$geo_accession, gsm = B$geo_accession, donor = don, disease = "not stated (series summary: previously healthy donors)", label = "beta", label_source = "GEO sample characteristic \"inferred_cell_type\" (authors)", stringsAsFactors = FALSE)
summ[[2]] <- data.frame(order = 2L, dataset = "GSE81547", source = "GEO ftp (GSE81547_RAW.tar, 176 MB, per-cell csv.gz)", decision = "attempt", label_source = "GEO sample characteristic \"inferred_cell_type\" (authors)", beta_cells = nrow(B), donors = length(unique(don)), reason = "", stringsAsFactors = FALSE)
## 3. E-MTAB-5061 (Segerstolpe 2016): 'inferred cell type' = beta cell in the SDRF; join key = Source Name; only cells with single cell quality OK carry an inferred type
S <- read.delim(file.path(DATA, "E-MTAB-5061", "E-MTAB-5061.sdrf.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", colClasses = "character", na.strings = character(0)); B <- S[S[["Characteristics [inferred cell type]"]] == "beta cell", ]
rows[[3]] <- data.frame(dataset = "E-MTAB-5061", cell_key = B[["Source Name"]], gsm = "", donor = B[["Characteristics [individual]"]], disease = ifelse(startsWith(B[["Characteristics [disease]"]], "type II"), "T2D", "non-diabetic"), label = "beta cell", label_source = "ArrayExpress SDRF \"inferred cell type\" (authors)", stringsAsFactors = FALSE)
summ[[3]] <- data.frame(order = 3L, dataset = "E-MTAB-5061", source = "ArrayExpress / BioStudies (processed.1.zip 159 MB holds the 723 MB rpkm+counts text; the zip is streamed, the 723 MB text is not downloaded)", decision = "attempt", label_source = "ArrayExpress SDRF \"inferred cell type\" (authors)", beta_cells = nrow(B), donors = length(unique(B[["Characteristics [individual]"]])), reason = "", stringsAsFactors = FALSE)
## 4. GSE85241 (Muraro 2016): the GEO record has 32 library-level samples and no per-cell annotation; clusters were found by the authors' StemID and are not in the GEO files
summ[[4]] <- data.frame(order = 4L, dataset = "GSE85241", source = "GEO ftp (counts csv.gz 17 MB, barcode list, demultiplexing readme)", decision = "skipped before expression download", label_source = "none in the public record (32 library-level GEO samples; supplementary files: counts, barcodes, readme)", beta_cells = 0L, donors = 0L,
                        reason = "no author-provided beta-cell labels in the GEO or repository files; protocol 4.3 requires author-provided labels and plain public download, and we do not assign cell types ourselves", stringsAsFactors = FALSE)
py_write_csv_repr(do.call(rbind, rows), file.path(OUT, "frozen_beta_cells.csv"), eol = "\r\n"); py_write_csv_repr(do.call(rbind, summ), file.path(OUT, "frozen_islet_datasets.csv"), eol = "\r\n"); cat("frozen lists written:", sum(vapply(rows, nrow, 1L)), "beta cells\n")
