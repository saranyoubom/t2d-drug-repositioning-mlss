## insilico.R - R port of reanalysis_17genes/insilico.py: in silico check of every primer pair of Table S1 (Additional file 1.docx) against its NCBI reference sequence (primer_lib.R, cached):
## match, mismatches (at most 2), reversed orientation, amplicon length against the table. Writes insilico_table.json into RA_OUT and prints one line per gene.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "docx_table.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE); source(file.path(HERE, "primer_lib.R"))
AF1 <- Sys.getenv("AF1_DOCX", normalizePath(file.path(R_PORT, "..", "..", "..", "Supplementary", "Additional file 1.docx")))
rows <- list()
for (c in p_table_rows(AF1)) { seqs <- ws_split(c[4]); if (length(seqs) < 2) next; rows[[length(rows) + 1]] <- list(strsplit(c[1], " ", fixed = TRUE)[[1]][1], c[2], seqs[1], seqs[2], as.integer(.allm("\\d+", c[5])[1])) }
out <- lapply(rows, function(r) p_check(r[[1]], r[[2]], r[[3]], r[[4]], r[[5]]))
write_text_bytes(py_json_typed(out, indent = NULL), file.path(OUT, "insilico_table.json"))
nn <- function(x) if (is.null(x)) "None" else as.character(x)
for (o in out) cat(sprintf("%-8s %-16s %-30s amp=%s table=%s | %s\n", o$gene, o$acc, o$result, nn(o$amplicon), nn(o$tablelen), substr(o$hdr, 1, 85)))
