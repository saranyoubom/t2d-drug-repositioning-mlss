## compare_wp5.R - WP5a: the R results (effect_sizes.R) against the archived Python results of wp5_reporting/.
## Usage: ES_OUT=<folder of the R run> ES_LABEL=effect_sizes.py Rscript effect_sizes.R ; ES_OUT=<same> Rscript compare_wp5.R
ES_HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(ES_HERE, "agree.R"))
ARC <- Sys.getenv("WP5_ARCHIVE", normalizePath(file.path(ES_HERE, "..", "wp5_reporting")))
OUT <- Sys.getenv("ES_OUT", file.path(ES_HERE, "results_R", "wp5"))
ok <- c(
  agree_csv(file.path(ARC, "effect_sizes.csv"), file.path(OUT, "effect_sizes.csv"))$ok,
  agree_json(file.path(ARC, "gate_reproduction.json"), file.path(OUT, "gate_reproduction.json"), tol = 2e-9)$ok)
## Markdown: the archived file has CRLF line ends and UTF-8 text; compare the text line by line and the raw bytes
a <- readLines(file.path(ARC, "effect_sizes.md"), encoding = "UTF-8", warn = FALSE); b <- readLines(file.path(OUT, "effect_sizes.md"), encoding = "UTF-8", warn = FALSE)
same_text <- identical(a, b); n_diff <- if (same_text) 0L else sum(a[seq_len(min(length(a), length(b)))] != b[seq_len(min(length(a), length(b)))]) + abs(length(a) - length(b))
cat(sprintf("%-6s %-42s %d lines (archive) / %d lines (R), %d lines differ; raw bytes identical: %s\n", if (same_text) "AGREE" else "DIFFER", "effect_sizes.md", length(a), length(b), n_diff,
            identical(readBin(file.path(ARC, "effect_sizes.md"), "raw", 1e7), readBin(file.path(OUT, "effect_sizes.md"), "raw", 1e7))))
if (!same_text) { k <- which(a[seq_len(min(length(a), length(b)))] != b[seq_len(min(length(a), length(b)))]); for (i in head(k, 6)) cat(sprintf("  line %d\n    Python: %s\n    R     : %s\n", i, substr(a[i], 1, 400), substr(b[i], 1, 400))) }
cat(sprintf("\nWP5a overall: %s\n", if (all(ok) && same_text) "ALL AGREE" else "DIFFERENCES FOUND"))
