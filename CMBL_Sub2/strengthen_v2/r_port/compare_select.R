## compare_select.R - select_datasets.R re-screening of the archived series: status, grouping rule, groups and the 12 sample accessions against the archived eligibility.csv and frozen_datasets.csv.
## Usage: Rscript compare_select.R   (RB_OUT = folder of the R run, default results_R/benchmark_route3_select; RB_ARCHIVE = archive folder). Only the series the R run screened are compared.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
ARC <- Sys.getenv("RB_ARCHIVE", normalizePath(file.path(HERE, "..", "..", "benchmark_route3"))); OUT <- Sys.getenv("RB_OUT", file.path(HERE, "results_R", "benchmark_route3_select"))
rd <- function(f) read.csv(f, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0))
ea <- rd(file.path(ARC, "eligibility.csv")); er <- rd(file.path(OUT, "eligibility.csv")); fa <- rd(file.path(ARC, "frozen_datasets.csv")); fr <- rd(file.path(OUT, "frozen_datasets.csv")); ok <- TRUE
key <- function(d) paste(d$platform, d$gse); m <- match(key(er), key(ea)); bad <- 0L
for (i in seq_len(nrow(er))) for (col in c("status", "how", "n_groups", "groups")) if (!identical(er[[col]][i], ea[[col]][m[i]])) { bad <- bad + 1L; cat(sprintf("DIFF eligibility %s %s %s:\n   Python: %s\n   R     : %s\n", er$platform[i], er$gse[i], col, ea[[col]][m[i]], er[[col]][i])) }
cat(sprintf("%-6s eligibility.csv: %d series screened, %d eligible in R / %d eligible in the archive among them, %d differing fields\n", if (bad == 0) "AGREE" else "DIFFER", nrow(er), sum(er$status == "eligible"), sum(ea$status[m] == "eligible"), bad))
fm <- match(key(fr), key(fa)); badf <- 0L
for (i in seq_len(nrow(fr))) for (col in c("how", "n_groups", "groups", "gsms")) if (is.na(fm[i]) || !identical(fr[[col]][i], fa[[col]][fm[i]])) { badf <- badf + 1L; cat(sprintf("DIFF frozen %s %s %s\n", fr$platform[i], fr$gse[i], col)) }
nfa <- sum(key(fa) %in% key(er)); cat(sprintf("%-6s frozen_datasets.csv: %d rows in R, %d archived rows among the screened series, %d differing fields\n", if (badf == 0 && nrow(fr) == nfa) "AGREE" else "DIFFER", nrow(fr), nfa, badf))
