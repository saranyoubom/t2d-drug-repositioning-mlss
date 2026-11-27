## test_pyformat.R - gate for the Python-style number formatting and parsing of pycompat.R (used by every writer of the port): py_float (= float()), .json_num (= repr()), py_fmt_g, py_fmt_f against reference/pyformat_reference.tsv
## (made once by make_pyformat_reference.py). Usage: Rscript rng/test_pyformat.R   (exit status 1 if any value differs)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "..", "pycompat.R")); suppressMessages(library(jsonlite))
R <- read.delim(file.path(HERE, "reference", "pyformat_reference.tsv"), colClasses = "character", na.strings = character(0), quote = "")
x <- py_float(R$text); got <- list(repr = vapply(x, .json_num, ""), g16 = py_fmt_g(x, 16L), g6 = py_fmt_g(x, 6L), f2 = py_fmt_f(x, 2L), f6 = py_fmt_f(x, 6L))
nbad <- 0L; for (k in names(got)) { bad <- which(got[[k]] != R[[k]]); nbad <- nbad + length(bad); if (length(bad)) cat(sprintf("%s differs for %d values, first: %s Python '%s' R '%s'\n", k, length(bad), R$text[bad[1]], R[[k]][bad[1]], got[[k]][bad[1]])) }
cat(sprintf("test_pyformat: %d values x %d formats, %d differences%s\n", nrow(R), length(got), nbad, if (nbad) "" else " - PASSED")); if (sys.nframe() == 0L) quit(status = if (nbad) 1L else 0L)
