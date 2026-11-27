## test_pyrandom.R - gate G0b: the R Mersenne Twister of pyrandom.R against Python's random.Random (reference/pyrandom_reference.csv made by make_pyrandom_reference.py). Usage: Rscript test_pyrandom.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "pyrandom.R"))
ref <- read.csv(file.path(HERE, "reference", "pyrandom_reference.csv")); ok <- TRUE
chk <- function(case, got) { want <- ref$value[ref$case == case]; good <- length(got) == length(want) && all(got == want); ok <<- ok && good; cat(sprintf("%s %-24s %d values\n", if (good) "OK  " else "FAIL", case, length(want))) }
g <- py_random(2026); chk("bits32_seed2026", vapply(1:2000, function(i) py_getrandbits(g, 32), 1))
g <- py_random(2^40 + 5); chk("bits32_seed2p40p5", vapply(1:700, function(i) py_getrandbits(g, 32), 1))
g <- py_random(2026); chk("bits7_seed2026", vapply(1:60, function(i) py_getrandbits(g, 7), 1))
g <- py_random(2026); chk("shuffle100_seed2026", py_shuffle(g, 0:99)); g <- py_random(2026); chk("shuffle5000_seed2026", py_shuffle(g, 0:4999)); g <- py_random(7); chk("shuffle37_seed7", py_shuffle(g, 0:36))
cat(if (ok) "PYRANDOM GATE PASSED\n" else "PYRANDOM GATE FAILED\n"); quit(status = if (ok) 0L else 1L)
