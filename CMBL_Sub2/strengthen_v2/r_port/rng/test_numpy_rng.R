## test_numpy_rng.R - gate G0: numpy_rng.R must reproduce numpy bit for bit.
## Reads reference/manifest.json and reference/*.f64 (made by make_numpy_reference.py) and replays the same operations in R.
## Run: Rscript test_numpy_rng.R            (set RNG_TEST_SEEDS=int2026,int7 to test only some seeds)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "numpy_rng.R"))
suppressMessages(library(jsonlite))
man <- fromJSON(file.path(here, "reference", "manifest.json"), simplifyVector = FALSE)
cat("numpy", man$numpy, "reference; R", R.version.string, "\n")
only <- strsplit(Sys.getenv("RNG_TEST_SEEDS", ""), ",")[[1]]

run_op <- function(g, op) {
  switch(op$op,
    random = np_random(g, op$n),
    random_matrix = np_random(g, op$rows * op$cols),
    integers = np_integers(g, op$low, op$high, op$n),
    integers_scalar = np_integers(g, 0, op$high, 1),
    standard_normal = np_standard_normal(g, op$n),
    normal = np_normal(g, op$loc, op$scale, op$n),
    normal_array = { loc <- unlist(op$loc); sc <- unlist(op$scale); z <- np_standard_normal(g, op$rows * length(loc)); rep(loc, op$rows) + rep(sc, op$rows) * z },
    permutation = np_permutation(g, op$n),
    permutation_array = np_permute(g, unlist(op$x)),
    choice_noreplace = np_choice_noreplace_idx(g, op$pop, op$k) * 10 + 1,
    choice_scalar = np_integers(g, 0, op$pop, 1) * 3,
    choice_replace = np_integers(g, 0, op$pop, op$n) * 3,
    stop("unknown op ", op$op))
}
describe <- function(op) paste0(op$op, "(", paste(unlist(op[setdiff(names(op), "op")])[1:min(3, length(op) - 1)], collapse = ","), ")")
compare <- function(a, b) {
  if (length(a) != length(b)) return(sprintf("length %d, expected %d", length(a), length(b)))
  bad <- which(!(a == b))
  if (length(bad) == 0) return(NULL)
  i <- bad[1]
  sprintf("%d of %d values differ; first at %d: got %.17g, expected %.17g", length(bad), length(a), i, a[i], b[i])
}

total_fail <- 0
for (s in man$seeds) {
  if (length(only) && !(s$name %in% only)) next
  f <- file.path(here, "reference", s$file)
  ref <- readBin(f, "double", n = file.size(f) / 8, size = 8)
  entropy <- unlist(s$value)
  nh <- s$n_head; nr <- s$n_raw
  cat(sprintf("\n== %s (entropy %s)\n", s$name, paste(entropy, collapse = ",")))
  t0 <- Sys.time()
  g <- np_rng(entropy)
  pool <- np_seedseq_pool(entropy)
  chk <- list(pool = compare(pool, ref[1:4]),
              words = compare(np_seedseq_words(pool), ref[5:12]),
              state = compare(g$state, ref[13:20]),
              inc = compare(g$inc, ref[21:28]))
  for (nm in names(chk)) { if (is.null(chk[[nm]])) cat("   ok   ", nm, "\n") else { cat("   FAIL ", nm, ":", chk[[nm]], "\n"); total_fail <- total_fail + 1 } }
  if (!is.null(chk$state) || !is.null(chk$inc)) { cat("   (seeding differs: later checks skipped)\n"); next }
  gr <- np_rng(entropy); w <- np_raw64(gr, nr)
  e1 <- compare(w$lo, ref[nh + 1:nr]); e2 <- compare(w$hi, ref[nh + nr + 1:nr])
  if (is.null(e1) && is.null(e2)) cat("   ok    raw 64-bit outputs:", nr, "\n") else { cat("   FAIL  raw outputs:", e1, e2, "\n"); total_fail <- total_fail + 1; next }
  pos <- nh + 2 * nr; npass <- 0; failed <- FALSE
  for (i in seq_along(s$ops)) {
    op <- s$ops[[i]]; len <- s$op_lengths[[i]]; expect <- ref[pos + seq_len(len)]; pos <- pos + len
    if (failed) next
    tm <- system.time(out <- tryCatch(run_op(g, op), error = function(e) conditionMessage(e)))[["elapsed"]]
    err <- if (is.character(out)) paste("error:", out) else compare(as.numeric(out), expect)
    if (is.null(err)) { npass <- npass + 1; if (tm > 1) cat(sprintf("   ok    op %2d %-28s %6.1f s\n", i, describe(op), tm)) }
    else { cat(sprintf("   FAIL  op %2d %s: %s\n", i, describe(op), err)); total_fail <- total_fail + 1; failed <- TRUE }
  }
  cat(sprintf("   operations passed: %d of %d%s   (%.1f s)\n", npass, length(s$ops), if (failed) " - later operations not checked" else "", as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("\nGATE G0: %s (%d failing checks)\n", if (total_fail == 0) "PASS" else "NOT YET", total_fail))
quit(status = if (total_fail == 0) 0 else 1)
