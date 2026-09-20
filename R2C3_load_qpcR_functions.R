## Load qpcR's R-level functions directly (bypassing the package's rgl Depends,
## which fails to install/load in some environments, e.g. no OpenGL/X11 dev
## headers available). qpcR is GPL>=2 licensed; its non-rgl dependencies
## (MASS, minpack.lm, robustbase, Matrix, methods) are ordinary CRAN packages.
## This gives the actual published Ritz & Spiess (2008) sliding-window
## efficiency algorithm instead of a hand-rolled reimplementation.
##
## If `library(qpcR)` already works in your environment, skip this script
## entirely and use that instead -- this workaround is only needed when qpcR's
## rgl dependency blocks a normal install/load.
##
## SRC_DIR must point to a local copy of qpcR's R/ source directory. Obtain
## one by either:
##   (a) downloading the qpcR source tarball from CRAN's archive
##       (https://cran.r-project.org/src/contrib/Archive/qpcR/) and extracting
##       its R/ subfolder, or
##   (b) if qpcR is already installed with source available,
##       file.path(find.package("qpcR"), "R") on some installations.
## Edit the path below to match.

library(MASS)
library(minpack.lm)
library(robustbase)
library(Matrix)
library(methods)

SRC_DIR <- "./qpcR_src/R"  # <-- edit to your local qpcR R/ source directory

if (!dir.exists(SRC_DIR)) {
  stop(sprintf(
    "SRC_DIR ('%s') not found. Set it to a local qpcR R/ source directory -- see this script's header comment.",
    SRC_DIR
  ))
}

rfiles <- list.files(SRC_DIR, pattern = "\\.[Rr]$", full.names = TRUE)
cat(sprintf("Sourcing %d files...\n", length(rfiles)))

ok <- c(); failed <- list()
for (f in rfiles) {
  res <- try(source(f), silent = TRUE)
  if (inherits(res, "try-error")) {
    failed[[basename(f)]] <- attr(res, "condition")$message
  } else {
    ok <- c(ok, basename(f))
  }
}
cat(sprintf("Sourced OK: %d\n", length(ok)))
if (length(failed) > 0) {
  cat("Failed on first pass:\n")
  for (n in names(failed)) cat(" -", n, ":", failed[[n]], "\n")
}

## load internal package data (model definitions etc.)
sysdata_path <- file.path(SRC_DIR, "sysdata.rda")
if (file.exists(sysdata_path)) {
  load(sysdata_path, envir = .GlobalEnv)
  cat("Loaded sysdata.rda\n")
}

## retry failed files once (order-of-definition issues resolve once everything else is loaded)
if (length(failed) > 0) {
  still_failed <- list()
  for (n in names(failed)) {
    f <- file.path(SRC_DIR, n)
    res <- try(source(f), silent = TRUE)
    if (inherits(res, "try-error")) {
      still_failed[[n]] <- attr(res, "condition")$message
    }
  }
  cat(sprintf("\nAfter retry, still failed: %d\n", length(still_failed)))
  for (n in names(still_failed)) cat(" -", n, ":", still_failed[[n]], "\n")
}

cat("\nDone. Testing pcrfit + sliwin availability:\n")
cat("pcrfit exists:", exists("pcrfit"), "\n")
cat("sliwin exists:", exists("sliwin"), "\n")
cat("efficiency exists:", exists("efficiency"), "\n")
cat("takeoff exists:", exists("takeoff"), "\n")
