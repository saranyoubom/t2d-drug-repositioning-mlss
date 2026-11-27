## agree.R - compare results of the R port with the saved Python results (the "agreement report" of R_PORT_PLAN.md).
## agree_csv(): a Python CSV (floats written with %.6g) against an R data frame or CSV; agree_json(): a Python JSON against an R list or JSON.
## Numbers agree when they are equal or |a - b| <= ABS_FLOOR + tol * max(|a|, |b|): a relative tolerance (default 2e-6 for CSV, 2e-9 for JSON) with an absolute floor of 1e-12 for rounding noise around 0; strings, logicals and NA must be equal.
## Changed 26 Nov 2026 (independent audit, audit_261126/AUDIT_REPORT.md O3): the rule was |a - b| <= tol * max(1, |b|), an ABSOLUTE tolerance below 1 (a P value of 1e-7 agreed with one of 1.5e-6).
## The largest relative difference printed is |a - b| / max(|a|, |b|) over the values that differ by more than the floor (0 when none does).
suppressMessages(library(jsonlite))
ABS_FLOOR <- 1e-12
.is_na_like <- function(x) is.null(x) || length(x) == 0 || (length(x) == 1 && (is.na(x) || identical(x, "")))
.as_logical_py <- function(x) if (is.character(x) && all(x %in% c("True", "False", NA, ""))) ifelse(x == "", NA, x == "True") else x
agree_csv <- function(py_file, r, tol = 2e-6, key = NULL, label = basename(py_file), quiet = FALSE) {
  py <- read.csv(py_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM")
  rd <- if (is.character(r)) read.csv(r, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM") else as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE)
  msgs <- character(0)
  if (!identical(names(py), names(rd))) msgs <- c(msgs, sprintf("columns differ: Python [%s] R [%s]", paste(names(py), collapse = ","), paste(names(rd), collapse = ",")))
  if (nrow(py) != nrow(rd)) msgs <- c(msgs, sprintf("rows: Python %d, R %d", nrow(py), nrow(rd)))
  worst <- 0; nbad <- 0; ncmp <- 0
  if (!length(msgs)) {
    if (!is.null(key)) { py <- py[do.call(order, py[key]), , drop = FALSE]; rd <- rd[do.call(order, rd[key]), , drop = FALSE] }
    for (cn in names(py)) {
      a <- .as_logical_py(py[[cn]]); b <- .as_logical_py(rd[[cn]])
      if (is.numeric(a) && (is.numeric(b) || is.logical(b))) {
        b <- as.numeric(b); both_na <- is.na(a) & is.na(b); one_na <- xor(is.na(a), is.na(b)); same <- !is.na(a) & !is.na(b) & a == b
        dif <- abs(a - b); m <- pmax(abs(a), abs(b)); d <- ifelse(!same & !is.na(dif) & is.finite(dif) & dif > ABS_FLOOR, dif / m, 0)
        bad <- which(!same & !both_na & (one_na | is.na(dif) | !is.finite(dif) | dif > ABS_FLOOR + tol * m)); ncmp <- ncmp + length(a); worst <- max(worst, suppressWarnings(max(d[is.finite(d)], 0)))
        if (length(bad)) { nbad <- nbad + length(bad); msgs <- c(msgs, sprintf("column %s: %d values differ, first at row %d (Python %s, R %s)", cn, length(bad), bad[1], format(a[bad[1]]), format(b[bad[1]]))) }
      } else {
        a <- as.character(a); b <- as.character(b); ncmp <- ncmp + length(a); bad <- which(!((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)))
        if (length(bad)) { nbad <- nbad + length(bad); msgs <- c(msgs, sprintf("column %s: %d entries differ, first at row %d (Python '%s', R '%s')", cn, length(bad), bad[1], a[bad[1]], b[bad[1]])) }
      }
    }
  }
  ok <- length(msgs) == 0
  if (!quiet) cat(sprintf("%-6s %-42s %d rows x %d columns, %d values compared, largest relative difference %.2g%s\n", if (ok) "AGREE" else "DIFFER", label, nrow(py), ncol(py), ncmp, worst, if (ok) "" else paste0("\n         ", paste(head(msgs, 4), collapse = "\n         "))))
  invisible(list(ok = ok, msgs = msgs, n = ncmp, worst = worst))
}
agree_json <- function(py_file, r, tol = 2e-9, label = basename(py_file), skip = character(0), quiet = FALSE) {
  py <- if (is.character(py_file)) fromJSON(py_file, simplifyVector = FALSE) else py_file; rr <- if (is.character(r)) fromJSON(r, simplifyVector = FALSE) else r   # py_file may be an already loaded list (give label)
  msgs <- character(0); n <- 0; worst <- 0
  walk <- function(a, b, path) {
    if (path %in% skip) return(invisible())
    if (is.null(a) || (length(a) == 1 && is.atomic(a) && is.na(a))) { n <<- n + 1; if (!.is_na_like(b) && !(is.numeric(b) && is.nan(b))) msgs <<- c(msgs, sprintf("%s: Python null, R %s", path, format(b)[1])); return(invisible()) }
    if (is.list(a)) {
      if (!is.list(b) && !is.atomic(b)) { msgs <<- c(msgs, sprintf("%s: Python list, R %s", path, class(b)[1])); return(invisible()) }
      if (is.atomic(b) && !is.list(b)) b <- as.list(b)
      if (!is.null(names(a))) {
        miss <- setdiff(names(a), names(b)); extra <- setdiff(names(b), names(a))
        if (length(miss)) msgs <<- c(msgs, sprintf("%s: keys missing in R: %s", path, paste(head(miss, 5), collapse = ",")))
        if (length(extra)) msgs <<- c(msgs, sprintf("%s: extra keys in R: %s", path, paste(head(extra, 5), collapse = ",")))
        for (k in intersect(names(a), names(b))) walk(a[[k]], b[[k]], paste0(path, "/", k))
      } else {
        if (length(a) != length(b)) { msgs <<- c(msgs, sprintf("%s: length Python %d, R %d", path, length(a), length(b))); return(invisible()) }
        for (i in seq_along(a)) walk(a[[i]], b[[i]], paste0(path, "[", i, "]"))
      }
    } else if (is.numeric(a)) {
      n <<- n + 1
      if (!(is.numeric(b) || is.logical(b)) || length(b) != 1) { msgs <<- c(msgs, sprintf("%s: Python number %s, R %s", path, format(a), paste(format(b), collapse = " "))); return(invisible()) }
      bv <- as.numeric(b); if (!is.na(bv) && isTRUE(a == bv)) return(invisible())
      dif <- abs(a - bv); m <- max(abs(a), abs(bv)); if (!is.na(dif) && is.finite(dif) && dif > ABS_FLOOR) worst <<- max(worst, dif / m)
      if (is.na(dif) || !is.finite(dif) || dif > ABS_FLOOR + tol * m) msgs <<- c(msgs, sprintf("%s: Python %s, R %s", path, format(a, digits = 15), format(bv, digits = 15)))
    } else { n <<- n + 1; if (!identical(as.character(a), as.character(b))) msgs <<- c(msgs, sprintf("%s: Python '%s', R '%s'", path, as.character(a)[1], as.character(b)[1])) }
  }
  walk(py, rr, "")
  ok <- length(msgs) == 0
  if (!quiet) cat(sprintf("%-6s %-42s %d values compared, largest relative difference %.2g%s\n", if (ok) "AGREE" else "DIFFER", label, n, worst, if (ok) "" else paste0("\n         ", paste(head(msgs, 6), collapse = "\n         "))))
  invisible(list(ok = ok, msgs = msgs, n = n, worst = worst))
}
