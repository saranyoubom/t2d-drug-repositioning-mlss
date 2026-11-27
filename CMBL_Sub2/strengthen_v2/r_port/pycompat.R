## pycompat.R - small helpers that let the R ports write the same files as the Python scripts (CSV in pandas style, ordered string sort, numpy-style helpers).
.fmt_py_num <- function(x, nd) { out <- rep("", length(x)); ok <- !is.na(x); if (any(ok)) out[ok] <- py_fmt_g(x[ok], nd + 2L); out }   # (py_fmt_g is defined further down in this file)
.csv_field <- function(s) ifelse(is.na(s), "", ifelse(grepl("[\",\n\r]", s), paste0("\"", gsub("\"", "\"\"", s), "\""), s))
## pandas DataFrame.to_csv(index=False, float_format='%.{nd+2}g'): doubles with %.6g, integers as integers, booleans True/False, NA empty, minimal quoting
py_write_csv <- function(df, path, nd = 4) {
  cols <- lapply(df, function(x) {
    if (is.logical(x)) ifelse(is.na(x), "", ifelse(x, "True", "False"))
    else if (is.integer(x)) ifelse(is.na(x), "", as.character(x))
    else if (is.numeric(x)) .fmt_py_num(x, nd)
    else .csv_field(as.character(x))
  })
  lines <- c(paste(.csv_field(names(df)), collapse = ","), if (nrow(df)) do.call(paste, c(unname(cols), sep = ",")))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "wb"); writeLines(lines, con, sep = "\n", useBytes = TRUE); close(con); invisible(path)
}
## Python sorted() on strings (code-point order) and the stable index sort of a numeric vector
py_sort <- function(x) sort(as.character(x), method = "radix")
py_order <- function(...) order(..., method = "radix")
## decimal text -> double, correctly rounded like Python's float(): R's as.numeric rounds twice (long double) and is one unit in the last place off for some decimals ("0.529113"); jsonlite parses with the C library's strtod
.parse_exact <- function(s) as.numeric(unlist(jsonlite::fromJSON(paste0("[", paste(s, collapse = ","), "]"), simplifyVector = FALSE), use.names = FALSE))      # s must be valid JSON numbers
py_float <- function(s) {      # Python float(text) for a character vector (leading +, ".5", "5.", inf, nan accepted; anything else is NA)
  t <- trimws(as.character(s)); out <- rep(NA_real_, length(t))
  inf <- grepl("^[+-]?inf(inity)?$", t, ignore.case = TRUE); out[inf] <- ifelse(startsWith(t[inf], "-"), -Inf, Inf); out[grepl("^[+-]?nan$", t, ignore.case = TRUE)] <- NaN
  num <- grepl("^[+-]?(([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?)$", t)
  if (any(num)) { u <- sub("^[+]", "", t[num]); u <- sub("^(-?)[.]", "\\10.", u); u <- sub("[.]($|[eE])", ".0\\1", u); u <- sub("^(-?)0+([0-9])", "\\1\\2", u); out[num] <- .parse_exact(u) }
  out
}
## json.dump(indent=1)-compatible writer for nested lists (floats with up to 17 significant digits, NaN and Inf as null like the `clean` helper of the scripts)
.json_num <- function(x) {
  if (is.na(x) || !is.finite(x)) {          # options(py_json_nan = TRUE): Python json.dump(allow_nan=True) style, NaN -> NaN, Inf -> Infinity (R NA stays null = Python None)
    if (isTRUE(getOption("py_json_nan"))) { if (is.nan(x)) return("NaN"); if (!is.na(x)) return(if (x > 0) "Infinity" else "-Infinity") }
    return("null")
  }
  if (x == 0) return(if (1 / x < 0) "-0.0" else "0.0")
  if (x == floor(x) && abs(x) < 1e15) return(paste0(sprintf("%.0f", x), ".0"))      # whole numbers (exact below 2^53)
  ## the shortest digits that round-trip (the correctly rounded 1 to 17 digit candidates from the exact binary value, exact ties to even, checked with a correctly rounded parser: R's sprintf and as.numeric are not), laid out like Python's repr(float):
  ## fixed notation for decimal exponents -4 to 15 (with ".0" for whole numbers), else d.ddde+XX
  neg <- x < 0; ed <- .exact_digits(abs(x)); cd <- vector("list", 17L)
  for (p in 1:17) { r <- .round_digits(ed$d, p); cd[[p]] <- list(d = if (r$carry) substr(r$d, 1, p) else r$d, e = ed$e + as.integer(r$carry)) }
  i <- which(.parse_exact(vapply(cd, function(z) paste0(substr(z$d, 1, 1), ".", if (nchar(z$d) > 1) substring(z$d, 2) else "0", "e", z$e), "")) == abs(x))[1]
  ex <- cd[[i]]$e; digs <- sub("0+$", "", cd[[i]]$d); if (!nzchar(digs)) digs <- "0"; nd <- nchar(digs)
  r <- if (ex >= -4 && ex < 16) {
    if (ex >= 0) paste0(if (nd <= ex + 1) paste0(digs, strrep("0", ex + 1 - nd)) else substr(digs, 1, ex + 1), ".", if (nd > ex + 1) substring(digs, ex + 2) else "0")
    else paste0("0.", strrep("0", -ex - 1), digs)
  } else paste0(substr(digs, 1, 1), if (nd > 1) paste0(".", substring(digs, 2)) else "", "e", if (ex < 0) "-" else "+", sprintf("%02d", abs(ex)))
  paste0(if (neg) "-" else "", r)
}
.json_str <- function(s) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", s)), "\"")
py_json <- function(o, indent = 1, level = 0) {
  pad <- strrep(" ", indent * (level + 1)); padc <- strrep(" ", indent * level)
  if (is.null(o)) return("null")
  if (is.list(o)) {
    if (length(o) == 0) return(if (is.null(names(o))) "[]" else "{}")
    if (!is.null(names(o))) return(paste0("{\n", paste0(pad, .json_str(names(o)), ": ", vapply(o, py_json, "", indent = indent, level = level + 1), collapse = ",\n"), "\n", padc, "}"))
    return(paste0("[\n", paste0(pad, vapply(o, py_json, "", indent = indent, level = level + 1), collapse = ",\n"), "\n", padc, "]"))
  }
  if (length(o) != 1) return(py_json(as.list(o), indent, level))
  if (is.logical(o)) return(if (is.na(o)) "null" else if (o) "true" else "false")
  if (is.numeric(o)) return(if (is.integer(o)) as.character(o) else if (!is.na(o) && is.finite(o) && o == round(o) && abs(o) < 1e15) sprintf("%.0f", o) else .json_num(o))
  .json_str(as.character(o))
}
py_write_json <- function(o, path) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); writeLines(py_json(o), path, useBytes = TRUE); invisible(path) }
## numpy-style helpers
np_percentile <- function(x, p) unname(quantile(x, p / 100, type = 7))
np_std <- function(x) sqrt(mean((x - mean(x))^2))
triu_rowmajor <- function(n) t(combn(n, 2))      # numpy.triu_indices(n, 1) as an (n(n-1)/2) x 2 matrix, 1-based, row-major order

## pandas DataFrame.to_csv(index=False) with the default float format (shortest round-trip repr, like Python's repr); NaN and NA empty, booleans True/False, minimal quoting, rows end with "\n"
py_write_csv_repr <- function(df, path, na = "", eol = "\n") {         # eol = "\r\n" for pandas on Windows (to_csv writes os.linesep); an integer column with a missing value is a float column in pandas ("18.0")
  cols <- lapply(df, function(x) {
    if (is.logical(x)) ifelse(is.na(x), na, ifelse(x, "True", "False"))
    else if (is.integer(x) && !anyNA(x)) as.character(x)
    else if (is.integer(x)) vapply(x, function(v) if (is.na(v)) na else .json_num(as.numeric(v)), "")
    else if (is.numeric(x)) vapply(x, function(v) if (is.na(v)) na else .json_num(v), "")
    else .csv_field(ifelse(is.na(x), na, as.character(x)))
  })
  lines <- c(paste(.csv_field(names(df)), collapse = ","), if (nrow(df)) do.call(paste, c(unname(cols), sep = ",")))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); con <- file(path, "wb"); writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = eol), eol))), con); close(con); invisible(path)
}

## ---------------------------------------------------------------- typed writers (moved here from effect_sizes_part4.R on 23 Nov 2026; used by several ports)
## Typed JSON: R integer = Python int, R double = Python float (whole floats keep ".0"), character = str. indent = 1 reproduces json.dump(indent=1); indent = NULL the compact one-line json.dump default
## (", " and ": " separators). Length-1 vectors are scalars: pass list(x) for a one-element Python list.
py_json_typed <- function(o, indent = 1, level = 0) {
  if (is.null(o)) return("null")
  if (is.list(o)) {
    if (length(o) == 0) return(if (is.null(names(o))) "[]" else "{}")
    items <- vapply(o, py_json_typed, "", indent = indent, level = level + 1); isd <- !is.null(names(o))
    kv <- if (isd) paste0(.json_str(names(o)), ": ", items) else items; open <- if (isd) "{" else "["; close <- if (isd) "}" else "]"
    if (is.null(indent)) return(paste0(open, paste(kv, collapse = ", "), close))
    return(paste0(open, "\n", paste0(strrep(" ", indent * (level + 1)), kv, collapse = ",\n"), "\n", strrep(" ", indent * level), close))
  }
  if (length(o) != 1) return(py_json_typed(as.list(o), indent, level))
  if (is.logical(o)) return(if (is.na(o)) "null" else if (o) "true" else "false")
  if (is.integer(o)) return(as.character(o))
  if (is.numeric(o)) return(.json_num(o))
  .json_str(as.character(o))
}
.NL <- if (.Platform$OS.type == "windows") "\r\n" else "\n"      # Python text-mode files (md, json) end lines with os.linesep; csv.DictWriter always writes CRLF
write_text_bytes <- function(txt, path) { dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); con <- file(path, "wb"); writeBin(charToRaw(enc2utf8(txt)), con); close(con); invisible(path) }
write_json_typed <- function(o, path, indent = 1) write_text_bytes(gsub("\n", .NL, py_json_typed(o, indent), fixed = TRUE), path)
## csv.DictWriter (default dialect) with the value rule of the scripts: floats '%.6g', everything else str(); None -> empty; fields with , " CR LF are quoted; lines end with CRLF
.dw_val <- function(v, fmt = "%.6g") {
  if (is.null(v) || length(v) == 0) return("")
  if (is.logical(v)) return(if (is.na(v)) "" else if (v) "True" else "False")
  if (is.integer(v)) return(as.character(v))
  if (is.numeric(v)) {
    if (is.nan(v)) return("nan"); if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
    if (identical(fmt, "repr")) return(.json_num(v))                                          # str(float): the shortest repr, as csv.DictWriter writes a float that was not formatted
    m <- regmatches(fmt, regexec("^%[.]([0-9]+)([fg])$", fmt))[[1]]                      # '%.Ng' / '%.Nf': Python's correctly rounded formatter (see py_fmt_g)
    if (length(m)) return(if (m[3] == "g") py_fmt_g(v, as.integer(m[2])) else py_fmt_f(v, as.integer(m[2])))
    return(sprintf(fmt, v))
  }
  as.character(v)
}
.dw_quote <- function(s) if (grepl("[\",\r\n]", s)) paste0("\"", gsub("\"", "\"\"", s, fixed = TRUE), "\"") else s
dictwriter_csv <- function(recs, cols, path, fmt = "%.6g", bom = FALSE) {
  lines <- c(paste(vapply(cols, .dw_quote, ""), collapse = ","), vapply(recs, function(r) paste(vapply(cols, function(k) .dw_quote(.dw_val(r[[k]], fmt)), ""), collapse = ","), ""))
  write_text_bytes(paste0(if (bom) "﻿" else "", paste(lines, collapse = "\r\n"), "\r\n"), path)
}

## numpy.percentile(x, 100 * p) with the default 'linear' method, operation for operation (numpy's _lerp: a + (b - a) * t, and b - (b - a) * (1 - t) for t >= 0.5); R's quantile(type = 7)
## computes (1 - h) * a + h * b, which differs in the last bit. NaN values must be removed by the caller.
np_quantile <- function(x, p) {
  x <- sort(x); n <- length(x); v <- (n - 1) * p; lo <- floor(v); t <- v - lo; a <- x[lo + 1]; b <- x[min(lo + 2, n)]; d <- b - a
  r <- a + d * t; if (t >= 0.5) r <- b - d * (1 - t); r
}
## statsmodels proportion_confint(method = 'wilson'), same operation order; z = scipy.stats.norm.isf(0.025) (R's qnorm(0.975) differs by 9e-16)
PY_Z975 <- 1.9599639845400545
py_wilson <- function(k, n, z = PY_Z975) {
  if (n == 0) return(c(NaN, NaN))
  q <- k / n; z2 <- z^2; den <- 1 + z2 / n; centre <- (q + z2 / (2 * n)) / den; dist <- z * sqrt(q * (1 - q) / n + z2 / (4 * n^2)); dist <- dist / den
  c(centre - dist, centre + dist)
}

## ---------------------------------------------------------------- correctly rounded printf (Python's '%.Nf' and '%.Ng')
## R on Windows (UCRT) rounds near-ties differently from Python (David Gay's dtoa: exact binary value, exact ties to even): sprintf("%.1f", 1.85) is "1.8" in R and "1.9" in Python
## (1.85 is 1.850000000000000088...), and sprintf("%.6g", 1234565) is "1.23457e+06" against Python's "1.23456e+06". UCRT's printed digits are not exact either for small numbers
## (6.5753507660701875e-11 is wrong from the 18th digit), so the exact decimal expansion of the double (m * 2^e, a finite decimal) is computed here with base-1e7 integer arithmetic
## and rounded here (half to even only when the remaining digits are exactly 5000...).
.PW <- new.env()
.big_norm <- function(r) { repeat { cr <- floor(r / 1e7); if (!any(cr > 0)) break; r <- c(r - cr * 1e7, 0) + c(0, cr) }; n <- length(r); while (n > 1L && r[n] == 0) n <- n - 1L; r[seq_len(n)] }      # carry propagation of little-endian base-1e7 limbs
.big_pow <- function(base, k) { key <- paste0("p", base); L <- .PW[[key]]; if (is.null(L)) L <- list(1); if (length(L) <= k) { while (length(L) <= k) L[[length(L) + 1L]] <- .big_norm(L[[length(L)]] * base); assign(key, L, envir = .PW) }; L[[k + 1L]] }
.exact_digits <- function(v) {        # v > 0 finite: list(d = all digits of the exact decimal expansion from the first non-zero digit, e = decimal exponent of that digit)
  e <- max(floor(log2(v)) - 52, -1074); m <- v / 2^e; while (m != floor(m)) { e <- e - 1; m <- m * 2 }      # v = m * 2^e with an integer m
  P <- if (e >= 0) .big_pow(2, e) else .big_pow(5, -e)      # m * 2^e, or m * 5^-e / 10^-e
  ml <- c(m %% 1e7, (m %/% 1e7) %% 1e7, m %/% 1e14); np <- length(P); r <- numeric(np + 2L); for (i in 1:3) r[i:(i + np - 1L)] <- r[i:(i + np - 1L)] + ml[i] * P
  r <- .big_norm(r); n <- length(r); d <- paste0(sprintf("%.0f", r[n]), if (n > 1L) paste(sprintf("%07.0f", r[(n - 1L):1L]), collapse = "") else ""); d <- sub("^0+", "", d)
  list(d = d, e = nchar(d) - 1L + as.integer(min(e, 0)))
}
.round_digits <- function(d, keep) {  # round the digit string d to `keep` digits (keep >= 0); returns list(d = rounded digits (keep digits, or keep + 1 on carry out), carry = TRUE if a digit was added)
  if (keep >= nchar(d)) return(list(d = paste0(d, strrep("0", keep - nchar(d))), carry = FALSE))
  head <- substr(d, 1, keep); rest <- substr(d, keep + 1, nchar(d)); first <- as.integer(substr(rest, 1, 1))
  up <- if (first > 5) TRUE else if (first < 5) FALSE else { if (grepl("[1-9]", substr(rest, 2, nchar(rest)))) TRUE else (keep > 0 && as.integer(substr(head, keep, keep)) %% 2 == 1) }
  if (!up) return(list(d = head, carry = FALSE))
  if (keep == 0) return(list(d = "1", carry = TRUE))
  ch <- as.integer(strsplit(head, "")[[1]]); i <- keep
  while (i >= 1 && ch[i] == 9L) { ch[i] <- 0L; i <- i - 1L }
  if (i == 0) return(list(d = paste0("1", paste(ch, collapse = "")), carry = TRUE))
  ch[i] <- ch[i] + 1L; list(d = paste(ch, collapse = ""), carry = FALSE)
}
py_fmt_f <- function(x, nd) {         # '%.{nd}f' % x
  vapply(x, function(v) {
    if (is.na(v)) return(if (is.nan(v)) "nan" else "None"); if (!is.finite(v)) return(if (v > 0) "inf" else "-inf")
    neg <- v < 0 || (v == 0 && 1 / v < 0); a <- abs(v)
    if (a == 0) return(paste0(if (neg) "-" else "", "0", if (nd > 0) paste0(".", strrep("0", nd)) else ""))
    ed <- .exact_digits(a); e <- ed$e; d <- ed$d
    # position of the rounding point inside d: digits d[1..(e + 1 + nd)] are kept (e + 1 = number of integer digits)
    keep <- e + 1L + nd
    if (keep < 0) { digs <- "0"; ip <- "0"; fp <- strrep("0", nd) }
    else {
      r <- .round_digits(d, keep); digs <- r$d; if (r$carry) e <- e + 1L
      # digs represents the number with (e + 1 + nd) digits, the last nd being decimals
      tot <- nchar(digs); ni <- tot - nd
      if (ni <= 0) { ip <- "0"; fp <- paste0(strrep("0", -ni), digs) } else { ip <- substr(digs, 1, ni); fp <- substr(digs, ni + 1, tot) }
      if (keep == 0 && !r$carry) { ip <- "0"; fp <- strrep("0", nd) }
      if (keep == 0 && r$carry) { ip <- "0"; fp <- paste0(strrep("0", max(nd - 1, 0)), "1"); if (nd == 0) ip <- "1" }
    }
    out <- if (nd > 0) paste0(ip, ".", substr(fp, 1, nd)) else ip
    if (neg && grepl("[1-9]", out)) out <- paste0("-", out) else if (neg) out <- paste0("-", out)
    out
  }, "", USE.NAMES = FALSE)
}
py_fmt_g <- function(x, p = 6L) {     # '%.{p}g' % x
  vapply(x, function(v) {
    if (is.na(v)) return(if (is.nan(v)) "nan" else "None"); if (!is.finite(v)) return(if (v > 0) "inf" else "-inf")
    neg <- v < 0 || (v == 0 && 1 / v < 0); a <- abs(v); if (a == 0) return(if (neg) "-0" else "0")
    ed <- .exact_digits(a); r <- .round_digits(ed$d, p); e <- ed$e; digs <- r$d
    if (r$carry) { e <- e + 1L; digs <- substr(digs, 1, p) }
    sign <- if (neg) "-" else ""
    if (e < -4L || e >= p) {
      m <- sub("0+$", "", substr(digs, 2, p)); ex <- sprintf("%s%02d", if (e < 0) "-" else "+", abs(e))
      paste0(sign, substr(digs, 1, 1), if (nzchar(m)) paste0(".", m) else "", "e", ex)
    } else if (e >= 0) {
      ip <- substr(digs, 1, e + 1); fp <- sub("0+$", "", substr(digs, e + 2, p)); paste0(sign, ip, if (nzchar(fp)) paste0(".", fp) else "")
    } else { fp <- sub("0+$", "", paste0(strrep("0", -e - 1L), digs)); paste0(sign, "0.", fp) }
  }, "", USE.NAMES = FALSE)
}

"%+%" <- function(a, b) paste0(a, b)      # string concatenation like Python's +
sha256_file <- function(path) { if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(path, algo = "sha256", file = TRUE)); h <- system2("certutil", c("-hashfile", shQuote(normalizePath(path)), "SHA256"), stdout = TRUE); tolower(gsub("\\s", "", h[2])) }      # hashlib.sha256 of a file
py_round <- function(x, nd) py_float(py_fmt_f(x, nd))      # Python round(x, nd): the correctly rounded decimal read back as a float (R's round() and as.numeric can differ in the last digit)
## Python's '%' operator for the specs the scripts use: %s %d %i %f %g %% with the flags - + space 0, a width and a precision; numbers go through py_fmt_f / py_fmt_g (correct rounding), %s of a double is repr(float), %d truncates.
## py_sprintf("%-6s %+.2f | %5d", "Wnt5a", 0.123, 12L) = "%-6s %+.2f | %5d" % ('Wnt5a', 0.123, 12). Scalars only (one value per spec).
py_sprintf <- function(fmt, ...) {
  args <- list(...); k <- 0L; pat <- "%([-+ 0#]*)([0-9]+)?(?:[.]([0-9]+))?([sdifg%])"; parts <- regmatches(fmt, gregexpr(pat, fmt, perl = TRUE))[[1]]; lit <- regmatches(fmt, gregexpr(pat, fmt, perl = TRUE), invert = TRUE)[[1]]; out <- lit[1]
  for (i in seq_along(parts)) {
    mm <- regmatches(parts[i], regexec(pat, parts[i], perl = TRUE))[[1]]; flags <- mm[2]; width <- if (nzchar(mm[3])) as.integer(mm[3]) else 0L; prec <- if (nzchar(mm[4])) as.integer(mm[4]) else NA_integer_; ty <- mm[5]
    if (ty == "%") s <- "%" else {
      k <- k + 1L; v <- args[[k]]
      s <- switch(ty,
        s = { x <- if (is.character(v)) v else if (is.logical(v)) (if (isTRUE(v)) "True" else "False") else if (is.integer(v)) as.character(v) else .json_num(v); if (!is.na(prec)) substr(x, 1, prec) else x },
        d = , i = format(trunc(as.numeric(v)), scientific = FALSE, trim = TRUE),
        f = py_fmt_f(as.numeric(v), if (is.na(prec)) 6L else prec),
        g = py_fmt_g(as.numeric(v), if (is.na(prec)) 6L else prec))
      if (ty != "s" && !startsWith(s, "-")) { if (grepl("+", flags, fixed = TRUE)) s <- paste0("+", s) else if (grepl(" ", flags, fixed = TRUE)) s <- paste0(" ", s) }
      if (nchar(s) < width) {
        pad <- width - nchar(s)
        s <- if (grepl("-", flags, fixed = TRUE)) paste0(s, strrep(" ", pad)) else if (grepl("0", flags, fixed = TRUE) && ty != "s") { sg <- sub("^([-+ ]?).*$", "\\1", s); paste0(sg, strrep("0", pad), substring(s, nchar(sg) + 1L)) } else paste0(strrep(" ", pad), s)
      }
    }
    out <- paste0(out, s, lit[i + 1L])
  }
  out
}

## numpy.add.reduce on a contiguous float64 vector (pairwise summation, numpy/core/src/umath/loops_utils.h.src): n < 8 sequential; n <= 128 eight accumulators combined as ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
## then the remainder added one by one; larger n: split at n/2 rounded down to a multiple of 8. R's sum() and mean() accumulate in long double and so can differ in the last bit.
np_sum <- function(a) {
  pw <- function(a) {
    n <- length(a)
    if (n < 8) { res <- 0; for (v in a) res <- res + v; return(res) }
    if (n <= 128) {
      r <- a[1:8]; i <- 8; lim <- n - (n %% 8)
      while (i < lim) { r <- r + a[i + 1:8]; i <- i + 8 }
      res <- ((r[1] + r[2]) + (r[3] + r[4])) + ((r[5] + r[6]) + (r[7] + r[8])); while (i < n) { i <- i + 1; res <- res + a[i] }
      return(res)
    }
    n2 <- n %/% 2; n2 <- n2 - n2 %% 8; pw(a[1:n2]) + pw(a[(n2 + 1):n])
  }
  0 + pw(as.numeric(a))
}
np_mean <- function(a) np_sum(a) / length(a)
np_median <- function(x) { x <- sort(x); n <- length(x); if (n %% 2) x[(n + 1) / 2] else (x[n / 2] + x[n / 2 + 1]) / 2 }      # numpy.median: the middle value, or the mean of the two middle ones
