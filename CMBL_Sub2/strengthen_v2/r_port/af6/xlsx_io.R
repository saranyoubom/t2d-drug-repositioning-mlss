## xlsx_io.R - pandas-like tables and an OOXML workbook writer for build_af6.R. The cells are written exactly as openpyxl writes them: inline strings, numbers formatted '%.16g', booleans as t="b", Arial 9 (bold in header rows), column widths and a frozen header row.
## Table: list(names, cols, dtype); cols = one list per column with one R value per row (integer, double, logical, character; NA = missing); dtype = "int" | "float" | "bool" | "obj" (the pandas dtype, which decides str() for the column widths).
suppressMessages({ library(xml2); library(zip) })
pint <- function(x) if (abs(x) <= 2147483647) as.integer(x) else structure(x, class = "bigint")
missing_v <- function(v) is.null(v) || (length(v) == 1L && is.atomic(v) && is.na(v))
pstr <- function(v) if (is.logical(v)) (if (v) "True" else "False") else if (is.integer(v)) as.character(v) else if (inherits(v, "bigint")) format(unclass(v), scientific = FALSE) else if (is.numeric(v)) .json_num(v) else v      # Python str() of a cell value
tbl <- function(names, cols, dtype) list(names = names, cols = cols, dtype = dtype)
nrow_t <- function(t) if (length(t$cols)) length(t$cols[[1]]) else 0L
## ---- pandas.read_csv
PD_NA <- c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null")
.E10 <- as.numeric(paste0("1e", 0:308))
pd_strtod <- function(s) vapply(s, function(z) {      # pandas precise_xstrtod (the default float parser of read_csv): at most 17 significant digits accumulated in a double, then one multiplication or division by a power of ten
  ch <- utf8ToInt(z); n <- length(ch); i <- 1L; neg <- FALSE; while (i <= n && ch[i] == 32L) i <- i + 1L
  if (i <= n && (ch[i] == 45L || ch[i] == 43L)) { neg <- ch[i] == 45L; i <- i + 1L }
  if (i <= n && (ch[i] == 73L || ch[i] == 105L)) return(if (neg) -Inf else Inf)
  num <- 0; ex <- 0L; nd <- 0L; ndec <- 0L
  while (i <= n && ch[i] >= 48L && ch[i] <= 57L) { if (nd < 17L) { num <- num * 10 + (ch[i] - 48L); nd <- nd + 1L } else ex <- ex + 1L; i <- i + 1L }
  if (i <= n && ch[i] == 46L) { i <- i + 1L; while (nd < 17L && i <= n && ch[i] >= 48L && ch[i] <= 57L) { num <- num * 10 + (ch[i] - 48L); i <- i + 1L; nd <- nd + 1L; ndec <- ndec + 1L }; if (nd >= 17L) while (i <= n && ch[i] >= 48L && ch[i] <= 57L) i <- i + 1L; ex <- ex - ndec }
  if (neg) num <- -num
  if (i <= n && (ch[i] == 101L || ch[i] == 69L)) { i <- i + 1L; eneg <- FALSE; if (i <= n && (ch[i] == 45L || ch[i] == 43L)) { eneg <- ch[i] == 45L; i <- i + 1L }; e <- 0L; while (i <= n && ch[i] >= 48L && ch[i] <= 57L) { e <- e * 10L + (ch[i] - 48L); i <- i + 1L }; ex <- ex + (if (eneg) -e else e) }
  if (ex > 308L) (if (neg) -Inf else Inf) else if (ex > 0L) num * .E10[ex + 1L] else if (ex < -308L) num / .E10[-308L - ex + 1L] / .E10[309L] else num / .E10[-ex + 1L]
}, 1, USE.NAMES = FALSE)
read_rows <- function(path, sep = ",") {      # all rows of a delimited text file as a character matrix-like data frame (header included), quoted fields may hold separators and line breaks
  L <- readLines(path, warn = FALSE, encoding = "UTF-8"); if (!length(L)) return(data.frame()); L[1] <- sub("^﻿", "", L[1])
  tc <- textConnection(L); nf <- count.fields(tc, sep = sep, quote = "\"", comment.char = ""); close(tc); n <- max(nf, na.rm = TRUE)
  read.table(text = L, sep = sep, quote = "\"", header = FALSE, colClasses = "character", na.strings = character(0), comment.char = "", fill = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, col.names = paste0("V", seq_len(n)), check.names = FALSE, stringsAsFactors = FALSE)
}
infer_col <- function(s, na_values) {      # pandas type inference of one column of strings
  na <- s %in% na_values; v <- s[!na]; vt <- trimws(v); out <- rep(list(NA), length(s))
  if (!length(v)) return(list(type = "float", v = out))
  if (all(grepl("^[+-]?[0-9]+$", vt))) { num <- as.numeric(vt); if (!any(na)) return(list(type = "int", v = lapply(num, pint))); out[!na] <- as.list(num); return(list(type = "float", v = out)) }
  if (all(grepl("^[+-]?(([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?|[Ii][Nn][Ff]([Ii][Nn][Ii][Tt][Yy])?)$", vt))) { out[!na] <- as.list(pd_strtod(vt)); return(list(type = "float", v = out)) }
  if (all(v %in% c("True", "False", "TRUE", "FALSE", "true", "false"))) { out[!na] <- as.list(tolower(v) == "true"); return(list(type = if (any(na)) "obj" else "bool", v = out)) }
  out[!na] <- as.list(v); list(type = "obj", v = out)
}
mangle <- function(h) { seen <- list(); for (i in seq_along(h)) { k <- h[i]; c0 <- if (is.null(seen[[k]])) 0L else seen[[k]]; if (c0 > 0L) h[i] <- paste0(k, ".", c0); seen[[k]] <- c0 + 1L }; h }
pd_read_csv <- function(path, sep = ",", na_values = PD_NA) {
  x <- read_rows(path, sep); hdr <- unlist(x[1, ], use.names = FALSE); blankh <- which(hdr == ""); hdr[blankh] <- paste0("Unnamed: ", blankh - 1L)      # pandas names an empty header 'Unnamed: n'
  x <- x[-1, , drop = FALSE]; res <- lapply(seq_along(hdr), function(j) infer_col(x[[j]], na_values))
  tbl(mangle(hdr), lapply(res, function(r) r$v), vapply(res, function(r) r$type, ""))
}
## ---- DataFrame(list of dicts), concat, filter, merge
kind_of <- function(x) if (missing_v(x)) "na" else if (is.logical(x)) "bool" else if (is.integer(x) || inherits(x, "bigint")) "int" else if (is.numeric(x)) "float" else "str"
unify <- function(v, kinds_dtype, any_missing) {      # pandas dtype of a column made of parts with the given dtypes (and missing values)
  ks <- unique(kinds_dtype)
  if (!length(ks)) "obj" else if (all(ks == "int")) (if (any_missing) "float" else "int") else if (all(ks %in% c("int", "float"))) "float" else if (all(ks == "bool")) (if (any_missing) "obj" else "bool") else "obj"
}
as_float_col <- function(v) lapply(v, function(x) if (missing_v(x)) NA else as.numeric(unclass(x)))
df_from_records <- function(recs) {
  keys <- unique(unlist(lapply(recs, names))); cols <- list(); dtype <- character(0)
  for (k in keys) { v <- lapply(recs, function(r) { x <- r[[k]]; if (is.null(x)) NA else x }); kinds <- vapply(v, kind_of, ""); ty <- unify(v, kinds[kinds != "na"], any(kinds == "na")); if (ty == "float") v <- as_float_col(v); cols[[length(cols) + 1L]] <- v; dtype <- c(dtype, ty) }
  tbl(keys, cols, dtype)
}
concat <- function(ts) {
  keys <- unique(unlist(lapply(ts, function(t) t$names))); cols <- list(); dtype <- character(0)
  for (k in keys) {
    parts <- lapply(ts, function(t) { j <- match(k, t$names); if (is.na(j)) list(v = rep(list(NA), nrow_t(t)), d = character(0), miss = nrow_t(t) > 0) else list(v = t$cols[[j]], d = t$dtype[j], miss = any(vapply(t$cols[[j]], missing_v, TRUE)) && t$dtype[j] != "float") })
    v <- unlist(lapply(parts, function(p) p$v), recursive = FALSE); d <- unlist(lapply(parts, function(p) p$d)); ty <- unify(v, d, any(vapply(parts, function(p) p$miss, TRUE))); if (ty == "float") v <- as_float_col(v)
    cols[[length(cols) + 1L]] <- v; dtype <- c(dtype, ty)
  }
  tbl(keys, cols, dtype)
}
filter_rows <- function(t, keep) tbl(t$names, lapply(t$cols, function(v) v[keep]), t$dtype)
col_of <- function(t, name) t$cols[[match(name, t$names)]]
chr_of <- function(t, name) vapply(col_of(t, name), function(x) if (missing_v(x)) NA_character_ else as.character(x), "")
drop_col <- function(t, name) { j <- match(name, t$names); tbl(t$names[-j], t$cols[-j], t$dtype[-j]) }
insert_col <- function(t, name, v, dtype, at = 1L) { n <- nrow_t(t); v <- if (length(v) == 1L && n != 1L) rep(list(v), n) else as.list(v); j <- at - 1L; tbl(append(t$names, name, j), append(t$cols, list(v), j), append(t$dtype, dtype, j)) }
merge_left <- function(left, right, on) {      # left join on the key columns, keeping the order of the left rows; the first matching right row
  key <- function(t) do.call(paste, c(lapply(on, function(k) chr_of(t, k)), sep = "\r")); m <- match(key(left), key(right))
  for (k in setdiff(right$names, on)) { v <- lapply(m, function(i) if (is.na(i)) NA else col_of(right, k)[[i]]); left <- tbl(c(left$names, k), c(left$cols, list(v)), c(left$dtype, right$dtype[match(k, right$names)])) }
  left
}
## ---- the workbook: a list of sheets (environments); cells are added in chunks (row, col, value, bold) so that large sheets stay fast
sheet_new <- function(name) { e <- new.env(); e$name <- name; e$chunks <- list(); e$widths <- numeric(0); e$wcol <- integer(0); e$freeze <- 1L; e }
sheet_cells <- function(s, r, c, v, bold = FALSE) { s$chunks[[length(s$chunks) + 1L]] <- list(r = r, c = c, v = v, b = rep_len(bold, length(r))); invisible(s) }
sheet_cell <- function(s, r, c, v, bold = FALSE) sheet_cells(s, r, c, list(v), bold)
sheet_width <- function(s, col, w) { i <- match(col, s$wcol); if (is.na(i)) { s$wcol <- c(s$wcol, col); s$widths <- c(s$widths, w) } else s$widths[i] <- w; invisible(s) }
col_letter <- function(j) vapply(j, function(k) { s <- ""; while (k > 0) { m <- (k - 1) %% 26; s <- paste0(LETTERS[m + 1], s); k <- (k - 1) %/% 26 }; s }, "")
xml_esc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE), fixed = TRUE)
cell_xml <- function(ref, v, style) {
  if (is.logical(v)) sprintf("<c r=\"%s\" s=\"%d\" t=\"b\"><v>%d</v></c>", ref, style, as.integer(v))
  else if (is.numeric(v) || inherits(v, "bigint")) sprintf("<c r=\"%s\" s=\"%d\" t=\"n\"><v>%s</v></c>", ref, style, { x <- as.numeric(unclass(v)); if (is.nan(x) || !is.finite(x)) "" else { s <- py_fmt_g(x, 16L); if (s == "-0") "0" else s } })      # openpyxl's later load / save turns "-0" into the integer 0
  else if (v == "") sprintf("<c r=\"%s\" s=\"%d\" t=\"inlineStr\"></c>", ref, style)      # openpyxl writes an empty text without an <is> element
  else { sp <- if (v != trimws(v) || grepl("\n", v, fixed = TRUE)) " xml:space=\"preserve\"" else ""; sprintf("<c r=\"%s\" s=\"%d\" t=\"inlineStr\"><is><t%s>%s</t></is></c>", ref, style, sp, xml_esc(v)) }
}
sheet_xml <- function(s) {
  R <- unlist(lapply(s$chunks, function(k) k$r)); C <- unlist(lapply(s$chunks, function(k) k$c)); V <- unlist(lapply(s$chunks, function(k) k$v), recursive = FALSE); B <- unlist(lapply(s$chunks, function(k) k$b))
  o <- order(R, C); r <- R[o]; cc <- C[o]; v <- V[o]; b <- B[o]; ref <- if (length(cc)) paste0(col_letter(seq_len(max(cc)))[cc], r) else character(0)
  cells <- mapply(cell_xml, ref, v, ifelse(b, 1L, 2L), SIMPLIFY = TRUE, USE.NAMES = FALSE)
  rows <- if (length(r)) vapply(split(cells, factor(r, levels = unique(r))), paste, "", collapse = "") else character(0)
  rowx <- paste0(sprintf("<row r=\"%d\">", as.integer(names(rows))), rows, "</row>", collapse = "")
  w <- order(s$wcol); colsx <- if (length(s$wcol)) paste0("<cols>", paste0(sprintf("<col width=\"%s\" customWidth=\"1\" min=\"%d\" max=\"%d\"/>", vapply(s$widths[w], function(x) py_fmt_g(x, 16L), ""), s$wcol[w], s$wcol[w]), collapse = ""), "</cols>") else ""
  top <- sprintf("A%d", s$freeze + 1L); pane <- sprintf("<sheetViews><sheetView workbookViewId=\"0\"><pane ySplit=\"%d\" topLeftCell=\"%s\" activePane=\"bottomLeft\" state=\"frozen\"/><selection pane=\"bottomLeft\" activeCell=\"A1\" sqref=\"A1\"/></sheetView></sheetViews>", s$freeze, top)
  dim <- if (length(r)) sprintf("A1:%s%d", col_letter(max(cc)), max(r)) else "A1"
  paste0("<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\"><sheetPr><outlinePr summaryBelow=\"1\" summaryRight=\"1\"/><pageSetUpPr/></sheetPr><dimension ref=\"", dim, "\"/>", pane,
         "<sheetFormatPr baseColWidth=\"8\" defaultRowHeight=\"15\"/>", colsx, "<sheetData>", rowx, "</sheetData><pageMargins left=\"0.75\" right=\"0.75\" top=\"1\" bottom=\"1\" header=\"0.5\" footer=\"0.5\"/></worksheet>")
}
STYLES_XML <- paste0("<styleSheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\"><numFmts count=\"0\"/><fonts count=\"3\"><font><name val=\"Calibri\"/><family val=\"2\"/><sz val=\"11\"/></font><font><name val=\"Arial\"/><b val=\"1\"/><sz val=\"9\"/></font><font><name val=\"Arial\"/><sz val=\"9\"/></font></fonts>",
                     "<fills count=\"2\"><fill><patternFill/></fill><fill><patternFill patternType=\"gray125\"/></fill></fills><borders count=\"1\"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\"/></cellStyleXfs>",
                     "<cellXfs count=\"3\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/><xf numFmtId=\"0\" fontId=\"1\" fillId=\"0\" borderId=\"0\" applyFont=\"1\" xfId=\"0\"/><xf numFmtId=\"0\" fontId=\"2\" fillId=\"0\" borderId=\"0\" applyFont=\"1\" xfId=\"0\"/></cellXfs><cellStyles count=\"1\"><cellStyle name=\"Normal\" xfId=\"0\" builtinId=\"0\"/></cellStyles></styleSheet>")
write_workbook <- function(sheets, path) {
  tmp <- tempfile("af6_"); dir.create(file.path(tmp, "_rels"), recursive = TRUE); dir.create(file.path(tmp, "xl", "_rels"), recursive = TRUE); dir.create(file.path(tmp, "xl", "worksheets")); n <- length(sheets)
  put <- function(rel, txt) { con <- file(file.path(tmp, rel), "wb"); writeBin(charToRaw(enc2utf8(paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n", txt))), con); close(con) }
  put("[Content_Types].xml", paste0("<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/>",
    "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/><Override PartName=\"/xl/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>",
    paste0(sprintf("<Override PartName=\"/xl/worksheets/sheet%d.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>", seq_len(n)), collapse = ""), "</Types>"))
  put("_rels/.rels", "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/></Relationships>")
  names <- vapply(sheets, function(s) s$name, "")
  put("xl/workbook.xml", paste0("<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"><bookViews><workbookView tabRatio=\"600\" activeTab=\"0\"/></bookViews><sheets>",
    paste0(sprintf("<sheet name=\"%s\" sheetId=\"%d\" state=\"visible\" r:id=\"rId%d\"/>", xml_esc(names), seq_len(n), seq_len(n)), collapse = ""), "</sheets></workbook>"))
  put("xl/_rels/workbook.xml.rels", paste0("<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">", paste0(sprintf("<Relationship Id=\"rId%d\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet%d.xml\"/>", seq_len(n), seq_len(n)), collapse = ""),
    sprintf("<Relationship Id=\"rId%d\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>", n + 1L), "</Relationships>"))
  put("xl/styles.xml", STYLES_XML)
  for (i in seq_len(n)) put(sprintf("xl/worksheets/sheet%d.xml", i), sheet_xml(sheets[[i]]))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE); if (file.exists(path)) file.remove(path); zpath <- file.path(normalizePath(dirname(path)), basename(path))      # absolute: zip() works inside `root`
  zip::zip(zipfile = zpath, files = c("[Content_Types].xml", "_rels/.rels", "xl/workbook.xml", "xl/_rels/workbook.xml.rels", "xl/styles.xml", sprintf("xl/worksheets/sheet%d.xml", seq_len(n))), root = tmp, mode = "mirror")
  unlink(tmp, recursive = TRUE); invisible(path)
}
