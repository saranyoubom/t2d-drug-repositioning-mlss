## docx_table.R - read the first table of a Word file through its XML (no docx package needed): the replacement for python-docx `docx.Document(path).tables[0].rows`.
## docx_table0_rows(path) returns a list of rows, each a character vector with the text of the row's cells, with python-docx's cell.text rules: runs of a paragraph joined, <w:br/> as "\n", <w:tab/> as "\t",
## paragraphs joined by "\n". Horizontally merged cells appear once (a single <w:tc>; python-docx repeats them, so a section row has length < the table width here where python-docx has c[0] == c[1]).
## Used by effsens.R, tmnn.R, tmcal.R, evidence.R, insilico.R.
.html_unescape <- function(s) {
  s <- gsub("&lt;", "<", s, fixed = TRUE); s <- gsub("&gt;", ">", s, fixed = TRUE); s <- gsub("&quot;", "\"", s, fixed = TRUE); s <- gsub("&apos;", "'", s, fixed = TRUE)
  m <- gregexpr("&#[0-9]+;", s); for (x in unique(unlist(regmatches(s, m)))) s <- gsub(x, intToUtf8(as.integer(sub("&#([0-9]+);", "\\1", x))), s, fixed = TRUE)
  gsub("&amp;", "&", s, fixed = TRUE)
}
.allm <- function(p, x) { r <- regmatches(x, gregexpr(p, x, perl = TRUE))[[1]]; if (length(r)) r else character(0) }
.docx_par_text <- function(p) {
  tok <- .allm("<w:t(?: [^>]*)?>[^<]*</w:t>|<w:br(?: [^>]*)?/>|<w:tab/>", p)
  paste(vapply(tok, function(t) if (startsWith(t, "<w:br")) "\n" else if (startsWith(t, "<w:tab")) "\t" else .html_unescape(sub("</w:t>$", "", sub("^<w:t(?: [^>]*)?>", "", t, perl = TRUE))), ""), collapse = "")
}
docx_table0_rows <- function(path) {
  tmp <- tempfile(); dir.create(tmp); unzip(path, files = "word/document.xml", exdir = tmp)
  x <- readChar(file.path(tmp, "word", "document.xml"), file.size(file.path(tmp, "word", "document.xml")), useBytes = TRUE); Encoding(x) <- "UTF-8"
  tbl <- regmatches(x, regexpr("(?s)<w:tbl>.*?</w:tbl>", x, perl = TRUE))
  lapply(.allm("(?s)<w:tr[ >].*?</w:tr>", tbl), function(tr)
    vapply(.allm("(?s)<w:tc>.*?</w:tc>", tr), function(tc) paste(vapply(.allm("(?s)<w:p[ >].*?</w:p>", tc), .docx_par_text, ""), collapse = "\n"), ""))
}
