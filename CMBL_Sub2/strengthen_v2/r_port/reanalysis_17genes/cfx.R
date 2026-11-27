## cfx.R - R port of reanalysis_17genes/cfx.py: reader for the Bio-Rad CFX Maestro exports (SBA/Rev2_250626/Export_qPCR/<run> -  <kind>.xlsx).
## Those workbooks have lower-case zip entry names with backslashes ("xl\sharedstrings.xml"), which openpyxl rejects (cfx.py repairs the names) and readxl cannot open; the sheets are read here with xml2.
## cfx_sheet_rows(run, kind) returns a named list: sheet title -> list of rows (row 1 first), each row a list of cells (NULL for empty), padded to the widest cell present, like openpyxl iter_rows(values_only=True).
suppressMessages(library(xml2))
CFX_EXP <- function() {
  cand <- c(Sys.getenv("CFX_EXPORT", ""), file.path(R_PORT, "..", "..", "..", "..", "SBA", "Rev2_250626", "Export_qPCR")); cand <- cand[nzchar(cand) & dir.exists(cand)]
  if (!length(cand)) stop("CFX export folder not found; set CFX_EXPORT"); normalizePath(cand[1])
}
.zip_read <- function(path, ent, name_lower) {           # bytes of the entry whose normalized (lower-case, '/') name is name_lower
  i <- which(tolower(gsub("\\", "/", ent$Name, fixed = TRUE)) == name_lower); if (!length(i)) return(NULL)
  con <- unz(path, ent$Name[i[1]], "rb"); on.exit(close(con)); readBin(con, "raw", ent$Length[i[1]])
}
.col_index <- function(ref) { l <- gsub("[0-9]", "", ref); vapply(strsplit(l, ""), function(ch) { v <- match(ch, LETTERS); sum(v * 26^(rev(seq_along(v)) - 1)) }, 1) }
cfx_load <- function(path) {
  ent <- unzip(path, list = TRUE)
  ss_raw <- .zip_read(path, ent, "xl/sharedstrings.xml"); sst <- character(0)
  if (!is.null(ss_raw)) { d <- xml_ns_strip(read_xml(ss_raw)); sst <- vapply(xml_find_all(d, "//si"), function(si) paste(xml_text(xml_find_all(si, ".//t")), collapse = ""), "") }
  wb <- rawToChar(.zip_read(path, ent, "xl/workbook.xml")); rels <- rawToChar(.zip_read(path, ent, "xl/_rels/workbook.xml.rels"))
  sheets <- regmatches(wb, gregexpr("<sheet [^>]*/>", wb))[[1]]
  attr <- function(s, a) { m <- regmatches(s, regexpr(sprintf("%s=\"[^\"]*\"", a), s)); sub(sprintf("^%s=\"", a), "", sub("\"$", "", m)) }
  rel <- regmatches(rels, gregexpr("<Relationship [^>]*/>", rels))[[1]]; rel_target <- setNames(vapply(rel, function(s) attr(s, "Target"), ""), vapply(rel, function(s) attr(s, "Id"), ""))
  out <- list()
  for (s in sheets) {
    title <- attr(s, "name"); title <- gsub("&amp;", "&", gsub("&lt;", "<", gsub("&gt;", ">", gsub("&quot;", "\"", title, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)
    rid <- regmatches(s, regexpr("[A-Za-z0-9]+:id=\"[^\"]*\"", s)); rid <- sub("^[^=]*=\"", "", sub("\"$", "", rid))      # the namespace prefix of the relationship id differs between exporters (r:id, d3p1:id)
    tgt <- tolower(gsub("\\", "/", rel_target[[rid]], fixed = TRUE)); tgt <- sub("^/", "", tgt); if (!startsWith(tgt, "xl/")) tgt <- paste0("xl/", tgt)
    d <- xml_ns_strip(read_xml(.zip_read(path, ent, tgt))); cells <- xml_find_all(d, "//sheetData/row/c")
    if (!length(cells)) { out[[title]] <- list(); next }
    ref <- xml_attr(cells, "r"); tt <- xml_attr(cells, "t"); v <- xml_text(xml_find_first(cells, "v")); isv <- xml_text(xml_find_first(cells, "is"))
    rown <- as.integer(gsub("[A-Z]", "", ref)); coln <- .col_index(ref); nr <- max(rown); nc <- max(coln)
    mat <- replicate(nr, replicate(nc, NULL, simplify = FALSE), simplify = FALSE)
    for (k in seq_along(cells)) {
      t <- tt[k]; val <- NULL
      if (!is.na(t) && t == "s") { if (!is.na(v[k])) val <- sst[as.integer(v[k]) + 1L] }
      else if (!is.na(t) && t == "inlineStr") { if (!is.na(isv[k])) val <- isv[k] }
      else if (!is.na(t) && t == "b") { if (!is.na(v[k])) val <- v[k] == "1" }
      else if (!is.na(t) && (t == "str" || t == "e")) { if (!is.na(v[k])) val <- v[k] }
      else if (!is.na(v[k]) && nzchar(v[k])) val <- as.numeric(v[k])
      if (!is.null(val)) mat[[rown[k]]][[coln[k]]] <- val
    }
    out[[title]] <- mat
  }
  out
}
cfx_runs <- function() py_sort(unique(sub(" -  .*", "", basename(list.files(CFX_EXP(), pattern = "\\.xlsx$", full.names = TRUE)))))
cfx_sheet_rows <- function(run, kind) cfx_load(file.path(CFX_EXP(), paste0(run, " -  ", kind, ".xlsx")))
