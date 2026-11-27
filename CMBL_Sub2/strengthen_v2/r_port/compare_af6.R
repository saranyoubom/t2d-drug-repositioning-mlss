## compare_af6.R - cell by cell comparison of the R build of Additional file 6 (results_R/af6/Additional file 6.xlsx, or AF6_OUT) with the submitted workbook (wp6a_assets/Additional file 6_JTM_working_201126b.xlsx, or AF6_REF):
## sheet names and order, the cells of every sheet (reference, type, value, bold or not), column widths and frozen rows. Numbers are compared as written (text of the <v> element); numbers that differ as text are also compared as values (relative difference).
## Usage: Rscript compare_af6.R [max lines per sheet, default 6]
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); suppressMessages(library(xml2))
REF <- Sys.getenv("AF6_REF", file.path(R_PORT, "..", "wp6a_assets", "Additional file 6_JTM_working_201126b.xlsx")); NEW <- Sys.getenv("AF6_OUT", file.path(R_PORT, "results_R", "af6", "Additional file 6.xlsx"))
a <- commandArgs(trailingOnly = TRUE); MAXL <- if (length(a)) as.integer(a[1]) else 6L
read_book <- function(path) {
  ent <- function(name) { con <- unz(path, name, "rb"); on.exit(close(con)); ch <- list(); repeat { b <- readBin(con, "raw", 4e6); if (!length(b)) break; ch[[length(ch) + 1L]] <- b }; s <- rawToChar(do.call(c, ch)); Encoding(s) <- "UTF-8"; s }
  wb <- read_xml(ent("xl/workbook.xml")); wb <- xml_ns_strip(wb); sh <- xml_find_all(wb, "//sheet"); names <- xml_attr(sh, "name")
  st <- xml_ns_strip(read_xml(ent("xl/styles.xml"))); xf <- xml_find_all(st, "//cellXfs/xf"); fonts <- xml_find_all(st, "//fonts/font"); bold_font <- vapply(fonts, function(f) length(xml_find_all(f, "./b")) > 0, TRUE); font_name <- vapply(fonts, function(f) { v <- xml_attr(xml_find_first(f, "./name"), "val"); paste(v, xml_attr(xml_find_first(f, "./sz"), "val")) }, "")
  sty <- list(bold = bold_font[as.integer(xml_attr(xf, "fontId")) + 1L], font = font_name[as.integer(xml_attr(xf, "fontId")) + 1L])
  sheets <- lapply(seq_along(names), function(i) { d <- xml_ns_strip(read_xml(ent(sprintf("xl/worksheets/sheet%d.xml", i)))); c <- xml_find_all(d, "//sheetData/row/c"); s <- as.integer(xml_attr(c, "s")); s[is.na(s)] <- 0L; t <- xml_attr(c, "t"); t[is.na(t)] <- "n"
    v <- ifelse(t == "inlineStr", xml_text(xml_find_first(c, "./is")), xml_text(xml_find_first(c, "./v")))
    cl <- xml_find_all(d, "//cols/col"); pane <- xml_find_first(d, "//pane")
    list(cells = data.frame(ref = xml_attr(c, "r"), t = t, v = v, bold = sty$bold[s + 1L], font = sty$font[s + 1L], stringsAsFactors = FALSE), widths = setNames(as.numeric(xml_attr(cl, "width")), xml_attr(cl, "min")), ysplit = xml_attr(pane, "ySplit")) })
  list(names = names, sheets = setNames(sheets, names))
}
A <- read_book(REF); B <- read_book(NEW)
cat(sprintf("reference: %s (%d sheets)\nR build:   %s (%d sheets)\n", basename(REF), length(A$names), basename(NEW), length(B$names)))
cat(if (identical(A$names, B$names)) "SAME   sheet names and order\n" else sprintf("DIFFER sheet names or order; only in reference: %s; only in R: %s\n", paste(setdiff(A$names, B$names), collapse = ", "), paste(setdiff(B$names, A$names), collapse = ", ")))
tot <- 0L; bad_sheets <- 0L; same_sheets <- 0L; numdiff_max <- 0
for (nm in intersect(A$names, B$names)) {
  x <- A$sheets[[nm]]; y <- B$sheets[[nm]]; mx <- merge(x$cells, y$cells, by = "ref", all = TRUE, suffixes = c(".ref", ".R")); tot <- tot + nrow(x$cells)
  miss <- is.na(mx$t.R); extra <- is.na(mx$t.ref); both <- !miss & !extra
  dt <- both & mx$t.ref != mx$t.R; dv <- both & !dt & (is.na(mx$v.ref) != is.na(mx$v.R) | (!is.na(mx$v.ref) & !is.na(mx$v.R) & mx$v.ref != mx$v.R)); db <- both & (mx$bold.ref != mx$bold.R | mx$font.ref != mx$font.R)
  numd <- dv & mx$t.ref == "n"; if (any(numd)) { u <- suppressWarnings(as.numeric(mx$v.ref[numd])); w <- suppressWarnings(as.numeric(mx$v.R[numd])); r <- abs(u - w) / pmax(1, abs(u)); numdiff_max <- max(numdiff_max, r[is.finite(r)], 0) }
  wd <- !identical(x$widths[order(as.integer(names(x$widths)))], y$widths[order(as.integer(names(y$widths)))]); fz <- !identical(x$ysplit, y$ysplit)
  probs <- c(if (any(miss)) sprintf("%d cells missing in R", sum(miss)), if (any(extra)) sprintf("%d extra cells in R", sum(extra)), if (any(dt)) sprintf("%d types differ", sum(dt)), if (any(dv)) sprintf("%d values differ", sum(dv)), if (any(db)) sprintf("%d fonts differ", sum(db)), if (wd) "column widths differ", if (fz) "frozen rows differ")
  if (length(probs)) {
    bad_sheets <- bad_sheets + 1L; cat(sprintf("DIFFER %-32s %s\n", nm, paste(probs, collapse = "; ")))
    k <- which(miss | extra | dt | dv | db); for (i in utils::head(k, MAXL)) cat(sprintf("         %-7s ref [%s %s] R [%s %s]\n", mx$ref[i], mx$t.ref[i], substr(mx$v.ref[i], 1, 70), mx$t.R[i], substr(mx$v.R[i], 1, 70)))
    if (wd) { ks <- union(names(x$widths), names(y$widths)); dk <- ks[vapply(ks, function(j) !identical(x$widths[j], y$widths[j]), TRUE)]; cat("         widths (column: ref, R):", paste(utils::head(sprintf("%s: %s, %s", dk, x$widths[dk], y$widths[dk]), MAXL), collapse = " | "), "\n") }
  } else same_sheets <- same_sheets + 1L
}
cat(sprintf("%d of %d sheets identical (cells, types, values, fonts, widths, frozen rows); %d differ; %d reference cells; largest relative difference among numbers written differently: %.2g\n", same_sheets, length(intersect(A$names, B$names)), bad_sheets, tot, numdiff_max))
cat(sprintf("%-6s %-42s %d sheets, %d values compared, largest relative difference %.2g\n", if (!bad_sheets && identical(A$names, B$names)) "AGREE" else "DIFFER", "Additional file 6 workbook", length(A$names), tot, numdiff_max))      # the line agreement_report.R reads
