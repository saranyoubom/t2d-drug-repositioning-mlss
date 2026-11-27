## wp11_lib.R - shared helpers of the R port of wp11_plate_provenance (plate provenance of the qPCR sheets, 8-9 Nov 2026): paths, a reader for the instrument workbooks (CFX Maestro 'Quantification ... Results.xlsx'), and the small functions
## that several check scripts share (norm, grp, cq_table, read_plate, cq_at, numpy-style arange / std). Source it after defining R_PORT. Read only: the raw exports are never written (the results go to results_R/wp11/).
## Folders (the same as in the Python scripts): P = project folder (CMBL_Sub2_291026), EXD = P/Exported_data (the user's exports of 8 Nov: Before = threshold as found, After = 150), SBA = P/../SBA/Rev2_250626/Export_qPCR (2025 re-exports),
## RT = P/../../Undergrad_Glucose uptake signaling/RT-PCR (26 Aug and 9 Nov exports; WP11_RT overrides), W1 = the matrix rows of WP1 (ins1_sheet_cq.csv, gene_run_map.csv; the archived ones, WP1_DIR overrides).
suppressMessages(library(xml2)); source(file.path(R_PORT, "pycompat.R"))
P <- normalizePath(file.path(R_PORT, "..", "..", "..")); SBA <- normalizePath(file.path(P, "..", "SBA", "Rev2_250626", "Export_qPCR")); EXD <- file.path(P, "Exported_data")
RT <- Sys.getenv("WP11_RT", normalizePath(file.path(P, "..", "..", "Undergrad_Glucose uptake signaling", "RT-PCR"))); W1 <- Sys.getenv("WP1_DIR", file.path(P, "CMBL_Sub2", "strengthen_v2", "wp1_ins1_structure"))
ARC <- Sys.getenv("WP11_ARCHIVE", file.path(P, "CMBL_Sub2", "strengthen_v2", "wp11_plate_provenance")); RES <- Sys.getenv("WP11_OUT", file.path(R_PORT, "results_R", "wp11")); dir.create(RES, showWarnings = FALSE, recursive = TRUE)

## ---- xlsx reader: cached cell values only (openpyxl data_only=True); the zip entries of these exports differ in case from the standard names, so entries are found without regard to case
.xml_text_of <- function(path, entry) { con <- unz(path, entry, "rb"); on.exit(close(con)); chunks <- list(); repeat { b <- readBin(con, "raw", 4e6); if (!length(b)) break; chunks[[length(chunks) + 1L]] <- b }; s <- rawToChar(do.call(c, chunks)); Encoding(s) <- "UTF-8"; s }
.tag_attr <- function(tag, a) { m <- regmatches(tag, regexec(paste0("(?:^|[\\s<])", a, "=\"([^\"]*)\""), tag, perl = TRUE))[[1]]; if (length(m)) m[2] else NA_character_ }
.unxml <- function(s) { s <- gsub("&lt;", "<", s, fixed = TRUE); s <- gsub("&gt;", ">", s, fixed = TRUE); s <- gsub("&quot;", "\"", s, fixed = TRUE); s <- gsub("&apos;", "'", s, fixed = TRUE); gsub("&amp;", "&", s, fixed = TRUE) }
xlsx_load <- function(path) {
  nm <- unzip(path, list = TRUE)$Name; ent <- setNames(nm, tolower(gsub("\\", "/", nm, fixed = TRUE))); get <- function(low) if (is.na(ent[low])) NULL else .xml_text_of(path, ent[[low]])      # the entry names of these files use backslashes ("xl\\workbook.xml") and lower case
  wbx <- get("xl/workbook.xml"); sheets <- regmatches(wbx, gregexpr("<sheet\\b[^>]*>", wbx, perl = TRUE))[[1]]; snames <- .unxml(vapply(sheets, .tag_attr, "", a = "name")); rids <- vapply(sheets, .tag_attr, "", a = "(?:[A-Za-z0-9_]+:)?id")      # the relationship attribute carries an arbitrary namespace prefix here (d3p1:id)
  rel <- get("xl/_rels/workbook.xml.rels"); rt <- regmatches(rel, gregexpr("<Relationship\\b[^>]*>", rel, perl = TRUE))[[1]]; rmap <- setNames(vapply(rt, .tag_attr, "", a = "Target"), vapply(rt, .tag_attr, "", a = "Id"))
  sst <- get("xl/sharedstrings.xml"); shared <- character(0)
  if (!is.null(sst)) { d <- xml_ns_strip(read_xml(sst)); shared <- vapply(xml_find_all(d, "//si"), function(si) paste(xml_text(xml_find_all(si, ".//t[not(ancestor::rPh)]")), collapse = ""), "") }
  list(path = path, get = get, sheetnames = unname(snames), targets = unname(vapply(rids, function(r) { t <- rmap[[r]]; tolower(if (startsWith(t, "/")) substring(t, 2L) else paste0("xl/", t)) }, "")), shared = shared)
}
xlsx_rows <- function(wb, sheet = 1L) {      # all rows of one sheet (name or position) as a list of rows; each row a list of length max_col with NA for an empty cell, numbers as double, text as character
  i <- if (is.character(sheet)) match(sheet, wb$sheetnames) else sheet; if (is.na(i)) stop("no sheet ", sheet); d <- xml_ns_strip(read_xml(wb$get(wb$targets[i])))
  cells <- xml_find_all(d, "//sheetData/row/c"); ref <- xml_attr(cells, "r"); typ <- xml_attr(cells, "t"); row <- as.integer(sub("^[A-Z]+", "", ref))
  col <- vapply(sub("[0-9]+$", "", ref), function(x) { v <- utf8ToInt(x) - 64L; sum(v * 26^rev(seq_along(v) - 1L)) }, 1, USE.NAMES = FALSE); val <- vector("list", length(cells))
  vtxt <- xml_text(xml_find_first(cells, "./v")); itxt <- xml_text(xml_find_first(cells, "./is"))      # one XPath call for all cells (a call per cell takes minutes on the amplification sheets)
  for (k in seq_along(cells)) { t <- typ[k]
    val[[k]] <- if (identical(t, "inlineStr")) (if (is.na(itxt[k])) "" else itxt[k]) else if (is.na(vtxt[k])) NA else if (identical(t, "s")) wb$shared[as.integer(vtxt[k]) + 1L] else if (identical(t, "str")) vtxt[k]
                else if (identical(t, "b")) vtxt[k] == "1" else if (identical(t, "e")) NA else as.numeric(vtxt[k]) }
  nr <- if (length(row)) max(row) else 0L; nc <- if (length(col)) max(col) else 0L; out <- lapply(seq_len(nr), function(r) as.list(rep(NA, nc)))
  for (k in seq_along(cells)) if (!is.null(val[[k]])) out[[row[k]]][[col[k]]] <- val[[k]]
  out
}
.WBC <- new.env(); load_wb <- function(path) { if (is.null(.WBC[[path]])) .WBC[[path]] <- xlsx_load(path); .WBC[[path]] }
.RWC <- new.env(); sheet0 <- function(path) { if (is.null(.RWC[[path]])) { wb <- load_wb(path); .RWC[[path]] <- xlsx_rows(wb, if ("0" %in% wb$sheetnames) "0" else 1L) }; .RWC[[path]] }      # the sheet '0' (or the first sheet) of a workbook, cached

## ---- the functions of check_threshold150.py / reconstruct_threshold.py
norm <- function(s) gsub("[^a-z0-9]", "", tolower(as.character(s)))
grp <- function(s) { s <- if (is.na(s)[1]) "none" else tolower(as.character(s)); if (grepl("ctrl", s, fixed = TRUE) || grepl("control", s, fixed = TRUE)) "CTRL" else if (grepl("dapt", s, fixed = TRUE)) "DAPT" else if (grepl("dkk", s, fixed = TRUE)) "DKK-1" else NA_character_ }
isnum <- function(x) is.numeric(x) && !is.na(x)
blank <- function(s) is.na(s)[1] || identical(s, "")
.target_sample <- function(t, s) { t <- as.character(t); if (grepl("_", t, fixed = TRUE) && blank(s)) { k <- regexpr("_[^_]*$", t); s <- substring(t, k + 1L); t <- substr(t, 1L, k - 1L) }; list(t = t, s = s) }      # '<gene>_<group>' target with an empty sample: split at the last underscore
cq_table <- function(path) {      # table(): list keyed 'gene\rgroup' of list(well, cq) from 'Quantification Cq Results' (sheet '0'), in file order
  rows <- sheet0(path); h <- vapply(rows[[1]], function(c) if (is.na(c)[1]) "None" else as.character(c), ""); wi <- match("Well", h); ti <- match("Target", h); si <- match("Sample", h); ci <- match("Cq", h); out <- list()
  for (r in rows[-1]) { if (is.na(r[[ti]]) || !isnum(r[[ci]])) next
    ts <- .target_sample(r[[ti]], r[[si]]); key <- paste(norm(sub("^[0-9]+[.]\\s*", "", ts$t)), grp(ts$s), sep = "\r"); if (is.null(out[[key]])) out[[key]] <- list(); out[[key]][[length(out[[key]]) + 1L]] <- list(well = as.character(r[[wi]]), cq = as.numeric(r[[ci]])) }
  out
}
"%+%" <- function(a, b) paste0(a, b)
plog_new <- function() { e <- new.env(); e$lines <- character(0); e }      # a logger that keeps the lines for the output file and prints them (named plog_*, not log: base log() must stay usable)
plog <- function(e, s = "") { e$lines <- c(e$lines, s); cat(s, "\n", sep = "") }
key_of <- function(gene, group) paste(gene, group, sep = "\r"); split_key <- function(key) strsplit(key, "\r", fixed = TRUE)[[1]]
cq_at <- function(cyc, y, T, start = 8) {      # first crossing of T by linear interpolation (cycles after `start`)
  n <- length(y); i <- which(y[-n] < T & y[-1] >= T & cyc[-1] > start); if (!length(i)) return(NA_real_); i <- i[1]
  cyc[i] + (T - y[i]) / (y[i + 1] - y[i]) * (cyc[i + 1] - cyc[i])
}
read_plate <- function(base) {      # base = path without the ' -  Quantification ...xlsx' suffix -> list(cyc, curves = named list of baseline-subtracted RFU curves by well, tab = data.frame(well, gene, group, cq))
  rows <- sheet0(paste0(base, " -  Quantification Amplification Results.xlsx")); head <- rows[[1]]; hc <- vapply(head, function(c) if (is.character(c)) c else "", ""); ci <- match("Cycle", hc)
  cyc <- vapply(rows[-1], function(r) if (isnum(r[[ci]])) r[[ci]] else NA_real_, 1); cyc <- cyc[!is.na(cyc)]; curves <- list()
  for (j in seq_along(head)) if (grepl("^[A-H][0-9]+$", hc[j])) { w <- sprintf("%s%02d", substr(hc[j], 1, 1), as.integer(substring(hc[j], 2))); curves[[w]] <- vapply(rows[-1][seq_along(cyc)], function(r) if (isnum(r[[j]])) r[[j]] else NA_real_, 1) }
  rr <- sheet0(paste0(base, " -  Quantification Cq Results.xlsx")); h <- vapply(rr[[1]], function(c) if (is.na(c)[1]) "" else as.character(c), ""); wi <- match("Well", h); ti <- match("Target", h); si <- match("Sample", h); qi <- match("Cq", h); tab <- list()
  for (r in rr[-1]) { if (is.na(r[[ti]]) || !isnum(r[[qi]])) next
    ts <- .target_sample(r[[ti]], r[[si]]); tab[[length(tab) + 1L]] <- data.frame(well = as.character(r[[wi]]), gene = norm(sub("^[0-9]+[.]\\s*", "", ts$t)), group = grp(ts$s), cq = as.numeric(r[[qi]]), stringsAsFactors = FALSE) }
  list(cyc = cyc, curves = curves, tab = do.call(rbind, tab))
}
## numpy.arange(start, stop, step) for floats: length ceil((stop - start) / step); element 1 = start + step, the others start + i * ((start + step) - start)
np_arange <- function(start, stop, step) { n <- ceiling((stop - start) / step); if (n <= 0) return(numeric(0)); if (n == 1) return(start); d <- (start + step) - start; v <- start + (0:(n - 1)) * d; v[2] <- start + step; v }
np_std1 <- function(x) sqrt(np_sum((x - np_mean(x))^2) / (length(x) - 1))      # numpy std(ddof = 1)
read_csv_chr <- function(p) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = character(0), fileEncoding = "UTF-8")      # csv.DictReader: every field a string
write_txt <- function(lines, path) write_text_bytes(gsub("\n", .NL, paste0(paste(lines, collapse = "\n"), "\n"), fixed = TRUE), path)      # text-mode file: lines end with os.linesep