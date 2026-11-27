## primer_lib.R - helpers of the primer checks, R port of the functions of reanalysis_17genes/insilico.py (get, rc, find, fetch, check) and tmnn.py (tm), plus the Table S1 reader.
## fetch(acc) asks NCBI E-utilities (efetch, FASTA) like the Python script (4 tries, 0.4 s pause) and keeps every sequence in a local cache folder (PRIMER_CACHE, default
## results_R/reanalysis_17genes/ncbi_cache/<accession>.fasta): later runs use the cache and need no network, so the numbers cannot drift with NCBI updates. Delete a file to refetch it.
## The User-Agent carries no e-mail address (the Python script has a placeholder).
suppressMessages(library(curl))
PRIMER_CACHE <- Sys.getenv("PRIMER_CACHE", file.path(OUT, "ncbi_cache"))
.fasta_cache <- new.env()
p_get <- function(u) {
  for (i in 1:4) {
    r <- tryCatch({ h <- new_handle(useragent = "primer-check/1.0 (R port)", timeout = 40); x <- curl_fetch_memory(u, h); if (x$status_code == 200) rawToChar(x$content) else NULL }, error = function(e) NULL)
    if (!is.null(r)) return(r); Sys.sleep(1.5)
  }
  ""
}
p_rc <- function(s) { ch <- rev(strsplit(s, "")[[1]]); m <- c(A = "T", C = "G", G = "C", T = "A", a = "T", c = "G", g = "C", t = "A"); toupper(paste(ifelse(ch %in% names(m), m[ch], ch), collapse = "")) }
p_find <- function(seq, p, maxmm = 2) {              # exact match, else the first window with the fewest mismatches (at most maxmm); pos is 0-based like Python
  p <- toupper(p); L <- nchar(p); n <- nchar(seq); i <- regexpr(p, seq, fixed = TRUE)[1]; if (i > 0) return(list(mm = 0L, pos = i - 1L)); if (n < L) return(list(mm = NA, pos = -1L))
  sc <- utf8ToInt(seq); pc <- utf8ToInt(p); nw <- n - L + 1; mism <- integer(nw); for (k in seq_len(L)) mism <- mism + (sc[k:(k + nw - 1)] != pc[k])
  best <- min(mism); if (best > maxmm) return(list(mm = NA, pos = -1L)); list(mm = best, pos = which(mism == best)[1] - 1L)
}
p_fetch <- function(acc) {
  if (!is.null(.fasta_cache[[acc]])) return(.fasta_cache[[acc]])
  f <- file.path(PRIMER_CACHE, paste0(acc, ".fasta"))
  if (file.exists(f)) fa <- readChar(f, file.size(f), useBytes = TRUE) else {
    fa <- p_get(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=%s", acc)); Sys.sleep(0.4)
    if (nzchar(fa)) { dir.create(PRIMER_CACHE, showWarnings = FALSE, recursive = TRUE); write_text_bytes(fa, f) }
  }
  ln <- strsplit(gsub("\r", "", fa, fixed = TRUE), "\n", fixed = TRUE)[[1]]; .fasta_cache[[acc]] <- list(hdr = if (length(ln)) ln[1] else "", seq = toupper(paste(ln[-1], collapse = "")))
}
p_check <- function(gene, acc, F, R, tablelen = NULL) {
  fz <- p_fetch(acc); hdr <- fz$hdr; seq <- fz$seq; f <- p_find(seq, F); r <- p_find(seq, p_rc(R)); h110 <- substr(hdr, 1, 110)
  if (is.na(f$mm) || is.na(r$mm)) {
    f2 <- p_find(seq, p_rc(F)); r2 <- p_find(seq, R)
    if (!is.na(f2$mm) && !is.na(r2$mm)) return(list(gene = gene, acc = acc, hdr = h110, result = "REVERSED ORIENTATION", F_mm = f2$mm, R_mm = r2$mm))
    return(list(gene = gene, acc = acc, hdr = h110, result = "NO MATCH", F_mm = f$mm, R_mm = r$mm, seqlen = nchar(seq)))
  }
  amp <- r$pos + nchar(R) - f$pos
  list(gene = gene, acc = acc, hdr = h110, result = if (f$mm == 0 && r$mm == 0) "MATCH" else sprintf("MATCH with %d+%d mismatches", f$mm, r$mm), F_mm = f$mm, R_mm = r$mm, amplicon = amp, tablelen = tablelen,
       len_ok = (is.null(tablelen) || amp == tablelen), F_pos = f$pos + 1L, R_pos = r$pos + nchar(R))
}
## nearest-neighbour melting temperature of tmnn.py (SantaLucia-type table, salt correction as in the script)
NN <- list(AA = c(-7.9, -22.2), TT = c(-7.9, -22.2), AT = c(-7.2, -20.4), TA = c(-7.2, -21.3), CA = c(-8.5, -22.7), TG = c(-8.5, -22.7), GT = c(-8.4, -22.4), AC = c(-8.4, -22.4),
           CT = c(-7.8, -21.0), AG = c(-7.8, -21.0), GA = c(-8.2, -22.2), TC = c(-8.2, -22.2), CG = c(-10.6, -27.2), GC = c(-9.8, -24.4), GG = c(-8.0, -19.9), CC = c(-8.0, -19.9))
p_tm <- function(seq, mono = 50e-3, div = 1.5e-3, dntp = 0.6e-3, conc = 50e-9) {
  s <- toupper(seq); ch <- strsplit(s, "")[[1]]; dH <- 0; dS <- 0
  for (k in seq_len(length(ch) - 1)) { v <- NN[[paste0(ch[k], ch[k + 1])]]; dH <- dH + v[1]; dS <- dS + v[2] }
  for (e in c(ch[1], ch[length(ch)])) if (e %in% c("G", "C")) { dH <- dH + 0.1; dS <- dS - 2.8 } else { dH <- dH + 2.3; dS <- dS + 4.1 }
  na <- (mono * 1000 + 120 * sqrt(max(0.0, div * 1000 - dntp * 1000))) / 1000; dS <- dS + 0.368 * (length(ch) - 1) * log(na)
  dH * 1000 / (dS + 1.987 * log(conc / 4)) - 273.15
}
p_pred <- function(a) { gc <- 100 * sum(strsplit(a, "")[[1]] %in% c("G", "C")) / nchar(a); 81.5 + 0.41 * gc - 675 / nchar(a) - 21.6 }   # predicted amplicon melt (0.05 M Na+ term of tmcal.py is 0)
## Table S1 of Additional file 1 (first table, rows after the two header rows; section rows skipped) as a data frame of the cells the scripts use; python-docx cell texts
p_table_rows <- function(path) {
  out <- list()
  for (c in docx_table0_rows(path)[-(1:2)]) { if (length(c) < 7 || c[1] == c[2]) next; out[[length(out) + 1]] <- c }
  out
}
ws_split <- function(s) { x <- strsplit(trimws(s), "[[:space:]]+")[[1]]; x[nzchar(x)] }          # Python str.split()
