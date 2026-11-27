## insilico2.R - R port of reanalysis_17genes/insilico2.py: searches NCBI (esearch + esummary) for the rat RefSeq mRNA records of Hes1, Hey1, Tcf7 and Ryr3 and checks the printed and the new primer pairs
## against every NM_/XM_ candidate. The esearch/esummary answers are cached like the sequences (PRIMER_CACHE), so a later run needs no network. Writes insilico2.json into RA_OUT; prints the candidates and matches.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "docx_table.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE); source(file.path(HERE, "primer_lib.R"))
cached_get <- function(u, key) { f <- file.path(PRIMER_CACHE, paste0(key, ".json")); if (file.exists(f)) return(readChar(f, file.size(f), useBytes = TRUE)); x <- p_get(u); Sys.sleep(0.4); if (nzchar(x)) { dir.create(PRIMER_CACHE, showWarnings = FALSE, recursive = TRUE); write_text_bytes(x, f) }; x }
esearch <- function(term, n = 8) {
  key <- paste0("esearch_", gsub("[^A-Za-z0-9]+", "_", term)); j <- fromJSON(cached_get(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmode=json&retmax=%d&term=%s", n, URLencode(term, reserved = TRUE)), key), simplifyVector = FALSE)
  ids <- unlist(j$esearchresult$idlist); if (!length(ids)) return(list())
  s <- fromJSON(cached_get(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id=", paste(ids, collapse = ",")), paste0("esummary_", gsub("[^A-Za-z0-9]+", "_", term))), simplifyVector = FALSE)$result
  lapply(ids, function(i) list(s[[i]]$caption, substr(s[[i]]$title, 1, 90), s[[i]]$slen, if (is.null(s[[i]]$statusflag)) "" else s[[i]]$statusflag))
}
P <- list(`Hes1_table(canine)` = c("GAGAAGGCGGACATTCTGGA", "ACCTCGTTCATACACTCGCTG"), Hes1_r_NEW = c("CACCGGACAAACCAAAGACAGC", "GAGGTGCTTCACTGTCATTTCC"), Hes1_r_old = c("CTGGAGAAGGCAGACATTCTGG", "GTCACCTCGTTCATGCACTCG"),
          `Hey1_table(canine)` = c("ACCTGAAAATGCTGCACACG", "GCTGGGAGGCGTAGTTGTTA"), Tcf7_table = c("CGTCCTTGATGCTGGGATCTG", "GATGATTGGCTTCTTCGCCTCC"), Tcf7_r_NEW = c("GAGCACACTTCGCAGAGACTTTTC", "GTAGTGCATGAGCAGATTGAAGGC"),
          Ryr3_table = c("CAACCTACTACCCCGAAGGAAC", "GAGTGGGCTCTGCTGTTAGG"))
cands <- list()
for (g in c("Hes1", "Hey1", "Tcf7", "Ryr3")) { r <- esearch(sprintf("%s[Gene Name] AND Rattus norvegicus[Organism] AND (biomol_mrna[PROP]) AND srcdb_refseq[PROP]", g)); cands[[g]] <- r; cat("##", g, "\n"); for (cc in r) cat("    ", paste(vapply(cc, function(x) if (is.null(x)) "None" else as.character(x), ""), collapse = " | "), "\n") }
res <- list(); nn <- function(x) if (is.null(x)) "None" else as.character(x)
for (name in names(P)) {
  g <- strsplit(name, "_", fixed = TRUE)[[1]][1]
  for (cc in cands[[g]]) {
    cap <- cc[[1]]; if (!(startsWith(cap, "NM_") || startsWith(cap, "XM_"))) next
    o <- p_check(name, cap, P[[name]][1], P[[name]][2]); res[[length(res) + 1]] <- list(primer = name, acc = cap, title = cc[[2]], result = o$result, amp = o$amplicon, F_mm = o$F_mm, R_mm = o$R_mm)
    cat(sprintf("%-20s vs %-16s %-28s amp=%s F_mm=%s R_mm=%s\n", name, cap, o$result, nn(o$amplicon), nn(o$F_mm), nn(o$R_mm)))
  }
}
write_text_bytes(py_json_typed(list(cands = cands, res = res), indent = NULL), file.path(OUT, "insilico2.json"))
