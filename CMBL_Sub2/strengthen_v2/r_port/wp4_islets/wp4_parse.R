## wp4_parse.R - R port of wp4_islets/parse_gse81608.py, parse_new_datasets.py and the loader of the earlier GSE84133 analysis (islet_GSE84133/load_beta_counts.py + wp4_core.load_baron): stream the expression files of the public islet datasets
## and keep ONLY the author-labelled beta-cells (of the frozen list frozen_beta_cells.csv, or the type 'beta' of the GEO sample records, or assigned_cluster == 'beta') as sparse matrices (cells x genes) with the number of detected genes and the
## library size over ALL rows. Data are untrusted: only read, never executed. Memory: row chunks (about 190 MB are free on the laptop). Needs wp4_core.R (new_cells). The parsed Cells are cached as <data>/parsed/<tag>.rds.
## Datasets (tags): gse81608 (RPKM, Entrez ids), gse86469 (RSEM expected counts, Ensembl ids), gse81547 (htseq counts per cell, RAW.tar), emtab5061 (rpkm + counts text inside a zip), gse84133_human and gse84133_mouse (inDrop UMI counts).
suppressMessages(library(data.table))
FROZEN <- function() read.csv(file.path(WP4_DIR, "frozen_beta_cells.csv"), stringsAsFactors = FALSE, colClasses = "character")
SPECIAL <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
chunk_columns <- function(path, sep, want_names, label_col = 1L, chunk = 1000L, skip_header_cols = NULL, check_header = NULL) {      # rows of a (gz) text table in chunks; keeps the wanted columns by header name; returns row labels, column totals over all rows, detected counts and the nonzero triplets
  con <- gzfile(path, "rt"); on.exit(close(con)); hdr <- strsplit(readLines(con, n = 1L, warn = FALSE), sep, fixed = TRUE)[[1]]; hdr <- gsub("^\"+|\"+$", "", hdr)
  if (!is.null(check_header)) check_header(hdr); pos <- match(want_names, hdr); if (anyNA(pos)) stop("columns missing from the file: ", paste(head(want_names[is.na(pos)], 5), collapse = ", "))
  nc <- length(pos); tot <- numeric(nc); ndet <- numeric(nc); ri <- list(); ci <- list(); xv <- list(); labs <- list(); nrow_all <- 0L
  repeat {
    ch <- readLines(con, n = chunk, warn = FALSE); if (!length(ch)) break
    sel <- sort(unique(c(label_col, pos))); part <- fread(text = paste(ch, collapse = "\n"), sep = sep, header = FALSE, select = sel, colClasses = list(character = label_col), data.table = FALSE, quote = "\"", fill = TRUE)
    lab <- gsub("^\"+|\"+$", "", part[[match(label_col, sel)]]); m <- as.matrix(part[, match(pos, sel), drop = FALSE]); storage.mode(m) <- "double"; m[is.na(m)] <- 0
    tot <- tot + colSums(m); ndet <- ndet + colSums(m > 0); nz <- which(m != 0, arr.ind = TRUE)
    if (nrow(nz)) { ri[[length(ri) + 1L]] <- nz[, 1] + nrow_all; ci[[length(ci) + 1L]] <- nz[, 2]; xv[[length(xv) + 1L]] <- m[nz] }
    labs[[length(labs) + 1L]] <- lab; nrow_all <- nrow_all + nrow(m)
  }
  list(labels = unlist(labs), tot = tot, ndet = ndet, i = unlist(ri), j = unlist(ci), x = unlist(xv), nrow = nrow_all, ncol = nc)
}
sparse_cells <- function(r) sparseMatrix(i = r$j, j = r$i, x = r$x, dims = c(r$ncol, r$nrow))             # cells x genes

parse_gse81608 <- function(D) {      # RPKM table of Xin et al. 2016 (39851 genes x 1600 cells): the 503 beta-cells (type 'beta' in the GEO sample records); depth covariate = detected genes over ALL rows
  con <- gzfile(file.path(D, "series_matrix.txt.gz"), "rt"); L <- readLines(con, warn = FALSE, encoding = "UTF-8"); close(con)
  fields <- function(l) gsub("^\"+|\"+$", "", strsplit(l, "\t", fixed = TRUE)[[1]][-1]); ch <- lapply(L[startsWith(L, "!Sample_characteristics_ch1")], fields); title <- fields(L[startsWith(L, "!Sample_title")][1])
  col <- function(i) sub("^.*?: ", "", ch[[i + 1]], perl = TRUE); sample <- paste0("Sample_", vapply(strsplit(trimws(title), "[[:space:]]+"), function(w) w[length(w)], ""))
  meta <- data.frame(sample = sample, donor = col(1), condition = col(2), type = col(6), stringsAsFactors = FALSE); b <- meta[meta$type == "beta", ]
  r <- chunk_columns(file.path(D, "human_islets_rpkm.txt.gz"), "\t", b$sample, check_header = function(h) stopifnot("column order differs from the metadata" = identical(h[-1], meta$sample)))
  cat(sprintf("GSE81608 beta cells %d genes %d nonzero %d | donors %d | T2D cells %d | median detected genes %d\n", r$ncol, r$nrow, length(r$x), length(unique(b$donor)), sum(b$condition == "T2D"), as.integer(median(r$ndet))))
  new_cells("GSE81608", "human", "rpkm", "entrez", r$labels, sparse_cells(r), b$donor, b$condition, r$tot, r$ndet)
}
parse_gse86469 <- function(D) {      # RSEM expected counts csv.gz: rows = Ensembl gene ids, columns = cell names; kind 'count', library size = column sum over all rows
  F <- FROZEN(); F <- F[F$dataset == "GSE86469", ]; r <- chunk_columns(file.path(D, "GSE86469_RSEM_expected_counts.csv.gz"), ",", F$cell_key, chunk = 2000L)
  nonens <- r$labels[!startsWith(r$labels, "ENSG")]; stopifnot("rows that are not Ensembl genes" = !length(nonens), anyDuplicated(r$labels) == 0)
  cat(sprintf("GSE86469 beta cells %d rows %d nonzero %d | donors %d | median detected genes %d | median library size %.0f\n", r$ncol, r$nrow, length(r$x), length(unique(F$donor)), as.integer(median(r$ndet)), median(r$tot)))
  new_cells("GSE86469", "human", "count", "ensembl", r$labels, sparse_cells(r), F$donor, F$disease, r$tot, r$ndet)
}
parse_gse81547 <- function(D) {      # RAW.tar of per-cell htseq count tables (gene symbol, count); ERCC and the 5 htseq special rows are dropped before library size / detected genes
  F <- FROZEN(); F <- F[F$dataset == "GSE81547", ]; want <- F$gsm; tarf <- file.path(D, "GSE81547_RAW.tar"); nm <- untar(tarf, list = TRUE); pre <- vapply(strsplit(nm, "_", fixed = TRUE), `[`, "", 1)
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE)); untar(tarf, files = nm[pre %in% want], exdir = tmp, tar = "internal"); cols <- list(); genes0 <- NULL
  for (f in nm[pre %in% want]) {
    g <- strsplit(f, "_", fixed = TRUE)[[1]][1]; d <- fread(file.path(tmp, f), sep = "\t", header = FALSE, colClasses = c("character", "numeric"), data.table = FALSE)
    if (is.null(genes0)) genes0 <- d[[1]]; stopifnot("row order differs" = identical(d[[1]], genes0)); cols[[g]] <- d[[2]]
  }
  stopifnot("frozen GSM missing from the tar" = all(want %in% names(cols))); keep <- !(startsWith(genes0, "ERCC") | genes0 %in% SPECIAL); A <- do.call(cbind, cols[want])[keep, , drop = FALSE]; genes <- genes0[keep]; stopifnot(anyDuplicated(genes) == 0)
  tot <- colSums(A); ndet <- colSums(A > 0); cat(sprintf("GSE81547 beta cells %d rows kept %d nonzero %d | donors %d | median detected genes %d | median library size %.0f\n", ncol(A), nrow(A), sum(A != 0), length(unique(F$donor)), as.integer(median(ndet)), median(tot)))
  new_cells("GSE81547", "human", "count", "symbol", genes, Matrix(t(A), sparse = TRUE), F$donor, F$disease, tot, ndet)
}
zip_line_reader <- function(zf, name, chunk = 4e6) {      # line batches of a zip member read as raw bytes (readLines does not work on unz connections: no seek); a batch ends at a newline, the rest is kept for the next one
  con <- unz(zf, name, "rb"); rest <- raw(0); eof <- FALSE
  list(next_lines = function() {
    repeat {
      if (!eof) { b <- readBin(con, "raw", chunk); if (!length(b)) eof <- TRUE else rest <<- c(rest, b) }
      nl <- which(rest == as.raw(10)); if (length(nl) || eof) break
    }
    if (!length(rest)) return(NULL)
    if (length(nl)) { last <- nl[length(nl)]; txt <- rest[seq_len(last)]; rest <<- if (last < length(rest)) rest[(last + 1L):length(rest)] else raw(0) } else { txt <- rest; rest <<- raw(0) }
    strsplit(sub("\r?\n$", "", rawToChar(txt)), "\r?\n")[[1]]
  }, close = function() close(con))
}
parse_emtab5061 <- function(D, log_dir = NULL) {      # count block of the authors' rpkm + counts text inside processed.1.zip (723 MB text streamed from the zip): header '#samples' + 3514 ids; rows = symbol, RefSeq id, 3514 RPKM, 3514 counts; rows 1-26179 RefSeq genes, then 92 ERCC, then eGFP
  F <- FROZEN(); F <- F[F$dataset == "E-MTAB-5061", ]; NS <- 3514L; NG <- 26179L
  S <- read.delim(file.path(D, "E-MTAB-5061.sdrf.txt"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", colClasses = "character", na.strings = character(0)); lab_of <- setNames(S[["Characteristics [inferred cell type]"]], S[["Source Name"]])
  lab_of[lab_of %in% c("", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null")] <- "nan"            # pandas missing values become 'nan' in the labels
  zf <- file.path(D, "E-MTAB-5061.processed.1.zip"); zl <- unzip(zf, list = TRUE); name <- zl$Name[1]; stopifnot(name == "pancreas_refseq_rpkms_counts_3514sc.txt")
  rd <- zip_line_reader(zf, name); on.exit(rd$close()); first <- rd$next_lines(); ids <- strsplit(first[1], "\t", fixed = TRUE)[[1]][-1]; stopifnot(length(ids) == NS); pend <- first[-1]
  pos <- match(F$cell_key, ids); stopifnot(!anyNA(pos)); cols <- 2L + NS + pos; genes <- character(NG); rows <- matrix(0, NG, length(cols)); ins <- NULL; gcg <- NULL; i <- 0L; done <- FALSE
  while (!done) {
    if (!length(pend)) { pend <- rd$next_lines(); if (is.null(pend)) break }
    for (line in pend) {
      p <- strsplit(line, "\t", fixed = TRUE)[[1]]; stopifnot("row has wrong number of fields" = length(p) == 2L + 2L * NS)
      if (i == NG) { stopifnot("row is not the first ERCC row" = startsWith(p[1], "ERCC")); done <- TRUE; break }
      i <- i + 1L; genes[i] <- p[1]; rows[i, ] <- as.numeric(p[cols]); if (p[1] == "INS") ins <- as.numeric(p[(3L + NS):(2L + 2L * NS)]); if (p[1] == "GCG") gcg <- as.numeric(p[(3L + NS):(2L + 2L * NS)])
    }
    pend <- character(0)
  }
  stopifnot(!is.null(ins), !is.null(gcg), i == NG)
  lab <- ifelse(ids %in% names(lab_of), lab_of[ids], "unlabelled"); chk <- list(); skip <- c("unlabelled", "nan", "co-expression cell", "unclassified endocrine cell", "not applicable")      # identity check (not a selection step)
  for (g in list(list("INS", ins, "beta cell"), list("GCG", gcg, "alpha cell"))) {
    vec <- g[[2]]; own <- g[[3]]; own_mean <- mean(vec[lab == own]); labs <- sort(unique(lab), method = "radix"); others <- list()
    for (l in labs) if (!(l %in% c(own, skip)) && sum(lab == l) >= 10) others[[l]] <- mean(vec[lab == l])
    chk[[g[[1]]]] <- list(cell_type = own, mean_count_own = own_mean, max_mean_count_other = max(unlist(others)), other_means = others, passed = own_mean > 2 * max(unlist(others)))
  }
  if (!is.null(log_dir)) write_json_typed(chk, file.path(log_dir, "identity_check_E-MTAB-5061.json"))
  cat("identity check", paste(vapply(names(chk), function(g) sprintf("%s: %.1f vs %.1f %s", g, chk[[g]]$mean_count_own, chk[[g]]$max_mean_count_other, chk[[g]]$passed), ""), collapse = "; "), "\n")
  stopifnot("identity check failed: the column alignment is wrong" = all(vapply(chk, function(v) v$passed, TRUE)))
  X <- rowsum(rows, group = genes, reorder = FALSE); ndup <- NG - nrow(X); tot <- colSums(X); ndet <- colSums(X > 0)
  cat(sprintf("E-MTAB-5061 beta cells %d rows kept %d (%d repeated symbols summed) nonzero %d | donors %d | median detected genes %d | median library size %.0f\n", ncol(X), nrow(X), ndup, sum(X != 0), length(unique(F$donor)), as.integer(median(ndet)), median(tot)))
  new_cells("E-MTAB-5061", "human", "count", "symbol", rownames(X), Matrix(t(X), sparse = TRUE), F$donor, F$disease, tot, ndet)
}
parse_baron_donor <- function(f) {      # one GSE84133 umifm_counts file: beta-cells as sparse counts, total UMI of every cell, INS / GCG counts for the identity check, cluster labels
  con <- gzfile(f, "rt"); on.exit(close(con)); header <- gsub("\"", "", strsplit(readLines(con, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1]], fixed = TRUE); stopifnot(header[2] == "barcode", header[3] == "assigned_cluster")
  genes <- header[-(1:3)]; gi <- setNames(seq_along(genes), genes); ins_cols <- unname(gi[intersect(c("INS", "Ins1", "Ins2"), names(gi))]); gcg_col <- if ("GCG" %in% names(gi)) gi[["GCG"]] else if ("Gcg" %in% names(gi)) gi[["Gcg"]] else NA
  cl <- character(0); tots <- numeric(0); ins <- numeric(0); gcg <- numeric(0); bi <- list(); bj <- list(); bx <- list(); bbc <- character(0); btot <- numeric(0); nb <- 0L
  repeat {
    ch <- readLines(con, n = 500L, warn = FALSE); if (!length(ch)) break
    for (line in ch) {
      p <- strsplit(line, ",", fixed = TRUE)[[1]]; x <- as.integer(p[-(1:3)]); stopifnot(length(x) == length(genes)); c1 <- gsub("\"", "", p[3], fixed = TRUE); tot <- sum(x)
      cl <- c(cl, c1); tots <- c(tots, tot); ins <- c(ins, if (length(ins_cols)) sum(x[ins_cols]) else 0); gcg <- c(gcg, if (!is.na(gcg_col)) x[gcg_col] else 0)
      if (c1 == "beta") { nb <- nb + 1L; nz <- which(x != 0L); bi[[nb]] <- rep(nb, length(nz)); bj[[nb]] <- nz; bx[[nb]] <- x[nz]; bbc <- c(bbc, gsub("\"", "", p[2], fixed = TRUE)); btot <- c(btot, tot) }
    }
  }
  B <- sparseMatrix(i = unlist(bi), j = unlist(bj), x = as.numeric(unlist(bx)), dims = c(nb, length(genes))); cp10k <- function(v) 1e4 * v / pmax(tots, 1)
  sanity <- lapply(setNames(sort(unique(cl), method = "radix"), sort(unique(cl), method = "radix")), function(cc) list(n = sum(cl == cc), INS_cp10k = mean(cp10k(ins)[cl == cc]), GCG_cp10k = mean(cp10k(gcg)[cl == cc])))
  cat(sprintf("%s cells %d genes %d beta %d | median UMI/cell %.0f\n", basename(f), length(cl), length(genes), nb, median(tots))); list(genes = genes, B = B, beta_total_umi = btot, clusters = sanity)
}
parse_baron <- function(D, species) {      # human1-4 (human4 has type 2 diabetes in the GEO records) or mouse1-2; identity check per donor as in wp4_core.load_baron
  tarf <- file.path(D, "GSE84133_RAW.tar"); nm <- untar(tarf, list = TRUE); donors <- if (species == "human") paste0("human", 1:4) else paste0("mouse", 1:2)
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE)); mats <- list(); genes0 <- NULL; tots <- numeric(0); donor <- character(0)
  for (d in donors) {
    f <- nm[grepl(paste0("_", d, "_"), nm, fixed = TRUE)][1]; untar(tarf, files = f, exdir = tmp); m <- parse_baron_donor(file.path(tmp, f)); cl <- m$clusters
    ins_b <- cl$beta$INS_cp10k; ins_o <- max(vapply(cl[names(cl) != "beta"], function(v) v$INS_cp10k, 1)); gcg_a <- cl$alpha$GCG_cp10k; gcg_o <- max(vapply(cl[names(cl) != "alpha"], function(v) v$GCG_cp10k, 1))
    if (!(ins_b > 2 * ins_o && gcg_a > 2 * gcg_o)) stop(sprintf("identity check failed in %s", d)); if (!is.null(genes0)) stopifnot("gene lists differ between donors" = identical(m$genes, genes0))
    genes0 <- m$genes; mats[[d]] <- m$B; tots <- c(tots, m$beta_total_umi); donor <- c(donor, rep(d, nrow(m$B))); unlink(file.path(tmp, f))
  }
  B <- do.call(rbind, mats); disease <- ifelse(donor == "human4", "T2D", "non-diabetic")
  new_cells(paste0("GSE84133_", species), species, "umi", "symbol", genes0, B, donor, disease, tots, as.numeric(Matrix::rowSums(B > 0)))
}
## load (parse and cache) the Cells of a tag
load_cells <- function(tag, data_dir, log_dir = NULL, refresh = FALSE) {
  cache <- file.path(data_dir, "parsed", paste0(tag, ".rds")); if (file.exists(cache) && !refresh) return(readRDS(cache)); dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
  cells <- switch(tag, gse81608 = parse_gse81608(file.path(data_dir, "GSE81608")), gse86469 = parse_gse86469(file.path(data_dir, "GSE86469")), gse81547 = parse_gse81547(file.path(data_dir, "GSE81547")),
                  emtab5061 = parse_emtab5061(file.path(data_dir, "E-MTAB-5061"), log_dir), gse84133_human = parse_baron(file.path(data_dir, "GSE84133"), "human"), gse84133_mouse = parse_baron(file.path(data_dir, "GSE84133"), "mouse"), stop("unknown tag ", tag))
  saveRDS(cells, cache); cells
}
