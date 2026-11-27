## metacell.R - R port of wp4_islets/metacell.py and metacell_genes.py: protocol 4.5 (metacell sensitivity), details fixed in PRESPEC_WP4_choices.md section C.
## Per dataset and metacell size (20 primary; 10 and 50 sensitivity): 20 random partitions (partition i uses numpy default_rng(5000 + i)); donors with >= 40 beta-cells; the cells of a donor are shuffled and cut into consecutive blocks of
## exactly `size` cells, the remainder is dropped, a donor needs >= 2 blocks. UMI and counts are summed (then log2(CP10K + 1) or log2(CPM + 1) with the summed library size), RPKM is averaged (log2(mean + 1)). Donor-centred, ranked,
## Spearman, no depth term. Evaluable gene = nonzero in >= 50% of the metacells. Matched null rebuilt on the metacell data (deciles of metacell detection among evaluable genes that are not one of the 17, 2000 draws, seed 2026) with the
## same excess statistic E as the cell-level analysis.
## Output (WP4_MC_OUT, default results_R/wp4/metacell): metacell_runs.csv (one row per dataset, size, partition), metacell_summary.csv (median and range over the partitions, compared with the cell-level 10% analysis of results/index.csv:
## WP4_RES), metacell_gene_evaluability.csv and metacell_pair_evaluability.csv (metacell_genes.py: which genes and pairs are evaluable in the partitions).
## Usage: Rscript metacell.R [tag ...]    (default: all six datasets); Rscript metacell.R genes   (the gene / pair evaluability tables from the partitions; needs metacell_runs.csv)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp4_core.R")); source(file.path(HERE, "wp4_parse.R")); options(py_json_nan = TRUE)
MC_EXACT <- Sys.getenv("MC_EXACT") == "1"
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); RES <- Sys.getenv("WP4_RES", normalizePath(file.path(R_PORT, "..", "wp4_islets", "results"))); OUT <- Sys.getenv("WP4_MC_OUT", file.path(R_PORT, "results_R", "wp4", if (MC_EXACT) "metacell_exact_rule" else "metacell")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
## detection of each gene = fraction of the (meta)cells in which it is nonzero. The Python code takes scipy's mean of a sparse boolean column: 1/n (rounded) added once per detecting cell, sequentially. A gene detected in exactly half of an even
## number of metacells can then come out below 0.5 (0.49999999999999961 for 124 metacells) and is NOT evaluable, although the pre-specified rule is "at least 50%". The default reproduces this (the archived results and the manuscript's
## metacell numbers were made so); MC_EXACT=1 applies the rule exactly (count / n) and writes to results_R/wp4/metacell_exact_rule.
det_cols <- function(S) { n <- nrow(S); k <- as.numeric(Matrix::colSums(S > 0)); if (MC_EXACT) return(k / n); s <- 0; tab <- numeric(n + 1L); for (i in seq_len(n)) { s <- s + 1 / n; tab[i + 1L] <- s }; tab[k + 1] }
TAGS <- c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061", "gse84133_mouse"); SIZES <- c(20L, 10L, 50L); NPART <- 20L; MIN_CELLS <- 40L; DET <- 0.50
RUN_COLS <- c("metacells", "donors_used", "eligible_genes", "genes_evaluable", "consensus_pairs_evaluable", "median_rho", "null_median_rho", "mean_z", "mean_z_null", "E", "SE_E", "tag", "size", "partition")
metacell_blocks <- function(cells, size, i) {      # the partition of partition i: blocks of cell positions (1-based) and their donors
  rng <- np_rng(5000L + i); blocks <- list(); bdon <- character(0)
  for (d in unique(cells$donor)) {
    idx <- which(cells$donor == d); if (length(idx) < MIN_CELLS) next; nb <- length(idx) %/% size; if (nb < 2) next
    perm <- np_permute(rng, idx - 1) + 1; for (b in seq_len(nb)) { blocks[[length(blocks) + 1]] <- perm[((b - 1) * size + 1):(b * size)]; bdon <- c(bdon, d) }
  }
  list(blocks = blocks, bdon = bdon)
}
block_matrix <- function(bl, size, n) sparseMatrix(i = rep(seq_along(bl$blocks), each = size), j = unlist(bl$blocks), x = 1, dims = c(length(bl$blocks), n))
one_run <- function(cells, size, i) {
  bl <- metacell_blocks(cells, size, i); nm <- length(bl$blocks); if (nm < 4) return(list(metacells = nm, donors_used = length(unique(bl$bdon))))
  M <- block_matrix(bl, size, nrow(cells$B)); S <- M %*% cells$B; tot <- as.numeric(M %*% cells$tot); det <- det_cols(S); elig <- which(det >= DET); D <- as.matrix(S[, elig, drop = FALSE])
  X <- if (cells$kind == "umi") log2(D / tot * 1e4 + 1) else if (cells$kind == "count") log2(D / tot * 1e6 + 1) else log2(D / size + 1)
  Xc <- center_by_donor(X, bl$bdon, "count"); Rk <- ranks_of(Xc); sym <- function(g) cells_label(cells, g); labels_e <- cells$genes[elig]; idx_e <- setNames(seq_along(labels_e), labels_e); det_e <- det[elig]
  ct <- consensus_tables(); own <- vapply(GENES, sym, ""); mp <- make_pool(det_e, labels_e, own); ev_genes <- GENES[vapply(GENES, function(g) sym(g) %in% labels_e, TRUE)]
  pairs <- list(); for (i1 in seq_along(GENES)) for (j1 in seq_along(GENES)) if (j1 > i1) { a <- GENES[i1]; cc <- GENES[j1]; if (ct$key(a, cc) %in% ct$consset && sym(a) %in% labels_e && sym(cc) %in% labels_e) pairs[[length(pairs) + 1]] <- c(a, cc) }
  out <- list(metacells = nm, donors_used = length(unique(bl$bdon)), eligible_genes = length(elig), genes_evaluable = length(ev_genes), consensus_pairs_evaluable = length(pairs)); if (!length(pairs)) return(out)
  pidx <- lapply(pairs, function(p) c(idx_e[[sym(p[1])]], idx_e[[sym(p[2])]])); rho <- vapply(pidx, function(p) sum(Rk[p[1], ] * Rk[p[2], ]), 1); z <- atanh(pmin(pmax(rho, -CLIP), CLIP))
  dr <- draws_main(mp$pool, mp$bins, pidx, NN, SEED_NULL); V <- null_values(Rk, dr$IA, dr$IC); zm <- rowMeans(atanh(pmin(pmax(V, -CLIP), CLIP)))
  c(out, list(median_rho = median(rho), null_median_rho = mean(apply(V, 1, median)), mean_z = mean(z), mean_z_null = mean(zm), E = mean(z) - mean(zm), SE_E = sd(zm)))
}
main <- function(tags) {
  rows <- list()
  for (tag in tags) {
    cells <- load_cells(tag, DATA); t0 <- Sys.time()
    for (size in SIZES) {
      first <- length(rows) + 1L; for (i in 0:(NPART - 1L)) { r <- one_run(cells, size, i); r$tag <- tag; r$size <- size; r$partition <- i; rows[[length(rows) + 1]] <- r }
      sub <- rows[first:length(rows)]; g <- function(k) { v <- vapply(sub, function(r) if (is.null(r[[k]])) NA_real_ else as.numeric(r[[k]]), 1); v }
      cat(sprintf("%-15s size %2d | metacells %s donors %s | evaluable genes %s pairs %s | median rho %s | E %s\n", tag, size, format(median(g("metacells"))), format(median(g("donors_used"))), format(median(g("genes_evaluable"), na.rm = TRUE)), format(median(g("consensus_pairs_evaluable"), na.rm = TRUE)),
                  if (any(!is.na(g("median_rho")))) sprintf("%.3f", median(g("median_rho"), na.rm = TRUE)) else "-", if (any(!is.na(g("E")))) sprintf("%.3f", median(g("E"), na.rm = TRUE)) else "-"))
    }
    cat(sprintf("%s done in %.0f s\n", tag, as.numeric(difftime(Sys.time(), t0, units = "secs")))); rm(cells); gc(FALSE)
    R <- do.call(rbind, lapply(rows, function(r) { x <- lapply(RUN_COLS, function(k) if (is.null(r[[k]])) NA else r[[k]]); names(x) <- RUN_COLS; as.data.frame(x, stringsAsFactors = FALSE) })); py_write_csv_repr(R, file.path(OUT, "metacell_runs.csv"), eol = "\r\n")
  }
  idx <- read.csv(file.path(RES, "index.csv"), stringsAsFactors = FALSE, check.names = FALSE); cell <- idx[idx$variant == "primary" & idx$mode == "unadjusted", ]; S <- list()
  for (tag in unique(R$tag)) for (size in unique(R$size[R$tag == tag])) {
    g <- R[R$tag == tag & R$size == size, ]; f3 <- function(col) if (col %in% names(g) && any(!is.na(g[[col]]))) c(median(g[[col]], na.rm = TRUE), min(g[[col]], na.rm = TRUE), max(g[[col]], na.rm = TRUE)) else c(NA_real_, NA_real_, NA_real_)
    row <- list(tag = tag, size = size, partitions = nrow(g), partitions_with_pairs = sum(!is.na(g$consensus_pairs_evaluable) & g$consensus_pairs_evaluable > 0), donors_used_median = median(g$donors_used), metacells_median = median(g$metacells))
    for (cc in c("genes_evaluable", "consensus_pairs_evaluable", "median_rho", "null_median_rho", "E", "SE_E")) { v <- f3(cc); row[[paste0(cc, "_med")]] <- v[1]; row[[paste0(cc, "_min")]] <- v[2]; row[[paste0(cc, "_max")]] <- v[3] }
    ci <- match(tag, cell$tag); row$cell_level_pairs <- if (!is.na(ci)) as.integer(cell$consensus_pairs_evaluable[ci]) else NA_integer_; row$cell_level_median_rho <- if (!is.na(ci)) cell$median_rho_consensus[ci] else NA_real_
    row$gain_pairs_median <- row$consensus_pairs_evaluable_med - row$cell_level_pairs; S[[length(S) + 1]] <- row
  }
  Sd <- do.call(rbind, lapply(S, as.data.frame, stringsAsFactors = FALSE)); py_write_csv_repr(Sd, file.path(OUT, "metacell_summary.csv"), eol = "\r\n"); print(Sd)
}
evaluable_genes <- function(cells, size, i) {
  bl <- metacell_blocks(cells, size, i); nm <- length(bl$blocks); if (nm < 4) return(list(ev = NULL, nm = nm))
  M <- block_matrix(bl, size, nrow(cells$B)); S <- M %*% cells$B; det <- det_cols(S)
  list(ev = vapply(GENES, function(g) { l <- cells_label(cells, g); l %in% cells$genes && det[cells$gi[[l]]] >= DET }, TRUE), nm = nm)
}
genes_main <- function() {
  runs <- read.csv(file.path(OUT, "metacell_runs.csv"), stringsAsFactors = FALSE); ct <- consensus_tables(); allpairs <- list(); for (i in seq_along(GENES)) for (j in seq_along(GENES)) if (j > i && ct$key(GENES[i], GENES[j]) %in% ct$consset) allpairs[[length(allpairs) + 1]] <- c(GENES[i], GENES[j])
  rows <- list(); prow <- list(); bad <- 0L
  for (tag in TAGS) {
    cells <- load_cells(tag, DATA); det0 <- det_cols(cells$B); cell_ev <- vapply(GENES, function(g) { l <- cells_label(cells, g); l %in% cells$genes && det0[cells$gi[[l]]] >= 0.10 }, TRUE)
    for (size in SIZES) {
      cnt <- setNames(integer(length(GENES)), GENES); pcnt <- integer(length(allpairs)); n_ok <- 0L
      for (i in 0:(NPART - 1L)) {
        r <- evaluable_genes(cells, size, i); if (is.null(r$ev)) next; n_ok <- n_ok + 1L; cnt <- cnt + r$ev; for (k in seq_along(allpairs)) pcnt[k] <- pcnt[k] + (r$ev[[allpairs[[k]][1]]] && r$ev[[allpairs[[k]][2]]])
        rr <- runs[runs$tag == tag & runs$size == size & runs$partition == i, ][1, ]; if (as.integer(rr$genes_evaluable) != sum(r$ev) || as.integer(rr$metacells) != r$nm) bad <- bad + 1L
      }
      for (g in GENES) rows[[length(rows) + 1]] <- list(tag = tag, size = size, gene = g, partitions_evaluable = as.integer(cnt[[g]]), partitions_analysed = n_ok, cell_level_evaluable = unname(cell_ev[[g]]))
      for (k in seq_along(allpairs)) prow[[length(prow) + 1]] <- list(tag = tag, size = size, gene1 = allpairs[[k]][1], gene2 = allpairs[[k]][2], partitions_evaluable = as.integer(pcnt[k]), partitions_analysed = n_ok)
    }
    cat(tag, "done\n"); rm(cells); gc(FALSE)
  }
  py_write_csv_repr(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)), file.path(OUT, "metacell_gene_evaluability.csv"), eol = "\r\n")
  py_write_csv_repr(do.call(rbind, lapply(prow, as.data.frame, stringsAsFactors = FALSE)), file.path(OUT, "metacell_pair_evaluability.csv"), eol = "\r\n"); cat("runs where the recomputed counts differ from metacell_runs.csv:", bad, "\n")
}
a <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L) { if (identical(a, "genes")) genes_main() else main(if (length(a)) a else TAGS) }
