## wp4_core.R - R port of wp4_islets/wp4_core.py: shared loaders and statistics of WP4 (public islet single-cell cross-check, version 2; protocol v2 section 4). No side effects at source time; define R_PORT (the r_port folder) first.
## Same statistics as the earlier scripts (log-normalised expression of beta-cells, donor-centring, rank correlation, genes evaluable when detected in >= 10% of beta-cells, expression-matched null of random gene pairs from the same
## detection deciles, 2000 draws seed 2026, pair-level seed 2027), with the additions of protocol 4.2 / 4.4: depth adjustment (partial Spearman: donor-centred expression regressed on donor-centred depth and depth^2, residuals
## ranked; df = cells - donors - 4), excess E = mean Fisher z of the evaluable consensus pairs minus the mean over the 2000 matched random pair sets, null rebuilt on the adjusted ranks with the SAME random draws.
## The random draws replicate the loop of the Python call for call (numpy PCG64 through rng/numpy_rng.R); the legacy numeric type of UMI data (float32: log2(CP10K + 1), donor centring with sequential float32 sums) is emulated
## by rounding to single precision after every operation (f32()); the C library's log2f is replaced by the correctly rounded float32 of log2().
suppressMessages({ library(Matrix); library(jsonlite) })
if (!exists("R_PORT")) stop("define R_PORT (the r_port folder) before sourcing wp4_core.R")
source(file.path(R_PORT, "pycompat.R")); source(file.path(R_PORT, "rng", "numpy_rng.R"))
WP4_DIR <- Sys.getenv("WP4_DIR", normalizePath(file.path(R_PORT, "..", "wp4_islets")))             # the folder with the data files that ship with the repository (ensembl_ids_17genes.csv, frozen_*.csv)
R17 <- Sys.getenv("R17_DIR", normalizePath(file.path(R_PORT, "..", "..", "R_17genes")))             # the authors' R files: Fig_4d_consensus_pairs.csv, Fig_6b_*.csv
NN <- 2000L; SEED_NULL <- 2026L; SEED_PAIR <- 2027L
GENES <- c("Hes1", "Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")      # order as in the earlier scripts
ENTREZ <- c(Hes1 = 3280, Wnt5a = 7474, Wnt5b = 81029, Wnt9a = 7483, Lef1 = 51176, Tcf7 = 6932, Tcf7l2 = 6934, Itpr1 = 3708, Ryr2 = 6262, Ryr3 = 6263, Cacna1c = 775, Cacna1d = 776, Kcnj11 = 3767, Ptbp1 = 5725, Glut2 = 6514, Ins1 = 3630, Glp1r = 2740)

f32 <- function(x) { d <- dim(x); r <- readBin(writeBin(as.numeric(x), raw(), size = 4L), "numeric", size = 4L, n = length(x)); if (!is.null(d)) dim(r) <- d; r }      # round to single precision
symbol_of <- function(species, g) if (species == "human") { if (g == "Glut2") "SLC2A2" else if (g == "Ins1") "INS" else toupper(g) } else { if (g == "Glut2") "Slc2a2" else g }
.ENS <- new.env()
ensembl_map <- function() { if (is.null(.ENS$m)) { d <- read.csv(file.path(WP4_DIR, "ensembl_ids_17genes.csv"), stringsAsFactors = FALSE); .ENS$m <- setNames(d$ensembl, d$gene) }; .ENS$m }
consensus_tables <- function() {      # the 21 consensus pairs (order of Fig_4d_consensus_pairs.csv), set keys, the in vitro edge r (Hes1-normalised) and the centrality degree
  cp <- read.csv(file.path(R17, "Fig_4d_consensus_pairs.csv"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)[[1]]; cons <- strsplit(cp, "_", fixed = TRUE)
  key <- function(a, b) paste(sort(c(a, b), method = "radix"), collapse = "|"); e <- read.csv(file.path(R17, "Fig_6b_Edge_Correlation_Statistics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  redge <- setNames(e[["Correlation (r)"]], mapply(key, e$Gene1, e$Gene2)); cen <- read.csv(file.path(R17, "Fig_6b_Gene_Centrality_Summary.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  list(cons = cons, consset = vapply(cons, function(p) key(p[1], p[2]), ""), redge = redge, cen = setNames(cen$Degree, cen$Gene), key = key)
}
## beta-cells of one dataset: B sparse (cells x genes) of raw counts (kind umi / count) or RPKM; id_type 'entrez', 'ensembl' or 'symbol'; gene labels are character
new_cells <- function(name, species, kind, id_type, genes, B, donor, disease, tot, ndet) {
  genes <- as.character(genes); B <- as(B, "CsparseMatrix"); stopifnot(nrow(B) == length(donor), ncol(B) == length(genes), length(tot) == length(donor), length(ndet) == length(donor), length(disease) == length(donor))
  list(name = name, species = species, kind = kind, id_type = id_type, genes = genes, gi = setNames(seq_along(genes), genes), B = B, donor = as.character(donor), disease = as.character(disease), tot = as.numeric(tot), ndet = as.numeric(ndet))
}
cells_n <- function(cells) nrow(cells$B); cells_donors <- function(cells) length(unique(cells$donor))
cells_depth <- function(cells) if (cells$kind %in% c("umi", "count")) log10(cells$tot) else log10(cells$ndet)        # protocol 4.2: log10 library size (UMI / read counts), log10 detected genes (RPKM)
cells_label <- function(cells, g) if (cells$id_type == "entrez") as.character(ENTREZ[[g]]) else if (cells$id_type == "ensembl") ensembl_map()[[g]] else symbol_of(cells$species, g)
cells_subset <- function(cells, mask, name = NULL) { m <- as.logical(mask); new_cells(if (is.null(name)) cells$name else name, cells$species, cells$kind, cells$id_type, cells$genes, cells$B[m, , drop = FALSE], cells$donor[m], cells$disease[m], cells$tot[m], cells$ndet[m]) }

## ------------------------------------------------------------------------------------------------------------------------------------------ expression frame
lognorm <- function(cells, elig) {      # dense log2-normalised expression of the eligible genes (cells x genes): UMI -> log2(CP10K + 1) in float32; read counts -> log2(CPM + 1); RPKM -> log2(RPKM + 1)
  D <- as.matrix(cells$B[, elig, drop = FALSE])
  if (cells$kind == "umi") { x <- f32(D / cells$tot); x <- f32(x * 1e4); x <- f32(x + 1); return(f32(log2(x))) }
  if (cells$kind == "count") return(log2(D / cells$tot * 1e6 + 1)); log2(D + 1)
}
center_by_donor <- function(X, donor, kind) {
  for (d in unique(donor)) {
    ix <- which(donor == d)
    if (kind == "umi") { s <- X[ix[1], ]; for (i in ix[-1]) s <- f32(s + X[i, ]); mu <- f32(s / length(ix)); X[ix, ] <- f32(sweep(X[ix, , drop = FALSE], 2, mu, "-")) }      # numpy float32 mean over axis 0: sequential float32 sum
    else X[ix, ] <- sweep(X[ix, , drop = FALSE], 2, colMeans(X[ix, , drop = FALSE]), "-")
  }
  X
}
rank_std <- function(part) {      # cells x genes chunk -> genes x cells matrix of mean-centred, unit-norm ranks (average ranks for ties); the dot product of two rows is the Spearman rho (the 1e-12 of the earlier scripts is kept)
  R <- t(apply(part, 2, rank)); R <- R - rowMeans(R); R / (sqrt(rowSums(R^2)) + 1e-12)
}
ranks_of <- function(Xc, chunk = 1500) { out <- matrix(0, ncol(Xc), nrow(Xc)); for (s in seq(1, ncol(Xc), by = chunk)) { e <- min(s + chunk - 1, ncol(Xc)); out[s:e, ] <- rank_std(Xc[, s:e, drop = FALSE]) }; out }
depth_design <- function(depth, donor) {      # Z (cells x 2): donor-centred depth and donor-centred squared depth (unit variance), the donor-fixed-effects model with common depth slopes
  d <- as.numeric(depth); for (g in unique(donor)) { m <- donor == g; d[m] <- d[m] - mean(d[m]) }
  d2 <- d^2; for (g in unique(donor)) { m <- donor == g; d2[m] <- d2[m] - mean(d2[m]) }
  cbind(d / np_std(d), d2 / np_std(d2))
}
adjusted_ranks <- function(Xc, Z, chunk = 1500) {      # ranks of the residuals of Xc (donor-centred) after regression on Z, genes x cells, standardised as rank_std
  sv <- svd(Z); Zp <- sv$v %*% diag(1 / sv$d, length(sv$d)) %*% t(sv$u); out <- matrix(0, ncol(Xc), nrow(Xc))
  for (s in seq(1, ncol(Xc), by = chunk)) { e <- min(s + chunk - 1, ncol(Xc)); part <- Xc[, s:e, drop = FALSE]; out[s:e, ] <- rank_std(part - Z %*% (Zp %*% part)) }
  out
}

## ------------------------------------------------------------------------------------------------------------------------------------------ matched null
make_pool <- function(det_e, labels_e, own) {      # detection deciles of the eligible genes (rank-based qcut) and, per decile, the positions of the eligible genes that are not among the 17 analysed genes
  rk <- rank(det_e, ties.method = "first"); edges <- vapply(seq(0, 1, length.out = 11), function(p) np_quantile(rk, p), 1)
  bins <- as.integer(cut(rk, breaks = edges, include.lowest = TRUE, labels = FALSE)) - 1L
  list(bins = bins, pool = lapply(0:9, function(k) which(bins == k & !(labels_e %in% own))))
}
.choice1 <- function(rng, v) v[np_integers(rng, 0, length(v), 1) + 1L]            # rng.choice(v) for a 1-D array
draws_main <- function(pool, bins, pair_idx, nn = NN, seed = SEED_NULL) {      # index draws of the main null in the order of the earlier loop: for each draw, for each evaluable consensus pair: ia, ic (ic redrawn while equal to ia)
  rng <- np_rng(seed); np <- length(pair_idx); IA <- matrix(0L, nn, np); IC <- matrix(0L, nn, np)
  for (s in seq_len(nn)) for (j in seq_len(np)) { a <- pair_idx[[j]][1]; cc <- pair_idx[[j]][2]; ia <- .choice1(rng, pool[[bins[a] + 1L]]); ic <- .choice1(rng, pool[[bins[cc] + 1L]]); while (ic == ia) ic <- .choice1(rng, pool[[bins[cc] + 1L]]); IA[s, j] <- ia; IC[s, j] <- ic }
  list(IA = IA, IC = IC)
}
draws_pair <- function(pool, bins, pair_idx, nn = NN, seed = SEED_PAIR) {      # index draws of the pair-level null (pair-major order of the earlier loop, separate generator)
  rng <- np_rng(seed); np <- length(pair_idx); IA <- matrix(0L, nn, np); IC <- matrix(0L, nn, np)
  for (j in seq_len(np)) { a <- pair_idx[[j]][1]; cc <- pair_idx[[j]][2]; for (k in seq_len(nn)) { ia <- .choice1(rng, pool[[bins[a] + 1L]]); ic <- .choice1(rng, pool[[bins[cc] + 1L]]); while (ic == ia) ic <- .choice1(rng, pool[[bins[cc] + 1L]]); IA[k, j] <- ia; IC[k, j] <- ic } }
  list(IA = IA, IC = IC)
}
null_values <- function(Rk, IA, IC, step = 250L) {      # Spearman rho of the drawn gene pairs (nn x pairs), in chunks of draws (memory)
  out <- matrix(0, nrow(IA), ncol(IA))
  for (j in seq_len(ncol(IA))) for (s in seq(1, nrow(IA), by = step)) { e <- min(s + step - 1, nrow(IA)); out[s:e, j] <- rowSums(Rk[IA[s:e, j], , drop = FALSE] * Rk[IC[s:e, j], , drop = FALSE]) }
  out
}
## scipy.stats.mannwhitneyu(x, y) with the defaults (two-sided, continuity correction, method 'auto': exact when one sample has at most 8 values and there are no ties, else asymptotic with the tie correction)
py_mannwhitney_p <- function(x, y) {
  n1 <- length(x); n2 <- length(y); r <- rank(c(x, y)); t <- as.numeric(table(r)); U1 <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2; U <- max(U1, n1 * n2 - U1); ties <- any(t > 1)
  p <- if (n1 > 8 && n2 > 8 || ties) { mu <- n1 * n2 / 2; n <- n1 + n2; s <- sqrt(n1 * n2 / 12 * ((n + 1) - sum(t^3 - t) / (n * (n - 1)))); z <- (U - mu - 0.5) / s; 2 * pnorm(z, lower.tail = FALSE) }
       else 2 * pwilcox(as.integer(U) - 1, min(n1, n2), max(n1, n2), lower.tail = FALSE)
  min(max(p, 0), 1)
}
bh_adjust <- function(p) p.adjust(p, method = "BH")

## ------------------------------------------------------------------------------------------------------------------------------------------ analysis of one dataset
analyse <- function(cells, thr = 0.10, modes = c("unadjusted", "adjusted"), nn = NN) {
  ct <- consensus_tables(); cons <- ct$cons; consset <- ct$consset; redge <- ct$redge; key <- ct$key; sym <- function(g) cells_label(cells, g); n <- cells_n(cells); nd <- cells_donors(cells)
  det <- as.numeric(Matrix::colMeans(cells$B > 0)); elig <- which(det >= thr); labels_e <- cells$genes[elig]; idx <- setNames(seq_along(labels_e), labels_e); det_e <- det[elig]
  own <- vapply(GENES, sym, ""); mp <- make_pool(det_e, labels_e, own); bins <- mp$bins; pool <- mp$pool; Xc <- center_by_donor(lognorm(cells, elig), cells$donor, cells$kind)
  syms <- vapply(GENES, sym, ""); ing <- syms %in% cells$genes; evg <- syms %in% labels_e                                                                          # detection of the 17 genes (all genes of the data, evaluable or not)
  D <- data.frame(gene_invitro = GENES, gene = unname(syms), in_data = unname(ing), beta_cells = rep(as.integer(n), 17), detected_pct = ifelse(ing, 100 * det[match(syms, cells$genes)], NA_real_), evaluable = unname(evg), stringsAsFactors = FALSE)
  pr_a <- unlist(lapply(seq_along(GENES), function(i) rep(GENES[i], length(GENES) - i))); pr_c <- unlist(lapply(seq_along(GENES), function(i) GENES[-seq_len(i)]))
  out <- list(); draws <- NULL
  for (mode in modes) {
    if (mode == "unadjusted") { Rk <- ranks_of(Xc); dfree <- n - nd - 2 } else if (mode == "adjusted") { Rk <- adjusted_ranks(Xc, depth_design(cells_depth(cells), cells$donor)); dfree <- n - nd - 4 } else stop("mode")
    ea <- unname(syms[pr_a]); ec <- unname(syms[pr_c]); evb <- (ea %in% labels_e) & (ec %in% labels_e)
    rho <- rep(NA_real_, length(ea)); for (k in which(evb)) rho[k] <- sum(Rk[idx[[ea[k]]], ] * Rk[idx[[ec[k]]], ])
    tt <- rho * sqrt(dfree / (1 - rho^2)); pv <- ifelse(evb, 2 * pt(abs(tt), dfree, lower.tail = FALSE), NA_real_); keys <- mapply(key, pr_a, pr_c)
    P <- data.frame(gene1 = pr_a, gene2 = pr_c, consensus_invitro = keys %in% consset, r_invitro_Hes1norm = unname(redge[keys]), evaluable = evb,
                    detected_pct_1 = ifelse(ea %in% cells$genes, 100 * det[match(ea, cells$genes)], NA_real_), detected_pct_2 = ifelse(ec %in% cells$genes, 100 * det[match(ec, cells$genes)], NA_real_), rho = rho, p = pv, stringsAsFactors = FALSE)
    P$q <- NA_real_; P$q[evb] <- bh_adjust(P$p[evb]); ev <- P[P$evaluable, ]; c_ev <- ev[ev$consensus_invitro, ]; o_ev <- ev[!ev$consensus_invitro, ]; rc <- c_ev$rho; ro <- o_ev$rho
    s <- list(dataset = cells$name, mode = mode, detection_threshold = thr, beta_cells = as.integer(n), donors = as.integer(nd), df_for_p = as.integer(dfree), genes_in_data = sum(D$in_data), genes_evaluable = sum(D$evaluable),
              genes_not_evaluable = as.list(D$gene_invitro[!D$evaluable]), eligible_genes_in_pool = length(elig), consensus_pairs_total = length(cons), consensus_pairs_evaluable = nrow(c_ev), other_pairs_evaluable = nrow(o_ev))
    s <- c(s, list(median_rho_consensus = if (length(rc)) median(rc) else NULL, median_rho_other = if (length(ro)) median(ro) else NULL, mannwhitney_p_consensus_vs_other = if (length(rc) > 1 && length(ro) > 1) py_mannwhitney_p(rc, ro) else NULL,
                   consensus_positive = sum(rc > 0), consensus_sig_positive_q05 = if (length(rc)) sum(rc > 0 & c_ev$q < 0.05) else 0L, consensus_sig_negative_q05 = if (length(rc)) sum(rc < 0 & c_ev$q < 0.05) else 0L, consensus_rho_ge_0_2 = sum(rc >= 0.2)))
    s$consensus_sig_positive_q05 <- as.integer(s$consensus_sig_positive_q05); s$consensus_sig_negative_q05 <- as.integer(s$consensus_sig_negative_q05); s$consensus_positive <- as.integer(s$consensus_positive); s$consensus_rho_ge_0_2 <- as.integer(s$consensus_rho_ge_0_2)
    s$genes_in_data <- as.integer(s$genes_in_data); s$genes_evaluable <- as.integer(s$genes_evaluable); s$eligible_genes_in_pool <- as.integer(s$eligible_genes_in_pool); s$consensus_pairs_total <- as.integer(s$consensus_pairs_total)
    s$consensus_pairs_evaluable <- as.integer(s$consensus_pairs_evaluable); s$other_pairs_evaluable <- as.integer(s$other_pairs_evaluable)
    cn <- c_ev[, c("gene1", "gene2", "r_invitro_Hes1norm", "rho", "p", "q")]; rownames(cn) <- NULL
    if (nrow(c_ev) >= 1) {
      pidx <- lapply(seq_len(nrow(c_ev)), function(i) c(idx[[sym(c_ev$gene1[i])]], idx[[sym(c_ev$gene2[i])]]))
      if (is.null(draws)) draws <- list(m = draws_main(pool, bins, pidx, nn, SEED_NULL), p = draws_pair(pool, bins, pidx, nn, SEED_PAIR))
      V <- null_values(Rk, draws$m$IA, draws$m$IC); meds <- apply(V, 1, median); means <- rowMeans(V); zm <- rowMeans(atanh(V)); obs_med <- median(rc); obs_mean <- mean(rc); z_obs <- mean(atanh(rc))
      s <- c(s, list(null_median_rho_mean = mean(meds), null_median_rho_95pct = np_quantile(meds, 0.95), observed_median_rho = obs_med, empirical_p_median_ge_observed = (1 + sum(meds >= obs_med)) / (nn + 1), empirical_p_mean_ge_observed = (1 + sum(means >= obs_mean)) / (nn + 1),
                     mean_z_observed = z_obs, mean_z_null = mean(zm), E = z_obs - mean(zm), SE_E = sd(zm), empirical_p_E = (1 + sum(zm >= z_obs)) / (nn + 1)))
      V2 <- null_values(Rk, draws$p$IA, draws$p$IC); cn$null_mean <- colMeans(V2); cn$null_95pct <- vapply(seq_len(ncol(V2)), function(j) np_quantile(V2[, j], 0.95), 1)
      cn$empirical_p_ge_observed <- vapply(seq_len(ncol(V2)), function(j) (1 + sum(V2[, j] >= cn$rho[j])) / (nn + 1), 1)
    }
    E_ <- ev[!is.na(ev$q) & ev$q < 0.05 & abs(ev$rho) >= 0.2, ]; hd <- vapply(GENES, function(g) as.integer(sum(E_$gene1 == g | E_$gene2 == g)), 1L)
    s$edges_q05_rho02 <- as.integer(nrow(E_)); s$invitro_hubs_degree <- as.list(hd[c("Glp1r", "Tcf7l2", "Cacna1d", "Ptbp1")])
    prow <- P[(P$gene1 == "Glp1r" & P$gene2 == "Cacna1d") | (P$gene1 == "Cacna1d" & P$gene2 == "Glp1r"), ][1, ]
    s$Cacna1d_Glp1r <- list(evaluable = isTRUE(prow$evaluable), rho = if (isTRUE(prow$evaluable)) prow$rho else NULL, p = if (isTRUE(prow$evaluable)) prow$p else NULL,
                            detected_pct_Glp1r = if (prow$gene1 == "Glp1r") prow$detected_pct_1 else prow$detected_pct_2, detected_pct_Cacna1d = if (prow$gene1 == "Glp1r") prow$detected_pct_2 else prow$detected_pct_1)
    out[[mode]] <- list(pairs = P, cons = cn, detection = D, summary = s); rm(Rk); gc(FALSE)
  }
  out
}

## ------------------------------------------------------------------------------------------------------------------------------------------ shared preparation of the add-on analyses (gene-graph null, positive control, metacells)
CLIP <- 0.999999
prep_adjusted <- function(cells, thr = 0.10) {      # eligible genes, detection deciles, donor-centred expression, depth design, adjusted ranks and the evaluable consensus pairs (GENES order) of one dataset
  ct <- consensus_tables(); consset <- ct$consset; key <- ct$key; sym <- function(g) cells_label(cells, g)
  det <- as.numeric(Matrix::colMeans(cells$B > 0)); elig <- which(det >= thr); labels_e <- cells$genes[elig]; idx <- setNames(seq_along(labels_e), labels_e); det_e <- det[elig]; own <- vapply(GENES, sym, "")
  mp <- make_pool(det_e, labels_e, own); Xc <- center_by_donor(lognorm(cells, elig), cells$donor, cells$kind); Z <- depth_design(cells_depth(cells), cells$donor); Rk <- adjusted_ranks(Xc, Z)
  pairs <- list(); for (i in seq_along(GENES)) for (j in seq_along(GENES)) if (j > i) { a <- GENES[i]; cc <- GENES[j]; if (key(a, cc) %in% consset && sym(a) %in% labels_e && sym(cc) %in% labels_e) pairs[[length(pairs) + 1]] <- c(a, cc) }
  list(elig = elig, labels_e = labels_e, idx = idx, det_e = det_e, own = unname(own), bins = mp$bins, pool = mp$pool, Xc = Xc, Z = Z, Rk = Rk, pairs = pairs, pidx = lapply(pairs, function(p) c(idx[[sym(p[1])]], idx[[sym(p[2])]])), sym = sym)
}
