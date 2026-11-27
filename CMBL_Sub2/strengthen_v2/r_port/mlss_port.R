## mlss_port.R - R port of wp3_mlss/mlss_port.py: the MLSS v4.0 as a function (the same scoring as R_17genes/Fig_6c_MLSS_v4.0.r, written as functions so that
## scenarios can be run), plus the upstream network step of Fig_4c/4d/4e/6b (four normalization schemes -> consensus pairs -> degree).
## MLSS = [0.45 C + 0.30 B + 0.10 V + 0.15 P] x 9 + Synergy - Antagonism (see mlss_port.py for the definitions and the details of the R script that it reproduces).
## Gene pairs are keyed "A|B" with the two names in code-point order; drug-gene potency keys are "DRUG\tGENE".

MLSS_W_STD <- c(0.45, 0.30, 0.10, 0.15); MLSS_BASE_MULT <- 9.0
MLSS_GENES17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
MLSS_GROUP <- rep(1:3, each = 4); MLSS_DENOMINATORS <- c("Hes1", "Kcnj11", "Ins1"); MLSS_KEY <- "CACNA1D|GLP1R"
mlss_pkey <- function(a, b) { ab <- sort(c(a, b), method = "radix"); paste(ab[1], ab[2], sep = "|") }
mlss_pkeys <- function(a, b) mapply(mlss_pkey, a, b, USE.NAMES = FALSE)
mlss_order_key <- function(s) order(tolower(s), s, method = "radix")       # Python r_sort_key: (casefold, s), the approximation of R's locale collation used by the port

## ---- inputs of Fig_6c ------------------------------------------------------------------------------------------------
mlss_read_drugcentral <- function(path) {
  dc <- read.csv(path, colClasses = "character", na.strings = "", check.names = FALSE, quote = "\"", comment.char = "")
  dc$ACT_VALUE <- suppressWarnings(as.numeric(dc$ACT_VALUE)); dc$G <- toupper(dc$GENE); dc
}
mlss_potency_from <- function(dc) {      # named vector "DRUG\tGENE" -> max ACT_VALUE over human IC50/Ki/EC50 rows with RELATION '='
  p <- dc[dc$ACT_TYPE %in% c("IC50", "Ki", "EC50") & !is.na(dc$ORGANISM) & dc$ORGANISM == "Homo sapiens" & !is.na(dc$RELATION) & dc$RELATION == "=" & !is.na(dc$ACT_VALUE), ]
  key <- paste(p$DRUG_NAME, p$G, sep = "\t"); v <- tapply(p$ACT_VALUE, key, max); setNames(as.numeric(v), names(v))
}
mlss_targets_from <- function(dc, hubs) {  # named list drug -> hub genes (unique, order of first appearance)
  hubs <- toupper(hubs); t <- dc[!is.na(dc$G) & dc$G %in% hubs & !is.na(dc$ORGANISM) & dc$ORGANISM == "Homo sapiens", c("DRUG_NAME", "GENE", "TARGET_NAME", "ACTION_TYPE", "G")]
  t <- t[!duplicated(t[, c("DRUG_NAME", "GENE", "TARGET_NAME", "ACTION_TYPE")]), ]; out <- list()
  for (i in seq_len(nrow(t))) { d <- t$DRUG_NAME[i]; g <- t$G[i]; if (is.null(out[[d]])) out[[d]] <- character(0); if (!(g %in% out[[d]])) out[[d]] <- c(out[[d]], g) }
  out
}
mlss_hubs_from_centrality <- function(path, min_degree = 5) { h <- read.csv(path, stringsAsFactors = FALSE); list(hubs = toupper(h$Gene[h$Degree >= min_degree]), table = h) }
mlss_corr_from_matrix_csv <- function(path) {   # Fig_4e_Hes1_significant_correlations_only.csv -> named vector "A|B" -> r (NA cells dropped)
  m <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE); genes <- toupper(m[[1]]); cols <- toupper(names(m)[-1]); out <- numeric(0)
  for (i in seq_along(genes)) for (j in seq_along(cols)) { a <- genes[i]; b <- cols[j]; v <- m[i, j + 1]; if (a != b && !is.na(v)) { k <- mlss_pkey(a, b); if (is.na(out[k])) out[k] <- as.numeric(v) } }
  out
}

## ---- the scoring (steps 7-14 of the R script) --------------------------------------------------------------------------
mlss_mean_pic50 <- function(drug, tg, pot) { v <- pot[paste(drug, unique(tg), sep = "\t")]; v <- v[!is.na(v)]; if (length(v)) mean(v) else 0 }
mlss_mult <- function(da, gi, db, gj, pot) {
  pa <- pot[paste(da, gi, sep = "\t")]; pb <- pot[paste(db, gj, sep = "\t")]
  if (is.na(pa) || is.na(pb)) 0 else min(max(pa, pb) / 9.0, 1.0)
}
mlss_table <- function(targets, pot, corr, weights = MLSS_W_STD, base_mult = MLSS_BASE_MULT, corr_v = NULL) {
  cv <- if (is.null(corr_v)) corr else corr_v; ds <- names(targets); ds <- ds[mlss_order_key(ds)]; rows <- list()
  for (j in seq_along(ds)) for (i in seq_len(j - 1)) {
    a <- ds[i]; b <- ds[j]; ta <- unique(targets[[a]]); tb <- unique(targets[[b]])
    inter <- length(intersect(ta, tb)); uni <- length(union(ta, tb)); C <- if (uni > 0) 1 - inter / uni else 0
    Bb <- if ((length(ta) + length(tb)) > 0) 1 - abs(length(ta) - length(tb)) / (length(ta) + length(tb)) else 0
    rs <- numeric(0); rsv <- numeric(0); syn <- 0; ant <- 0
    for (gi in ta) for (gj in tb) if (gi != gj) {
      k <- mlss_pkey(gi, gj); kc <- corr[k]; kv <- cv[k]
      if (!is.na(kv)) rsv <- c(rsv, kv)
      if (!is.na(kc)) { rs <- c(rs, kc); m <- mlss_mult(a, gi, b, gj, pot); if (kc > 0) syn <- syn + kc * m * 1 else if (kc < 0) ant <- ant + abs(kc) * m * 1 }
    }
    pos <- rs[rs > 0]; neg <- rs[rs < 0]; mabs <- if (length(rsv)) mean(abs(rsv)) else 0
    mA <- mlss_mean_pic50(a, ta, pot); mB <- mlss_mean_pic50(b, tb, pot); mp <- (mA + mB) / 2.0; P <- min(max(mp - 4, 0), 5) / 5
    rows[[length(rows) + 1]] <- list(Drug_A = a, Drug_B = b, Targets_A = paste(ta, collapse = "|"), Targets_B = paste(tb, collapse = "|"), n_targets_A = length(ta), n_targets_B = length(tb), n_targets_combined = length(ta) + length(tb),
      Mean_pIC50_A = mA, Mean_pIC50_B = mB, Mean_pIC50 = mp, Complementarity_C = C, Balance_B = Bb, Mean_Abs_r = mabs, Mean_Pos_r = if (length(pos)) mean(pos) else 0, Mean_Neg_r = if (length(neg)) mean(neg) else 0,
      Potency_P = P, Synergy_Bonus = syn, Antagonism_Penalty = ant)
  }
  t <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)); if (is.null(t) || nrow(t) == 0) return(data.frame())
  lo <- min(t$Mean_Abs_r); hi <- max(t$Mean_Abs_r); t$Network_Coverage_V <- if (hi > lo) (t$Mean_Abs_r - lo) / (hi - lo) else 0
  t$Base_Score <- (weights[1] * t$Complementarity_C + weights[2] * t$Balance_B + weights[3] * t$Network_Coverage_V + weights[4] * t$Potency_P) * base_mult
  t[["MLSS_v4.0"]] <- t$Base_Score + t$Synergy_Bonus - t$Antagonism_Penalty; t$Combination <- paste(t$Drug_A, t$Drug_B, sep = " + ")
  t <- t[order(-t[["MLSS_v4.0"]]), ]; rownames(t) <- NULL; t$Rank <- seq_len(nrow(t)); t
}
mlss_run_reported <- function(dc, hubs, corr) mlss_table(mlss_targets_from(dc, hubs), mlss_potency_from(dc), corr)

## ---- upstream: INS-1 matrix -> four schemes -> consensus -> degree --------------------------------------------------------
mlss_gapdh_variant <- function(X, which) {          # C0 / C1 recalculation of the 26 August DAPT values of five genes (constants of reanalysis_17genes/gapdh_exact.py)
  if (which == "reported") return(X)
  hk <- c(22.5685093010354, 22.7982726003435, 22.5685093010354, 22.7982726003435); raw <- c(A10 = 24.9873996418614, B10 = 24.1318270358876, C10 = 22.5685093010354, D10 = 22.7982726003435)
  ref <- if (which == "C0") mean(c(raw[["C10"]], raw[["D10"]])) else mean(raw); d <- X
  for (g in c("Tcf7l2", "Kcnj11", "Cacna1d", "Cacna1c", "Ptbp1")) d[5:8, g] <- X[5:8, g] * 2.0^(ref - hk)
  d
}
mlss_pearson_p <- function(M, df_resid = NULL) {
  n <- nrow(M); r <- suppressWarnings(cor(M)); df <- if (is.null(df_resid)) n - 2 else df_resid
  t <- r * sqrt(df / (1 - r^2)); p <- 2 * pt(-abs(t), df); diag(p) <- NA_real_; list(r = r, p = p)
}
mlss_four_schemes <- function(X) {
  out <- list(Raw = X)
  for (d in MLSS_DENOMINATORS) out[[paste0(d, "-normalized")]] <- X[, setdiff(colnames(X), d), drop = FALSE] / X[, d]
  out
}
mlss_network_from_matrix <- function(X) {
  sch <- mlss_four_schemes(X); res <- list(); sig <- list()
  for (k in names(sch)) {
    m <- sch[[k]]; rp <- mlss_pearson_p(m); g <- colnames(m); ij <- t(combn(length(g), 2)); keys <- mlss_pkeys(g[ij[, 1]], g[ij[, 2]])
    res[[k]] <- list(keys = keys, r = rp$r[ij], p = rp$p[ij]); sig[[k]] <- keys[rp$p[ij] < 0.05]
  }
  cons <- Reduce(intersect, sig); deg <- list()
  for (pr in cons) for (g in strsplit(pr, "|", fixed = TRUE)[[1]]) deg[[g]] <- (if (is.null(deg[[g]])) 0L else deg[[g]]) + 1L
  list(res = res, sig_counts = vapply(sig, length, 1L), consensus = py_sort(cons), degree = deg)
}
mlss_hes1_lookup <- function(net, alpha = 0.05, source = "Hes1-normalized", signif_only = TRUE) {
  s <- net$res[[source]]; keep <- if (signif_only) s$p < alpha else rep(TRUE, length(s$p)); keys <- toupper(s$keys[keep])
  keys <- vapply(strsplit(keys, "|", fixed = TRUE), function(p) mlss_pkey(p[1], p[2]), ""); setNames(s$r[keep], keys)
}
mlss_adjusted_lookup <- function(X, log = FALSE, df_resid = 8, alpha = 0.05, signif_only = TRUE) {
  M <- if (log) base::log(X) else X; R <- M
  for (g in 1:3) R[MLSS_GROUP == g, ] <- M[MLSS_GROUP == g, ] - matrix(colMeans(M[MLSS_GROUP == g, , drop = FALSE]), 4, ncol(M), byrow = TRUE)
  rp <- mlss_pearson_p(R, df_resid); cols <- toupper(colnames(X)); ij <- t(combn(length(cols), 2)); keep <- if (signif_only) rp$p[ij] < alpha else rep(TRUE, nrow(ij))
  setNames(rp$r[ij][keep], mlss_pkeys(cols[ij[keep, 1]], cols[ij[keep, 2]]))
}
