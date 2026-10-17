# Treatment-adjusted correlation reanalysis of the Fig. 4 co-expression network,
# with the CORRECT degrees of freedom, plus a group-preserving permutation null.
# Base R only. R implementation of analysis/corr_reanalysis_df8_and_null.py,
# written so the two can be cross-checked.
#
# Inputs (Fig. 4c exports of the original pipeline, unchanged; 12 rows =
# 4 technical replicates x 3 groups; row order CTRL x4, DAPT x4, DKK-1 x4,
# confirmed row by row against the Treatment column of Supplementary Data 1):
#   Fig_4c_Raw_expression_data.csv, Fig_4c_Hes1_normalized_data.csv,
#   Fig_4c_Kcnj11_normalized_data.csv, Fig_4c_Ins1_normalized_data.csv
# Ryr1 is dropped (excluded everywhere else in the study).
#
# Definitions
#   pooled     : Pearson r across all 12 samples, df = 12 - 2 = 10
#   adjusted   : Pearson r of residuals after removing group means
#                (partial correlation given treatment), df = 12 - 2 - 2 = 8
#   BH         : Benjamini-Hochberg FDR within each normalization scheme,
#                over all gene pairs of that scheme (190 for Raw, 171 otherwise)
#   consensus  : pair significant in all four schemes
#   null       : each gene's group means are kept and its within-group
#                residuals (log scale) are permuted independently within group,
#                which preserves treatment structure and removes within-group
#                co-regulation; the four schemes are rebuilt from the permuted
#                raw matrix, so the shared-denominator structure is preserved.
#
# Usage: Rscript R_corr_reanalysis_df8_and_null.R [data_dir] [n_perm] [out_dir]

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "input_data"
n_perm   <- if (length(args) >= 2) as.integer(args[2]) else 2000L
out_dir  <- if (length(args) >= 3) args[3] else "."
SEED     <- 20261015
set.seed(SEED)
invisible(Sys.setlocale("LC_COLLATE", "C"))

files <- c(Raw    = "Fig_4c_Raw_expression_data.csv",
           Hes1   = "Fig_4c_Hes1_normalized_data.csv",
           Kcnj11 = "Fig_4c_Kcnj11_normalized_data.csv",
           Ins1   = "Fig_4c_Ins1_normalized_data.csv")
S <- lapply(files, function(f) {
  d <- read.csv(file.path(data_dir, f), check.names = FALSE)
  d[, setdiff(names(d), "Ryr1"), drop = FALSE]
})
stopifnot(all(vapply(S, nrow, 1L) == 12))
all_genes <- names(S$Raw)                                # 20 genes (Ryr1 dropped)
grp <- factor(rep(c("CTRL", "DAPT", "DKK-1"), each = 4), levels = c("CTRL", "DAPT", "DKK-1"))

resid_grp <- function(M) apply(M, 2, function(x) x - ave(x, grp))

pair_code <- function(g1, g2) {                          # order-free integer key
  i <- match(g1, all_genes); j <- match(g2, all_genes)
  pmin(i, j) * 100L + pmax(i, j)
}

corr_p <- function(M, df) {
  r  <- cor(M)
  ut <- upper.tri(r)
  ix <- which(ut, arr.ind = TRUE)
  rv <- r[ut]
  t  <- rv * sqrt(df / pmax(1 - rv^2, 1e-12))
  list(g1 = colnames(M)[ix[, 1]], g2 = colnames(M)[ix[, 2]], r = rv,
       p = 2 * pt(-abs(t), df))
}

analyse <- function(sch, adjusted) {
  lapply(sch, function(d) {
    M  <- as.matrix(d)
    df <- if (adjusted) 12 - 2 - 2 else 12 - 2
    if (adjusted) M <- resid_grp(M)
    o  <- corr_p(M, df)
    data.frame(g1 = o$g1, g2 = o$g2, r = o$r, p = o$p,
               q = p.adjust(o$p, method = "BH"), stringsAsFactors = FALSE)
  })
}

sig_codes  <- function(d, col, thr = 0.05) pair_code(d$g1[d[[col]] < thr], d$g2[d[[col]] < thr])
consensus  <- function(res, col) Reduce(intersect, lapply(res, sig_codes, col = col))

summarize <- function(res, label) {
  cat(sprintf("\n== %s\n", label))
  for (cc in c("p", "q")) {
    nm <- if (cc == "p") "uncorrected p<0.05" else "BH q<0.05"
    counts <- vapply(res, function(d) sum(d[[cc]] < 0.05), 1L)
    cat(sprintf("  %-20s per-scheme {%s} | consensus(all 4) = %d\n", nm,
                paste(sprintf("'%s': %d", names(counts), counts), collapse = ", "),
                length(consensus(res, cc))))
  }
}

cat("N rows = 12 | pooled df = 10 | treatment-adjusted (partial) df = 8\n")
pooled <- analyse(S, FALSE)
adj8   <- analyse(S, TRUE)
summarize(pooled, "POOLED across 3 treatment groups")
summarize(adj8,   "TREATMENT-ADJUSTED (residuals on group), correct df = 8")

# sign consistency of the pooled p<0.05 consensus (the manuscript's stated criterion)
cons <- sort(consensus(pooled, "p"))
sign_of <- function(res, code) vapply(res, function(d) {
  k <- pair_code(d$g1, d$g2); sign(d$r[match(code, k)])
}, 1)
incons <- cons[vapply(cons, function(cd) length(unique(sign_of(pooled, cd))) > 1, TRUE)]
dec <- function(cd) paste(all_genes[cd %/% 100L], all_genes[cd %% 100L], sep = "-")
cat(sprintf("\nPooled p<0.05 consensus = %d; sign-INCONSISTENT across schemes = %d -> sign-consistent consensus = %d\n",
            length(cons), length(incons), length(cons) - length(incons)))
cat("  inconsistent:", paste(sort(vapply(incons, dec, "")), collapse = ", "), "\n")

# ---- permutation null --------------------------------------------------------
raw <- as.matrix(S$Raw)
L   <- log(raw)
fitted <- L - resid_grp(L)                    # group means (log scale)
E      <- resid_grp(L)

schemes_from_raw <- function(rawp) {
  out <- list(Raw = rawp)
  for (den in c("Hes1", "Kcnj11", "Ins1"))
    out[[den]] <- rawp[, setdiff(colnames(rawp), den), drop = FALSE] / rawp[, den]
  out
}

idx_by_grp <- split(seq_len(12), grp)
cons_n <- integer(n_perm); cons_signed_n <- integer(n_perm)
per_scheme_n <- matrix(0L, n_perm, 4, dimnames = list(NULL, names(files)))
for (b in seq_len(n_perm)) {
  Ep <- E
  for (idx in idx_by_grp) Ep[idx, ] <- apply(E[idx, , drop = FALSE], 2, sample)
  sch <- schemes_from_raw(exp(fitted + Ep))
  res <- lapply(sch, function(d) { o <- corr_p(as.matrix(d), 10); data.frame(g1 = o$g1, g2 = o$g2, r = o$r, p = o$p, stringsAsFactors = FALSE) })
  per_scheme_n[b, ] <- vapply(res[names(files)], function(d) sum(d$p < 0.05), 1L)
  cc <- Reduce(intersect, lapply(res, sig_codes, col = "p"))
  cons_n[b] <- length(cc)
  cons_signed_n[b] <- sum(vapply(cc, function(cd) length(unique(sign_of(res, cd))) == 1, TRUE))
}
write.csv(data.frame(rep = seq_len(n_perm), per_scheme_n, consensus = cons_n, consensus_sign_consistent = cons_signed_n),
          file.path(out_dir, "R_null_counts.csv"), row.names = FALSE)
q95 <- function(x) sprintf("%.0f-%.0f", quantile(x, 0.025, names = FALSE), quantile(x, 0.975, names = FALSE))
cat(sprintf("\n== PERMUTATION NULL (%d reps; independent within-group residuals; group means retained; seed %d)\n", n_perm, SEED))
for (k in names(files))
  cat(sprintf("  %-7s pooled p<0.05 pairs: mean %.1f (95%% %s)\n", k, mean(per_scheme_n[, k]), q95(per_scheme_n[, k])))
n_obs <- length(cons); n_obs_s <- length(cons) - length(incons)
cat(sprintf("  4-scheme consensus pairs: mean %.1f (95%% %s); P(null >= %d) = %.3f\n", mean(cons_n), q95(cons_n), n_obs, mean(cons_n >= n_obs)))
cat(sprintf("  ... sign-consistent consensus: mean %.1f (95%% %s); P(null >= %d) = %.3f\n", mean(cons_signed_n), q95(cons_signed_n), n_obs_s, mean(cons_signed_n >= n_obs_s)))

# ---- smallest |r| that can pass BH at df = 8 ----------------------------------
# With BH, the single best pair passes only if p <= q * 1/m. This is the minimum
# |r| a pair needs even in the most favourable case (it is the only discovery).
cat("\n== Minimum |r| to pass BH at df = 8 (single best pair; m = number of pairs tested)\n")
for (m in c(190, 171)) for (qq in c(0.05, 0.10)) {
  tcrit <- qt(1 - (qq / m) / 2, 8)
  cat(sprintf("  m = %d, q < %.2f: p <= %.2e -> |r| >= %.3f\n", m, qq, qq / m, tcrit / sqrt(tcrit^2 + 8)))
}

# ---- per-pair tables ------------------------------------------------------------
all_rows <- do.call(rbind, lapply(names(files), function(k) {
  a <- pooled[[k]]; b <- adj8[[k]]
  data.frame(scheme = k, gene1 = a$g1, gene2 = a$g2,
             r_pooled = a$r, p_pooled_df10 = a$p, q_pooled_BH = a$q,
             r_adjusted = b$r, p_adjusted_df8 = b$p, q_adjusted_BH = b$q,
             stringsAsFactors = FALSE)
}))
write.csv(all_rows, file.path(out_dir, "R_all_pairs_statistics.csv"), row.names = FALSE)

cons_tab <- do.call(rbind, lapply(cons, function(cd) {
  row <- data.frame(gene1 = all_genes[cd %/% 100L], gene2 = all_genes[cd %% 100L], stringsAsFactors = FALSE)
  for (k in names(files)) {
    a <- pooled[[k]]; b <- adj8[[k]]; i <- match(cd, pair_code(a$g1, a$g2))
    row[[paste0("r_pooled_", k)]]   <- a$r[i]
    row[[paste0("r_adj_", k)]]      <- b$r[i]
    row[[paste0("p_adj_", k)]]      <- b$p[i]
    row[[paste0("q_adj_", k)]]      <- b$q[i]
  }
  row$sign_consistent_pooled <- length(unique(sign_of(pooled, cd))) == 1
  row$n_schemes_p_adj_lt_0.05 <- sum(vapply(adj8, function(d) d$p[match(cd, pair_code(d$g1, d$g2))] < 0.05, TRUE))
  row$n_schemes_q_adj_lt_0.05 <- sum(vapply(adj8, function(d) d$q[match(cd, pair_code(d$g1, d$g2))] < 0.05, TRUE))
  row
}))
write.csv(cons_tab, file.path(out_dir, "R_SuppTable3_v2_consensus_pairs.csv"), row.names = FALSE)
cat(sprintf("\nWrote R_all_pairs_statistics.csv (%d rows) and R_SuppTable3_v2_consensus_pairs.csv (%d rows)\n", nrow(all_rows), nrow(cons_tab)))
cat("\nsessionInfo():\n"); print(R.version.string)
