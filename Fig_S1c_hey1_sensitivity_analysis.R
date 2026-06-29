# ==============================================================================
# Fig_S1c_hey1_sensitivity_analysis.R
# Hey1-for-Hes1 normalization sensitivity reanalysis (Supplementary Table 2)
#
# Tests whether the 27-pair consensus co-expression network is robust to
# substituting Hey1 for Hes1 as the Notch-target normalization denominator.
#
# A gene pair (g1, g2) is analytically validated if |r| >= 0.70 with consistent
# sign across every normalization condition in which NEITHER gene serves as the
# denominator. This matches the recurrence criterion used for the 27-pair
# consensus network reported in the manuscript.
#
# Two schemes are evaluated:
#   Original    : raw + /Hes1 + /Kcnj11 + /Ins1
#   Alternative : raw + /Hey1 + /Kcnj11 + /Ins1
#
# Input
# -----
#   Supplementary_Data_1.xlsx — sheet "Supplementary Data 1"
#   Rows 5+: gene name (col A) + 12 replicate dCt values
#   Columns: CTRL (cols G-J), DAPT (cols M-P), DKK-1 (cols S-V)
#   Set DATA_FILE below to the path of this file.
#
# Output
# ------
#   Console: validated pairs per scheme, overlap summary
#   Fig_S1c_hey1_retention_summary.csv: pair-level retention table
#
# Dependencies: readxl (install.packages("readxl"))
# ==============================================================================

library(readxl)

# ── User-configurable path ────────────────────────────────────────────────────
DATA_FILE <- "Supplementary_Data_1.xlsx"
OUT_CSV   <- "Fig_S1c_hey1_retention_summary.csv"
# ─────────────────────────────────────────────────────────────────────────────

THRESHOLD <- 0.70

# ── Load data ─────────────────────────────────────────────────────────────────
raw_sheet <- read_excel(DATA_FILE, sheet = "Supplementary Data 1",
                        col_names = FALSE, skip = 4)

# Extract gene names and replicate columns (G-J = cols 7-10, M-P = 13-16, S-V = 19-22)
# Columns are 1-indexed; skip = 4 means row 5 is now row 1 of raw_sheet
CTRL_COLS  <- 7:10
DAPT_COLS  <- 13:16
DKK1_COLS  <- 19:22

EXCLUDE_STARTS <- c("Abbreviations", "Trend", "CTRL", "Hey1,", "†", "Ins2")

raw_data <- list()
for (i in seq_len(nrow(raw_sheet))) {
  gene <- raw_sheet[[i, 1]]
  if (is.na(gene) || !is.character(gene)) next
  if (any(startsWith(gene, EXCLUDE_STARTS)))  next

  vals <- unlist(raw_sheet[i, c(CTRL_COLS, DAPT_COLS, DKK1_COLS)])
  if (any(is.na(vals)) || !all(sapply(vals, is.numeric))) next
  raw_data[[gene]] <- as.numeric(vals)
}

all_genes <- sort(names(raw_data))
cat(sprintf("Genes loaded (%d): %s\n\n", length(all_genes), paste(all_genes, collapse = ", ")))

# ── Validated-pair computation ────────────────────────────────────────────────
compute_validated_pairs <- function(norm_gene_list) {
  # Build expression matrix for each condition
  conditions <- list(raw = do.call(rbind, lapply(raw_data, identity)))
  rownames(conditions$raw) <- all_genes

  for (ng in norm_gene_list) {
    denom <- raw_data[[ng]]
    mat   <- do.call(rbind, lapply(raw_data, function(x) x / denom))
    rownames(mat) <- all_genes
    conditions[[ng]] <- mat
  }

  validated <- list()
  gene_pairs <- combn(all_genes, 2, simplify = FALSE)

  for (pair in gene_pairs) {
    g1 <- pair[1]; g2 <- pair[2]

    # Conditions where both genes are present as targets (not denominator)
    applicable <- names(conditions)[
      names(conditions) == "raw" | (names(conditions) != g1 & names(conditions) != g2)
    ]

    rs <- c()
    valid <- TRUE
    for (cname in applicable) {
      mat <- conditions[[cname]]
      x <- mat[g1, ]; y <- mat[g2, ]
      if (sd(x) < 1e-10 || sd(y) < 1e-10) { valid <- FALSE; break }
      r <- cor(x, y, method = "pearson")
      if (is.na(r) || abs(r) < THRESHOLD) { valid <- FALSE; break }
      rs <- c(rs, r)
    }

    if (valid && length(rs) > 0) {
      if (all(rs > 0) || all(rs < 0)) {
        key <- paste(sort(c(g1, g2)), collapse = "__")
        validated[[key]] <- list(
          gene1      = sort(c(g1, g2))[1],
          gene2      = sort(c(g1, g2))[2],
          rs         = rs,
          conditions = applicable,
          mean_r     = mean(rs)
        )
      }
    }
  }
  validated
}

# ── Run both schemes ──────────────────────────────────────────────────────────
cat("=== ORIGINAL: raw + /Hes1 + /Kcnj11 + /Ins1 ===\n")
orig <- compute_validated_pairs(c("Hes1", "Kcnj11", "Ins1"))
cat(sprintf("Validated pairs: %d\n", length(orig)))
for (v in orig[order(sapply(orig, function(x) -abs(x$mean_r)))]) {
  cat(sprintf("  %-12s — %-12s  [%s]  mean r = %.3f\n",
              v$gene1, v$gene2,
              paste(v$conditions, collapse = ", "),
              v$mean_r))
}

cat("\n=== ALTERNATIVE: raw + /Hey1 + /Kcnj11 + /Ins1 ===\n")
alt <- compute_validated_pairs(c("Hey1", "Kcnj11", "Ins1"))
cat(sprintf("Validated pairs: %d\n", length(alt)))
for (v in alt[order(sapply(alt, function(x) -abs(x$mean_r)))]) {
  cat(sprintf("  %-12s — %-12s  [%s]  mean r = %.3f\n",
              v$gene1, v$gene2,
              paste(v$conditions, collapse = ", "),
              v$mean_r))
}

# ── Overlap summary ───────────────────────────────────────────────────────────
o_keys <- names(orig); a_keys <- names(alt)
shared  <- intersect(o_keys, a_keys)
only_o  <- setdiff(o_keys, a_keys)
only_a  <- setdiff(a_keys, o_keys)

cat(sprintf("\n=== OVERLAP SUMMARY ===\n"))
cat(sprintf("Original pairs:             %d\n", length(o_keys)))
cat(sprintf("Alternative (Hey1) pairs:   %d\n", length(a_keys)))
cat(sprintf("Shared (retained):          %d  (%.1f%% of original)\n",
            length(shared), length(shared) / length(o_keys) * 100))
cat(sprintf("Lost (original only):       %d\n", length(only_o)))
cat(sprintf("Gained (alternative only):  %d\n", length(only_a)))

if (length(only_o) > 0) {
  cat("\n  Pairs lost when Hes1 is replaced by Hey1:\n")
  for (k in only_o) cat(sprintf("    %s — %s\n", orig[[k]]$gene1, orig[[k]]$gene2))
}
if (length(only_a) > 0) {
  cat("\n  New pairs appearing only under Hey1 normalization:\n")
  for (k in only_a) cat(sprintf("    %s — %s\n", alt[[k]]$gene1, alt[[k]]$gene2))
}

# ── Export retention table ────────────────────────────────────────────────────
all_keys  <- union(o_keys, a_keys)
out_rows  <- lapply(all_keys, function(k) {
  if (k %in% o_keys) {
    v <- orig[[k]]
    data.frame(
      Gene1      = v$gene1,
      Gene2      = v$gene2,
      In_Original  = TRUE,
      In_Alternative = k %in% a_keys,
      Mean_r_original = round(v$mean_r, 3),
      stringsAsFactors = FALSE
    )
  } else {
    v <- alt[[k]]
    data.frame(
      Gene1      = v$gene1,
      Gene2      = v$gene2,
      In_Original  = FALSE,
      In_Alternative = TRUE,
      Mean_r_original = NA_real_,
      stringsAsFactors = FALSE
    )
  }
})
out_df <- do.call(rbind, out_rows)
out_df <- out_df[order(!out_df$In_Original, !out_df$In_Alternative), ]
write.csv(out_df, OUT_CSV, row.names = FALSE)
cat(sprintf("\nRetention table saved: %s\n", OUT_CSV))

cat(sprintf("\nIns2 in panel: %s\n", "Ins2" %in% all_genes))
cat("(Ins2 substitution for Ins1 would require new RT-qPCR experiments.)\n")
