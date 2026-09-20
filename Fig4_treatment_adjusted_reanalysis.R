# Reanalysis of the Fig.4 consensus co-expression network for npjSBA Rev3.
# Addresses Reviewer 3 Major Concern 2: treatment-pooled correlations without
# multiple-testing correction, and the 28/27/19 consensus-pair-count discrepancy.
#
# Inputs (already exported by the original pipeline, unchanged):
#   Figures/Fig_4c_Hes1_normalized_data.csv
#   Figures/Fig_4c_Kcnj11_normalized_data.csv
#   Figures/Fig_4c_Ins1_normalized_data.csv
#   Figures/Fig_4c_Raw_expression_data.csv
#
# Row order (12 rows = 4 reps x 3 groups) inferred from Hes1 expression pattern
# (Notch target; DAPT is a Notch inhibitor so Hes1 should be lowest in DAPT):
#   rows 1-4   = CTRL
#   rows 5-8   = DAPT
#   rows 9-12  = DKK-1

# Run from a directory containing the R_scripts_Rev1 pipeline's Fig_4c_*.csv
# outputs (see README.md, Input Data Files Required).
BASE <- "."

files <- list(
  "Hes1-normalized"   = "Fig_4c_Hes1_normalized_data.csv",
  "Kcnj11-normalized" = "Fig_4c_Kcnj11_normalized_data.csv",
  "Ins1-normalized"   = "Fig_4c_Ins1_normalized_data.csv",
  "Raw"               = "Fig_4c_Raw_expression_data.csv"
)

group <- factor(rep(c("CTRL", "DAPT", "DKK1"), each = 4), levels = c("CTRL", "DAPT", "DKK1"))

sanity_check_group_order <- function(raw_df) {
  # Hes1 is not in the normalized files (used as denominator), but IS in Raw.
  if ("Hes1" %in% names(raw_df)) {
    means <- tapply(raw_df$Hes1, group, mean)
    cat("Hes1 mean by assumed group (expect DAPT lowest):\n")
    print(means)
    cat("\n")
  }
}

pooled_correlations <- function(df) {
  genes <- names(df)
  pairs <- combn(genes, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    ct <- cor.test(df[[pr[1]]], df[[pr[2]]])
    data.frame(gene1 = pr[1], gene2 = pr[2], r = unname(ct$estimate), p_raw = ct$p.value)
  })
  out <- do.call(rbind, rows)
  out$p_bh <- p.adjust(out$p_raw, method = "BH")
  out
}

treatment_adjusted_correlations <- function(df) {
  # residualize each gene on group (treatment contrasts, CTRL as reference) via lm(),
  # then correlate residuals
  resid_df <- as.data.frame(lapply(df, function(y) resid(lm(y ~ group))))
  names(resid_df) <- names(df)
  genes <- names(df)
  pairs <- combn(genes, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pr) {
    ct <- cor.test(resid_df[[pr[1]]], resid_df[[pr[2]]])
    data.frame(gene1 = pr[1], gene2 = pr[2], r = unname(ct$estimate), p_raw = ct$p.value)
  })
  out <- do.call(rbind, rows)
  out$p_bh <- p.adjust(out$p_raw, method = "BH")
  out
}

# CORRECTION (14-Sep-2026): Fig_4c_*.csv include a Ryr1 column, but the manuscript's
# 20-gene correlation panel deliberately excludes Ryr1 (Methods: "excluding Ryr1,
# undetected across all conditions"; also absent from Fig 3e's gene list and from
# Fig_3g_..._exclude_Ryr1.svg). Dropping Ryr1 here reproduces the manuscript's existing
# per-scheme significant-pair counts (169/118/98/80) and 27-pair consensus exactly --
# confirming Ryr1's presence in the raw CSVs is a pipeline artefact, not intended input.
DROP_GENES <- c("Ryr1")

results <- list()
for (label in names(files)) {
  df <- read.csv(file.path(BASE, files[[label]]))
  df <- df[, !(names(df) %in% DROP_GENES), drop = FALSE]
  if (label == "Raw") {
    sanity_check_group_order(df)
  }
  results[[label]] <- list(
    pooled   = pooled_correlations(df),
    adjusted = treatment_adjusted_correlations(df)
  )
}

sig_pairs <- function(df, pcol, thresh = 0.05) {
  sub <- df[df[[pcol]] < thresh, ]
  pairs <- mapply(function(a, b) paste(sort(c(a, b)), collapse = "||"),
                   sub$gene1, sub$gene2)
  unique(pairs)
}

cat(strrep("=", 70), "\n", sep = "")
cat("PER-METHOD SIGNIFICANT PAIR COUNTS\n")
cat(strrep("=", 70), "\n", sep = "")
for (scheme in c("pooled", "adjusted")) {
  for (pc in list(list(col = "p_raw", label = "uncorrected p<0.05"),
                   list(col = "p_bh", label = "BH-corrected q<0.05"))) {
    counts <- sapply(names(files), function(label) length(sig_pairs(results[[label]][[scheme]], pc$col)))
    cat(sprintf("%-10s | %-22s | %s\n", scheme, pc$label,
                paste(sprintf("%s=%d", names(counts), counts), collapse = ", ")))
  }
}

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("CONSENSUS PAIR COUNTS (intersection across all 4 normalization schemes)\n")
cat(strrep("=", 70), "\n", sep = "")

consensus <- list()
for (scheme in c("pooled", "adjusted")) {
  for (pc in list(list(col = "p_raw", label = "uncorrected"),
                   list(col = "p_bh", label = "BH-corrected"))) {
    sets <- lapply(names(files), function(label) sig_pairs(results[[label]][[scheme]], pc$col))
    cons <- Reduce(intersect, sets)
    key <- paste(scheme, pc$label, sep = "||")
    consensus[[key]] <- cons
    cat(sprintf("%-10s x %-13s -> %d consensus pairs\n", scheme, pc$label, length(cons)))
  }
}

pair_to_tuple_str <- function(p) {
  parts <- strsplit(p, "\\|\\|")[[1]]
  sprintf("('%s', '%s')", parts[1], parts[2])
}

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("REPRODUCE CURRENT MANUSCRIPT METHOD (pooled, uncorrected) -- sanity check vs code's '28'\n")
cat(strrep("=", 70), "\n", sep = "")
old_pairs <- sort(consensus[["pooled||uncorrected"]])
cat("[", paste(sapply(old_pairs, pair_to_tuple_str), collapse = ", "), "]\n", sep = "")

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("CORRECTLY-SPECIFIED CONSENSUS (treatment-adjusted, BH-corrected)\n")
cat(strrep("=", 70), "\n", sep = "")
final <- consensus[["adjusted||BH-corrected"]]
cat(sprintf("N = %d\n", length(final)))
for (p in sort(final)) {
  parts <- strsplit(p, "\\|\\|")[[1]]
  cat(sprintf("  ('%s', '%s')\n", parts[1], parts[2]))
}

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("OVERLAP: how many of the current 28 survive proper treatment-adjustment + BH?\n")
cat(strrep("=", 70), "\n", sep = "")
old <- consensus[["pooled||uncorrected"]]
new <- consensus[["adjusted||BH-corrected"]]
cat(sprintf("Old (pooled, uncorrected): %d\n", length(old)))
cat(sprintf("New (adjusted, BH): %d\n", length(new)))
cat(sprintf("Overlap: %d\n", length(intersect(old, new))))
cat(sprintf("Lost (in old, not new): %d\n", length(setdiff(old, new))))
cat(sprintf("Gained (in new, not old): %d\n", length(setdiff(new, old))))
