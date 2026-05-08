################################################################################
# FIGURE 4e - HES1-NORMALIZED CORRELATION MATRIX
# SIGNIFICANT PAIRS ONLY + ANNOTATE 28 CONSENSUS PAIRS WITH ASTERISKS
# Revised: Hmisc removed; Ryr1 excluded (QC failure); cor/cor_pmat used only
################################################################################

rm(list = ls())

library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

################################################################################
# LOAD DATA
################################################################################

cat("\nLoading Hes1-normalized data...\n")
hes1_norm <- read.csv("Fig_4c_Hes1_normalized_data.csv")

# Safety: drop Ryr1 if present
hes1_norm <- hes1_norm[, colnames(hes1_norm) != "Ryr1"]

cat("Genes:", ncol(hes1_norm), "\n")
cat("Samples:", nrow(hes1_norm), "\n")

cat("\nLoading consensus pairs...\n")
consensus_df    <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
consensus_pairs <- consensus_df$Consensus_Pairs
cat("Consensus pairs loaded:", length(consensus_pairs), "\n")
print(consensus_pairs)

################################################################################
# CALCULATE CORRELATIONS
################################################################################

cat("\nCalculating correlation matrix...\n")

corr_e <- cor(hes1_norm, use = "pairwise.complete.obs")
p_mat_e <- cor_pmat(hes1_norm)

cat("Correlation matrix dimensions:", dim(corr_e), "\n")

################################################################################
# FILTER TO SIGNIFICANT CORRELATIONS ONLY
################################################################################

cat("\nFiltering significant correlations...\n")

sig_corr_e <- corr_e
sig_corr_e[p_mat_e >= 0.05] <- NA

n_sig <- sum(p_mat_e < 0.05, na.rm = TRUE) / 2
cat("Significant pairs (p < 0.05):", n_sig, "\n")

################################################################################
# CREATE CONSENSUS PAIR ANNOTATION MATRIX
################################################################################

cat("\nCreating consensus pair annotations...\n")

genes <- rownames(corr_e)
consensus_labels <- matrix("", nrow = nrow(corr_e), ncol = ncol(corr_e))
rownames(consensus_labels) <- genes
colnames(consensus_labels) <- genes

for (pair in consensus_pairs) {
  pair_genes <- strsplit(pair, "_")[[1]]
  gene1 <- pair_genes[1]
  gene2 <- pair_genes[2]
  idx1  <- which(genes == gene1)
  idx2  <- which(genes == gene2)
  if (length(idx1) > 0 && length(idx2) > 0) {
    consensus_labels[idx1, idx2] <- "*"
    consensus_labels[idx2, idx1] <- "*"
  }
}

n_consensus_annotated <- sum(consensus_labels == "*") / 2
cat("Consensus pairs annotated:", n_consensus_annotated, "\n")

################################################################################
# EXPORT DATA
################################################################################

cat("\nExporting data...\n")

write.csv(as.data.frame(sig_corr_e),
          "Fig_4e_Hes1_significant_correlations_only.csv", row.names = TRUE)

consensus_matrix_df <- as.data.frame(consensus_labels)
write.csv(consensus_matrix_df,
          "Fig_4e_consensus_annotations.csv", row.names = TRUE)

cat("Exported: Fig_4e_Hes1_significant_correlations_only.csv\n")
cat("Exported: Fig_4e_consensus_annotations.csv\n")

################################################################################
# VERIFY KEY PAIRS
################################################################################

cat("\n=== VERIFICATION: Sample correlation values ===\n")
test_pairs <- c("Cacna1d_Wnt5a", "Cacna1d_Ptbp1", "Ptbp1_Wnt5a")
for (pair in test_pairs) {
  pair_genes <- strsplit(pair, "_")[[1]]
  if (pair_genes[1] %in% genes && pair_genes[2] %in% genes) {
    val   <- corr_e[pair_genes[1], pair_genes[2]]
    p_val <- p_mat_e[pair_genes[1], pair_genes[2]]
    cat(sprintf("%s: r = %.10f (p = %.4e)\n", pair, val, p_val))
  }
}
cat("===============================================\n\n")

################################################################################
# PREPARE PLOT DATA
################################################################################

cat("\nPreparing plot data...\n")

corr_e_display          <- round(corr_e, 2)
sig_corr_e_display      <- corr_e_display
sig_corr_e_display[p_mat_e >= 0.05] <- NA

corr_long     <- melt(sig_corr_e_display)
consensus_long <- melt(consensus_labels)

plot_data <- data.frame(
  Gene1       = corr_long$Var1,
  Gene2       = corr_long$Var2,
  Correlation = corr_long$value,
  Consensus   = consensus_long$value
)
plot_data <- plot_data[!is.na(plot_data$Correlation), ]

cat("Total data points to plot:", nrow(plot_data), "\n")

################################################################################
# CREATE PLOT
################################################################################

cat("\nGenerating Figure 4e...\n")

panel_e <- ggplot(plot_data, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = Correlation), color = "white", linewidth = 0.2) +
  geom_text(aes(label = Consensus),
            size = 3, color = "black", fontface = "bold", family = "Arial") +
  scale_fill_gradient2(
    low      = npg_colors[4],
    mid      = "white",
    high     = npg_colors[1],
    midpoint = 0,
    limits   = c(-1, 1),
    na.value = "grey85",
    name     = "Pearson correlation\ncoefficient (r)"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    legend.position      = "top",
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.spacing   = margin(5),
    legend.justification = c("right", "top"),
    legend.text          = element_text(size = 8, family = "Arial"),
    legend.title         = element_text(size = 8, colour = "black", angle = 0,
                                        vjust = 0.95, hjust = 1, family = "Arial"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    axis.text.x          = element_text(size = 8, colour = "black", face = "italic",
                                        angle = 90, hjust = 1, vjust = 0.5,
                                        family = "Arial"),
    axis.text.y          = element_text(size = 8, colour = "black", face = "italic",
                                        family = "Arial"),
    panel.grid           = element_blank()
  )

################################################################################
# SAVE FIGURE
################################################################################

cat("\nSaving figure...\n")

ggsave("Fig_4e_Significant_pairs_consensus_annotated.png", plot = panel_e,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)
ggsave("Fig_4e_Significant_pairs_consensus_annotated.svg", plot = panel_e,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

################################################################################
# SUMMARY
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("FIGURE 4e COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")
cat("Data source    : Hes1-normalized gene expression\n")
cat("Genes analysed :", ncol(hes1_norm), "(Ryr1 excluded - QC failure)\n")
cat("Total pairs    :", nrow(corr_e) * (nrow(corr_e) - 1) / 2, "\n")
cat("Significant    :", n_sig, "(p < 0.05)\n")
cat("Consensus      :", n_consensus_annotated, "(all 4 methods)\n")
cat("% consensus/sig:",
    round(100 * n_consensus_annotated / n_sig, 1), "%\n\n")
cat("Exports:\n")
cat("  Fig_4e_Hes1_significant_correlations_only.csv (full precision)\n")
cat("  Fig_4e_consensus_annotations.csv\n")
cat("  Fig_4e_Significant_pairs_consensus_annotated.png (600 DPI)\n")
cat("  Fig_4e_Significant_pairs_consensus_annotated.svg\n\n")
cat("NOTE: CSV exported at full precision; figure displays rounded (2 dp).\n")
cat(rep("=", 60), "\n\n", sep = "")

################################################################################
# END OF SCRIPT
################################################################################

