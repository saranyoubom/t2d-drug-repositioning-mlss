################################################################################
# FIGURE 4e - HES1-NORMALIZED CORRELATION MATRIX
# SIGNIFICANT PAIRS ONLY + ANNOTATE 28 CONSENSUS PAIRS WITH ASTERISKS
# Author: Revised for publication
# Date: January 2026
# Color Palette: NPG (Nature Publishing Group)
# Font: Arial
# Dependencies: Uses Hes1-normalized data CSV from Fig_4c_script.R
#               Uses consensus pairs CSV from Fig_4d_script.R
################################################################################

rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")

# Load required libraries
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(ggsci)

# Extract NPG color palette
npg_colors <- pal_npg("nrc")(10)

################################################################################
# LOAD HES1-NORMALIZED DATA FROM CSV
################################################################################

cat("\nLoading Hes1-normalized data from CSV file...\n")
hes1_norm <- read.csv("Fig_4c_Hes1_normalized_data.csv")

cat("Data loaded successfully!\n")
cat("Genes:", ncol(hes1_norm), "\n")
cat("Samples:", nrow(hes1_norm), "\n")

################################################################################
# LOAD CONSENSUS PAIRS FROM CSV
################################################################################

cat("\nLoading consensus pairs from CSV file...\n")
consensus_df <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
consensus_pairs <- consensus_df$Consensus_Pairs

cat("Consensus pairs loaded:", length(consensus_pairs), "\n")
print(consensus_pairs)

################################################################################
# CALCULATE CORRELATIONS
################################################################################

cat("\nCalculating correlation matrix...\n")

# Calculate correlation matrix WITH MORE PRECISION
corr_e <- cor(hes1_norm, use = "pairwise.complete.obs")

# Compute p-values
p_mat_e <- cor_pmat(hes1_norm)

cat("Correlation matrix dimensions:", dim(corr_e), "\n")

################################################################################
# FILTER TO SHOW ONLY SIGNIFICANT CORRELATIONS
################################################################################

cat("\nFiltering significant correlations...\n")

# Create filtered correlation matrix - only significant pairs
sig_corr_e <- corr_e
sig_corr_e[p_mat_e >= 0.05] <- NA  # Remove non-significant

# Count significant correlations
n_sig <- sum(p_mat_e < 0.05, na.rm = TRUE) / 2
cat("Significant pairs (p < 0.05):", n_sig, "\n")

################################################################################
# CREATE CONSENSUS PAIR ANNOTATIONS (ASTERISKS)
################################################################################

cat("\nCreating consensus pair annotations...\n")

# Create annotation matrix for consensus pairs only
consensus_labels <- matrix("", nrow = nrow(corr_e), ncol = ncol(corr_e))
genes <- rownames(corr_e)

# Mark consensus pairs with asterisk
for (pair in consensus_pairs) {
  pair_genes <- strsplit(pair, "_")[[1]]
  gene1 <- pair_genes[1]
  gene2 <- pair_genes[2]
  
  # Find indices
  idx1 <- which(genes == gene1)
  idx2 <- which(genes == gene2)
  
  if (length(idx1) > 0 && length(idx2) > 0) {
    # Mark both upper and lower triangle
    consensus_labels[idx1, idx2] <- "*"
    consensus_labels[idx2, idx1] <- "*"
  }
}

# Count how many consensus pairs were annotated
n_consensus_annotated <- sum(consensus_labels == "*") / 2
cat("Consensus pairs annotated:", n_consensus_annotated, "\n")

################################################################################
# EXPORT DATA - CRITICAL FOR DOWNSTREAM ANALYSIS
################################################################################

cat("\nExporting data...\n")

# Export significant correlations with FULL PRECISION (not rounded to 2 decimals)
sig_corr_df <- as.data.frame(sig_corr_e)
write.csv(sig_corr_df, "Fig_4e_Hes1_significant_correlations_only.csv", row.names = TRUE)

# Export consensus annotations
consensus_matrix_df <- as.data.frame(consensus_labels)
rownames(consensus_matrix_df) <- rownames(corr_e)
colnames(consensus_matrix_df) <- colnames(corr_e)
write.csv(consensus_matrix_df, "Fig_4e_consensus_annotations.csv", row.names = TRUE)

cat("✓ Exported: Fig_4e_Hes1_significant_correlations_only.csv (FULL PRECISION)\n")
cat("✓ Exported: Fig_4e_consensus_annotations.csv\n")

################################################################################
# VERIFY EXPORTED VALUES
################################################################################

cat("\n=== VERIFICATION: Sample correlation values ===\n")
test_pairs <- c("Cacna1d_Wnt5a", "Cacna1d_Ptbp1", "Ptbp1_Wnt5a")
for (pair in test_pairs) {
  pair_genes <- strsplit(pair, "_")[[1]]
  if (pair_genes[1] %in% genes && pair_genes[2] %in% genes) {
    val <- corr_e[pair_genes[1], pair_genes[2]]
    p_val <- p_mat_e[pair_genes[1], pair_genes[2]]
    cat(sprintf("%s: r = %.10f (p = %.4e)\n", pair, val, p_val))
  }
}
cat("===============================================\n\n")

################################################################################
# CREATE PLOT - SIGNIFICANT PAIRS ONLY + CONSENSUS ASTERISKS
################################################################################

cat("\nGenerating Panel (e): Significant Pairs with Consensus Annotations...\n")

# Round for DISPLAY ONLY (not for export)
corr_e_display <- round(corr_e, 2)
sig_corr_e_display <- corr_e_display
sig_corr_e_display[p_mat_e >= 0.05] <- NA

# Prepare data for ggplot - ONLY SIGNIFICANT CORRELATIONS
corr_long <- melt(sig_corr_e_display)
consensus_long <- melt(consensus_labels)

plot_data <- data.frame(
  Gene1 = corr_long$Var1,
  Gene2 = corr_long$Var2,
  Correlation = corr_long$value,
  Consensus = consensus_long$value
)

# Remove NA rows (non-significant correlations)
plot_data <- plot_data[!is.na(plot_data$Correlation), ]

cat("Total data points to plot:", nrow(plot_data), "\n")

# Create the plot with NPG colors
# Only significant correlations shown, asterisks mark consensus pairs
panel_e <- ggplot(plot_data, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = Correlation), color = "white", size = 0.5) +
  # FIXED: Use 'Consensus' variable to show asterisks
  geom_text(aes(label = Consensus),
            size = 3, color = "black", fontface = "bold", family = "Arial") +
  scale_fill_gradient2(
    low = npg_colors[4],      # Blue for negative
    mid = "white",            # White for zero
    high = npg_colors[1],     # Red/Orange for positive
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",      # Grey for non-significant
    name = expression(atop("Pearson correlation", "coefficient (" * italic(r) * ")"))
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = margin(5),
        legend.justification = c("right", "top"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8, family = "Arial"),
        legend.title = element_text(size = 8, colour = "black", angle = 0, 
                                    vjust = 0.95, hjust = 1, family = "Arial"),
        axis.text.x = element_text(size = 8, colour = "black", face = "italic", 
                                   angle = 90, hjust = 1, vjust = 0.5, family = "Arial"),
        axis.text.y = element_text(size = 8, colour = "black", face = "italic", 
                                   family = "Arial"),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

panel_e

################################################################################
# SAVE FIGURE
################################################################################

cat("\nSaving figure...\n")

ggsave("Fig_4e_Significant_pairs_consensus_annotated.png", plot = panel_e,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)
ggsave("Fig_4e_Significant_pairs_consensus_annotated.svg", plot = panel_e,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

cat("✓ Fig_4e_Significant_pairs_consensus_annotated.png (600 DPI)\n")
cat("✓ Fig_4e_Significant_pairs_consensus_annotated.svg\n")

################################################################################
# SUMMARY
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("PANEL 4e COMPLETE!\n")
cat(rep("=", 80), "\n\n", sep="")
cat("Data source: Hes1-normalized gene expression\n")
cat("Display: ONLY significant correlations (p < 0.05) shown with color\n")
cat("Annotations: Asterisks (*) mark the 28 consensus pairs\n\n")
cat("Consensus pairs definition:\n")
cat("  Pairs that are significant across ALL 4 normalization methods:\n")
cat("  - Hes1-normalized\n")
cat("  - Kcnj11-normalized\n")
cat("  - Ins1-normalized\n")
cat("  - Raw expression\n\n")
cat("Color gradient:\n")
cat(sprintf("  Blue (%s): Negative correlations (r < 0)\n", npg_colors[4]))
cat("  White: Zero/No correlation (r ≈ 0)\n")
cat(sprintf("  Red/Orange (%s): Positive correlations (r > 0)\n", npg_colors[1]))
cat("  Grey: Non-significant (p ≥ 0.05) - not displayed\n\n")
cat("Exported CSV files:\n")
cat(" ✓ Fig_4e_Hes1_significant_correlations_only.csv (FULL PRECISION)\n")
cat(" ✓ Fig_4e_consensus_annotations.csv\n\n")
cat("Generated figures:\n")
cat(" ✓ Fig_4e_Significant_pairs_consensus_annotated.png (600 DPI)\n")
cat(" ✓ Fig_4e_Significant_pairs_consensus_annotated.svg\n\n")
cat("Statistics:\n")
cat(" Total gene pairs:", nrow(corr_e) * (nrow(corr_e) - 1) / 2, "\n")
cat(" Significant pairs (p < 0.05):", n_sig, "\n")
cat(" Consensus pairs annotated:", n_consensus_annotated, "\n")
cat(" Percentage consensus of significant:", 
    round(100 * n_consensus_annotated / n_sig, 1), "%\n")
cat("\n")

cat("IMPORTANT: CSV exported with FULL PRECISION for downstream analysis!\n")
cat("Figure display uses rounded values (2 decimals) for readability.\n\n")

################################################################################
# END OF SCRIPT
################################################################################