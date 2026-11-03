################################################################################
# FIGURE 4c - NORMALIZATION METHODS COMPARISON BAR CHART
# Author: Revised for publication
# Date: January 2026
# Color Palette: NPG (Nature Publishing Group)
# Font: Arial
################################################################################

rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")

# Load required libraries
library(ggplot2)
library(readxl)
library(reshape2)
library(Hmisc)
library(ggsci)

# Extract NPG color palette
npg_colors <- pal_npg("nrc")(10)

################################################################################
# LOAD DATA
################################################################################

cat("\nLoading gene expression data...\n")
gene_expr_data <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx", sheet = 1)

################################################################################
# CORRELATION ANALYSIS FUNCTION
################################################################################

perform_correlation_analysis <- function(data, method_name) {
  data_clean <- na.omit(data)
  cor_results <- rcorr(as.matrix(data_clean), type = "pearson")
  cor_matrix <- cor_results$r
  p_matrix <- cor_results$P
  
  genes <- colnames(data_clean)
  n_genes <- length(genes)
  
  sig_count <- 0
  strong_count <- 0
  very_strong_count <- 0
  
  for (i in 1:(n_genes-1)) {
    for (j in (i+1):n_genes) {
      p_val <- p_matrix[i, j]
      r_val <- abs(cor_matrix[i, j])
      
      if (!is.na(p_val) && p_val <= 0.05) {
        sig_count <- sig_count + 1
        if (r_val >= 0.7) strong_count <- strong_count + 1
        if (r_val >= 0.8) very_strong_count <- very_strong_count + 1
      }
    }
  }
  
  return(data.frame(
    Method = method_name,
    Significant = sig_count,
    Strong = strong_count,
    VeryStrong = very_strong_count
  ))
}

################################################################################
# PERFORM NORMALIZATIONS AND ANALYSES
################################################################################

cat("\nPerforming normalization analyses...\n")

results_list <- list()

# 1. Hes1-normalized
cat("  - Hes1-normalized...\n")
hes1_norm <- gene_expr_data
for (col in colnames(hes1_norm)) {
  if (col != "Hes1") hes1_norm[[col]] <- hes1_norm[[col]] / hes1_norm$Hes1
}
hes1_norm <- hes1_norm[, colnames(hes1_norm) != "Hes1"]
results_list[[1]] <- perform_correlation_analysis(hes1_norm, "Hes1-\nnormalized")
write.csv(hes1_norm, "Fig_4c_Hes1_normalized_data.csv", row.names = FALSE)

# 2. Kcnj11-normalized
cat("  - Kcnj11-normalized...\n")
kcnj11_norm <- gene_expr_data
for (col in colnames(kcnj11_norm)) {
  if (col != "Kcnj11") kcnj11_norm[[col]] <- kcnj11_norm[[col]] / kcnj11_norm$Kcnj11
}
kcnj11_norm <- kcnj11_norm[, colnames(kcnj11_norm) != "Kcnj11"]
results_list[[2]] <- perform_correlation_analysis(kcnj11_norm, "Kcnj11-\nNormalized")
write.csv(kcnj11_norm, "Fig_4c_Kcnj11_normalized_data.csv", row.names = FALSE)

# 3. Ins1-normalized
cat("  - Ins1-normalized...\n")
ins1_norm <- gene_expr_data
for (col in colnames(ins1_norm)) {
  if (col != "Ins1") ins1_norm[[col]] <- ins1_norm[[col]] / ins1_norm$Ins1
}
ins1_norm <- ins1_norm[, colnames(ins1_norm) != "Ins1"]
results_list[[3]] <- perform_correlation_analysis(ins1_norm, "Ins1-\nnormalized")
write.csv(ins1_norm, "Fig_4c_Ins1_normalized_data.csv", row.names = FALSE)

# 4. Raw expression
cat("  - Raw expression...\n")
results_list[[4]] <- perform_correlation_analysis(gene_expr_data, "Raw\nexpression")
write.csv(gene_expr_data, "Fig_4c_Raw_expression_data.csv", row.names = FALSE)

# Combine results
comparison_data <- do.call(rbind, results_list)

################################################################################
# EXPORT COMPARISON DATA
################################################################################

cat("\nExporting comparison data...\n")

write.csv(comparison_data, "Fig_4c_normalization_comparison_data.csv", row.names = FALSE)

################################################################################
# CREATE BAR CHART
################################################################################

cat("\nGenerating Panel (c): Bar Chart...\n")

# Reshape for plotting
comparison_long <- melt(comparison_data, id.vars = "Method")
comparison_long$Method <- factor(comparison_long$Method, 
                                 levels = comparison_data$Method)

# Create bar chart with NPG colors
panel_c <- ggplot(comparison_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = value), 
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 2.5, family = "Arial") +
  scale_fill_manual(
    values = c(npg_colors[4], npg_colors[2], npg_colors[3]),
    labels = expression(atop("Significant", "(" * italic(p) * " ≤ 0.05)"),
                        atop("Strong", "(|" * italic(r) * "| ≥ 0.7)"),
                        atop("Very Strong", "(|" * italic(r) * "| ≥ 0.8)"))
  ) +
  # x labels in the order of the factor levels: Hes1, Kcnj11, Ins1, raw (gene symbols in italic)
  scale_x_discrete(labels = expression(atop(italic("Hes1") * "-", "normalized"),
                                       atop(italic("Kcnj11") * "-", "normalized"),
                                       atop(italic("Ins1") * "-", "normalized"),
                                       atop("Raw", "expression"))) +
  labs(x = NULL, y = "Value", fill = NULL) +  # legend title "Variable" dropped: the legend was wider than the panel and "Very Strong" was cut off
  theme_minimal() +
  theme(legend.position = "top",
    axis.text.x = element_text(size = 8, color = "black", hjust = 0.5, family = "Arial"),
    axis.text.y = element_text(size = 8, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 8, face = "bold", family = "Arial"),
    legend.key.width  = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.margin=margin(0,0,0,0),
    legend.box.spacing = margin(5),
    legend.text = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    #legend.key.size = unit(0.4, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

################################################################################
# SAVE FIGURE
################################################################################

cat("\nSaving figure...\n")

ggsave("Fig_4c_Normalization_comparison.png", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)
ggsave("Fig_4c_Normalization_comparison.svg", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)

################################################################################
# SUMMARY
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("PANEL 4c COMPLETE!\n")
cat(rep("=", 80), "\n\n", sep="")
cat("Exported CSV files:\n")
cat(" - Fig_4c_normalization_comparison_data.csv (summary statistics)\n")
cat(" - Fig_4c_Hes1_normalized_data.csv\n")
cat(" - Fig_4c_Kcnj11_normalized_data.csv\n")
cat(" - Fig_4c_Ins1_normalized_data.csv\n")
cat(" - Fig_4c_Raw_expression_data.csv\n\n")
cat("Generated figures:\n")
cat(" - Fig_4c_Normalization_comparison.png\n")
cat(" - Fig_4c_Normalization_comparison.svg\n\n")
cat("Summary of correlations:\n")
print(comparison_data)
cat("\n")

################################################################################
# END OF SCRIPT
################################################################################
