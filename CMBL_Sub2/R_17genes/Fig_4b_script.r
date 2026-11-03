################################################################################
# FIGURE 4b - CORRELATION MATRIX (SIGNIFICANT ONLY)
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
library(ggcorrplot)
library(dplyr)
library(ggsci)

# Extract NPG color palette
npg_colors <- pal_npg("nrc")(10)
cat("NPG Color Palette:\n")
print(npg_colors)

################################################################################
# LOAD DATA
################################################################################

cat("\nLoading gene expression data...\n")
dat <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx", sheet = 1)

################################################################################
# CALCULATE CORRELATIONS
################################################################################

cat("\nCalculating correlation matrix...\n")

# Calculate correlation matrix
corr <- round(cor(dat, use = "pairwise.complete.obs"), 4)

# Compute p-values
p.mat <- round(cor_pmat(dat), 4)

# Extract significant correlations
sig_corr <- corr
sig_corr[p.mat > 0.05] <- NA

################################################################################
# EXPORT DATA
################################################################################

cat("\nExporting data...\n")

# Export full correlation matrix
write.csv(corr, "Fig_4b_full_correlation_matrix.csv", row.names = TRUE)

# Export p-value matrix
write.csv(p.mat, "Fig_4b_pvalue_matrix.csv", row.names = TRUE)

# Export significant correlations only
sig_corr_text <- as.data.frame(sig_corr)
write.csv(sig_corr_text, "Fig_4b_significant_correlations.csv", row.names = TRUE)

################################################################################
# CREATE PLOT
################################################################################

cat("\nGenerating Panel (b): Correlation Matrix...\n")

# Plot correlation matrix with NPG colors
panel_b <- ggcorrplot(sig_corr,
                      hc.order = FALSE,
                      outline.col = "grey",
                      type = "full",
                      show.legend = TRUE) +
  scale_fill_gradientn(colors = c(npg_colors[4], "white", npg_colors[1]),
                       values = scales::rescale(c(-1, 0, 1)),
                       limits = c(-1, 1),
                       name = expression(atop("Pearson correlation", "coefficient (" * italic(r) * ")")),
                       na.value = "grey85") +
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

################################################################################
# SAVE FIGURE
################################################################################

cat("\nSaving figure...\n")

ggsave("Fig_4b_Correlation_matrix.png", plot = panel_b,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)
ggsave("Fig_4b_Correlation_matrix.svg", plot = panel_b,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

################################################################################
# SUMMARY
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("PANEL 4b COMPLETE!\n")
cat(rep("=", 80), "\n\n", sep="")
cat("Exported CSV files:\n")
cat(" - Fig_4b_full_correlation_matrix.csv\n")
cat(" - Fig_4b_pvalue_matrix.csv\n")
cat(" - Fig_4b_significant_correlations.csv\n\n")
cat("Generated figures:\n")
cat(" - Fig_4b_Correlation_matrix.png\n")
cat(" - Fig_4b_Correlation_matrix.svg\n")
cat("\n")

################################################################################
# END OF SCRIPT
################################################################################
