# Gene Correlation Matrix - Figure 4b
 
rm(list = ls())

library(ggplot2)
library(readxl)
library(ggcorrplot)
library(dplyr)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Load data
dat <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx", sheet = 1)

# Calculate correlations
corr <- round(cor(dat, use = "pairwise.complete.obs"), 2)
p.mat <- round(cor_pmat(dat), 2)

# Extract significant correlations only
sig_corr <- corr
sig_corr[p.mat > 0.05] <- NA

# Export
write.csv(corr, "Fig_4b_full_correlation_matrix.csv", row.names = TRUE)
write.csv(p.mat, "Fig_4b_pvalue_matrix.csv", row.names = TRUE)
write.csv(sig_corr, "Fig_4b_significant_correlations.csv", row.names = TRUE)

# Plot
panel_b <- ggcorrplot(sig_corr,
                      hc.order = FALSE,
                      outline.col = "white",
                      type = "full",
                      show.legend = TRUE) +
  scale_fill_gradientn(colors = c(npg_colors[4], "white", npg_colors[1]),
                       values = scales::rescale(c(-1, 0, 1)),
                       limits = c(-1, 1), 
                       name = "Pearson correlation\ncoefficient (r)",
                       na.value = "grey85") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
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

ggsave("Fig_4b_Correlation_matrix.png", plot = panel_b,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

ggsave("Fig_4b_Correlation_matrix.svg", plot = panel_b,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

