# =============================================================================
# Normalization Methods Comparison - Figure 4c
# Revised: Hmisc replaced with base R cor.test(); Ryr1 excluded (QC failure)
# =============================================================================

rm(list = ls())

library(ggplot2)
library(readxl)
library(reshape2)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Load data ──────────────────────────────────────────────────────────────────
gene_expr_data <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx", sheet = 1)

# ── Exclude Ryr1 (non-specific melt curve; QC failure) ────────────────────────
gene_expr_data <- gene_expr_data[, colnames(gene_expr_data) != "Ryr1"]

# ── Correlation analysis function (base R, no Hmisc) ──────────────────────────
perform_correlation_analysis <- function(data, method_name) {
  data_clean <- as.matrix(na.omit(data))
  n_genes <- ncol(data_clean)
  sig_count <- 0
  strong_count <- 0
  very_strong_count <- 0
  
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      test  <- cor.test(data_clean[, i], data_clean[, j], method = "pearson")
      p_val <- test$p.value
      r_val <- abs(test$estimate)
      if (!is.na(p_val) && p_val < 0.05) {
        sig_count <- sig_count + 1
        if (r_val > 0.7) strong_count <- strong_count + 1
        if (r_val > 0.8) very_strong_count <- very_strong_count + 1
      }
    }
  }
  
  return(data.frame(
    Method    = method_name,
    Significant  = sig_count,
    Strong       = strong_count,
    VeryStrong   = very_strong_count
  ))
}

# ── Hes1-normalized ───────────────────────────────────────────────────────────
hes1_norm <- gene_expr_data
for (col in colnames(hes1_norm)) {
  if (col != "Hes1") hes1_norm[[col]] <- hes1_norm[[col]] / hes1_norm$Hes1
}
hes1_norm <- hes1_norm[, colnames(hes1_norm) != "Hes1"]
results_list <- list()
results_list[[1]] <- perform_correlation_analysis(hes1_norm, "Hes1-\nnormalized")
write.csv(hes1_norm, "Fig_4c_Hes1_normalized_data.csv", row.names = FALSE)

# ── Kcnj11-normalized ─────────────────────────────────────────────────────────
kcnj11_norm <- gene_expr_data
for (col in colnames(kcnj11_norm)) {
  if (col != "Kcnj11") kcnj11_norm[[col]] <- kcnj11_norm[[col]] / kcnj11_norm$Kcnj11
}
kcnj11_norm <- kcnj11_norm[, colnames(kcnj11_norm) != "Kcnj11"]
results_list[[2]] <- perform_correlation_analysis(kcnj11_norm, "Kcnj11-\nnormalized")
write.csv(kcnj11_norm, "Fig_4c_Kcnj11_normalized_data.csv", row.names = FALSE)

# ── Ins1-normalized ───────────────────────────────────────────────────────────
ins1_norm <- gene_expr_data
for (col in colnames(ins1_norm)) {
  if (col != "Ins1") ins1_norm[[col]] <- ins1_norm[[col]] / ins1_norm$Ins1
}
ins1_norm <- ins1_norm[, colnames(ins1_norm) != "Ins1"]
results_list[[3]] <- perform_correlation_analysis(ins1_norm, "Ins1-\nnormalized")
write.csv(ins1_norm, "Fig_4c_Ins1_normalized_data.csv", row.names = FALSE)

# ── Raw expression ────────────────────────────────────────────────────────────
results_list[[4]] <- perform_correlation_analysis(gene_expr_data, "Raw\nexpression")
write.csv(gene_expr_data, "Fig_4c_Raw_expression_data.csv", row.names = FALSE)

# ── Combine and export results ────────────────────────────────────────────────
comparison_data <- do.call(rbind, results_list)
write.csv(comparison_data, "Fig_4c_normalization_comparison_data.csv", row.names = FALSE)

# ── Reshape for plotting ──────────────────────────────────────────────────────
comparison_long <- melt(comparison_data, id.vars = "Method")
comparison_long$Method <- factor(comparison_long$Method, levels = comparison_data$Method)

# ── Plot ──────────────────────────────────────────────────────────────────────
panel_c <- ggplot(comparison_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = value),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 2.5, family = "Arial") +
  scale_fill_manual(
    values = c(npg_colors[4], npg_colors[2], npg_colors[3]),
    labels = c("Significant\n(p<0.05)", "Strong\n(|r|>0.7)", "Very Strong\n(|r|>0.8)")
  ) +
  labs(x = NULL, y = "Count", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position  = "top",
    axis.text.x      = element_text(size = 8, color = "black", family = "Arial"),
    axis.text.y      = element_text(size = 8, color = "black", family = "Arial"),
    axis.title.y     = element_text(size = 8, face = "bold", family = "Arial"),
    legend.text      = element_text(size = 8, family = "Arial"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# ── Save outputs ──────────────────────────────────────────────────────────────
ggsave("Fig_4c_Normalization_comparison.png", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)

ggsave("Fig_4c_Normalization_comparison.svg", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)

