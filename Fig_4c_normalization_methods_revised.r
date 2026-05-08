# Normalization Methods Comparison - Figure 4c

rm(list = ls())

library(ggplot2)
library(readxl)
library(reshape2)
library(Hmisc)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

gene_expr_data <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx",
                             sheet = 1)

# ── Helper: count significant correlations (vectorised) ─────────────────────

perform_correlation_analysis <- function(data, method_name) {
  data_clean <- na.omit(data)
  cor_res    <- rcorr(as.matrix(data_clean), type = "pearson")
  r_upper    <- cor_res$r[upper.tri(cor_res$r)]
  p_upper    <- cor_res$P[upper.tri(cor_res$P)]
  
  sig_mask <- !is.na(p_upper) & p_upper < 0.05
  r_abs    <- abs(r_upper)
  
  data.frame(
    Method      = method_name,
    Significant = sum(sig_mask),
    Strong      = sum(sig_mask & r_abs > 0.7),
    VeryStrong  = sum(sig_mask & r_abs > 0.8)
  )
}

# ── Normalization helper ─────────────────────────────────────────────────────

normalize_by <- function(data, ref_gene, label) {
  norm        <- data
  norm[]      <- lapply(names(data), function(col) {
    if (col == ref_gene) data[[col]] else data[[col]] / data[[ref_gene]]
  })
  norm[[ref_gene]] <- NULL
  list(data = norm, label = label)
}

# ── Run all four methods ─────────────────────────────────────────────────────

hes1   <- normalize_by(gene_expr_data, "Hes1",   "Hes1-\nnormalized")
kcnj11 <- normalize_by(gene_expr_data, "Kcnj11", "Kcnj11-\nnormalized")
ins1   <- normalize_by(gene_expr_data, "Ins1",   "Ins1-\nnormalized")

write.csv(hes1$data,   "Fig_4c_Hes1_normalized_data.csv",   row.names = FALSE)
write.csv(kcnj11$data, "Fig_4c_Kcnj11_normalized_data.csv", row.names = FALSE)
write.csv(ins1$data,   "Fig_4c_Ins1_normalized_data.csv",   row.names = FALSE)
write.csv(gene_expr_data, "Fig_4c_Raw_expression_data.csv", row.names = FALSE)

results_list <- list(
  perform_correlation_analysis(hes1$data,       hes1$label),
  perform_correlation_analysis(kcnj11$data,     kcnj11$label),
  perform_correlation_analysis(ins1$data,       ins1$label),
  perform_correlation_analysis(gene_expr_data,  "Raw\nexpression")
)

comparison_data <- do.call(rbind, results_list)
write.csv(comparison_data, "Fig_4c_normalization_comparison_data.csv", row.names = FALSE)

# ── Permutation null: shuffle each column INDEPENDENTLY ─────────────────────
# FIX: row-shuffle is order-invariant for Pearson correlation (SD = 0 bug).
# Each gene must be shuffled separately to break all gene-gene relationships.

set.seed(42)
n_perm    <- 1000
n_samples <- nrow(gene_expr_data)
perm_sig  <- numeric(n_perm)

cat("Running permutation null...\n")
for (perm in seq_len(n_perm)) {
  if (perm %% 100 == 0) cat(sprintf("  %d / %d\n", perm, n_perm))
  
  # Independent column-wise shuffle — the only valid approach
  shuffled    <- as.data.frame(lapply(gene_expr_data, sample))
  perm_result <- perform_correlation_analysis(shuffled, "perm")
  perm_sig[perm] <- perm_result$Significant
}

null_mean <- round(mean(perm_sig), 1)
null_q95  <- round(quantile(perm_sig, 0.95), 1)

# ── Plot ─────────────────────────────────────────────────────────────────────

comparison_long        <- melt(comparison_data, id.vars = "Method")
comparison_long$Method <- factor(comparison_long$Method,
                                 levels = comparison_data$Method)

panel_c <- ggplot(comparison_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = value),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 2.5, family = "Arial") +
  geom_hline(yintercept = null_mean,
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("text", x = 0.6, y = null_mean + 1.5,
           label = paste0("Permutation null\n(mean = ", null_mean,
                          ", 95th = ", null_q95, ")"),   # FIX: \n not \\n
           hjust = 0, size = 2.2, color = "grey40", family = "Arial") +
  scale_fill_manual(
    values = c(npg_colors[4], npg_colors[2], npg_colors[3]),
    labels = c("Significant\n(p<0.05)",        # FIX: \n not \\n
               "Strong\n(|r|>0.7)",
               "Very Strong\n(|r|>0.8)")
  ) +
  labs(x = NULL, y = "Count", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position      = "top",
    axis.text.x          = element_text(size = 8, color = "black", family = "Arial"),
    axis.text.y          = element_text(size = 8, color = "black", family = "Arial"),
    axis.title.y         = element_text(size = 8, face = "bold", family = "Arial"),
    legend.text          = element_text(size = 8, family = "Arial"),
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("Fig_4c_Normalization_comparison.png", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)
ggsave("Fig_4c_Normalization_comparison.svg", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)

# ── Export permutation null stats ────────────────────────────────────────────

perm_null_stats <- data.frame(
  Metric = c("Mean_significant_under_null", "Median_significant_under_null",
             "Q95_significant_under_null",  "SD_significant_under_null",
             "N_permutations"),
  Value  = c(null_mean,
             round(median(perm_sig), 1),
             null_q95,
             round(sd(perm_sig), 2),
             n_perm)
)
write.csv(perm_null_stats, "Fig_4c_permutation_null_stats.csv", row.names = FALSE)

cat("\n=== Permutation Null Summary ===\n")
print(perm_null_stats)

