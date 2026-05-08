# MLSS v4.0 Validation - Sensitivity and Robustness Analysis

rm(list = ls())
set.seed(42)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(reshape2)
library(gridExtra)

npg_colors <- pal_npg("nrc")(10)

# ── Load data ─────────────────────────────────────────────────────────────────
mlss_data <- read.csv("Fig_6c_mlss_all_combinations.csv", stringsAsFactors = FALSE)
top_15     <- mlss_data %>% head(15)

# ── Weight scenarios ──────────────────────────────────────────────────────────
weight_scenarios <- data.frame(
  Scenario = c("MLSS (Standard)", "C-Heavy (60:25:15)", "C-Heavy (55:20:25)",
               "B-Heavy (35:50:15)", "V-Heavy (35:30:35)", "P-Heavy (30:20:10:40)",
               "Balanced (40:40:10:10)", "Network-Focus (30:20:40:10)"),
  C_Weight = c(0.45, 0.60, 0.55, 0.35, 0.35, 0.30, 0.40, 0.30),
  B_Weight = c(0.30, 0.25, 0.20, 0.50, 0.30, 0.20, 0.40, 0.20),
  V_Weight = c(0.10, 0.10, 0.15, 0.10, 0.30, 0.10, 0.10, 0.40),
  P_Weight = c(0.15, 0.05, 0.10, 0.05, 0.05, 0.40, 0.10, 0.10),
  NPG_Color = npg_colors[1:8],
  stringsAsFactors = FALSE
)

# ── Recalculate scores for each scenario (top 15) ────────────────────────────
top_15_data    <- mlss_data %>% head(15)
scenario_results <- data.frame(Combination  = top_15_data$Combination,
                               Original_MLSS = top_15_data$MLSS_v4.0)

for (i in 1:nrow(weight_scenarios)) {
  scenario_name <- weight_scenarios$Scenario[i]
  c_w <- weight_scenarios$C_Weight[i]
  b_w <- weight_scenarios$B_Weight[i]
  v_w <- weight_scenarios$V_Weight[i]
  p_w <- weight_scenarios$P_Weight[i]
  
  new_base_score <- (c_w * top_15_data$Complementarity_C +
                       b_w * top_15_data$Balance_B +
                       v_w * top_15_data$Network_Coverage_V +
                       p_w * top_15_data$Potency_P) * 9
  
  new_mlss <- new_base_score + top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty
  scenario_results[[scenario_name]] <- new_mlss
}

# ── Ranking correlation analysis ──────────────────────────────────────────────
rank_data          <- apply(scenario_results[, -1], 2, rank)
rank_correlation   <- cor(rank_data, method = "spearman")
original_correlations <- rank_correlation[1, -1]   # correlations vs. original

# ── Top-5 stability across scenarios ─────────────────────────────────────────
top_5_variations <- list()

for (i in 2:ncol(scenario_results)) {
  scenario_name   <- colnames(scenario_results)[i]
  top5_scenario   <- which(rank(scenario_results[[i]],             ties.method = "first") <= 5)
  top5_original   <- which(rank(scenario_results$Original_MLSS,    ties.method = "first") <= 5)
  top_5_variations[[scenario_name]] <- length(intersect(top5_scenario, top5_original))
}

# ── Ablation study ────────────────────────────────────────────────────────────
ablation_results <- data.frame(Combination  = top_15_data$Combination,
                               Original_MLSS = top_15_data$MLSS_v4.0,
                               stringsAsFactors = FALSE)

ablation_results$Without_C <- ((0.30/0.55 * top_15_data$Balance_B +
                                  0.10/0.55 * top_15_data$Network_Coverage_V +
                                  0.15/0.55 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_B <- ((0.45/0.70 * top_15_data$Complementarity_C +
                                  0.10/0.70 * top_15_data$Network_Coverage_V +
                                  0.15/0.70 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_V <- ((0.45/0.90 * top_15_data$Complementarity_C +
                                  0.30/0.90 * top_15_data$Balance_B +
                                  0.15/0.90 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_P <- ((0.45/0.85 * top_15_data$Complementarity_C +
                                  0.30/0.85 * top_15_data$Balance_B +
                                  0.10/0.85 * top_15_data$Network_Coverage_V) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Base_Only  <- (0.45 * top_15_data$Complementarity_C +
                                  0.30 * top_15_data$Balance_B +
                                  0.10 * top_15_data$Network_Coverage_V +
                                  0.15 * top_15_data$Potency_P) * 9

ablation_impact <- data.frame(
  Component     = c("Complementarity (C)", "Balance (B)", "Coverage (V)",
                    "Potency (P)", "Synergy/Antagonism"),
  Avg_Rank_Shift = NA,
  Avg_Score_Drop = NA,
  Top5_Stability = NA,
  NPG_Color      = npg_colors[1:5],
  stringsAsFactors = FALSE
)

ablation_cols <- c("Without_C", "Without_B", "Without_V", "Without_P", "Base_Only")

for (i in seq_along(ablation_cols)) {
  ablation_col <- ablation_cols[i]
  rank_diff    <- abs(rank(ablation_results$Original_MLSS) - rank(ablation_results[[ablation_col]]))
  ablation_impact$Avg_Rank_Shift[i] <- mean(rank_diff, na.rm = TRUE)
  ablation_impact$Avg_Score_Drop[i] <- mean(ablation_results$Original_MLSS - ablation_results[[ablation_col]], na.rm = TRUE)
  top5_orig    <- which(rank(ablation_results$Original_MLSS,          ties.method = "first") <= 5)
  top5_abl     <- which(rank(ablation_results[[ablation_col]],         ties.method = "first") <= 5)
  ablation_impact$Top5_Stability[i] <- length(intersect(top5_orig, top5_abl)) / 5
}

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Sensitivity heatmap  →  Fig_S2a
# ═══════════════════════════════════════════════════════════════════════════════
heatmap_data <- scenario_results %>%
  select(-Original_MLSS) %>%
  arrange(match(Combination, top_15$Combination))

rownames(heatmap_data) <- heatmap_data$Combination
heatmap_data           <- heatmap_data[, -1]

heatmap_normalized <- apply(heatmap_data, 2, function(x)
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

heatmap_long <- reshape2::melt(as.matrix(heatmap_normalized))
colnames(heatmap_long) <- c("Combination", "Scenario", "Normalized_Score")

p_heatmap <- ggplot(heatmap_long, aes(x = Scenario, y = Combination, fill = Normalized_Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = npg_colors[6], high = npg_colors[1],
                      name = "Normalized\nMLSS Score") +
  labs(title = "Weight Scenario Sensitivity (Top 15)",
       x = "Weight Scenario", y = "Drug Combination") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.title   = element_text(size = 12, color = "black", face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "right",
    panel.grid   = element_blank()
  )

ggsave("Fig_S2a_mlss_sensitivity_top15.png", p_heatmap,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_S2a_mlss_sensitivity_top15.svg", p_heatmap,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Weight robustness bar chart  →  Fig_6d
# ═══════════════════════════════════════════════════════════════════════════════
top_5_variations_filtered <- top_5_variations[names(top_5_variations) != "Original_MLSS"]
alternative_correlations  <- original_correlations[names(original_correlations) != "MLSS (Standard)"]

robustness_data <- data.frame(
  Scenario       = names(top_5_variations_filtered),
  Top5_Stability = unlist(top_5_variations_filtered) / 5,
  stringsAsFactors = FALSE
)

robustness_data$Rank_Correlation <- c(1.0, alternative_correlations)

robustness_data$NPG_Color <- sapply(robustness_data$Scenario, function(s) {
  idx <- which(weight_scenarios$Scenario == s)
  if (length(idx) > 0) weight_scenarios$NPG_Color[idx] else npg_colors[1]
})

robustness_long <- robustness_data %>%
  pivot_longer(
    cols      = c(Top5_Stability, Rank_Correlation),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = factor(Metric,
                         levels = c("Top5_Stability", "Rank_Correlation"),
                         labels = c("Top 5 Stability", "Rank Correlation")))

p_robustness <- ggplot(robustness_long,
                       aes(x = reorder(Scenario, -Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Top 5 Stability" = npg_colors[1], "Rank Correlation" = npg_colors[2]),
    name   = "Robustness Metric"
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed",  color = npg_colors[6], linewidth = 0.5) +
  geom_hline(yintercept = 0.85, linetype = "dotted",  color = npg_colors[6], linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "MLSS Robustness to Weight Changes (Top 15)",
    x     = "Weight Scenario",
    y     = "Robustness Score (0–1)"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title         = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    plot.title         = element_text(size = 16, face = "bold", hjust = 1),
    legend.position    = "bottom",
    legend.title       = element_text(size = 12, color = "black", face = "bold"),
    legend.text        = element_text(size = 12, color = "black"),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("Fig_6d_Weight_Robustness_Top15.png", plot = p_robustness,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_6d_Weight_Robustness_Top15.svg", plot = p_robustness,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")

