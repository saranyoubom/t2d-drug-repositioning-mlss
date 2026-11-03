###############################################################################
#
# Fig 6d VALIDATION SCRIPT - REVISED FOR MLSS v4.0 WITH NPG PALETTE
# USING TOP 15 COMBINATIONS FOR VALIDATION
#
# Validates Fig_6c_MLSS_v4.0.r output using sensitivity & ablation analysis
# with Nature Portfolio Group (NPG) professional color scheme
#
# INPUT:
#   - Fig_6c_MLSS_ALL_Combinations.csv (from Fig_6c_MLSS_v4.0.r)
#
# OUTPUT (Publication-Ready, 600 DPI):
#   - Fig_6d_MLSS_Sensitivity_Top15.png (top 15 combinations)
#   - Fig_6d_MLSS_Sensitivity_Top15.svg
#   - Fig_6d_Weight_Robustness_Top15.png
#   - Fig_6d_Weight_Robustness_Top15.svg
#   - Fig_S2b_Component_Importance_Top15.png
#   - Fig_S2b_Component_Importance_Top15.svg
#   - Fig_6d_Validation_Report_Top15.txt
#
# MLSS v4.0 FORMULA (validated):
#   MLSS_v4.0 = [C×0.45 + B×0.30 + V×0.10 + P×0.15] × 9
#             + Σ(r_pos × Potency_Mult) × 1
#             - Σ(|r_neg| × Potency_Mult) × 1
#
# Validation Tests (focused on Top 15):
#   1. Weight sensitivity (±20% variation on top 15)
#   2. Component ablation (leave-one-out on top 15)
#   3. Ranking correlation across weight schemes (top 15)
#   4. Robustness metrics (top 15 stability and overlap)
#
# v1.2 - January 2026 - NPG Palette - TOP 15 FOCUS
#
###############################################################################

rm(list = ls())

# ============================================================================
# 0. SETUP & CONFIGURATION
# ============================================================================

setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")
set.seed(42)

required_packages <- c("dplyr", "tidyr", "ggplot2", "ggsci", "reshape2", "gridExtra")

invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

npg_colors <- pal_npg("nrc")(10)

cat("\n")
cat("████████████████████████████████████████████████████████████████████\n")
cat(" Fig 6d VALIDATION - MLSS v4.0 WITH NPG PALETTE\n")
cat(" Sensitivity & Robustness Analysis - TOP 15 COMBINATIONS\n")
cat("████████████████████████████████████████████████████████████████████\n\n")

# ============================================================================
# 1. LOAD AND VALIDATE DATA
# ============================================================================

cat("STEP 1: Loading MLSS v4.0 output from Fig_6c_MLSS_v4.0.r...\n\n")

mlss_data <- read.csv("Fig_6c_MLSS_ALL_Combinations.csv", stringsAsFactors = FALSE)

cat(sprintf(" ✓ Total combinations: %d\n", nrow(mlss_data)))
cat(sprintf(" ✓ Columns: %d\n", ncol(mlss_data)))
cat(sprintf(" ✓ MLSS score range: %.2f - %.2f\n",
            min(mlss_data$MLSS_v4.0, na.rm = TRUE),
            max(mlss_data$MLSS_v4.0, na.rm = TRUE)))

# Verify required columns
required_cols <- c("Complementarity_C", "Balance_B", "Network_Coverage_V", 
                   "Potency_P", "Mean_Pos_r", "Mean_Neg_r", "Synergy_Bonus",
                   "Antagonism_Penalty", "MLSS_v4.0", "Combination")

missing_cols <- setdiff(required_cols, colnames(mlss_data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
}

cat(sprintf(" ✓ All required columns present\n\n"))

# Extract top 15 for detailed analysis
top_15 <- mlss_data %>% head(15)

cat(sprintf(" ✓ Top 15 combos for validation\n"))
cat(sprintf(" ✓ Mean MLSS (top 15): %.3f\n", mean(top_15$MLSS_v4.0, na.rm = TRUE)))
cat(sprintf(" ✓ SD MLSS (top 15): %.3f\n\n", sd(top_15$MLSS_v4.0, na.rm = TRUE)))

# Print top 15 list
cat("TOP 15 DRUG COMBINATIONS:\n")
cat(sprintf("%2d. %s (MLSS=%.2f)\n", 
            1:15, 
            top_15$Combination, 
            top_15$MLSS_v4.0))
cat("\n")

# ============================================================================
# 2. DEFINE WEIGHT SCENARIOS FOR VALIDATION
# ============================================================================

cat("STEP 2: Defining weight scenarios for sensitivity analysis...\n\n")

weight_scenarios <- data.frame(
  Scenario = c(
    "MLSS (Standard)",
    "C-Heavy (60:25:15)",
    "C-Heavy (55:20:25)",
    "B-Heavy (35:50:15)",
    "V-Heavy (35:30:35)",
    "P-Heavy (30:20:10:40)",
    "Balanced (40:40:10:10)",
    "Network-Focus (30:20:40:10)"
  ),
  C_Weight = c(0.45, 0.60, 0.55, 0.35, 0.35, 0.30, 0.40, 0.30),
  B_Weight = c(0.30, 0.25, 0.20, 0.50, 0.30, 0.20, 0.40, 0.20),
  V_Weight = c(0.10, 0.10, 0.15, 0.10, 0.30, 0.10, 0.10, 0.40),
  P_Weight = c(0.15, 0.05, 0.10, 0.05, 0.05, 0.40, 0.10, 0.10),
  NPG_Color = c(npg_colors[1], npg_colors[2], npg_colors[3], npg_colors[4],
                npg_colors[5], npg_colors[6], npg_colors[7], npg_colors[8]),
  stringsAsFactors = FALSE
)

print(weight_scenarios[, 1:6])

cat("\n ✓ Weight scenarios defined\n\n")

# ============================================================================
# 3. RECALCULATE SCORES FOR EACH SCENARIO (TOP 15 ONLY)
# ============================================================================

cat("STEP 3: Recalculating MLSS scores for top 15 across weight scenarios...\n\n")

# Work with top 15 data
top_15_data <- mlss_data %>% head(15)

scenario_results <- data.frame(
  Combination = top_15_data$Combination,
  Original_MLSS = top_15_data$MLSS_v4.0
)

for (i in 1:nrow(weight_scenarios)) {
  scenario_name <- weight_scenarios$Scenario[i]
  c_w <- weight_scenarios$C_Weight[i]
  b_w <- weight_scenarios$B_Weight[i]
  v_w <- weight_scenarios$V_Weight[i]
  p_w <- weight_scenarios$P_Weight[i]
  
  # Recalculate base score for top 15
  new_base_score <- (c_w * top_15_data$Complementarity_C +
                     b_w * top_15_data$Balance_B +
                     v_w * top_15_data$Network_Coverage_V +
                     p_w * top_15_data$Potency_P) * 9
  
  # Final MLSS v4.0 (keep synergy/antagonism scaling = 1)
  new_mlss <- new_base_score + 
              top_15_data$Synergy_Bonus - 
              top_15_data$Antagonism_Penalty
  
  scenario_results[[scenario_name]] <- new_mlss
  
  cat(sprintf(" ✓ %s: mean MLSS = %.3f (SD=%.3f)\n",
              scenario_name, mean(new_mlss, na.rm = TRUE), sd(new_mlss, na.rm = TRUE)))
}

cat("\n")

# ============================================================================
# 4. RANKING CORRELATION ANALYSIS (TOP 15)
# ============================================================================

cat("STEP 4: Computing ranking correlations for top 15 across scenarios...\n\n")

# Convert to ranks
rank_data <- apply(scenario_results[, -1], 2, rank)

# Correlation matrix
rank_correlation <- cor(rank_data, method = "spearman")

# Extract correlations with original
original_correlations <- rank_correlation[1, -1]

cat(" Spearman Rank Correlation with Original MLSS (Top 15):\n\n")
for (i in 2:ncol(rank_correlation)) {
  scenario_name <- colnames(rank_correlation)[i]
  corr <- rank_correlation[1, i]
  cat(sprintf(" %-30s: %.4f\n", scenario_name, corr))
}

cat("\n")

# Compute robustness: how much do top 5 rankings change? (subset of top 15)
top_5_variations <- list()

for (i in 2:ncol(scenario_results)) {
  scenario_name <- colnames(scenario_results)[i]
  top_5_scenario <- which(rank(scenario_results[[i]], ties.method = "first") <= 5)
  top_5_original <- which(rank(scenario_results$Original_MLSS, ties.method = "first") <= 5)
  overlap <- length(intersect(top_5_scenario, top_5_original))
  
  top_5_variations[[scenario_name]] <- overlap
  
  cat(sprintf(" %s: %d/5 overlap with original (top 5)\n", scenario_name, overlap))
}

cat("\n✓ Ranking correlation analysis complete\n\n")

# ============================================================================
# 5. COMPONENT CONTRIBUTION ANALYSIS - TOP 15 ONLY
# ============================================================================

cat("STEP 5: Ablation study on top 15 - Component importance...\n\n")

# Leave-one-out: remove each component for top 15
ablation_results <- data.frame(
  Combination = top_15_data$Combination,
  Original_MLSS = top_15_data$MLSS_v4.0,
  stringsAsFactors = FALSE
)

# WITHOUT C (Complementarity)
ablation_results$Without_C <- ((0 * top_15_data$Complementarity_C +
                               0.30/0.55 * top_15_data$Balance_B +
                               0.10/0.55 * top_15_data$Network_Coverage_V +
                               0.15/0.55 * top_15_data$Potency_P) * 9) +
                              top_15_data$Synergy_Bonus - 
                              top_15_data$Antagonism_Penalty

# WITHOUT B (Balance)
ablation_results$Without_B <- ((0.45/0.70 * top_15_data$Complementarity_C +
                               0 * top_15_data$Balance_B +
                               0.10/0.70 * top_15_data$Network_Coverage_V +
                               0.15/0.70 * top_15_data$Potency_P) * 9) +
                              top_15_data$Synergy_Bonus - 
                              top_15_data$Antagonism_Penalty

# WITHOUT V (Coverage)
ablation_results$Without_V <- ((0.45/0.90 * top_15_data$Complementarity_C +
                               0.30/0.90 * top_15_data$Balance_B +
                               0 * top_15_data$Network_Coverage_V +
                               0.15/0.90 * top_15_data$Potency_P) * 9) +
                              top_15_data$Synergy_Bonus - 
                              top_15_data$Antagonism_Penalty

# WITHOUT P (Potency)
ablation_results$Without_P <- ((0.45/0.85 * top_15_data$Complementarity_C +
                               0.30/0.85 * top_15_data$Balance_B +
                               0.10/0.85 * top_15_data$Network_Coverage_V +
                               0 * top_15_data$Potency_P) * 9) +
                              top_15_data$Synergy_Bonus - 
                              top_15_data$Antagonism_Penalty

# WITHOUT Synergy/Antagonism
ablation_results$Base_Only <- ((0.45 * top_15_data$Complementarity_C +
                               0.30 * top_15_data$Balance_B +
                               0.10 * top_15_data$Network_Coverage_V +
                               0.15 * top_15_data$Potency_P) * 9)

# Calculate impact for top 15
ablation_impact <- data.frame(
  Component = c("Complementarity (C)", "Balance (B)", "Coverage (V)", "Potency (P)", "Synergy/Antagonism"),
  Avg_Rank_Shift = NA,
  Avg_Score_Drop = NA,
  Top5_Stability = NA,
  NPG_Color = npg_colors[1:5],
  stringsAsFactors = FALSE
)

ablation_cols <- c("Without_C", "Without_B", "Without_V", "Without_P", "Base_Only")

for (i in 1:length(ablation_cols)) {
  ablation_col <- ablation_cols[i]
  
  # Rank difference
  rank_diff <- abs(rank(ablation_results$Original_MLSS) - rank(ablation_results[[ablation_col]]))
  ablation_impact$Avg_Rank_Shift[i] <- mean(rank_diff, na.rm = TRUE)
  
  # Score drop
  score_drop <- ablation_results$Original_MLSS - ablation_results[[ablation_col]]
  ablation_impact$Avg_Score_Drop[i] <- mean(score_drop, na.rm = TRUE)
  
  # Top 5 stability (within top 15)
  top5_original <- which(rank(ablation_results$Original_MLSS, ties.method = "first") <= 5)
  top5_ablated <- which(rank(ablation_results[[ablation_col]], ties.method = "first") <= 5)
  overlap <- length(intersect(top5_original, top5_ablated))
  ablation_impact$Top5_Stability[i] <- overlap / 5
}

print(ablation_impact)

cat("\n✓ Ablation study complete (top 15)\n\n")

# ============================================================================
# 6. VISUALIZATION 1: SENSITIVITY HEATMAP (TOP 15)
# ============================================================================

cat("STEP 6: Creating sensitivity analysis visualization (top 15)...\n\n")

# Use all top 15 for heatmap
heatmap_data <- scenario_results %>%
  select(-Original_MLSS) %>%
  arrange(match(Combination, top_15$Combination))

rownames(heatmap_data) <- heatmap_data$Combination
heatmap_data <- heatmap_data[, -1]

# Normalize to 0-1 for visualization
heatmap_normalized <- apply(heatmap_data, 2, function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

heatmap_long <- reshape2::melt(as.matrix(heatmap_normalized))
colnames(heatmap_long) <- c("Combination", "Scenario", "Normalized_Score")

p_heatmap <- ggplot(heatmap_long, aes(x = Scenario, y = Combination, fill = Normalized_Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = npg_colors[6],
    high = npg_colors[1],
    name = "Normalized\nMLSS Score"
  ) +
  labs(
    title = "Weight Scenario Sensitivity (Top 15)",
    #subtitle = "MLSS scores across different weight parameterizations (normalized)",
    x = "Weight Scenario",
    y = "Drug Combination"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0, color = "grey40"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )
p_heatmap
ggsave(
  "Fig_S2a_MLSS_Sensitivity_Top15.png",
  plot = p_heatmap,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Fig_S2a_MLSS_Sensitivity_Top15.svg",
  plot = p_heatmap,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

cat(" ✓ Saved: Fig_S2a_MLSS_Sensitivity_Top15.png (600 DPI)\n")
cat(" ✓ Saved: Fig_S2a_MLSS_Sensitivity_Top15.svg (600 DPI)\n\n")

# ============================================================================
# 7. VISUALIZATION 2: WEIGHT ROBUSTNESS BAR CHART (TOP 15)
# ============================================================================

cat("STEP 7: Creating weight robustness visualization (top 15)...\n\n")

# Filter and prepare robustness data
top_5_variations_filtered <- top_5_variations[names(top_5_variations) != "Original_MLSS"]

alternative_correlations <- original_correlations[names(original_correlations) != "MLSS (Standard)"]

robustness_data <- data.frame(
  Scenario = names(top_5_variations_filtered),
  Top5_Stability = unlist(top_5_variations_filtered) / 5,
  stringsAsFactors = FALSE
)

robustness_data$Rank_Correlation <- c(1.0, alternative_correlations)

robustness_data$NPG_Color <- sapply(robustness_data$Scenario, function(s) {
  idx <- which(weight_scenarios$Scenario == s)
  if (length(idx) > 0) {
    return(weight_scenarios$NPG_Color[idx])
  } else {
    return(npg_colors[1])
  }
})

robustness_long <- robustness_data %>%
  pivot_longer(
    cols = c(Top5_Stability, Rank_Correlation),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("Top5_Stability", "Rank_Correlation"),
                    labels = c("Top 5 Stability", "Rank Correlation"))
  )

p_robustness <- ggplot(robustness_long, aes(x = reorder(Scenario, -Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Top 5 Stability" = npg_colors[1], "Rank Correlation" = npg_colors[2]),
    name = "Robustness Metric"
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = npg_colors[6], linewidth = 0.5) +
  geom_hline(yintercept = 0.85, linetype = "dotted", color = npg_colors[6], linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "MLSS Robustness to Weight Changes (Top 15)",
    #subtitle = "Stability of top 5 rankings within top 15 and correlation with original MLSS v4.0",
    x = "Weight Scenario",
    y = "Robustness Score (0-1)",
    #caption = "Dashed line: 0.95 (excellent); Dotted line: 0.85 (good)"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, color = "black",  face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 1),
    plot.subtitle = element_text(size = 12, hjust = 0, color = "grey40"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, color = "black",  face = "bold"),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  )
p_robustness
ggsave(
  "Fig_6d_Weight_Robustness_Top15.png",
  plot = p_robustness,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Fig_6d_Weight_Robustness_Top15.svg",
  plot = p_robustness,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

cat(" ✓ Saved: Fig_6d_Weight_Robustness_Top15.png (600 DPI)\n")
cat(" ✓ Saved: Fig_6d_Weight_Robustness_Top15.svg (600 DPI)\n\n")

# ============================================================================
# 8. VISUALIZATION 3: COMPONENT IMPORTANCE (TOP 15 ABLATION)
# ============================================================================

cat("STEP 8: Creating component importance visualization (top 15)...\n\n")

importance_plot_data <- ablation_impact %>%
  pivot_longer(
    cols = c(Avg_Rank_Shift, Avg_Score_Drop),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("Avg_Rank_Shift", "Avg_Score_Drop"),
                    labels = c("Avg Rank Shift", "Avg Score Drop")),
    Component = factor(Component, levels = ablation_impact$Component)
  )

p_importance <- ggplot(importance_plot_data, aes(x = reorder(Component, -abs(Value)), y = abs(Value), fill = NPG_Color)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3, width = 0.7) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_identity() +
  labs(
    title = "Component Importance in MLSS (Top 15)",
    #subtitle = "Ablation study: impact of removing each component from model",
    x = "Component",
    y = "Impact Metric (Left: Rank positions; Right: Score units)",
    #caption = "Higher values = more important components"
  ) +
  coord_flip() +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, hjust = 0, color = "grey40"),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )
p_importance
ggsave(
  "Fig_S2b_Component_Importance_Top15.png",
  plot = p_importance,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Fig_S2b_Component_Importance_Top15.svg",
  plot = p_importance,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

cat(" ✓ Saved: Fig_S2b_Component_Importance_Top15.png (600 DPI)\n")
cat(" ✓ Saved: Fig_S2b_Component_Importance_Top15.svg (600 DPI)\n\n")

# End of script
