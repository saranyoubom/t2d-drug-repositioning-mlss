# =============================================================================
# FigS1 — Observed vs. Permutation Null + Pearson r CI Validation
# =============================================================================
# OUTPUTS:
#   FigS1a_Validation_Summary.csv     → Supplementary source data for Fig. S1
#   FigS1a_Validation_BarChart.png   → Fig. S1a manuscript (600 dpi)
#   FigS1a_Validation_BarChart.svg   → Fig. S1a journal submission (vector)
#   FigS1b_Pearson_r_CI_Table.csv           → Supplementary source data (CI sensitivity)
#   FigS1b_Pearson_r_CI_Plot.png     → Fig. S1b manuscript (600 dpi)
#   FigS1b_Pearson_r_CI_Plot.svg     → Fig. S1b journal submission (vector)
#   Console: CI values               → Discussion + Reviewer Response
#   Console: r_min                   → Methods justification for |r| >= r_thresh
# =============================================================================

rm(list = ls())
set.seed(42)

library(ggplot2)
library(dplyr)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# =============================================================================
# HELPERS
# =============================================================================

col_shuffle <- function(df) as.data.frame(apply(df, 2, sample))

# count_sig_pairs — add p < 0.05 filter
count_sig_pairs <- function(mat, threshold = 0.7) {
  mat     <- as.matrix(na.omit(mat))
  n       <- nrow(mat)
  cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  upper   <- cor_mat[upper.tri(cor_mat)]
  t_stat  <- upper * sqrt(n - 2) / sqrt(1 - upper^2)
  p_val   <- 2 * pt(-abs(t_stat), df = n - 2)
  sum(abs(upper) >= threshold & p_val < 0.05, na.rm = TRUE)
}

# get_sig_pairs — add p < 0.05 filter
get_sig_pairs <- function(mat, threshold = 0.7) {
  mat     <- as.matrix(na.omit(mat))
  n       <- nrow(mat)
  cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  cn      <- colnames(cor_mat)
  pairs   <- character(0)
  for (i in seq_len(ncol(cor_mat) - 1))
    for (j in (i + 1):ncol(cor_mat)) {
      r     <- cor_mat[i, j]
      if (is.na(r)) next
      t_stat <- r * sqrt(n - 2) / sqrt(1 - r^2)
      p_val  <- 2 * pt(-abs(t_stat), df = n - 2)
      if (abs(r) >= threshold && p_val < 0.05)
        pairs <- c(pairs, paste(sort(c(cn[i], cn[j])), collapse = "-"))
    }
  unique(pairs)
}
make_normalized <- function(expr, ref_gene) {
  expr[, colnames(expr) != ref_gene, drop = FALSE] / expr[[ref_gene]]
}

pearson_r_ci <- function(r, n, alpha = 0.05) {
  z       <- atanh(r)
  se      <- 1 / sqrt(n - 3)
  z_crit  <- qnorm(1 - alpha / 2)
  ci_low  <- tanh(z - z_crit * se)
  ci_high <- tanh(z + z_crit * se)
  return(c(lower = round(ci_low, 3), upper = round(ci_high, 3)))
}

# Shared theme for Fig. S1a and Fig. S1b
theme_figs1 <- theme_minimal(base_size = 10) +
  theme(
    legend.position    = "none",
    axis.title         = element_text(size = 10, color = "black", family = "Arial"),
    axis.text          = element_text(size = 9,  color = "black", family = "Arial"),
    plot.title         = element_text(size = 11, color = "black", family = "Arial"),
    plot.caption       = element_text(size = 7,  color = "grey40", family = "Arial",
                                      lineheight = 1.3, hjust = 0),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.line.x        = element_line(color = "grey60", linewidth = 0.3)
  )

# Shared palette for Fig. S1a and Fig. S1b
fill_colors <- c(
  sig       = npg_colors[4],   # teal  — significant single methods
  ns        = "grey72",        # grey  — non-significant
  consensus = npg_colors[1]    # red   — consensus highlight
)

# =============================================================================
# DATA
# =============================================================================

expr_raw <- read.csv("Fig_4c_Raw_expression_data.csv", stringsAsFactors = FALSE)
expr_raw <- expr_raw[, colnames(expr_raw) != "Ryr1"]

n_genes  <- ncol(expr_raw)
n_pairs  <- n_genes * (n_genes - 1) / 2   # 190 pairs
n_obs    <- nrow(expr_raw)                 # derived — not hardcoded
r_thresh <- 0.7
n_perm   <- 1000

# =============================================================================
# OBSERVED COUNTS
# → Values cited in Fig. S1a legend and Discussion
# =============================================================================

obs <- c(
  Raw       = count_sig_pairs(expr_raw,                          r_thresh),
  Ins1      = count_sig_pairs(make_normalized(expr_raw, "Ins1"),   r_thresh),
  Kcnj11    = count_sig_pairs(make_normalized(expr_raw, "Kcnj11"), r_thresh),
  Hes1      = count_sig_pairs(make_normalized(expr_raw, "Hes1"),   r_thresh),
  Consensus = length(Reduce(intersect, list(
    get_sig_pairs(expr_raw,                          r_thresh),
    get_sig_pairs(make_normalized(expr_raw, "Ins1"),   r_thresh),
    get_sig_pairs(make_normalized(expr_raw, "Kcnj11"), r_thresh),
    get_sig_pairs(make_normalized(expr_raw, "Hes1"),   r_thresh)
  )))
)

cat("Observed pairs:", obs, "\n")

# =============================================================================
# PERMUTATION NULL
# → Null mean, SD, 95th pct cited in:
#     - Fig. S1a bar chart
#     - Fig. S1a legend
#     - Discussion (SNR, empirical p)
#     - Reviewer Response (R1 sample size comment)
# =============================================================================

null_mat <- matrix(NA_real_, nrow = n_perm, ncol = 5,
                   dimnames = list(NULL, c("Raw","Ins1","Kcnj11","Hes1","Consensus")))

for (i in seq_len(n_perm)) {
  p_raw    <- col_shuffle(expr_raw)
  p_ins1   <- col_shuffle(expr_raw)
  p_kcnj11 <- col_shuffle(expr_raw)
  p_hes1   <- col_shuffle(expr_raw)
  
  null_mat[i, "Raw"]    <- count_sig_pairs(p_raw, r_thresh)
  null_mat[i, "Ins1"]   <- count_sig_pairs(make_normalized(p_ins1,   "Ins1"),   r_thresh)
  null_mat[i, "Kcnj11"] <- count_sig_pairs(make_normalized(p_kcnj11, "Kcnj11"), r_thresh)
  null_mat[i, "Hes1"]   <- count_sig_pairs(make_normalized(p_hes1,   "Hes1"),   r_thresh)
  null_mat[i, "Consensus"] <- length(Reduce(intersect, list(
    get_sig_pairs(p_raw, r_thresh),
    get_sig_pairs(make_normalized(p_ins1,   "Ins1"),   r_thresh),
    get_sig_pairs(make_normalized(p_kcnj11, "Kcnj11"), r_thresh),
    get_sig_pairs(make_normalized(p_hes1,   "Hes1"),   r_thresh)
  )))
  
  if (i %% 100 == 0) cat("Permutation", i, "/", n_perm, "\n")
}

# =============================================================================
# SUMMARY TABLE
# → OUTPUT: FigS1a_Validation_Summary.csv
#     Supplementary source data for Fig. S1a; available to reviewers on request
# =============================================================================

method_labels <- c(
  Raw       = "Raw\nExpression",
  Ins1      = "Ins1-\nNormalized",
  Kcnj11    = "Kcnj11-\nNormalized",
  Hes1      = "Hes1-\nNormalized",
  Consensus = "4-Method\nConsensus"
)

summary_df <- data.frame(
  Method_key = names(obs),
  Method     = factor(method_labels[names(obs)], levels = method_labels),
  Observed   = as.numeric(obs),
  Null_Mean  = colMeans(null_mat),
  Null_SD    = apply(null_mat, 2, sd),
  Null_95    = apply(null_mat, 2, quantile, 0.95),
  Emp_p      = mapply(function(k, o) mean(null_mat[, k] >= o),
                      names(obs), as.numeric(obs))
) %>%
  mutate(
    FP_Rate_Pct = round(100 * Null_Mean / n_pairs, 2),
    SNR         = round(Observed / pmax(Null_Mean, 0.01), 1),
    sig_label   = case_when(
      Emp_p == 0   ~ "***",
      Emp_p < 0.01 ~ "**",
      Emp_p < 0.05 ~ "*",
      TRUE         ~ "n.s."
    ),
    Fill_group = case_when(
      Method_key == "Consensus" ~ "consensus",
      Emp_p >= 0.05             ~ "ns",
      TRUE                      ~ "sig"
    )
  )

write.csv(summary_df %>% select(-Fill_group),
          "FigS1a_Validation_Summary.csv", row.names = FALSE)
# OUTPUT: FigS1a_Validation_Summary.csv — source data for Fig. S1a

cat("\n=== Validation Summary ===\n")
print(summary_df[, c("Method_key","Observed","Null_Mean","Null_SD",
                     "FP_Rate_Pct","SNR","Emp_p","sig_label")])

# =============================================================================
# FIG. S1a — Observed vs. Permutation Null Bar Chart
# → OUTPUT: FigS1a_Validation_BarChart.png   manuscript (600 dpi)
# → OUTPUT: FigS1a_Validation_BarChart.svg   journal submission (vector)
# =============================================================================

summary_df <- summary_df %>%
  mutate(label_y = Observed + max(Observed) * 0.04)

p_s1a <- ggplot(summary_df, aes(x = Method, y = Observed, fill = Fill_group)) +
  
  geom_col(width = 0.62, alpha = 0.7, color = NA) +
  
  geom_errorbar(
    aes(ymin = Null_Mean - Null_SD,
        ymax = Null_Mean + Null_SD),
    width = 0.22, linewidth = 0.7, color = "grey25"
  ) +
  
  geom_point(
    aes(y = Null_Mean),
    shape = 18, size = 2.8, color = "grey25"
  ) +
  
  # ── Bar count labels (inside, near top of each bar) ──────────────────
  geom_text(
    aes(y = Observed * 0.92, label = Observed),
    size = 3.0, family = "Arial", fontface = "bold",
    color = "white", vjust = 1
  ) +
  
  # ── Significance labels (above bar top) ──────────────────────────────
  geom_text(
    aes(y = label_y, label = sig_label),
    size = 3.2, family = "Arial", color = "grey20"
  ) +
  
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  
  labs(x = NULL, y = "Significant Gene Pairs") +
  theme_figs1
ggsave("FigS1a_Validation_BarChart.png", plot = p_s1a,
       scale = 1, width = 12, height = 8, dpi = 600, units = "cm")
# OUTPUT: FigS1a_Validation_BarChart.png — Fig. S1a manuscript

ggsave("FigS1a_Validation_BarChart.svg", plot = p_s1a,
       scale = 1, width = 12, height = 8, dpi = 600, units = "cm")
# OUTPUT: FigS1a_Validation_BarChart.svg — Fig. S1a journal submission

cat("\nFig. S1a saved (PNG + SVG).\n")

# =============================================================================
# PEARSON r CI — Fisher z-transformation
# → VALUES CITED IN:
#     - Discussion: "r = r_thresh, n = n_obs: 95% CI approximately 0.21 to 0.91"
#     - Reviewer Response (R1 sample size comment)
# =============================================================================

cat(sprintf("\nr threshold : %.2f (from r_thresh)\n", r_thresh))
cat(sprintf("n           : %d   (nrow of expr_raw)\n", n_obs))

ci <- pearson_r_ci(r = r_thresh, n = n_obs)
cat(sprintf(
  "95%% CI for r = %.2f, n = %d (Fisher z-transformation): %.3f to %.3f\n",
  r_thresh, n_obs, ci["lower"], ci["upper"]
))
# Cited in manuscript as: "approximately 0.21 to 0.91"

# Minimum r for significance at n = n_obs, p < 0.05 two-tailed
t_crit <- qt(0.975, df = n_obs - 2)
r_min  <- t_crit / sqrt(t_crit^2 + n_obs - 2)
cat(sprintf(
  "Minimum |r| for p < 0.05 at n = %d (df = %d): %.3f\n",
  n_obs, n_obs - 2, r_min
))
cat(sprintf(
  "r_thresh (%.2f) > r_min (%.3f): threshold is more stringent than significance floor\n",
  r_thresh, r_min
))

# =============================================================================
# CI SENSITIVITY TABLE
# → OUTPUT: FigS1b_Pearson_r_CI_Table.csv — source data for Fig. S1b
# =============================================================================

r_values <- seq(0.60, 0.90, by = 0.10)
ci_table <- do.call(rbind, lapply(r_values, function(r) {
  ci_r <- pearson_r_ci(r, n_obs)
  data.frame(
    r        = r,
    n        = n_obs,
    CI_low   = ci_r["lower"],
    CI_high  = ci_r["upper"],
    CI_width = round(ci_r["upper"] - ci_r["lower"], 3)
  )
}))
rownames(ci_table) <- NULL

cat("\n=== 95% CI Sensitivity Table (n =", n_obs, ") ===\n")
print(ci_table)

write.csv(ci_table, "FigS1b_Pearson_r_CI_Table.csv", row.names = FALSE)
# OUTPUT: FigS1b_Pearson_r_CI_Table.csv — source data for Fig. S1b

# =============================================================================
# FIG. S1b — 95% CI Width for Pearson r at n = n_obs
# → OUTPUT: FigS1b_Pearson_r_CI_Plot.png   manuscript (600 dpi)
# → OUTPUT: FigS1b_Pearson_r_CI_Plot.svg   journal submission (vector)
# =============================================================================

ci_plot_df       <- ci_table
ci_plot_df$color <- ifelse(ci_plot_df$r == r_thresh,
                           npg_colors[1], npg_colors[4])

p_s1b <- ggplot(ci_plot_df, aes(x = r, y = r)) +
  
  geom_linerange(
    aes(ymin = CI_low, ymax = CI_high, color = color),
    linewidth = 2, show.legend = FALSE, alpha = 0.7
  ) +
  
  geom_point(
    aes(color = color), size = 3, show.legend = FALSE, shape = 21
  ) +
  
  # CI bound labels for r_thresh only
  geom_text(
    data = subset(ci_plot_df, r == r_thresh),
    aes(y = CI_low,  label = sprintf("%.3f", CI_low)),
    hjust = 1.3, size = 2.5, color = npg_colors[1]
  ) +
  geom_text(
    data = subset(ci_plot_df, r == r_thresh),
    aes(y = CI_high, label = sprintf("%.3f", CI_high)),
    hjust = 1.3, size = 2.5, color = npg_colors[1]
  ) +
  
  # CI width labels right of each bar
  geom_text(
    aes(y = r, label = sprintf("w = %.3f", CI_width), color = color),
    hjust = -0.25, size = 2.5, show.legend = FALSE
  ) +
  
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey70", linewidth = 0.5) +
  
  geom_vline(xintercept = r_thresh, linetype = "dashed",
             color = npg_colors[1], linewidth = 0.5, alpha = 0.5) +
  
  geom_vline(xintercept = r_min, linetype = "dashed",
             color = "grey30", linewidth = 0.5, alpha = 0.5) +
  
  scale_color_identity() +
  scale_x_continuous(limits = c(0.5, 1.0),
                     breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  scale_y_continuous(limits = c(-0.2, 1.2),
                     breaks = seq(0, 1, by = 0.2)) +
  
  annotate("text", x = r_thresh, y = 1.04,
           label = paste("|r| threshold =", r_thresh),
           color = npg_colors[1], size = 2.5, hjust = 0.5) +
  annotate("text", x = r_min, y = -0.07,
           label = paste0("r_min = ", round(r_min, 3), "\n(p < 0.05 floor)"),
           color = "grey40", size = 2.5, hjust = 0.5, lineheight = 1.1) +
  
  labs(
    x     = "Observed r",
    y     = "r and 95% CI bounds",
    title = paste0("95% CI Widths for Pearson r at n = ", n_obs)
  ) +
  theme_figs1 +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line.x        = element_blank()
  )

ggsave("FigS1b_Pearson_r_CI_Plot.png", plot = p_s1b,
       scale = 1, width = 12, height = 8, dpi = 600, units = "cm")
# OUTPUT: FigS1b_Pearson_r_CI_Plot.png — Fig. S1b manuscript

ggsave("FigS1b_Pearson_r_CI_Plot.svg", plot = p_s1b,
       scale = 1, width = 12, height = 8, dpi = 600, units = "cm")
# OUTPUT: FigS1b_Pearson_r_CI_Plot.svg — Fig. S1b journal submission

cat("\nFig. S1b saved (PNG + SVG).\n")
cat("\nAll outputs saved.\n")

# =============================================================================
# OUTPUT SUMMARY
# -----------------------------------------------------------------------------
# FigS1a_Validation_Summary.csv      Supplementary source data for Fig. S1a
# FigS1a_Validation_BarChart.png    Fig. S1a — manuscript (600 dpi)
# FigS1a_Validation_BarChart.svg    Fig. S1a — journal submission (vector)
# FigS1b_Pearson_r_CI_Table.csv            Supplementary source data for Fig. S1b
# FigS1b_Pearson_r_CI_Plot.png      Fig. S1b — manuscript (600 dpi)
# FigS1b_Pearson_r_CI_Plot.svg      Fig. S1b — journal submission (vector)
# Console: CI (0.211–0.909)         Discussion + Reviewer Response
# Console: r_min (0.576)            Methods justification for |r| >= r_thresh
# =============================================================================
