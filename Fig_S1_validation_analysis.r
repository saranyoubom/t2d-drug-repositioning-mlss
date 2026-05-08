# =============================================================================
# Fig. S1 — Multi-Normalization Consensus Validation (Simplified, 3 panels)
#
# DATA DEPENDENCIES:
#   Fig_4c_Raw_expression_data.csv
#   Fig_4c_Ins1_normalized_data.csv
#   Fig_4c_Kcnj11_normalized_data.csv
#   Fig_4c_Hes1_normalized_data.csv
#   Fig_4d_consensus_pairs.csv
#
# Run order: Fig_4c → Fig_4d → this script
# =============================================================================

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_Fig_S1 <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 9),
    axis.line          = element_line(linewidth = 0.4, colour = "grey30"),
    axis.ticks         = element_line(linewidth = 0.3, colour = "grey30"),
    axis.text          = element_text(colour = "grey20", family = "Arial"),
    axis.title         = element_text(colour = "grey20", family = "Arial"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.key.size    = unit(0.35, "cm"),
    legend.text        = element_text(size = 8, family = "Arial"),
    legend.title       = element_text(size = 8, face = "bold", family = "Arial"),
    plot.title         = element_text(face = "bold", size = 9,
                                      colour = "grey15", family = "Arial"),
    plot.subtitle      = element_text(size = 7.5, colour = "grey40",
                                      family = "Arial")
  )


# =============================================================================
# 0. LOAD DATA
# =============================================================================

raw_data  <- read.csv("Fig_4c_Raw_expression_data.csv",    check.names = FALSE)
ins1_norm <- read.csv("Fig_4c_Ins1_normalized_data.csv",   check.names = FALSE)
kcnj_norm <- read.csv("Fig_4c_Kcnj11_normalized_data.csv", check.names = FALSE)
hes1_norm <- read.csv("Fig_4c_Hes1_normalized_data.csv",   check.names = FALSE)

drop_ryr1 <- function(df) df[, colnames(df) != "Ryr1", drop = FALSE]
raw_data  <- drop_ryr1(raw_data)
ins1_norm <- drop_ryr1(ins1_norm)
kcnj_norm <- drop_ryr1(kcnj_norm)
hes1_norm <- drop_ryr1(hes1_norm)

consensus_df  <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
consensus_ids <- consensus_df[[1]]

n_obs <- nrow(raw_data)   # 12


# =============================================================================
# 1. HELPERS
# =============================================================================

compute_pairs <- function(df, method_label) {
  df_mat  <- as.matrix(na.omit(df))
  genes   <- colnames(df_mat)
  n_genes <- length(genes)
  out     <- vector("list", choose(n_genes, 2))
  k <- 0L
  for (i in seq_len(n_genes - 1)) {
    for (j in seq(i + 1, n_genes)) {
      ct  <- cor.test(df_mat[, i], df_mat[, j], method = "pearson")
      k   <- k + 1L
      out[[k]] <- data.frame(
        Gene1   = genes[i], Gene2 = genes[j],
        pair_id = paste(genes[i], genes[j], sep = "_"),
        r       = as.numeric(ct$estimate),
        p_val   = ct$p.value,
        CI_low  = ct$conf.int[1],
        CI_high = ct$conf.int[2],
        Method  = method_label,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

fisher_pool <- function(r_vec, n_per_method = 12, alpha = 0.05) {
  r_vec  <- pmin(pmax(r_vec, -0.9999), 0.9999)
  zs     <- atanh(r_vec)
  ws     <- rep(n_per_method - 3, length(r_vec))
  z_bar  <- weighted.mean(zs, ws)
  se     <- 1 / sqrt(sum(ws))
  z_crit <- qnorm(1 - alpha / 2)
  list(r_pool = tanh(z_bar),
       ci_lo  = tanh(z_bar - z_crit * se),
       ci_hi  = tanh(z_bar + z_crit * se))
}


# =============================================================================
# 2. COMPUTE CORRELATIONS
# =============================================================================

method_levels <- c("Raw\nexpression", "Ins1-\nnormalized",
                   "Kcnj11-\nnormalized", "Hes1-\nnormalized")

all_pairs <- bind_rows(
  compute_pairs(raw_data,  "Raw\nexpression"),
  compute_pairs(ins1_norm, "Ins1-\nnormalized"),
  compute_pairs(kcnj_norm, "Kcnj11-\nnormalized"),
  compute_pairs(hes1_norm, "Hes1-\nnormalized")
) %>%
  mutate(
    Method    = factor(Method, levels = method_levels),
    sig       = p_val < 0.05,
    consensus = pair_id %in% consensus_ids
  )


# =============================================================================
# 3. PERMUTATION NULL (1,000 permutations)
# =============================================================================

set.seed(2025)
run_permutation <- function(df_expr, n_perm = 1000, alpha = 0.05) {
  df_mat  <- as.matrix(na.omit(df_expr))
  n       <- nrow(df_mat)
  n_genes <- ncol(df_mat)
  null_counts <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    perm_mat <- apply(df_mat, 2, function(x) x[sample(n)])
    cnt <- 0L
    for (i in seq_len(n_genes - 1))
      for (j in seq(i + 1, n_genes)) {
        ct <- cor.test(perm_mat[, i], perm_mat[, j], method = "pearson")
        if (!is.na(ct$p.value) && ct$p.value < alpha) cnt <- cnt + 1L
      }
    null_counts[p] <- cnt
  }
  list(mean = mean(null_counts), sd = sd(null_counts),
       q95  = quantile(null_counts, 0.95))
}

cat("Running permutations...\n")
nulls <- list(
  "Raw\nexpression"     = run_permutation(raw_data),
  "Ins1-\nnormalized"   = run_permutation(ins1_norm),
  "Kcnj11-\nnormalized" = run_permutation(kcnj_norm),
  "Hes1-\nnormalized"   = run_permutation(hes1_norm)
)


# =============================================================================
# 4. PANEL a — Observed vs Permutation Null
# =============================================================================

obs_df <- all_pairs %>%
  filter(sig) %>%
  count(Method, name = "Observed") %>%
  mutate(
    Null_Mean  = sapply(as.character(Method), function(m) nulls[[m]]$mean),
    Null_SD    = sapply(as.character(Method), function(m) nulls[[m]]$sd),
    Emp_p      = pnorm(Observed, mean = Null_Mean, sd = Null_SD,
                       lower.tail = FALSE),
    sig_label  = case_when(
      Emp_p < 0.001 ~ "***",
      Emp_p < 0.01  ~ "**",
      Emp_p < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    Fill_group = "Single",
    label_y    = pmax(Observed, Null_Mean + Null_SD) + max(Observed) * 0.04
  )

# Add consensus row
consensus_row <- tibble(
  Method     = factor("4-Method\nConsensus",
                      levels = c(method_levels, "4-Method\nConsensus")),
  Observed   = length(consensus_ids),
  Null_Mean  = 0.05^4 * choose(ncol(raw_data), 2),
  Null_SD    = NA_real_,
  Emp_p      = NA_real_,
  sig_label  = "",
  Fill_group = "Consensus",
  label_y    = length(consensus_ids) + max(obs_df$Observed) * 0.04
)

summary_df <- bind_rows(obs_df, consensus_row) %>%
  mutate(Method = factor(Method,
                         levels = c(method_levels, "4-Method\nConsensus")))

fill_colors <- c("Single" = npg_colors[4], "Consensus" = npg_colors[1])

p_s1a <- ggplot(summary_df,
                aes(x = Method, y = Observed, fill = Fill_group)) +
  
  geom_col(width = 0.62, alpha = 0.7, color = NA) +
  
  geom_errorbar(
    data = filter(summary_df, !is.na(Null_SD)),
    aes(x = Method,
        ymin = Null_Mean - Null_SD,
        ymax = Null_Mean + Null_SD),
    width = 0.22, linewidth = 0.7, color = "grey25",
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = filter(summary_df, !is.na(Null_SD)),
    aes(x = Method, y = Null_Mean),
    shape = 18, size = 2.8, color = "grey25",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    aes(y = Observed * 0.92, label = Observed),
    size = 3.0, family = "Arial", fontface = "bold",
    color = "white", vjust = 1
  ) +
  
  geom_text(
    aes(y = label_y, label = sig_label),
    size = 3.5, family = "Arial", color = "grey20"
  ) +
  
  scale_fill_manual(values = fill_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title    = "(a)  Observed vs Permutation Null",
    subtitle = "\u25c6 = null mean \u00b1 SD (1,000 permutations, \u03b1 = 0.05)",
    x = NULL, y = "Significant Gene Pairs"
  ) +
  theme_Fig_S1 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))


# =============================================================================
# 5. PANEL b — Fisher-z pooled CI: single vs consensus (sign-stable pairs)
# =============================================================================

hey1_pairs <- c("Hey1_Cacna1c", "Tcf7l2_Hey1", "Hey1_Cacna1d",
                "Wnt5b_Hey1",   "Hey1_Ptbp1",  "Wnt5a_Hey1")

# Sign-stable consensus pairs (21)
ci_df <- all_pairs %>%
  filter(consensus, !pair_id %in% hey1_pairs) %>%
  group_by(pair_id) %>%
  summarise(
    pool    = list(fisher_pool(r, n_per_method = n_obs)),
    r_pool  = pool[[1]]$r_pool,
    lo_pool = pool[[1]]$ci_lo,
    hi_pool = pool[[1]]$ci_hi,
    r_raw   = r[Method == "Raw\nexpression"],
    lo_raw  = CI_low[Method  == "Raw\nexpression"],
    hi_raw  = CI_high[Method == "Raw\nexpression"],
    .groups = "drop"
  ) %>%
  mutate(
    width_single     = hi_raw  - lo_raw,
    width_pool       = hi_pool - lo_pool,
    CI_reduction_pct = (width_single - width_pool) / width_single * 100,
    pair_label       = gsub("_", "\u2013", pair_id),
    pair_label       = reorder(pair_label, r_pool)
  ) %>%
  select(-pool)

mean_reduction <- mean(ci_df$CI_reduction_pct, na.rm = TRUE)

ci_long <- ci_df %>%
  select(pair_label, width_single, width_pool) %>%
  pivot_longer(c(width_single, width_pool),
               names_to = "Approach", values_to = "CI_width") %>%
  mutate(Approach = recode(Approach,
                           width_single = "Single method (Raw)",
                           width_pool   = "4-method pooled"))

p_s1b <- ggplot(ci_long,
                aes(x = pair_label, y = CI_width,
                    colour = Approach, shape = Approach)) +
  geom_line(aes(group = pair_label), colour = "grey75", linewidth = 0.35) +
  geom_point(size = 2.0, alpha = 0.88) +
  scale_colour_manual(
    values = c("Single method (Raw)" = npg_colors[1],
               "4-method pooled"     = npg_colors[4])
  ) +
  scale_shape_manual(
    values = c("Single method (Raw)" = 16,
               "4-method pooled"     = 18)
  ) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.08))) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Mean CI reduction: %.1f%%", mean_reduction),
           hjust = 1.05, vjust = 1.6,
           size = 2.8, family = "Arial",
           colour = npg_colors[4], fontface = "bold") +
  coord_flip() +
  labs(
    title    = "(b)  95% CI Width: Single vs 4-Method Pooled",
    subtitle = "21 sign-stable consensus pairs | Fisher z meta-analytic pooling",
    x = NULL, y = "95% CI width",
    colour = NULL, shape = NULL
  ) +
  theme_Fig_S1 +
  theme(
    axis.text.y     = element_text(size = 6.5, face = "italic"),
    legend.position = "bottom"
  )


# =============================================================================
# 6. PANEL c — Hey1 sign-flip: self-correcting property of consensus filter
# =============================================================================

flip_df <- all_pairs %>%
  filter(consensus, pair_id %in% hey1_pairs) %>%
  group_by(pair_id) %>%
  summarise(
    pool    = list(fisher_pool(r, n_per_method = n_obs)),
    r_pool  = pool[[1]]$r_pool,
    lo_pool = pool[[1]]$ci_lo,
    hi_pool = pool[[1]]$ci_hi,
    .groups = "drop"
  ) %>%
  select(-pool) %>%
  left_join(
    all_pairs %>%
      filter(pair_id %in% hey1_pairs) %>%
      select(pair_id, Method, r) %>%
      pivot_wider(names_from = Method, values_from = r,
                  names_glue = "r_{Method}"),
    by = "pair_id"
  ) %>%
  mutate(
    pair_label = gsub("_", "\u2013", pair_id),
    pair_label = reorder(pair_label, r_pool)
  )

# Reshape for dot plot: per-method r values + pooled
flip_long <- all_pairs %>%
  filter(pair_id %in% hey1_pairs) %>%
  mutate(pair_label = gsub("_", "\u2013", pair_id)) %>%
  bind_rows(
    flip_df %>%
      transmute(pair_id, pair_label,
                Method = factor("4-Method\nPooled"),
                r = r_pool,
                CI_low  = lo_pool,
                CI_high = hi_pool)
  ) %>%
  mutate(
    Method     = factor(Method,
                        levels = c(method_levels, "4-Method\nPooled")),
    pair_label = reorder(pair_label, ifelse(Method == "Raw\nexpression", r, NA),
                         FUN = function(x) mean(x, na.rm = TRUE))
  )

method_dot_colors <- c(
  "Raw\nexpression"     = npg_colors[7],
  "Ins1-\nnormalized"   = npg_colors[3],
  "Kcnj11-\nnormalized" = npg_colors[2],
  "Hes1-\nnormalized"   = npg_colors[4],
  "4-Method\nPooled"    = "grey30"
)

p_s1c <- ggplot(flip_long,
                aes(x = pair_label, y = r,
                    colour = Method, shape = Method)) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  # ── FIXED: single aes() with x, ymin, ymax ──────────────────────────────
  geom_errorbar(
    data = filter(flip_long, Method == "4-Method\nPooled"),
    aes(x = pair_label, ymin = CI_low, ymax = CI_high),
    width = 0.25, linewidth = 0.6, colour = "grey30",
    inherit.aes = FALSE
  ) +
  
  geom_point(size = 2.2, alpha = 0.90) +
  
  scale_colour_manual(values = method_dot_colors) +
  scale_shape_manual(
    values = c(16, 16, 16, 16, 18)
  ) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.5)) +
  coord_flip() +
  labs(
    title    = "(c)  Hey1 Pairs: Sign Reversal Across Normalization Methods",
    subtitle = "Pooled CI crosses zero \u2192 correctly excluded by consensus filter",
    x = NULL, y = "Pearson r",
    colour = NULL, shape = NULL
  ) +
  theme_Fig_S1 +
  theme(
    axis.text.y     = element_text(size = 7.5, face = "italic"),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
  )


# =============================================================================
# 7. CSV EXPORTS
# =============================================================================

write.csv(
  summary_df %>%
    mutate(Method = gsub("\n", " ", as.character(Method))) %>%
    select(Method, Observed, Null_Mean, Null_SD, Emp_p, sig_label),
  "Fig_S1a_Observed_vs_Null.csv", row.names = FALSE, na = ""
)

write.csv(
  ci_df %>%
    mutate(pair_label = as.character(pair_label)) %>%
    select(pair_label, r_pool, lo_pool, hi_pool,
           r_raw, lo_raw, hi_raw,
           width_single, width_pool, CI_reduction_pct) %>%
    arrange(desc(r_pool)),
  "Fig_S1b_CI_width_pooled_vs_single.csv", row.names = FALSE, na = ""
)

write.csv(
  flip_df %>%
    mutate(pair_label = as.character(pair_label)) %>%
    select(pair_label, r_pool, lo_pool, hi_pool) %>%
    left_join(
      all_pairs %>%
        filter(pair_id %in% hey1_pairs) %>%
        mutate(Method = gsub("\n", " ", as.character(Method))) %>%
        select(pair_id, Method, r) %>%
        pivot_wider(names_from = Method, values_from = r),
      by = c("pair_label" = "pair_id")
    ),
  "Fig_S1c_Hey1_signflip.csv", row.names = FALSE, na = ""
)

cat("CSVs exported.\n")


# =============================================================================
# 8. ASSEMBLE — 3 rows, 15 cm × 20 cm
# =============================================================================

p_Fig_S1 <- p_s1a / p_s1b / p_s1c +
  plot_layout(heights = c(1, 1.2, 1)) +
  plot_annotation(
    #title    = "Fig. S1 \u2502 Multi-Normalization Consensus: Statistical Validation",
    subtitle = paste0(
      #"n = ", n_obs,
      #" (4 biological replicates \u00d7 3 conditions) \u2502 ",
      #length(consensus_ids),
      #" consensus pairs significant across all 4 normalization methods"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 10,
                                   family = "Arial", colour = "grey10"),
      plot.subtitle = element_text(size = 8, family = "Arial",
                                   colour = "grey40")
    )
  )


# =============================================================================
# 9. SAVE
# =============================================================================

ggsave("Fig_S1_Validation.svg",
       plot   = p_Fig_S1,
       width  = 15, height = 20, units = "cm",
       dpi    = 300)

ggsave("Fig_S1_Validation.png",
       plot   = p_Fig_S1,
       width  = 15, height = 20, units = "cm",
       dpi    = 300)

cat(sprintf(
  "\nFig. S1 saved (3 rows, 15 \u00d7 20 cm).\n  n=%d | consensus=%d | mean CI reduction=%.1f%%\n",
  n_obs, length(consensus_ids), mean_reduction
))

