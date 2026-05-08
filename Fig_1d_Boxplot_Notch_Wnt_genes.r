# Figure 1d — Boxplots: Notch targets & Wnt pathway gene expression
# Theme: matched exactly to Fig_2cdef | Layout: 5 columns × 2 rows | 16.2 × 12.8 cm

rm(list = ls())

library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(rstatix)
library(patchwork)
library(rlang)

npg_colors <- pal_npg("nrc")(10)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
dat <- read_excel("Fig_1d_Fold_change_mRNA_expression.xlsx", sheet = 1)

if (!"Gene" %in% colnames(dat)) {
  dat <- dat %>%
    pivot_longer(cols = -c(Treatment, Sample),
                 names_to  = "Gene",
                 values_to = "Fold_change")
}

gene_order   <- c("Hes1", "Hey1", "Wnt2", "Wnt2b", "Wnt5a",
                  "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2")
treat_levels <- c("CTRL", "DAPT", "DKK-1")

dat <- dat %>%
  filter(Gene %in% gene_order) %>%
  mutate(Gene        = factor(Gene,      levels = gene_order),
         Treatment   = factor(Treatment, levels = treat_levels),
         Fold_change = as.numeric(Fold_change))

# ── HELPER: STATISTICS + CSV EXPORT ──────────────────────────────────────────
run_stats <- function(df, gene_name, parametric = TRUE) {
  
  x   <- df$Fold_change
  grp <- df$Treatment
  
  sw <- tryCatch(shapiro.test(x),
                 error   = function(e) list(statistic = NA, p.value = NA),
                 warning = function(w) list(statistic = NA, p.value = NA))
  
  if (parametric) {
    fit           <- lm(Fold_change ~ Treatment, data = df)
    fit_anv       <- anova(fit)
    omni_stat     <- round(as.numeric(fit_anv[["F value"]][1]), 4)
    omni_p        <- round(as.numeric(fit_anv[["Pr(>F)"]][1]), 4)
    omni_name     <- "One-way ANOVA (F)"
    pw            <- pairwise.t.test(x, grp, p.adjust.method = "bonferroni")
    pw_raw        <- pairwise.t.test(x, grp, p.adjust.method = "none")
    post_hoc_name <- "Pairwise t-test (Bonferroni)"
  } else {
    omni          <- kruskal.test(Fold_change ~ Treatment, data = df)
    omni_stat     <- round(as.numeric(omni$statistic), 4)
    omni_p        <- round(as.numeric(omni$p.value), 4)
    omni_name     <- "Kruskal-Wallis (H)"
    pw            <- pairwise.wilcox.test(x, grp, p.adjust.method = "bonferroni")
    pw_raw        <- pairwise.wilcox.test(x, grp, p.adjust.method = "none")
    post_hoc_name <- "Pairwise Wilcoxon (Bonferroni)"
  }
  
  # Bonferroni-corrected p-values
  pw_df <- as.data.frame(as.table(pw$p.value))
  colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
  pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                               ifelse(pw_df$p_adjusted < 0.01,  "**",
                                      ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))
  
  # Raw (uncorrected) p-values
  pw_raw_df <- as.data.frame(as.table(pw_raw$p.value))
  colnames(pw_raw_df) <- c("Group1", "Group2", "p_raw")
  pw_raw_df <- pw_raw_df[!is.na(pw_raw_df$p_raw), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  
  # Merge raw p into pw_df and flag trends
  pw_df <- merge(pw_df, pw_raw_df[, c("Group1", "Group2", "p_raw")],
                 by = c("Group1", "Group2"), all.x = TRUE)
  pw_df$Trend <- ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted >= 0.05, "#",
                        ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted < 0.05,
                               pw_df$Significance, "ns"))
  
  desc <- df %>%
    group_by(Treatment) %>%
    summarise(Mean = mean(Fold_change), SEM = sd(Fold_change) / sqrt(n()),
              N = n(), .groups = "drop")
  
  sec1 <- data.frame(
    Section = c(paste0("Descriptive Statistics - ", gene_name), rep("", nrow(desc))),
    Col1    = c("Treatment",  as.character(desc$Treatment)),
    Col2    = c("Mean (AU)",  as.character(round(desc$Mean, 4))),
    Col3    = c("SEM",        as.character(round(desc$SEM,  4))),
    Col4    = c("N",          as.character(desc$N)),
    Col5    = NA,
    Col6    = NA
  )
  
  stat_df <- data.frame(
    Test = c(paste0("Shapiro-Wilk (", gene_name, ")"), omni_name),
    Statistic = c(round(as.numeric(sw$statistic), 4), omni_stat),
    p_value   = c(round(as.numeric(sw$p.value),   4), omni_p),
    Interpretation = c(
      ifelse(is.na(sw$p.value), "Undetermined",
             ifelse(sw$p.value > 0.05, "Normal", "Non-normal")),
      ifelse(is.na(omni_p), "Undetermined",
             ifelse(omni_p < 0.05, "Significant", "Not significant"))
    )
  )
  
  sec2 <- data.frame(
    Section = c(paste0("Normality & Omnibus - ", gene_name), rep("", nrow(stat_df))),
    Col1    = c("Test",           as.character(stat_df$Test)),
    Col2    = c("Statistic",      as.character(stat_df$Statistic)),
    Col3    = c("p_value",        as.character(stat_df$p_value)),
    Col4    = c("Interpretation", as.character(stat_df$Interpretation)),
    Col5    = NA,
    Col6    = NA
  )
  
  sec3 <- data.frame(
    Section = c(paste0("Pairwise - ", gene_name, " (", post_hoc_name, ")"),
                rep("", nrow(pw_df))),
    Col1    = c("Group1",       as.character(pw_df$Group1)),
    Col2    = c("Group2",       as.character(pw_df$Group2)),
    Col3    = c("p_raw",        as.character(round(pw_df$p_raw,      6))),
    Col4    = c("p_adjusted",   as.character(round(pw_df$p_adjusted, 6))),
    Col5    = c("Significance", as.character(pw_df$Significance)),
    Col6    = c("Trend_flag",   as.character(pw_df$Trend))
  )
  
  spacer <- data.frame(Section = "---", Col1 = NA, Col2 = NA, Col3 = NA,
                       Col4 = NA, Col5 = NA, Col6 = NA)
  out    <- rbind(sec1, spacer, sec2, spacer, sec3, spacer)
  
  list(csv = out, pw = pw_df)
}

# ── NORMALITY SCREEN (non-CTRL only) ─────────────────────────────────────────
safe_shapiro_p <- function(x) {
  x <- x[!is.na(x)]
  if (length(unique(x)) < 3 || length(x) < 3) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value,
           error   = function(e) NA_real_,
           warning = function(w) NA_real_)
}

parametric_map <- sapply(gene_order, function(g) {
  pvals <- dat %>%
    filter(Gene == g, Treatment != "CTRL") %>%
    group_by(Treatment) %>%
    summarise(p = safe_shapiro_p(Fold_change), .groups = "drop") %>%
    pull(p)
  all(pvals > 0.05, na.rm = TRUE)
})

cat("Normality screen (non-CTRL groups):\n")
print(data.frame(
  Gene       = gene_order,
  Parametric = unname(parametric_map),
  Test       = ifelse(unname(parametric_map),
                      "One-way ANOVA + Bonferroni",
                      "Kruskal-Wallis + Wilcoxon Bonferroni")
))

# ── RUN STATISTICS PER GENE ──────────────────────────────────────────────────
gene_dfs <- lapply(gene_order, function(g) dat %>% filter(Gene == g))
names(gene_dfs) <- gene_order

res <- lapply(gene_order, function(g) {
  cat("Processing:", g, "| parametric =", unname(parametric_map[g]), "\n")
  run_stats(gene_dfs[[g]], g, unname(parametric_map[g]))
})
names(res) <- gene_order

combined_csv <- do.call(rbind, lapply(res, `[[`, "csv"))
write.csv(combined_csv, "Fig_1d_stats.csv", row.names = FALSE)

# ── HELPER: build pw_plot with ns# embedded ───────────────────────────────────
make_pw_plot <- function(pw_df, y_positions) {
  pw_df %>%
    rename(group1 = Group1, group2 = Group2) %>%
    mutate(
      p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                            "ns#",
                            Significance),
      y.position   = y_positions
    )
}

# ── HELPER: boxplot builder ───────────────────────────────────────────────────
make_boxplot <- function(df, pw_plot, gene_lab, gene_key,
                         omnibus_method, y_lab_pos, y_lim_top,
                         show_y_title = TRUE,
                         show_y_text  = TRUE,
                         summary_y    = 0.25) {
  labs_vec <- setNames(gene_lab, gene_key)
  
  p <- ggboxplot(df,
                 x       = "Treatment",
                 y       = "Fold_change",
                 palette = "npg",
                 color   = "Treatment",
                 fill    = "Treatment",
                 alpha   = 0.5,
                 add     = "jitter",
                 ylim    = c(-0.1, y_lim_top)) +
    labs(y = if (show_y_title) "Fold change in\nmRNA expression (AU)" else NULL) +
    facet_wrap(vars(Gene), labeller = labeller(Gene = labs_vec)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text      = element_text(size = 8, colour = "black", face = "bold.italic"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = 8, colour = "black",
                                     angle = 45, hjust = 1, vjust = 1),
      axis.text.y     = if (show_y_text)
        element_text(size = 8, colour = "black", angle = 90, hjust = 0.5)
      else
        element_blank(),
      axis.title.y    = if (show_y_title) element_text(size = 8) else element_blank(),
      plot.margin     = margin(0, 0, 0, 0)
    ) +
    stat_compare_means(
      label.y     = y_lab_pos,
      label.x.npc = "center",
      size        = 2.5,
      method      = omnibus_method,
      hjust       = 0.5
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          y     = summary_y,
          label = paste0(round(mean(x), 2), "\n\u00b1 ",
                         round(sd(x) / sqrt(length(x)), 2))
        )
      },
      geom  = "text",
      vjust = 1,
      size  = 2.5,
      angle = 0
    )
  
  if (nrow(pw_plot) > 0) {
    p <- p + stat_pvalue_manual(pw_plot,
                                label      = "p.adj.signif",
                                tip.length = 0.03,
                                size       = 2.5,
                                step.increase = 0.001)
  }
  
  return(p)
}

save_plot <- function(plot, filename) {
  ggsave(paste0(filename, ".png"), plot = plot, scale = 1,
         width = 8, height = 5.3, units = "cm", dpi = 600, limitsize = TRUE)
  ggsave(paste0(filename, ".svg"), plot = plot, scale = 1,
         width = 8, height = 5.3, units = "cm", dpi = 600, limitsize = TRUE)
}

# ── GENE LABELS ──────────────────────────────────────────────────────────────
gene_labels <- c(
  "Hes1"   = "Hes1 (Notch target)",
  "Hey1"   = "Hey1 (Notch target)",
  "Wnt2"   = "Wnt2 (Wnt ligand)",
  "Wnt2b"  = "Wnt2b (Wnt ligand)",
  "Wnt5a"  = "Wnt5a (NC-Wnt)",
  "Wnt5b"  = "Wnt5b (NC-Wnt)",
  "Wnt9a"  = "Wnt9a (Wnt ligand)",
  "Lef1"   = "Lef1 (Wnt TF)",
  "Tcf7"   = "Tcf7 (Wnt TF)",
  "Tcf7l2" = "Tcf7l2 (Wnt TF)"
)

y_title_genes <- c("Hes1", "Wnt5b")
y_text_genes  <- c("Hes1", "Wnt5b")

# ── PASS 1: per-gene y_lim_top candidates ────────────────────────────────────
gene_y_info <- lapply(gene_order, function(g) {
  df      <- gene_dfs[[g]]
  pw_df   <- res[[g]]$pw
  n_pairs <- nrow(pw_df)
  y_max   <- max(df$Fold_change, na.rm = TRUE)
  
  if (n_pairs > 0) {
    y_pos     <- seq(y_max * 1.12, y_max * 1.50, length.out = n_pairs)
    y_lim_top <- max(y_pos) * 1.15
  } else {
    y_pos     <- numeric(0)
    y_lim_top <- y_max * 1.65
  }
  list(y_lim_top = y_lim_top)
})
names(gene_y_info) <- gene_order

# ── PASS 2: shared ceiling & bracket positions ────────────────────────────────
global_y_lim_top <- max(sapply(gene_y_info, `[[`, "y_lim_top"))
max_n_pairs      <- max(sapply(gene_order, function(g) nrow(res[[g]]$pw)))

global_y_pos_all <- seq(global_y_lim_top * 0.62,
                        global_y_lim_top * 0.87,
                        length.out = max(max_n_pairs, 1))
global_y_lab_pos <- global_y_lim_top * 0.97

cat("Global y_lim_top:", round(global_y_lim_top, 3),
    "| bracket range:", round(global_y_pos_all[1], 3), "-",
    round(global_y_pos_all[length(global_y_pos_all)], 3),
    "| label pos:", round(global_y_lab_pos, 3), "\n")

# ── BUILD PLOTS ──────────────────────────────────────────────────────────────
plot_list <- list()

for (g in gene_order) {
  df       <- gene_dfs[[g]]
  pw_df    <- res[[g]]$pw
  is_param <- unname(parametric_map[g])
  omni_m   <- ifelse(is_param, "anova", "kruskal.test")
  n_pairs  <- nrow(pw_df)
  y_pos_g  <- if (n_pairs > 0) global_y_pos_all[seq_len(n_pairs)] else numeric(0)
  
  pw_plot <- make_pw_plot(pw_df, y_pos_g)
  
  plot_list[[g]] <- make_boxplot(
    df             = df,
    pw_plot        = pw_plot,
    gene_lab       = gene_labels[g],
    gene_key       = g,
    omnibus_method = omni_m,
    y_lab_pos      = global_y_lab_pos,
    y_lim_top      = global_y_lim_top,
    show_y_title   = g %in% y_title_genes,
    show_y_text    = g %in% y_text_genes,
    summary_y      = global_y_lim_top * 0.07
  )
  
  cat("Plot built:", g,
      "| n_brackets =", n_pairs,
      "| bracket y =", paste(round(y_pos_g, 2), collapse = ", "), "\n")
}

# ── SAVE INDIVIDUAL PANELS ───────────────────────────────────────────────────
panel_files <- c(
  "Hes1"   = "Fig_1d_1_Boxplot_Hes1",
  "Hey1"   = "Fig_1d_2_Boxplot_Hey1",
  "Wnt2"   = "Fig_1d_3_Boxplot_Wnt2",
  "Wnt2b"  = "Fig_1d_4_Boxplot_Wnt2b",
  "Wnt5a"  = "Fig_1d_5_Boxplot_Wnt5a",
  "Wnt5b"  = "Fig_1d_6_Boxplot_Wnt5b",
  "Wnt9a"  = "Fig_1d_7_Boxplot_Wnt9a",
  "Lef1"   = "Fig_1d_8_Boxplot_Lef1",
  "Tcf7"   = "Fig_1d_9_Boxplot_Tcf7",
  "Tcf7l2" = "Fig_1d_10_Boxplot_Tcf7l2"
)

for (g in gene_order) save_plot(plot_list[[g]], panel_files[g])

# ── ASSEMBLE 5×2 COMBINED FIGURE ─────────────────────────────────────────────
combined <- wrap_plots(plot_list, ncol = 5, nrow = 2)

ggsave("Fig_1d_Boxplot_combined_expression.png", plot = combined, scale = 1,
       width = 16.6, height = 12.8, units = "cm", dpi = 600, limitsize = TRUE)
ggsave("Fig_1d_Boxplot_combined_expression.svg", plot = combined, scale = 1,
       width = 16.6, height = 12.8, units = "cm", dpi = 600, limitsize = TRUE)

cat("\nAll Fig 1d files saved.\nStats: Fig_1d_stats.csv\n")

