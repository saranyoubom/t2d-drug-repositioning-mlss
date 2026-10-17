# Supplementary Fig. S2a redrawn in R with one label format for all eight weight scenarios (C:B:V:P, %).
# The earlier figure (Fig_6d_Fig_S2_MLSS_validation.r) mixed 3-part labels (C:B:V+P lumped) and 4-part labels, e.g. "C-Heavy (60:25:15)".
# Same data (input_data/Fig_6c_MLSS_ALL_Combinations.csv), same calculation (base score x 9 + synergy bonus - antagonism penalty,
# top 15 combinations, min-max normalization within each scenario), same theme (ggplot2 theme_minimal(), ggsci NPG, Arial, bold axis titles).
# Usage: Rscript R_SuppFig2a_labels_npg.R [input_csv] [out_dir]
suppressPackageStartupMessages({ library(ggplot2); library(ggsci); library(reshape2) })
args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "input_data/Fig_6c_MLSS_ALL_Combinations.csv"
OUT <- if (length(args) >= 2) args[2] else "../Figures/Supplementary_figures"
npg <- pal_npg("nrc")(10)

d <- read.csv(IN, stringsAsFactors = FALSE)
top <- d[order(d$Rank), ][1:15, ]
sc <- data.frame(
  Scenario = c("Standard (45:30:10:15)", "C-Heavy 1 (60:25:10:05)", "C-Heavy 2 (55:20:15:10)", "B-Heavy (35:50:10:05)",
               "V-Heavy (35:30:30:05)", "P-Heavy (30:20:10:40)", "Balanced (40:40:10:10)", "Network-Focus (30:20:40:10)"),
  C = c(.45, .60, .55, .35, .35, .30, .40, .30), B = c(.30, .25, .20, .50, .30, .20, .40, .20),
  V = c(.10, .10, .15, .10, .30, .10, .10, .40), P = c(.15, .05, .10, .05, .05, .40, .10, .10), stringsAsFactors = FALSE)
stopifnot(all(abs(rowSums(sc[, c("C", "B", "V", "P")]) - 1) < 1e-12))

mat <- sapply(seq_len(nrow(sc)), function(i)
  (sc$C[i] * top$Complementarity_C + sc$B[i] * top$Balance_B + sc$V[i] * top$Network_Coverage_V + sc$P[i] * top$Potency_P) * 9 +
    top$Synergy_Bonus - top$Antagonism_Penalty)
colnames(mat) <- sc$Scenario; rownames(mat) <- top$Combination
stopifnot(max(abs(mat[, 1] - top$MLSS_v4.0)) < 1e-9)              # standard scheme reproduces the reported MLSS
rk <- apply(mat, 2, rank); rho <- cor(rk, method = "spearman")
cat("pairwise Spearman rho between the 8 schemes: min =", format(min(rho[upper.tri(rho)]), digits = 4),
    "| pairs =", sum(upper.tri(rho)), "\n")
stopifnot(sum(upper.tri(rho)) == 28, min(rho[upper.tri(rho)]) > 0.9999)

norm <- apply(mat, 2, function(x) (x - min(x)) / (max(x) - min(x)))
long <- melt(norm); colnames(long) <- c("Combination", "Scenario", "Score")
long$Scenario    <- factor(long$Scenario, levels = sc$Scenario)
long$Combination <- factor(long$Combination, levels = top$Combination)   # rank 1 at the bottom, as in the earlier figure

p <- ggplot(long, aes(Scenario, Combination, fill = Score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = npg[6], high = npg[1], name = "Normalized\nMLSS score", limits = c(0, 1)) +
  labs(title = "Weight scenario sensitivity (top 15)", x = "Weight scenario (C:B:V:P, %)", y = "Drug combination") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black", family = "Arial"),
        axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
        axis.title = element_text(size = 12, colour = "black", face = "bold", family = "Arial"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
        legend.text = element_text(size = 12, colour = "black", family = "Arial"),
        legend.title = element_text(size = 12, colour = "black", face = "bold", family = "Arial"),
        panel.grid = element_blank())
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "Fig_S2a_mlss_sensitivity_top15.png"), p, width = 17.7, height = 13.8, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_S2a_mlss_sensitivity_top15.svg"), p, width = 17.7, height = 13.8, units = "cm", bg = "white")
cat("written to", normalizePath(OUT), "\n")
