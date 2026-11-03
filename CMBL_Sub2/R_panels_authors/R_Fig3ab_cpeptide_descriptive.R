# Fig. 3a and 3b redrawn from the raw ELISA workbook: four culture wells per group and glucose concentration (as stated by the authors), raw signal instead of the LOESS curve.
# No LOESS smoothing and no confidence bands. Panel (b) carries the pairwise bars of the other panels (same procedure as the authors' Fig_3b script); the four wells
# per group are technical replicates from one preparation, so the tests are descriptive and within-preparation like all others in the paper.
# Same theme as the other figures: ggplot2 theme_minimal(), NPG treatment colours, 8-pt text.
#   (a) raw ELISA signal (dOD = OD450 - OD590) at 0, 2.8, 5.5, 22 and 44 mM; points = the four readings, line = their mean;
#       dashed lines mark the blank (0.026) and the lowest standard (100 pM, 0.044)
#   (b) at 5.5 mM: raw signal (left) and value per ug of group-mean DNA in workbook units (right); both axes start at zero
# Input: Cpeptide_raw_elisa_tidy.csv from elisa_raw_audit.py.   Usage: Rscript R_Fig3ab_cpeptide_descriptive.R <tidy.csv> <out_dir>
# Output: Fig_3a_*.png/.svg, Fig_3b_*.png/.svg and Fig_3b_pairwise_tests.csv (the numbers behind the bars of panel b)
suppressPackageStartupMessages({ library(ggplot2); library(ggsci) })
args <- commandArgs(trailingOnly = TRUE); tidy <- read.csv(args[1]); OUT <- args[2]
npg <- pal_npg("nrc")(10); cols <- c("CTRL" = npg[1], "DAPT" = npg[2], "DKK-1" = npg[3])
tidy$group <- factor(tidy$group, levels = names(cols)); tidy$Assay <- "Glucose-stimulated C-peptide secretion (GSCS)"
lev <- c(0, 2.8, 5.5, 22, 44); tidy$glc <- factor(tidy$condition_mM, levels = lev, labels = c("0", "2.8", "5.5", "22", "44"))
th <- theme_minimal() +
  theme(legend.position = "top", legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.5, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = margin(5),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8),
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5))

# ---- (a) ----------------------------------------------------------------------------------------------------------------------------------
mn <- aggregate(dOD ~ group + glc, tidy, mean); mn$Assay <- "Glucose-stimulated C-peptide secretion (GSCS)"
refs <- data.frame(y = c(0.026, 0.044), lab = c("blank", "lowest standard"))
pa <- ggplot(tidy, aes(glc, dOD, colour = group, group = group)) +
  geom_hline(data = refs, aes(yintercept = y), linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_line(data = mn, linewidth = 0.8) +
  geom_point(size = 1.5, alpha = 0.5) +
  facet_grid(cols = vars(Assay)) +
  scale_x_discrete() +
  scale_colour_manual(values = cols, name = "Treatment") +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.04)) +
  labs(x = "Glucose concentration (mM)", y = "\nELISA signal (ΔOD, 450−590 nm)") + th
ggsave(file.path(OUT, "Fig_3a_ELISA_signal_by_glucose.png"), pa, width = 8, height = 6.5, units = "cm", dpi = 600)
ggsave(file.path(OUT, "Fig_3a_ELISA_signal_by_glucose.svg"), pa, width = 8, height = 6.5, units = "cm")

# ---- (b) ----------------------------------------------------------------------------------------------------------------------------------
d55 <- tidy[tidy$condition_mM == 5.5, ]
long <- rbind(data.frame(group = d55$group, y = d55$dOD, panel = "Signal (ΔOD)"),
              data.frame(group = d55$group, y = d55$ng_per_ug_groupmeanDNA, panel = "Per µg DNA"))
long$panel <- factor(long$panel, levels = c("Signal (ΔOD)", "Per µg DNA"))
sm <- aggregate(y ~ group + panel, long, mean)
# Pairwise bars exactly as in the other panels and in the authors' Fig_3b_C-peptide_5.5mm_revised.r: a Shapiro-Wilk test of the twelve values guides the omnibus
# test (p < 0.05: Kruskal-Wallis with Bonferroni-corrected pairwise Wilcoxon rank-sum tests; otherwise one-way ANOVA with Bonferroni-corrected pairwise t-tests);
# "ns#" = uncorrected p < 0.05 but Bonferroni-adjusted p >= 0.05.
# Note: the four wells per group are technical replicates from one preparation (descriptive, within-preparation tests).
sig <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
pairs_df <- data.frame(g1 = c("CTRL", "CTRL", "DAPT"), g2 = c("DAPT", "DKK-1", "DKK-1"), x1 = c(1, 1, 2), x2 = c(2, 3, 3), lv = c(1.12, 1.27, 1.42))
res <- lapply(levels(long$panel), function(pn) {
  v <- long[long$panel == pn, ]; sw <- shapiro.test(v$y)
  if (sw$p.value < 0.05) {
    g <- kruskal.test(y ~ group, v); meth <- "Kruskal-Wallis"; gp <- g$p.value; gs <- unname(g$statistic)
    pa_ <- suppressWarnings(pairwise.wilcox.test(v$y, v$group, p.adjust.method = "bonferroni"))$p.value
    pr_ <- suppressWarnings(pairwise.wilcox.test(v$y, v$group, p.adjust.method = "none"))$p.value
  } else {
    a <- summary(aov(y ~ group, v))[[1]]; meth <- "Anova"; gp <- a[["Pr(>F)"]][1]; gs <- a[["F value"]][1]
    pa_ <- pairwise.t.test(v$y, v$group, p.adjust.method = "bonferroni")$p.value
    pr_ <- pairwise.t.test(v$y, v$group, p.adjust.method = "none")$p.value
  }
  dmax <- max(v$y); top <- dmax * 1.78
  b <- pairs_df; b$panel <- pn
  b$p_raw <- mapply(function(a1, a2) pr_[a2, a1], b$g1, b$g2); b$p_adj <- mapply(function(a1, a2) pa_[a2, a1], b$g1, b$g2)
  b$lab <- ifelse(b$p_raw < 0.05 & b$p_adj >= 0.05, "ns#", sig(b$p_adj))
  b$y <- dmax * b$lv; b$top <- top; b$tip <- 0.03 * top; b$gy <- dmax * 1.64
  b$method <- meth; b$omnibus_stat <- gs; b$omnibus_p <- gp; b$shapiro_W <- unname(sw$statistic); b$shapiro_p <- sw$p.value
  b$gl <- paste0('"', meth, ', "*italic(p)*" = ', format.pval(gp, digits = 2), '"')
  b
})
br <- do.call(rbind, res); br$panel <- factor(br$panel, levels = levels(long$panel))
write.csv(br[, c("panel", "g1", "g2", "method", "shapiro_W", "shapiro_p", "omnibus_stat", "omnibus_p", "p_raw", "p_adj", "lab")],
          file.path(OUT, "Fig_3b_pairwise_tests.csv"), row.names = FALSE)
gl <- unique(br[, c("panel", "gl", "gy")])
pb <- ggplot(long, aes(group, y, colour = group)) +
  geom_point(position = position_jitter(width = 0.1, height = 0, seed = 1), size = 1.5, alpha = 0.6) +
  geom_crossbar(data = sm, aes(ymin = y, ymax = y), width = 0.5, linewidth = 0.3) +
  geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, colour = "black", linewidth = 0.3) +
  geom_segment(data = br, aes(x = x1, xend = x1, y = y, yend = y - tip), inherit.aes = FALSE, colour = "black", linewidth = 0.3) +
  geom_segment(data = br, aes(x = x2, xend = x2, y = y, yend = y - tip), inherit.aes = FALSE, colour = "black", linewidth = 0.3) +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = y + 0.01 * top, label = lab), inherit.aes = FALSE, colour = "black", size = 2.5, vjust = 0) +
  geom_text(data = gl, aes(x = 2, y = gy, label = gl), inherit.aes = FALSE, colour = "black", size = 2.5, vjust = 0, parse = TRUE) +
  geom_blank(data = br, aes(x = x1, y = top), inherit.aes = FALSE) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = cols) + expand_limits(y = 0) +
  labs(x = NULL, y = "\n5.5 mM glucose") + th + theme(legend.position = "none", axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5))
# 8 x 5.46 cm: the height of the slot that panel (b) occupies in the composite of Fig. 3 (no legend row, no x-axis title)
ggsave(file.path(OUT, "Fig_3b_C-peptide_5.5mM.png"), pb, width = 8, height = 5.46, units = "cm", dpi = 600)
ggsave(file.path(OUT, "Fig_3b_C-peptide_5.5mM.svg"), pb, width = 8, height = 5.46, units = "cm")
cat("done\n")
