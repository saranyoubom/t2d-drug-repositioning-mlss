# Supplementary Fig. S1 in the theme of the main figures (ggplot2 theme_minimal(), ggsci NPG palette, Arial, black text, bold axis titles).
#   (a) observed pooled counts against the group-preserving permutation null (R_null_counts.csv, 2000 permutations)
#   (b) significant pairs before and after treatment adjustment (R_all_pairs_statistics.csv)
#   (c) pooled correlations of the six Hey1 consensus pairs under the four schemes (R_SuppTable3_v2_consensus_pairs.csv)
# Same data and checks as build_supp_fig1.py (matplotlib version, replaced because it did not match the theme of the other figures).
# Usage: Rscript R_SuppFig1_theme_npg.R [analysis_dir] [out_dir]
suppressPackageStartupMessages({ library(ggplot2); library(ggsci); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
A   <- if (length(args) >= 1) args[1] else "."
OUT <- if (length(args) >= 2) args[2] else "."
npg <- pal_npg("nrc")(10)       # 1 red, 2 cyan, 3 green, 4 dark blue, 5 salmon, 6 grey-blue
schemes <- c("Raw", "Hes1", "Kcnj11", "Ins1")
lab_s   <- c(Raw = "Raw\nexpression", Hes1 = "Hes1-\nnormalized", Kcnj11 = "Kcnj11-\nnormalized", Ins1 = "Ins1-\nnormalized")

null  <- read.csv(file.path(A, "R_null_counts.csv"))
pairs <- read.csv(file.path(A, "R_all_pairs_statistics.csv"), stringsAsFactors = FALSE)
sc    <- read.csv(file.path(A, "R_SuppTable3_v2_consensus_pairs.csv"), stringsAsFactors = FALSE)

cnt <- function(col, thr = 0.05) vapply(schemes, function(s) sum(pairs[pairs$scheme == s, col] < thr), 1L)
obs_p <- cnt("p_pooled_df10"); obs_q <- cnt("q_pooled_BH"); adj_p <- cnt("p_adjusted_df8"); adj_q <- cnt("q_adjusted_BH")
key <- apply(cbind(pairs$gene1, pairs$gene2), 1, function(x) paste(sort(x), collapse = "|"))
full <- names(which(table(key) == 4))
consensus <- function(col) sum(vapply(full, function(k) all(pairs[key == k, col] < 0.05), TRUE))
cons_p <- consensus("p_pooled_df10"); cons_q <- consensus("q_pooled_BH"); cons_ap <- consensus("p_adjusted_df8"); cons_aq <- consensus("q_adjusted_BH")
obs_cons <- nrow(sc); obs_sign <- sum(sc$sign_consistent_pooled %in% c(TRUE, "TRUE"))
stopifnot(identical(unname(c(obs_p, cons_p, obs_cons, obs_sign)), c(80L, 169L, 118L, 98L, 27L, 27L, 21L)))

base <- theme_minimal() +
  theme(text = element_text(family = "Arial"),
        axis.text = element_text(size = 8, colour = "black", family = "Arial"),
        axis.title = element_text(size = 8, face = "bold", colour = "black", family = "Arial"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5, colour = "black", family = "Arial"),
        legend.text = element_text(size = 8, colour = "black", family = "Arial"),
        legend.title = element_blank(), legend.position = "top",
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.4, "cm"),
        plot.tag = element_text(size = 10, face = "bold", family = "Arial"), plot.tag.position = c(0.005, 0.995))

# ---- null summary helper: thin line 2.5-97.5 %, thick bar IQR, diamond = observed -------------------------------
summ <- function(df) do.call(rbind, lapply(split(df, df$grp), function(d) {
  q <- quantile(d$n, c(0.025, 0.25, 0.75, 0.975)); data.frame(grp = d$grp[1], lo = q[1], q1 = q[2], q3 = q[3], hi = q[4], obs = d$obs[1]) }))
violin_panel <- function(df, xlab_levels, title, ylab) {
  df$grp <- factor(df$grp, levels = xlab_levels); sm <- summ(df); sm$grp <- factor(sm$grp, levels = xlab_levels)
  ggplot(df, aes(grp, n)) +
    geom_violin(fill = npg[6], alpha = 0.40, colour = NA, scale = "width", width = 0.85) +
    geom_segment(data = sm, aes(x = grp, xend = grp, y = lo, yend = hi), colour = npg[4], linewidth = 0.4, inherit.aes = FALSE) +
    geom_segment(data = sm, aes(x = grp, xend = grp, y = q1, yend = q3), colour = npg[4], linewidth = 1.6, inherit.aes = FALSE) +
    geom_point(data = sm, aes(grp, obs), shape = 23, fill = "white", colour = "black", size = 2.6, stroke = 0.7, inherit.aes = FALSE) +
    geom_text(data = sm, aes(grp, obs, label = obs), nudge_x = 0.33, size = 2.8, family = "Arial", inherit.aes = FALSE) +
    labs(title = title, x = NULL, y = ylab) + base + theme(legend.position = "none")
}
va <- do.call(rbind, lapply(schemes, function(s) data.frame(grp = s, n = null[[s]], obs = obs_p[[s]])))
pa1 <- violin_panel(va, schemes, "Observed counts vs group-preserving null", expression(bold(paste("Pairs with pooled ", bolditalic(p), " < 0.05")))) +
  scale_x_discrete(labels = lab_s) + labs(tag = "(a)")
vc <- rbind(data.frame(grp = "All four\nschemes", n = null$consensus, obs = obs_cons), data.frame(grp = "+ same\nsign", n = null$consensus_sign_consistent, obs = obs_sign))
pa2 <- violin_panel(vc, c("All four\nschemes", "+ same\nsign"), "Consensus counts", "Consensus pairs")

# ---- (b) pooled vs adjusted counts ---------------------------------------------------------------------------------
cats <- c(unname(lab_s[schemes]), "All four\nschemes")
sers <- c("pooled p", "pooled q", "adjusted p", "adjusted q")
bdf <- data.frame(cat = factor(rep(cats, 4), levels = cats), series = factor(rep(sers, each = 5), levels = sers),
                  n = c(obs_p, cons_p, obs_q, cons_q, adj_p, cons_ap, adj_q, cons_aq))
pb <- ggplot(bdf, aes(cat, n, fill = series)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76) +
  geom_text(aes(label = n), position = position_dodge(width = 0.82), hjust = -0.25, angle = 90, size = 2.5, family = "Arial") +
  scale_fill_manual(values = c(npg[4], npg[2], npg[5], npg[1]),
                    labels = c(expression(paste("pooled ", italic(p), " < 0.05")), expression(paste("pooled BH ", italic(q), " < 0.05")),
                               expression(paste("adjusted ", italic(p), " < 0.05")), expression(paste("adjusted BH ", italic(q), " < 0.05")))) +
  scale_y_continuous(limits = c(0, 215), expand = expansion(mult = c(0, 0.02))) +
  labs(title = expression(bold(paste("Significant pairs before and after adjustment for treatment (", bolditalic(df), " = 8); only |", bolditalic(r), "| ≥ 0.91 can pass BH at ", bolditalic(q), " < 0.05"))),
       x = NULL, y = "Significant pairs", tag = "(b)") + base

# ---- (c) Hey1 pairs ---------------------------------------------------------------------------------------------------
hey <- sc[sc$gene1 == "Hey1" | sc$gene2 == "Hey1", ]
hey <- hey[order(hey$r_pooled_Raw), ]
hey$other <- ifelse(hey$gene2 == "Hey1", hey$gene1, hey$gene2)
hl <- do.call(rbind, lapply(seq_len(nrow(hey)), function(i) data.frame(pair = hey$other[i], scheme = schemes, r = unlist(hey[i, paste0("r_pooled_", schemes)], use.names = FALSE))))
hl$pair <- factor(hl$pair, levels = hey$other); hl$scheme <- factor(hl$scheme, levels = schemes)
hl$sign <- factor(ifelse(hl$r > 0, "positive", "negative"), levels = c("positive", "negative"))
rng <- aggregate(r ~ pair, hl, range); rng <- data.frame(pair = rng$pair, lo = rng$r[, 1], hi = rng$r[, 2])
ylab <- parse(text = sprintf("italic('Hey1')*' – '*italic('%s')", levels(hl$pair)))
pc <- ggplot(hl, aes(r, pair)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = rng, aes(x = lo, xend = hi, y = pair, yend = pair), colour = "grey70", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(aes(shape = scheme, fill = sign), size = 2.6, colour = "black", stroke = 0.35) +
  scale_shape_manual(values = c(Raw = 21, Hes1 = 22, Kcnj11 = 24, Ins1 = 23), labels = c("Raw expression", "Hes1-normalized", "Kcnj11-normalized", "Ins1-normalized")) +
  scale_fill_manual(values = c(positive = npg[1], negative = npg[4]), labels = c(expression(italic(r) > 0), expression(italic(r) < 0))) +
  scale_y_discrete(labels = ylab) + scale_x_continuous(limits = c(-1.05, 1.05), breaks = seq(-1, 1, 0.5)) +
  guides(fill = guide_legend(order = 2, nrow = 1, override.aes = list(shape = 21, size = 2.6)), shape = guide_legend(order = 1, nrow = 1, override.aes = list(fill = "grey85"))) +
  labs(title = expression(bold(paste(bolditalic("Hey1"), " pairs change sign between schemes"))), x = expression(bold(paste("Pooled Pearson ", bolditalic(r), " (12 observations)"))), y = NULL, tag = "(c)") + base +
  theme(axis.text.y = element_text(size = 8, colour = "black", family = "Arial"), legend.box = "vertical", legend.spacing.y = unit(0.05, "cm"))

fig <- ((pa1 | pa2) + plot_layout(widths = c(2.2, 1))) / pb / pc + plot_layout(heights = c(1.05, 0.95, 0.9))
png_path <- file.path(OUT, "Fig_S1_theme_npg.png")
ggsave(png_path, fig, device = ragg::agg_png, width = 17, height = 20, units = "cm", dpi = 600)
cat("wrote", png_path, "\nobserved pooled p<0.05:", obs_p, "| consensus", cons_p, "| adjusted p:", adj_p, "| adjusted q:", adj_q, "| cons adj p/q:", cons_ap, cons_aq, "\n")
