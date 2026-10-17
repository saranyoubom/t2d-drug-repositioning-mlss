# Supplementary Fig. S3 (melt curves of the 21 RT-qPCR assays) redrawn from the raw CFX Maestro 'Melt Curve Derivative Results'
# exports in the theme of the main figures: ggplot2 theme_minimal(), ggsci NPG palette, Arial, black text, treatment colours as in
# Figs 1-3 (CTRL = NPG red, DAPT = NPG cyan, DKK-1 = NPG green). Every well of the exported melt data is drawn.
# Inputs (made by extract_melt_curves.py): input_data/melt_curves_long.csv, input_data/melt_peaks.csv
# Outputs: Fig_S3_melt_curves_npg.png, melt_qc_per_well.csv (main peak, plateau, secondary low-temperature signal for every well)
# Usage: Rscript R_SuppFig3_melt_curves_npg.R [analysis_dir] [out_dir]
suppressPackageStartupMessages({ library(ggplot2); library(ggsci) })
args <- commandArgs(trailingOnly = TRUE)
A   <- if (length(args) >= 1) args[1] else "."
OUT <- if (length(args) >= 2) args[2] else "."
npg <- pal_npg("nrc")(10)
cols <- c("CTRL" = npg[1], "DAPT" = npg[2], "DKK-1" = npg[3])

cv <- read.csv(file.path(A, "input_data", "melt_curves_long.csv"), stringsAsFactors = FALSE)
pk <- read.csv(file.path(A, "input_data", "melt_peaks.csv"), stringsAsFactors = FALSE)
cv$id <- paste(cv$run, cv$well_label, cv$gene)

# ---- per-well quality control on the curves (all genes incl. Ryr1) --------------------------------------------
qc <- do.call(rbind, lapply(split(cv, cv$id), function(d) {
  d <- d[order(d$temperature), ]
  hi <- d$temperature >= 78
  main <- max(d$dRFU_dT[hi]); tm_peak <- d$temperature[hi][which.max(d$dRFU_dT[hi])]
  plateau <- median(d$dRFU_dT[d$temperature >= 71 & d$temperature <= 74])
  seg <- d$temperature >= 72 & d$temperature <= tm_peak - 4
  sec <- if (any(seg)) max(d$dRFU_dT[seg]) else NA
  data.frame(run = d$run[1], well = d$well_label[1], gene = d$gene[1], treatment = d$treatment[1], main_peak = main, Tm_curve = tm_peak,
             plateau = plateau, secondary_max = sec, secondary_above_plateau = sec - plateau,
             secondary_fraction = (sec - plateau) / main, no_main_peak = main < 20)
}))
qc$secondary_signal <- !qc$no_main_peak & qc$secondary_fraction >= 0.20      # secondary signal >= 20% of the main peak height
write.csv(qc, file.path(OUT, "melt_qc_per_well.csv"), row.names = FALSE)

order_genes <- c("Glp1r", "Glut2", "Ins1", "Hey1", "Hes1", "Wnt2", "Wnt2b", "Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2",
                 "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Gapdh")
cv <- cv[cv$gene %in% order_genes & cv$temperature >= 65 & cv$temperature <= 95, ]
pk <- pk[pk$gene %in% order_genes, ]
cv$gene <- factor(cv$gene, levels = order_genes)
cv$treatment <- factor(cv$treatment, levels = names(cols))

fmt <- function(x) formatC(x, format = "f", digits = 1)
ann <- do.call(rbind, lapply(order_genes, function(g) {
  d <- pk[pk$gene == g & !is.na(pk$Tm), ]; rng <- range(d$Tm); q <- qc[qc$gene == g, ]
  lab <- c(paste0("Tm ", if (rng[1] == rng[2]) fmt(rng[1]) else paste0(fmt(rng[1]), "–", fmt(rng[2])), " °C"), paste0("n = ", nrow(q), " wells"))
  if (any(q$secondary_signal)) lab <- c(lab, sprintf("secondary signal: %d", sum(q$secondary_signal)))
  if (any(q$no_main_peak)) lab <- c(lab, sprintf("no main peak: %d", sum(q$no_main_peak)))
  data.frame(gene = factor(g, levels = order_genes), label = paste(lab, collapse = "
"))
}))
ymax <- tapply(cv$dRFU_dT, cv$gene, max)[as.character(ann$gene)]; ymin <- tapply(cv$dRFU_dT, cv$gene, min)[as.character(ann$gene)]
nl <- lengths(strsplit(ann$label, "
"))
ann$y <- as.numeric(ymax + (0.10 + 0.115 * nl) * (ymax - ymin))      # top of the text band above the curves; geom_blank keeps it inside the panel

p <- ggplot(cv, aes(temperature, dRFU_dT, group = id, colour = treatment)) +
  geom_line(linewidth = 0.25, alpha = 0.75) +
  geom_blank(data = ann, aes(x = 65.5, y = y), inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = 65.5, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.2, lineheight = 0.95, colour = "black") +
  facet_wrap(~gene, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = cols, name = "Treatment") +
  scale_x_continuous(limits = c(65, 95), breaks = seq(65, 95, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.04))) +
  labs(x = "Temperature (°C)", y = "−d(RFU)/dT") +
  theme_minimal() +
  theme(legend.position = "top",
        text = element_text(family = "Arial"),
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic", family = "Arial"),
        axis.text.x = element_text(size = 7, colour = "black", family = "Arial"),
        axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
        axis.title = element_text(size = 8, face = "bold", colour = "black", family = "Arial"),
        legend.text = element_text(size = 8, colour = "black", family = "Arial"),
        legend.title = element_text(size = 8, face = "bold", colour = "black", family = "Arial"),
        panel.spacing = unit(0.35, "cm"),
        plot.margin = margin(4, 6, 4, 4)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1)))

png_path <- file.path(OUT, "Fig_S3_melt_curves_npg.png")
ggsave(png_path, p, device = ragg::agg_png, width = 17, height = 21, units = "cm", dpi = 600)
cat("wrote", png_path, "\n")
cat("\nQC summary (secondary signal >= 20% of main peak above the 71-74 C plateau; no main peak = max above 78 C < 20 RFU/dT)\n")
s <- aggregate(cbind(secondary_signal, no_main_peak) ~ gene, qc, sum); s$wells <- as.numeric(table(qc$gene)[s$gene])
print(s[order(-s$secondary_signal - s$no_main_peak), ], row.names = FALSE)
print(qc[qc$secondary_signal | (qc$no_main_peak & qc$gene != "Ryr1"), c("run", "well", "gene", "treatment", "main_peak", "secondary_above_plateau", "secondary_fraction", "no_main_peak")], row.names = FALSE)
