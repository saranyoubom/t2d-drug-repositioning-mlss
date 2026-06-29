# ==============================================================================
# Fig_S4_linregpcr_plot.R
# Supplementary Figure 4 — LinRegPCR amplification efficiency plots
#
# Reads the three CSV files produced by Fig_S4_linregpcr_extract.py and
# produces a 4 × 6 panel figure (21 genes) showing one representative
# amplification curve per treatment per gene.  The 5-cycle exponential-phase
# regression window is overlaid as a bold segment; per-condition E values
# (E = 10^slope) are annotated in treatment colour in the bottom-right corner
# of each panel.
#
# Gene panel order matches Supplementary Figure 3.
#
# Dependencies: ggplot2, dplyr (install.packages(c("ggplot2", "dplyr")))
# Usage:
#   1. Run Fig_S4_linregpcr_extract.py to generate the three CSV files.
#   2. Set DATA_DIR below to the directory containing those CSVs.
#   3. Set OUT_PNG to the desired output file path.
#   4. source("Fig_S4_linregpcr_plot.R")
# ==============================================================================

rm(list = ls())
library(ggplot2)
library(dplyr)

# ── User-configurable paths ───────────────────────────────────────────────────
DATA_DIR <- "."           # directory containing the three linregpcr_*.csv files
OUT_PNG  <- "Supplementary_Figure_S4_LinRegPCR.png"
# ─────────────────────────────────────────────────────────────────────────────

curves   <- read.csv(file.path(DATA_DIR, "linregpcr_curves.csv"),   stringsAsFactors = FALSE)
windows  <- read.csv(file.path(DATA_DIR, "linregpcr_windows.csv"),  stringsAsFactors = FALSE)
repannot <- read.csv(file.path(DATA_DIR, "linregpcr_repannot.csv"), stringsAsFactors = FALSE)

# ── Factor levels — order matches Supplementary Figure 3 ─────────────────────
gene_order <- c(
  "Glp1r",   "Glut2",   "Ins1",    "Hey1",
  "Hes1",    "Wnt2",    "Wnt2b",   "Wnt5a",
  "Wnt5b",   "Wnt9a",   "Lef1",    "Tcf7",
  "Tcf7l2",  "Itpr1",   "Ryr2",    "Ryr3",
  "Cacna1c", "Cacna1d", "Kcnj11",  "Ptbp1",
  "Gapdh"
)
cond_order  <- c("CTRL", "DAPT", "DKK-1")
cond_colors <- c(CTRL = "#2166ac", DAPT = "#d6604d", `DKK-1` = "#4dac26")

curves$gene      <- factor(curves$gene,      levels = gene_order)
windows$gene     <- factor(windows$gene,     levels = gene_order)
repannot$gene    <- factor(repannot$gene,    levels = gene_order)
curves$condition   <- factor(curves$condition,   levels = cond_order)
windows$condition  <- factor(windows$condition,  levels = cond_order)
repannot$condition <- factor(repannot$condition, levels = cond_order)

# ── Annotation labels (E value per representative well, per condition) ────────
# Stack three condition labels in bottom-right; vjust offsets upward from -Inf.
repannot <- repannot %>%
  mutate(
    label   = paste0(condition, "  E = ", formatC(E_rep, digits = 3, format = "f")),
    vjust_n = as.integer(condition)   # CTRL=1, DAPT=2, DKK-1=3
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  geom_line(
    data    = curves,
    mapping = aes(x = cycle, y = log10_rfu, colour = condition),
    linewidth = 0.5, alpha = 0.8
  ) +
  geom_line(
    data    = windows,
    mapping = aes(x = cycle, y = fit_log10, colour = condition),
    linewidth = 1.7, alpha = 1.0
  ) +
  geom_text(
    data    = repannot,
    mapping = aes(label = label, colour = condition, vjust = -vjust_n * 1.55),
    x = 39.5, y = -Inf, hjust = 1.0,
    size = 2.05, fontface = "plain", inherit.aes = FALSE
  ) +
  scale_colour_manual(
    name   = "Treatment",
    values = cond_colors,
    breaks = cond_order
  ) +
  facet_wrap(
    vars(gene),
    ncol   = 4,
    nrow   = 6,
    scales = "free_y"
  ) +
  labs(
    x       = "Cycle",
    y       = expression(log[10](RFU)),
    title   = "Supplementary Figure 4. Representative LinRegPCR amplification efficiency plots",
    caption = paste0(
      "Bold segment = 5-cycle exponential-phase window identified by LinRegPCR (R² ≥ 0.99, E ≥ 1.10); ",
      "representative well selected by highest R² in the exponential phase.\n",
      "E = 10^slope; 100% efficiency corresponds to E = 2.000."
    )
  ) +
  theme_grey() +
  theme(
    plot.title       = element_text(face = "bold", size = 9, hjust = 0),
    plot.caption     = element_text(size = 6.5, hjust = 0, colour = "grey30", lineheight = 1.3),
    strip.text       = element_text(face = "italic", size = 9, colour = "black"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    axis.text        = element_text(size = 7, colour = "black"),
    axis.title       = element_text(size = 8),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.8, "lines"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(6, 6, 4, 4)
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))

ggsave(
  filename = OUT_PNG,
  plot     = p,
  width    = 22,
  height   = 26,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)
cat("Saved:", OUT_PNG, "\n")
