# =============================================================================
# Four-way Venn Diagram - Figure 4d
# Revised: Hmisc replaced with base R cor.test(); Ryr1 excluded (QC failure)
# =============================================================================

rm(list = ls())

library(VennDiagram)
library(grid)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Load normalized data (already Ryr1-excluded from Fig 4c script) ───────────
hes1_norm   <- read.csv("Fig_4c_Hes1_normalized_data.csv")
kcnj11_norm <- read.csv("Fig_4c_Kcnj11_normalized_data.csv")
ins1_norm   <- read.csv("Fig_4c_Ins1_normalized_data.csv")
raw_data    <- read.csv("Fig_4c_Raw_expression_data.csv")

# ── Safety check: drop Ryr1 if still present ──────────────────────────────────
drop_ryr1 <- function(df) df[, colnames(df) != "Ryr1"]
hes1_norm   <- drop_ryr1(hes1_norm)
kcnj11_norm <- drop_ryr1(kcnj11_norm)
ins1_norm   <- drop_ryr1(ins1_norm)
raw_data    <- drop_ryr1(raw_data)

# ── Function to extract significant pairs (base R, no Hmisc) ─────────────────
get_pairs <- function(data) {
  data_clean <- as.matrix(na.omit(data))
  n_genes    <- ncol(data_clean)
  genes      <- colnames(data_clean)
  pairs      <- character()
  
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      test  <- cor.test(data_clean[, i], data_clean[, j], method = "pearson")
      p_val <- test$p.value
      if (!is.na(p_val) && p_val < 0.05) {
        pairs <- c(pairs, paste(genes[i], genes[j], sep = "_"))
      }
    }
  }
  return(pairs)
}

# ── Extract pairs from each method ────────────────────────────────────────────
pairs_list <- list(
  "Hes1-\nnormalized"   = get_pairs(hes1_norm),
  "Kcnj11-\nnormalized" = get_pairs(kcnj11_norm),
  "Ins1-\nnormalized"   = get_pairs(ins1_norm),
  "Raw\nexpression"     = get_pairs(raw_data)
)

# ── Identify consensus pairs (present in all 4 methods) ───────────────────────
consensus_pairs <- Reduce(intersect, pairs_list)

# ── Export pair lists and consensus ───────────────────────────────────────────
write.csv(data.frame(Pairs = pairs_list[[1]]),
          "Fig_4d_Hes1_normalized_pairs.csv",   row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[2]]),
          "Fig_4d_Kcnj11_normalized_pairs.csv", row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[3]]),
          "Fig_4d_Ins1_normalized_pairs.csv",   row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[4]]),
          "Fig_4d_Raw_expression_pairs.csv",    row.names = FALSE)
write.csv(data.frame(Consensus_Pairs = consensus_pairs),
          "Fig_4d_consensus_pairs.csv",         row.names = FALSE)

# ── Build Venn diagram object ─────────────────────────────────────────────────
venn.plot <- venn.diagram(
  x              = pairs_list,
  category.names = names(pairs_list),
  filename       = NULL,
  fill           = c(npg_colors[4], npg_colors[2], npg_colors[3], npg_colors[7]),
  col            = c(npg_colors[4], npg_colors[2], npg_colors[3], npg_colors[7]),
  lwd            = 1.5,
  alpha          = 0.5,
  cex            = 0.8,
  fontface       = "bold",
  fontfamily     = "Arial",
  cat.cex        = 0.7,
  cat.fontface   = "bold",
  cat.fontfamily = "Arial",
  cat.default.pos = "outer",
  cat.dist       = c(0.08, 0.08, 0.08, 0.08)
)

# ── Save PNG ──────────────────────────────────────────────────────────────────
png("Fig_4d_Venn_diagram.png",
    width = 9.16, height = 9.16, units = "cm", res = 600)
grid.draw(venn.plot)
dev.off()

# ── Save SVG ──────────────────────────────────────────────────────────────────
svg("Fig_4d_Venn_diagram.svg",
    width = 9.16 / 2.54, height = 9.16 / 2.54)
grid.draw(venn.plot)
dev.off()

