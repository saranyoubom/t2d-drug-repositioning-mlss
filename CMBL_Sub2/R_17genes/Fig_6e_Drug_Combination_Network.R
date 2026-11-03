###############################################################################
# FIGURE 6e – DRUG COMBINATION NETWORK (MLSS, NPG PALETTE)
#
# Data source:
#   - Fig_6c_MLSS_ALL_Combinations.csv   (MLSS output)
#
# Key design choices:
#   - Nodes = individual drugs
#   - Node color = primary pathway (Calcium / Incretin / Wnt / Metabolic / Other)
#   - Node size = number of appearances in Top N combinations (hub drugs bigger)
#   - Edges = drug pairs (Top N MLSS combinations)
#   - Edge width = MLSS_v4.0 score (stronger combinations = thicker)
#   - Edge linetype = Cross‑Pathway vs Same‑Pathway
#   - All edges colored grey70 (to keep focus on node categories)
#   - NPG palette for node fill/outline (ggsci::pal_npg)
###############################################################################

rm(list = ls())
set.seed(42)

#==============================================================================
# 0. SETUP & PACKAGES
#==============================================================================

# Adjust to your figure directory
setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")

required_packages <- c(
  "tidyverse",
  "igraph",
  "ggraph",
  "scales",
  "ggsci"
)

invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

# NPG palette
npg_colors <- pal_npg("nrc")(10)
cat("\n✓ NPG Color Palette loaded for Fig 6e\n")
cat(sprintf(" Colors: %s\n\n", paste(npg_colors[1:6], collapse = " | ")))

# Simple export helper (SVG + PNG, 600 DPI, SBA sizing)
export_figure <- function(filename_base, plot_object,
                          width_cm = 7.7, height_cm = 6) {
  ggsave(
    paste0(filename_base, ".svg"),
    plot_object,
    scale  = 2.3,
    width  = width_cm,
    height = height_cm,
    units  = "cm",
    dpi    = 600,
    bg     = "white"
  )
  ggsave(
    paste0(filename_base, ".png"),
    plot_object,
    scale  = 2.3,
    width  = width_cm,
    height = height_cm,
    units  = "cm",
    dpi    = 600,
    bg     = "white"
  )
  cat(sprintf("✓ Exported: %s (.svg & .png)\n", filename_base))
}

cat("████████████████████████████████████████████████████████████████████\n")
cat(" FIGURE 6e: DRUG COMBINATION NETWORK (MLSS, NPG PALETTE)\n")
cat("████████████████████████████████████████████████████████████████████\n\n")

#==============================================================================
# 1. LOAD MLSS DATA (ALL COMBINATIONS)
#==============================================================================

cat("STEP 1: Loading MLSS combination data...\n\n")

mlss_all <- read.csv("Fig_6c_MLSS_ALL_Combinations.csv",
                     stringsAsFactors = FALSE)

cat(sprintf("  ✓ Total combinations: %d\n", nrow(mlss_all)))
cat(sprintf("  ✓ Columns: %s\n\n",
            paste(colnames(mlss_all)[1:10], collapse = ", ")))

# Sanity: ensure required columns exist
required_cols <- c(
  "Rank", "Combination", "Drug_A", "Drug_B",
  "Pathway_A", "Pathway_B",
  "Cross_Pathway", "MLSS_v4.0"
)
missing_cols <- setdiff(required_cols, colnames(mlss_all))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing columns in MLSS file: %s",
               paste(missing_cols, collapse = ", ")))
}

#==============================================================================
# 2. SELECT TOP N COMBINATIONS FOR NETWORK
#==============================================================================

top_n <- 15  # change to 30/50 if needed

cat(sprintf("STEP 2: Selecting Top %d MLSS combinations...\n\n", top_n))

top_combos <- mlss_all %>%
  arrange(Rank) %>%
  slice(1:min(top_n, nrow(.))) %>%
  mutate(
    Rank = row_number(),
    Is_CrossPathway = ifelse(Cross_Pathway == "Yes",
                             "Cross-Pathway", "Same-Pathway")
  )

cat(sprintf("  ✓ Top %d combinations extracted\n", nrow(top_combos)))
cat(sprintf("  ✓ MLSS range (Top %d): %.2f – %.2f\n\n",
            nrow(top_combos),
            min(top_combos$MLSS_v4.0, na.rm = TRUE),
            max(top_combos$MLSS_v4.0, na.rm = TRUE)))

#==============================================================================
# 3. BUILD EDGE LIST AND NODE METADATA
#==============================================================================

cat("STEP 3: Building network edges and nodes...\n\n")

# Edge list for igraph
edge_list <- top_combos %>%
  transmute(
    from         = Drug_A,
    to           = Drug_B,
    weight       = MLSS_v4.0,
    relationship = Is_CrossPathway,
    path_A       = Pathway_A,
    path_B       = Pathway_B,
    combo_label  = Combination,
    rank         = Rank
  )

# All unique drugs in Top N
all_drugs <- sort(unique(c(top_combos$Drug_A, top_combos$Drug_B)))

# Pathway assignment for each drug:
#   if a drug appears in multiple combos with conflicting pathways,
#   choose the most frequent pathway; if tie → "Multi-pathway".
drug_pathways <- bind_rows(
  top_combos %>%
    select(Drug = Drug_A, Pathway = Pathway_A),
  top_combos %>%
    select(Drug = Drug_B, Pathway = Pathway_B)
) %>%
  group_by(Drug, Pathway) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Drug) %>%
  arrange(Drug, desc(n)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  group_by(Drug) %>%
  summarise(
    Primary_Pathway = if (n_distinct(Pathway) == 1) {
      first(Pathway)
    } else if (n_distinct(Pathway) > 1 && first(n) > nth(n, 2)) {
      # strict majority
      first(Pathway)
    } else {
      "Multi-pathway"
    },
    .groups = "drop"
  )

# Node table
node_list <- tibble(name = all_drugs) %>%
  left_join(drug_pathways,
            by = c("name" = "Drug")) %>%
  mutate(
    Primary_Pathway = ifelse(is.na(Primary_Pathway),
                             "Unknown", Primary_Pathway),
    appearance_count = sapply(name, function(d) {
      sum(top_combos$Drug_A == d | top_combos$Drug_B == d)
    })
  )

cat(sprintf("  ✓ Nodes (unique drugs): %d\n", nrow(node_list)))
cat(sprintf("  ✓ Average appearances per drug: %.2f\n\n",
            mean(node_list$appearance_count)))

#==============================================================================
# 4. DEFINE NPG COLOR MAPPING FOR PATHWAYS
#==============================================================================

cat("STEP 4: Defining NPG color mapping for pathways...\n\n")

# Align with Figure 6b categorical mapping as much as possible
pathway_colors_npg <- c(
  "Calcium Signaling"     = npg_colors[4], # dark blue
  "Incretin Signaling"    = npg_colors[1], # red/orange
  "Wnt Signaling"         = npg_colors[2], # green
  "Metabolic Regulation"  = npg_colors[3], # light blue
  "Multi-pathway"         = npg_colors[6], # purple
  "Unknown"               = "grey70"
)

cat("NPG Pathway Colors (Fig 6e):\n")
for (nm in names(pathway_colors_npg)) {
  cat(sprintf("  %-20s : %s\n", nm, pathway_colors_npg[[nm]]))
}
cat("\n")

#==============================================================================
# 5. CREATE iGRAPH OBJECT
#==============================================================================

cat("STEP 5: Creating igraph network object...\n\n")

g_drug <- graph_from_data_frame(
  d = edge_list,
  vertices = node_list,
  directed = FALSE
)

# Node attributes
V(g_drug)$appearance_count <- node_list$appearance_count
V(g_drug)$Primary_Pathway  <- node_list$Primary_Pathway

# Edge attributes
E(g_drug)$weight       <- edge_list$weight
E(g_drug)$relationship <- edge_list$relationship

hub_threshold <- 4
V(g_drug)$is_hub <- ifelse(V(g_drug)$appearance_count >= hub_threshold,
                           "Hub", "Non-hub")

cat(sprintf("  ✓ Network created: %d nodes, %d edges\n",
            vcount(g_drug), ecount(g_drug)))
cat(sprintf("  ✓ Hub threshold: ≥%d appearances\n\n", hub_threshold))

#==============================================================================
# 6. VISUALIZATION – NPG-PALETTE NETWORK
#==============================================================================

cat("STEP 6: Rendering NPG-style drug combination network...\n\n")

mlss_min <- min(E(g_drug)$weight, na.rm = TRUE)
mlss_max <- max(E(g_drug)$weight, na.rm = TRUE)

p_fig6e <- ggraph(g_drug, layout = "circle") +
  
  # Edges: width = MLSS, linetype = cross vs same pathway, color = grey70
  geom_edge_link(
    aes(
      width    = weight,
      linetype = relationship
    ),
    color = "grey70",
    alpha = 0.6,
    show.legend = TRUE
  ) +
  scale_edge_linetype_manual(
    name   = "Combination Type",
    values = c("Cross-Pathway" = "solid",
               "Same-Pathway"  = "dotted"),
    guide  = guide_legend(
      order = 3,
      override.aes = list(color = "black", width = 0.8)
    )
  ) +
  scale_edge_width_continuous(
    name   = "MLSS Score",
    range  = c(0.5, 3.5),
    breaks = c(mlss_min,
               (mlss_min + mlss_max) / 2,
               mlss_max),
    labels = c("Lower", "Intermediate", "Higher"),
    guide  = guide_legend(order = 4,
                          override.aes = list(color = "black"))
  ) +
  
  # Nodes: NPG colors by pathway, size by appearance_count
  geom_node_point(
    aes(
      size  = appearance_count,
      color = Primary_Pathway,
      fill  = Primary_Pathway
    ),
    shape  = 21,
    stroke = 1.2,
    alpha  = 0.95
  ) +
  scale_color_manual(
    name   = "Functional Category",
    values = pathway_colors_npg,
    guide  = guide_legend(order = 1)
  ) +
  scale_fill_manual(
    name   = "Functional Category",
    values = pathway_colors_npg,
    guide  = "none"
  ) +
  scale_size_continuous(
    name   = "Network Appearances",
    range  = c(5, 14),
    breaks = c(1, hub_threshold, max(node_list$appearance_count)),
    labels = c("1", paste0(hub_threshold, " (hub)"),
               as.character(max(node_list$appearance_count))),
    guide  = guide_legend(order = 2)
  ) +
  
  # Labels: bold for hub drugs
  geom_node_text(
    aes(
      label    = name,
      fontface = ifelse(appearance_count >= hub_threshold,
                        "bold", "plain")
    ),
    repel        = TRUE,
    size         = 4.5,
    max.overlaps = 30
  ) +
  
  theme_void(base_size = 12) +
  theme(
    plot.title      = element_text(
      face  = "bold",
      size  = 16,
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(0.45, "cm")
  ) +
  labs(
    title = sprintf("Drug Combination Network (Top %d)",
                    nrow(top_combos))
  )

print(p_fig6e)

# Export
export_figure(
  filename_base = "Fig_6e_Drug_Combination_Network",
  plot_object   = p_fig6e,
  width_cm      = 7.7,
  height_cm     = 6
)

#==============================================================================
# 7. EXPORT SIMPLE NODE/EDGE STATISTICS
#==============================================================================

node_stats <- node_list %>%
  mutate(
    Hub_Status = ifelse(appearance_count >= hub_threshold,
                        "Hub", "Non-hub")
  ) %>%
  arrange(desc(appearance_count))

edge_stats <- edge_list %>%
  arrange(desc(weight))

write.csv(node_stats,
          "Fig_6e_Drug_Combination_Network_NodeStats.csv",
          row.names = FALSE)
write.csv(edge_stats,
          "Fig_6e_Drug_Combination_Network_EdgeStats.csv",
          row.names = FALSE)

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat(" ✓ FIGURE 6e COMPLETE – NPG PALETTE, MLSS-BASED DRUG NETWORK\n")
cat("   Output:\n")
cat("    • Fig_6e_Drug_Combination_Network_NPG.svg / .png\n")
cat("    • Fig_6e_Drug_Combination_Network_NodeStats.csv\n")
cat("    • Fig_6e_Drug_Combination_Network_EdgeStats.csv\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# End of script
