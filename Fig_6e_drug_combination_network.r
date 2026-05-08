# Drug Combination Network - Figure 6e

rm(list = ls())
set.seed(42)

library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Load MLSS data ─────────────────────────────────────────────────────────────
mlss_all <- read.csv("Fig_6c_mlss_all_combinations.csv", stringsAsFactors = FALSE)

# Normalise column names (trim whitespace, handle R mangling)
colnames(mlss_all) <- trimws(colnames(mlss_all))

# ── Select top N combinations ──────────────────────────────────────────────────
top_n <- 15

# Auto-detect cross-pathway column or derive from Pathway_A vs Pathway_B
if ("Cross_Pathway" %in% colnames(mlss_all)) {
  cross_col <- "Cross_Pathway"
} else if ("Cross.Pathway" %in% colnames(mlss_all)) {
  cross_col <- "Cross.Pathway"
} else {
  cross_col <- NULL
}

top_combos <- mlss_all %>%
  arrange(Rank) %>%
  slice(1:min(top_n, nrow(.))) %>%
  mutate(
    Rank = row_number(),
    Is_CrossPathway = if (!is.null(cross_col)) {
      ifelse(.data[[cross_col]] == "Yes", "Cross-Pathway", "Same-Pathway")
    } else {
      ifelse(Pathway_A != Pathway_B, "Cross-Pathway", "Same-Pathway")
    }
  )

# ── Build edge list ────────────────────────────────────────────────────────────
edge_list <- top_combos %>%
  transmute(
    from         = Drug_A,
    to           = Drug_B,
    weight       = MLSS_v4.0,
    relationship = Is_CrossPathway,
    rank         = Rank
  )

# ── All unique drugs ───────────────────────────────────────────────────────────
all_drugs <- sort(unique(c(top_combos$Drug_A, top_combos$Drug_B)))

# ── Assign primary pathway to each drug ───────────────────────────────────────
drug_pathways <- bind_rows(
  top_combos %>% select(Drug = Drug_A, Pathway = Pathway_A),
  top_combos %>% select(Drug = Drug_B, Pathway = Pathway_B)
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
      first(Pathway)
    } else {
      "Multi-pathway"
    },
    .groups = "drop"
  )

# ── Node list ──────────────────────────────────────────────────────────────────
node_list <- tibble(name = all_drugs) %>%
  left_join(drug_pathways, by = c("name" = "Drug")) %>%
  mutate(
    Primary_Pathway  = ifelse(is.na(Primary_Pathway), "Unknown", Primary_Pathway),
    appearance_count = sapply(name, function(d) {
      sum(top_combos$Drug_A == d | top_combos$Drug_B == d)
    })
  )

# ── NPG colors for pathways ────────────────────────────────────────────────────
pathway_colors_npg <- c(
  "Calcium Signaling"   = npg_colors[4],
  "Incretin Signaling"  = npg_colors[1],
  "Wnt Signaling"       = npg_colors[2],
  "Metabolic Regulation"= npg_colors[3],
  "Multi-pathway"       = npg_colors[6],
  "Unknown"             = "grey70"
)

# ── Create igraph object ───────────────────────────────────────────────────────
g_drug <- graph_from_data_frame(
  d        = edge_list,
  vertices = node_list,
  directed = FALSE
)

V(g_drug)$appearance_count <- node_list$appearance_count
V(g_drug)$Primary_Pathway  <- node_list$Primary_Pathway
E(g_drug)$weight           <- edge_list$weight
E(g_drug)$relationship     <- edge_list$relationship

hub_threshold       <- 4
V(g_drug)$is_hub    <- ifelse(V(g_drug)$appearance_count >= hub_threshold, "Hub", "Non-hub")

# ── MLSS edge width range ──────────────────────────────────────────────────────
mlss_min <- min(E(g_drug)$weight, na.rm = TRUE)
mlss_max <- max(E(g_drug)$weight, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT — Drug Combination Network  →  Fig_6e
# ══════════════════════════════════════════════════════════════════════════════
p_Fig_6e <- ggraph(g_drug, layout = "circle") +
  geom_edge_link(
    aes(width = weight, linetype = relationship),
    color       = "grey70",
    alpha       = 0.6,
    show.legend = TRUE
  ) +
  scale_edge_linetype_manual(
    name   = "Combination Type",
    values = c("Cross-Pathway" = "solid", "Same-Pathway" = "dotted"),
    guide  = guide_legend(order = 3,
                          override.aes = list(color = "black", width = 0.8))
  ) +
  scale_edge_width_continuous(
    name   = "MLSS Score",
    range  = c(0.5, 3.5),
    breaks = c(mlss_min, (mlss_min + mlss_max) / 2, mlss_max),
    labels = c("Lower", "Intermediate", "Higher"),
    guide  = guide_legend(order = 4, override.aes = list(color = "black"))
  ) +
  geom_node_point(
    aes(size = appearance_count, color = Primary_Pathway, fill = Primary_Pathway),
    shape  = 21,
    stroke = 1.2,
    alpha  = 0.95
  ) +
  scale_color_manual(
    name  = "Functional Category",
    values = pathway_colors_npg,
    guide  = guide_legend(order = 1)
  ) +
  scale_fill_manual(
    name  = "Functional Category",
    values = pathway_colors_npg,
    guide  = "none"
  ) +
  scale_size_continuous(
    name   = "Network Appearances",
    range  = c(5, 14),
    breaks = c(1, hub_threshold, max(node_list$appearance_count)),
    labels = c("1",
               paste0(hub_threshold, " (hub)"),
               as.character(max(node_list$appearance_count))),
    guide  = guide_legend(order = 2)
  ) +
  geom_node_text(
    aes(label    = name,
        fontface = ifelse(appearance_count >= hub_threshold, "bold", "plain")),
    repel        = TRUE,
    size         = 4.5,
    max.overlaps = 30
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.45, "cm")
  ) +
  labs(title = sprintf("Drug Combination Network (Top %d)", nrow(top_combos)))

ggsave("Fig_6e_drug_combination_network.png", p_Fig_6e,
       scale = 2.3, width = 7.7, height = 6,
       units = "cm", dpi = 600, bg = "white")
ggsave("Fig_6e_drug_combination_network.svg", p_Fig_6e,
       scale = 2.3, width = 7.7, height = 6,
       units = "cm", dpi = 600, bg = "white")

# ── Export statistics ──────────────────────────────────────────────────────────
node_stats <- node_list %>%
  mutate(Hub_Status = ifelse(appearance_count >= hub_threshold, "Hub", "Non-hub")) %>%
  arrange(desc(appearance_count))

# Explicit namespace prevents tidyverse masking igraph's as_data_frame()
edge_stats <- igraph::as_data_frame(g_drug, what = "edges") %>%
  arrange(desc(weight)) %>%
  rename(
    Gene1              = from,
    Gene2              = to,
    `MLSS Score`       = weight,
    `Combination Type` = relationship
  )

write.csv(node_stats, "Fig_6e_drug_node_statistics.csv", row.names = FALSE)
write.csv(edge_stats, "Fig_6e_drug_edge_statistics.csv", row.names = FALSE)

