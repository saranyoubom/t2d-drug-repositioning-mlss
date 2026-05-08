# Regulatory Approval & Mechanistic Pathway Sankey - Figure 6f

rm(list = ls())
set.seed(42)

library(tidyverse)
library(ggalluvial)
library(scales)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Load data
mlss_top15 <- read.csv("Fig_6c_mlss_top15.csv", stringsAsFactors = FALSE)
gene_centrality <- read.csv("Fig_6b_gene_centrality_summary.csv", stringsAsFactors = FALSE)

# Load regulatory databases
fda_approved <- read.csv("FDA_Approved.csv", stringsAsFactors = FALSE)
ema_approved <- read.csv("EMA_Approved.csv", stringsAsFactors = FALSE)
pmda_approved <- read.csv("PMDA_Approved.csv", stringsAsFactors = FALSE)
triple_approved <- read.csv("FDA-EMA-PMDA_Approved.csv", stringsAsFactors = FALSE)

# Extract hub drugs from top 15
top_15 <- mlss_top15 %>% arrange(Rank) %>% slice(1:15)

drug_freq_top15 <- bind_rows(
  top_15 %>% select(Drug = Drug_A),
  top_15 %>% select(Drug = Drug_B)
) %>%
  count(Drug, name = "Frequency") %>%
  arrange(desc(Frequency))

hub_threshold <- 4
hub_drugs <- drug_freq_top15 %>%
  filter(Frequency >= hub_threshold) %>%
  pull(Drug)

# Standardize drug names for regulatory lookup
fda_drugs <- tolower(fda_approved[,2])
ema_drugs <- tolower(ema_approved[,2])
pmda_drugs <- tolower(pmda_approved[,2])
triple_drugs <- tolower(triple_approved[,2])

# Assign regulatory approval tiers
drug_approval <- drug_freq_top15 %>%
  filter(Drug %in% hub_drugs) %>%
  mutate(
    Drug_lower = tolower(Drug),
    FDA = Drug_lower %in% fda_drugs,
    EMA = Drug_lower %in% ema_drugs,
    PMDA = Drug_lower %in% pmda_drugs,
    Approval_Tier = case_when(
      Drug_lower %in% triple_drugs ~ "FDA-EMA-PMDA\n(Tier 1)",
      FDA & EMA ~ "FDA-EMA\n(Tier 2)",
      FDA ~ "FDA-Only\n(Tier 3)",
      TRUE ~ "Not Approved"
    ),
    Confidence = case_when(
      Approval_Tier == "FDA-EMA-PMDA\n(Tier 1)" ~ 3.0,
      Approval_Tier == "FDA-EMA\n(Tier 2)" ~ 2.0,
      Approval_Tier == "FDA-Only\n(Tier 3)" ~ 1.0,
      TRUE ~ 0.3
    )
  ) %>%
  select(Drug, Frequency, Approval_Tier, Confidence, FDA, EMA, PMDA)

# Map drugs to pathways
drug_pathway_map <- bind_rows(
  top_15 %>% select(Drug = Drug_A, Pathway = Pathway_A),
  top_15 %>% select(Drug = Drug_B, Pathway = Pathway_B)
) %>%
  filter(Drug %in% hub_drugs) %>%
  group_by(Drug) %>%
  summarise(Primary_Pathway = first(Pathway), .groups = "drop")

# Define action types
drug_action_map <- tribble(
  ~Drug, ~Action_Type,
  "isradipine", "Antagonist",
  "(S)-nitrendipine", "Antagonist",
  "nitrendipine", "Antagonist",
  "albiglutide", "Agonist",
  "liraglutide", "Agonist"
)

# Map pathways to hub genes
pathway_gene_map <- tribble(
  ~Primary_Pathway, ~Hub_Gene, ~Connection_Strength,
  "Calcium Signaling", "Cacna1d", 1.0,
  "Calcium Signaling", "Cacna1c", 0.9,
  "Incretin Signaling", "Glp1r", 1.0,
  "Incretin Signaling", "Tcf7", 0.85,
  "Incretin Signaling", "Ptbp1", 0.75
)

# Map genes/pathways to clinical outcomes
pathway_outcome_map <- tribble(
  ~Primary_Pathway, ~Clinical_Outcome,
  "Calcium Signaling", "Enhanced insulin\nsecretion",
  "Calcium Signaling", "β-cell\npreservation",
  "Incretin Signaling", "Enhanced insulin\nsecretion",
  "Incretin Signaling", "β-cell\npreservation",
  "Incretin Signaling", "Metabolic\nhomeostasis"
)

# Build 6-layer alluvial dataset
alluvial_data <- drug_approval %>%
  left_join(drug_pathway_map, by = "Drug") %>%
  left_join(drug_action_map, by = "Drug") %>%
  left_join(pathway_gene_map, by = "Primary_Pathway") %>%
  left_join(pathway_outcome_map, by = "Primary_Pathway") %>%
  mutate(Flow = Frequency * Connection_Strength * Confidence) %>%
  select(`Approval Status` = Approval_Tier, `Hub Drugs` = Drug,
         `Action Type` = Action_Type, `Hub Genes` = Hub_Gene,
         `Pathways` = Primary_Pathway, `Clinical Outcomes` = Clinical_Outcome,
         Flow) %>%
  group_by(`Approval Status`, `Hub Drugs`, `Action Type`, `Hub Genes`, `Pathways`, `Clinical Outcomes`) %>%
  summarise(Flow = sum(Flow), .groups = "drop") %>%
  filter(Flow > 0) %>%
  arrange(desc(Flow))

# NPG colors
approval_colors <- c(
  "FDA-EMA-PMDA\n(Tier 1)" = npg_colors[1],
  "FDA-EMA\n(Tier 2)" = npg_colors[3],
  "FDA-Only\n(Tier 3)" = npg_colors[5],
  "Not Approved" = "#A6A6A6"
)

pathway_colors <- c(
  "Calcium Signaling" = npg_colors[4],
  "Incretin Signaling" = npg_colors[1]
)

action_colors <- c(
  "Agonist" = npg_colors[2],
  "Antagonist" = npg_colors[4]
)

outcome_colors <- c(
  "Enhanced insulin\nsecretion" = npg_colors[2],
  "β-cell\npreservation" = npg_colors[6],
  "Metabolic\nhomeostasis" = npg_colors[3]
)

# Create alluvial diagram
p_sankey <- ggplot(
  alluvial_data,
  aes(axis1 = `Approval Status`, axis2 = `Hub Drugs`, axis3 = `Action Type`,
      axis4 = `Hub Genes`, axis5 = `Pathways`, axis6 = `Clinical Outcomes`, y = Flow)
) +
  geom_alluvium(aes(fill = `Pathways`), width = 3/8, alpha = 0.7, curve_type = "sigmoid") +
  geom_stratum(width = 3/8, fill = "white", color = "white", size = 0.6, alpha = 0.8) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(name = "Signaling Pathway", values = pathway_colors,
                    guide = guide_legend(override.aes = list(alpha = 0.9, size = 5), order = 1)) +
  scale_x_discrete(limits = c("Approval\nStatus", "Hub\nDrugs", "Action\nType",
                               "Hub\nGenes", "Pathways", "Clinical\nOutcomes"),
                   expand = c(0.08, 0.08)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "6-Layer Regulatory-Mechanistic Network for Type 2 Diabetes Drug Repositioning")

ggsave("Fig_6f_regulatory_mechanistic_sankey.png", p_sankey,
       width = 17, height = 4.47, units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_6f_regulatory_mechanistic_sankey.svg", p_sankey,
       width = 17, height = 4.47, units = "cm", scale = 2.3, dpi = 600, bg = "white")

# Export supplementary tables
write.csv(drug_approval, "Fig_6f_hub_drug_approval_profile.csv", row.names = FALSE)
write.csv(alluvial_data, "Fig_6f_six_layer_complete_mapping.csv", row.names = FALSE)
write.csv(pathway_gene_map, "Fig_6f_pathway_gene_connectivity.csv", row.names = FALSE)

