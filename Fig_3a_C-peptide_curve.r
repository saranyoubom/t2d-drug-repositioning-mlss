rm(list = ls())

library(ggplot2)
library(readr)
library(ggsci)
library(readxl)
library(dplyr)

# ── LOAD DATA ────────────────────────────────────────────────────────────────

dat2 <- read_excel("Fig_3a_Scatterplot_C_peptide.xlsx", sheet = 1)

dat2$Glucose <- as.numeric(dat2$Glucose)
dat2 <- filter(dat2, Glucose != 44)

# ── EXPORT LOESS SMOOTH + CI PER TREATMENT ───────────────────────────────────

# Generate predicted values across the glucose range for each Treatment
glucose_seq <- seq(min(dat2$Glucose), max(dat2$Glucose), length.out = 200)

loess_export <- lapply(unique(dat2$Treatment), function(trt) {
  sub   <- dat2 %>% filter(Treatment == trt)
  model <- loess(C_peptide ~ Glucose, data = sub, span = 0.75)  # default ggplot span
  pred  <- predict(model, newdata = data.frame(Glucose = glucose_seq), se = TRUE)
  
  data.frame(
    Treatment  = trt,
    Assay      = unique(sub$Assay)[1],
    Glucose    = glucose_seq,
    Fitted     = pred$fit,
    SE         = pred$se.fit,
    CI_lower   = pred$fit - 1.96 * pred$se.fit,
    CI_upper   = pred$fit + 1.96 * pred$se.fit
  )
})

loess_df <- bind_rows(loess_export)

# Clip CI to plot limits (y >= 0) to match the plot
loess_df <- loess_df %>%
  mutate(CI_lower = pmax(CI_lower, 0),
         CI_upper = pmin(CI_upper, 320))

write.csv(loess_df, "Fig_3a_LOESS_smooth_CI.csv", row.names = FALSE)

# ── PLOT ─────────────────────────────────────────────────────────────────────

plot_loess <- ggplot(dat2, aes(x = Glucose, y = C_peptide, color = Treatment)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  geom_point(size = 1.5, alpha = 0.5) +
  facet_grid(cols = vars(Assay)) +
  scale_color_npg() +
  labs(x = "Glucose concentration (mM)",
       y = "\nRelative C-peptide release (AU)") +
  lims(y = c(0, 320)) +
  theme_minimal() +
  theme(legend.position      = "top",
        legend.key.width     = unit(0.5, "cm"),
        legend.key.height    = unit(0.5, "cm"),
        strip.text           = element_text(size = 8, colour = "black", face = "bold"),
        legend.margin        = margin(0, 0, 0, 0),
        legend.box.spacing   = margin(5),
        legend.text          = element_text(size = 8),
        legend.title         = element_text(size = 8),
        axis.title.x         = element_text(size = 8),
        axis.text.x          = element_text(size = 8, colour = "black",
                                            angle = 0, vjust = 0.5),
        axis.text.y          = element_text(size = 8, colour = "black",
                                            angle = 90, hjust = 0.5),
        axis.title.y         = element_text(size = 8))

plot_loess

# ── SAVE ─────────────────────────────────────────────────────────────────────

ggsave("Fig_3a_Scatterplot_C_peptide.png", plot = plot_loess,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3a_Scatterplot_C_peptide.svg", plot = plot_loess,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

