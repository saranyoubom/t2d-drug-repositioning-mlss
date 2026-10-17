# Sensitivity of the 2-NBDG (glucose uptake) result to the undocumented analysis choices in the source workbook
# 2021-04-07_Gluclose uptake assay(.xlsx / _Original.xlsx): CTRL replicate 2 (184 raw vs 331 in the working blocks),
# exclusion of the high CTRL well (3382), background (negative-control wells), read (1st/2nd) and cell density (25K/50K).
# Input : input_data/nbdg_blocks.csv (cell addresses included)
# Output: 2NBDG_sensitivity_table.csv
# Tests (as in the published script): Kruskal-Wallis + pairwise Wilcoxon rank-sum, Bonferroni. Also one-way ANOVA with
# Bonferroni pooled-SD t-tests on log10 intensity, because a rank test with n = 3-4 per group cannot reach Bonferroni
# p < 0.05 and a log-scale parametric test is the usual alternative for fluorescence data.
args <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(args) >= 1) args[1] else "input_data/nbdg_blocks.csv"
outp <- if (length(args) >= 2) args[2] else "2NBDG_sensitivity_table.csv"
raw  <- read.csv(inp, stringsAsFactors = FALSE, colClasses = "character")
pubv <- as.numeric(strsplit(raw$CTRL[raw$block == "published_CTRL_vector"], ",")[[1]])
b <- raw[raw$block != "published_CTRL_vector", ]
for (k in c("rep", "blank", "CTRL", "DKK1", "DAPT")) b[[k]] <- as.numeric(b[[k]])
blk <- function(name) b[b$block == name, ]

sem <- function(x) sd(x) / sqrt(length(x))
run <- function(label, ctrl, dapt, dkk, note = "") {
  v <- c(ctrl, dapt, dkk)
  g <- factor(rep(c("CTRL", "DAPT", "DKK-1"), c(length(ctrl), length(dapt), length(dkk))), levels = c("CTRL", "DAPT", "DKK-1"))
  fold <- v / mean(ctrl)
  m <- tapply(fold, g, mean); s <- tapply(fold, g, sem)
  med <- tapply(v, g, median) / median(ctrl)
  kw <- kruskal.test(v ~ g)$p.value
  pw <- suppressWarnings(pairwise.wilcox.test(v, g, p.adjust.method = "bonferroni")$p.value)
  lg <- log10(v)
  an <- summary(aov(lg ~ g))[[1]][["Pr(>F)"]][1]
  pt <- pairwise.t.test(lg, g, p.adjust.method = "bonferroni")$p.value
  data.frame(scenario = label, n_CTRL = length(ctrl), n_DAPT = length(dapt), n_DKK1 = length(dkk),
             CTRL_mean = m["CTRL"], CTRL_sem = s["CTRL"], DAPT_mean = m["DAPT"], DAPT_sem = s["DAPT"], DKK1_mean = m["DKK-1"], DKK1_sem = s["DKK-1"],
             DAPT_median_fold = med["DAPT"], DKK1_median_fold = med["DKK-1"],
             KW_p = kw, Wilcox_bonf_CTRL_DAPT = pw["DAPT", "CTRL"], Wilcox_bonf_CTRL_DKK1 = pw["DKK-1", "CTRL"], Wilcox_bonf_DAPT_DKK1 = pw["DKK-1", "DAPT"],
             logANOVA_p = an, logT_bonf_CTRL_DAPT = pt["DAPT", "CTRL"], logT_bonf_CTRL_DKK1 = pt["DKK-1", "CTRL"], logT_bonf_DAPT_DKK1 = pt["DKK-1", "DAPT"],
             note = note, row.names = NULL, stringsAsFactors = FALSE)
}

r2 <- blk("25K_read2"); r1 <- blk("25K_read1"); h1 <- blk("50K_read1"); h2 <- blk("50K_read2")
res <- list()
# --- 25K, read 2 (the block used for the published figure) ---
res[[1]] <- run("S1 published: 25K read 2, CTRL rep2 = 331, CTRL 3382 excluded (CTRL only), no blank correction", pubv[1:3], r2$DAPT, r2$DKK1, "reproduces Supp Data 2 Fig_2b_stats")
res[[2]] <- run("S2 as S1 but CTRL rep2 = 184 (raw reading block)", c(r2$CTRL[1:3]), r2$DAPT, r2$DKK1)
res[[3]] <- run("S3 25K read 2, CTRL rep2 = 331, all four CTRL kept", pubv, r2$DAPT, r2$DKK1)
res[[4]] <- run("S4 25K read 2, CTRL rep2 = 184, all four CTRL kept", r2$CTRL, r2$DAPT, r2$DKK1)
res[[5]] <- run("S5 25K read 2, rep2 = 184, 3382 excluded, blank-subtracted by row", (r2$CTRL - r2$blank)[1:3], r2$DAPT - r2$blank, r2$DKK1 - r2$blank)
res[[6]] <- run("S6 25K read 2, rep2 = 184, all four kept, blank-subtracted by row", r2$CTRL - r2$blank, r2$DAPT - r2$blank, r2$DKK1 - r2$blank)
mb <- mean(r2$blank)
res[[7]] <- run("S7 25K read 2, rep2 = 184, all four kept, mean blank subtracted", r2$CTRL - mb, r2$DAPT - mb, r2$DKK1 - mb)
# --- 25K, read 1 ---
res[[8]] <- run("S8 25K read 1, all four CTRL kept", r1$CTRL, r1$DAPT, r1$DKK1)
res[[9]] <- run("S9 25K read 1, high CTRL (3435) excluded", r1$CTRL[1:3], r1$DAPT, r1$DKK1)
# --- 50K ---
res[[10]] <- run("S10 50K read 1, all four kept (no outlier)", h1$CTRL, h1$DAPT, h1$DKK1)
res[[11]] <- run("S11 50K read 2, all four kept (no outlier)", h2$CTRL, h2$DAPT, h2$DKK1)
res[[12]] <- run("S12 50K read 2, blank-subtracted by row", h2$CTRL - h2$blank, h2$DAPT - h2$blank, h2$DKK1 - h2$blank)
out <- do.call(rbind, res)
write.csv(out, outp, row.names = FALSE)

# --- paired-by-row ratio used in the original workbook (% 2-NBDG uptake = treated / CTRL of the same row x 100) ---
paired <- function(label, ctrl, x) { r <- 100 * x / ctrl; c(label = label, mean = round(mean(r), 1), sem = round(sem(r), 1)) }
pr <- rbind(
  c(g = "DAPT",  paired("25K read 2, raw (rep2 = 184)", r2$CTRL, r2$DAPT)), c(g = "DKK-1", paired("25K read 2, raw (rep2 = 184)", r2$CTRL, r2$DKK1)),
  c(g = "DAPT",  paired("25K read 2, rep2 = 331",       pubv,    r2$DAPT)), c(g = "DKK-1", paired("25K read 2, rep2 = 331",       pubv,    r2$DKK1)),
  c(g = "DAPT",  paired("50K read 2",                   h2$CTRL, h2$DAPT)), c(g = "DKK-1", paired("50K read 2",                   h2$CTRL, h2$DKK1)))
cat("\n=== Main table (fold change relative to the mean of the CTRL values used; mean +- SEM) ===\n")
f <- function(x) sprintf("%.3f", x)
for (i in seq_len(nrow(out))) {
  o <- out[i, ]
  cat(sprintf("%s\n   n %d/%d/%d | CTRL %s+-%s | DAPT %s+-%s (median %s) | DKK-1 %s+-%s (median %s)\n   KW p %s | Wilcoxon-Bonf C-D %s C-K %s D-K %s | log-ANOVA p %s | log-t Bonf C-D %s C-K %s D-K %s\n",
              o$scenario, o$n_CTRL, o$n_DAPT, o$n_DKK1, f(o$CTRL_mean), f(o$CTRL_sem), f(o$DAPT_mean), f(o$DAPT_sem), f(o$DAPT_median_fold),
              f(o$DKK1_mean), f(o$DKK1_sem), f(o$DKK1_median_fold), f(o$KW_p), f(o$Wilcox_bonf_CTRL_DAPT), f(o$Wilcox_bonf_CTRL_DKK1),
              f(o$Wilcox_bonf_DAPT_DKK1), f(o$logANOVA_p), f(o$logT_bonf_CTRL_DAPT), f(o$logT_bonf_CTRL_DKK1), f(o$logT_bonf_DAPT_DKK1)))
}
cat("\n=== Paired-by-row % of CTRL (original workbook analysis) ===\n"); print(as.data.frame(pr), row.names = FALSE)
cat("\nBlank (negative-control) wells, 25K read 2:", r2$blank, "| treated-well signals DAPT:", r2$DAPT, "DKK-1:", r2$DKK1, "\n")
