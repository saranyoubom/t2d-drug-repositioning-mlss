# R2C3: Post hoc amplification-efficiency estimation for every well of the
# actual RT-qPCR experimental runs, using the published LinRegPCR-equivalent
# sliding-window log-linear regression method (sliwin(), qpcR package;
# Ritz & Spiess 2008, Bioinformatics).
#
# Run after R2C3_extract_all_wells.py has produced all_wells_wide.csv.
#
# qpcR's normal library(qpcR) load path is blocked here by a broken `rgl`
# Depends; R2C3_load_qpcR_functions.R sources qpcR's R/ functions directly
# from an existing qpcR install, bypassing that. If qpcR loads normally in
# your environment, `library(qpcR)` is sufficient and this source() line can
# be replaced.
source("./R2C3_load_qpcR_functions.R")

PLOT_DIR <- "./sliwin_plots"
dir.create(PLOT_DIR, showWarnings = FALSE)

df <- read.csv("./all_wells_wide.csv", check.names = FALSE)
cyc_cols <- grep("^c[0-9]+$", names(df), value = TRUE)
cyc_cols <- cyc_cols[order(as.integer(sub("c", "", cyc_cols)))]
cycles_full <- as.integer(sub("c", "", cyc_cols))

results <- data.frame(Run=character(), Well=character(), Gene=character(),
                       fit_status=character(), eff=numeric(), Eff_pct=numeric(),
                       rsq=numeric(), win_lo=numeric(), win_hi=numeric(),
                       stringsAsFactors = FALSE)

plotted_genes <- character()

cat(sprintf("Processing %d wells...\n", nrow(df)))
for (i in seq_len(nrow(df))) {
  run <- df$Run[i]; well <- df$Well[i]; gene <- df$Gene[i]
  fluo_raw <- as.numeric(df[i, cyc_cols])
  valid <- !is.na(fluo_raw)
  cyc <- cycles_full[valid]; fluo <- fluo_raw[valid]
  if (length(cyc) < 10) next

  if (min(fluo, na.rm = TRUE) <= 0) fluo <- fluo - min(fluo, na.rm = TRUE) + 1
  d <- data.frame(cyc = cyc, fluo = fluo)

  fit <- try(pcrfit(d, cyc = "cyc", fluo = "fluo", model = b4, verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    results[nrow(results)+1,] <- list(run, well, gene, "FIT_FAILED", NA, NA, NA, NA, NA)
    next
  }
  sw <- try(sliwin(fit, plot = FALSE, verbose = FALSE), silent = TRUE)
  if (inherits(sw, "try-error")) {
    results[nrow(results)+1,] <- list(run, well, gene, "SLIWIN_FAILED", NA, NA, NA, NA, NA)
    next
  }
  results[nrow(results)+1,] <- list(run, well, gene, "OK", sw$eff, (sw$eff-1)*100, sw$rsq,
                                     sw$window[1], sw$window[2])

  # Save one representative converged plot per gene
  if (!(gene %in% plotted_genes)) {
    png_path <- file.path(PLOT_DIR, paste0(gene, "_", run, "_", gsub(":","-",well), ".png"))
    tryCatch({
      png(png_path, width = 900, height = 500)
      par(mfrow = c(1,2))
      plot(fit, main = paste(gene, run, well, "- sigmoidal fit"))
      sliwin(fit, plot = TRUE, verbose = FALSE)
      title(main = paste("Sliding-window log-linear region, eff=", round(sw$eff,3)))
      dev.off()
      plotted_genes <<- c(plotted_genes, gene)
    }, error = function(e) { try(dev.off(), silent = TRUE) })
  }
}

write.csv(results, "./sliwin_all_results.csv", row.names = FALSE)

cat("\n=== Per-gene convergence summary ===\n")
for (g in sort(unique(results$Gene))) {
  sub <- results[results$Gene == g,]
  n_ok <- sum(sub$fit_status == "OK")
  n_tot <- nrow(sub)
  eff_ok <- sub$Eff_pct[sub$fit_status == "OK"]
  if (n_ok > 0) {
    cat(sprintf("%-10s %d/%d converged. eff%% range [%.1f, %.1f], mean %.1f\n",
                g, n_ok, n_tot, min(eff_ok), max(eff_ok), mean(eff_ok)))
  } else {
    cat(sprintf("%-10s %d/%d converged. NO valid estimates.\n", g, n_ok, n_tot))
  }
}
cat(sprintf("\nPlots saved for %d genes to %s\n", length(plotted_genes), PLOT_DIR))
cat("Genes with a saved plot:", paste(plotted_genes, collapse=", "), "\n")
