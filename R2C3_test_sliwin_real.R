source("R2C3_load_qpcR_functions.R")

library(readxl)

# Point this at a local folder containing the raw CFX Maestro well-level
# fluorescence exports (.xlsx, baseline-subtracted). Not included in this
# deposit -- raw experimental data is provided separately per journal policy
# (see R_scripts_Rev1/README.md, "Input Data Files Required").
FIXED <- "./raw_cfx_exports"

run_well <- function(path, sheet, well_col_name) {
  df <- as.data.frame(read_excel(path, sheet = sheet))
  cyc <- df[[2]]  # column B = Cycle
  fluo <- df[[well_col_name]]
  # qpcR's sigmoidal models (fit in linear RFU space, not log) need positive
  # fluorescence; CFX Maestro's baseline-subtracted export can dip negative.
  # A constant vertical shift doesn't distort the sigmoid's shape/timing (unlike
  # a shift before a log-transform, which was the earlier retracted script's bug),
  # so this is safe here.
  if (min(fluo, na.rm = TRUE) <= 0) fluo <- fluo - min(fluo, na.rm = TRUE) + 1
  d <- data.frame(cyc = cyc, fluo = fluo)
  fit <- try(pcrfit(d, cyc = "cyc", fluo = "fluo", model = b4, verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat(well_col_name, "FIT FAILED:", attr(fit,"condition")$message, "\n")
    return(NULL)
  }
  sw <- try(sliwin(fit, plot = FALSE, verbose = FALSE), silent = TRUE)
  if (inherits(sw, "try-error")) {
    cat(well_col_name, "SLIWIN FAILED:", attr(sw,"condition")$message, "\n")
    return(NULL)
  }
  cat(sprintf("%s: eff=%.4f (%.1f%%) rsq=%.4f window=[%s]\n",
              well_col_name, sw$eff, (sw$eff - 1) * 100, sw$rsq,
              paste(round(sw$window,1), collapse=",")))
  return(sw)
}

cat("=== Kcnj11 (Aug26 run) -- compare against published Suppl Fig: DKK-1=1.699, DAPT=1.682, CTRL=1.641 ===\n")
df_hdr <- as.data.frame(read_excel(file.path(FIXED, "gr_aug26.xlsx"), sheet = "Kcnj11", n_max = 0))
wells <- names(df_hdr)[3:ncol(df_hdr)]
for (w in wells) {
  run_well(file.path(FIXED, "gr_aug26.xlsx"), "Kcnj11", w)
}

cat("\n=== Gapdh well A10 (Aug26 run) -- I manually sanity-checked this earlier, found eff~62% in the true exponential region ===\n")
run_well(file.path(FIXED, "gr_aug26.xlsx"), "Gapdh", names(as.data.frame(read_excel(file.path(FIXED,"gr_aug26.xlsx"), sheet="Gapdh", n_max=0)))[3])
