## wp4_meta_analysis.R - R port of wp4_islets/meta_analysis.py: protocol 4.4 (primary endpoint of WP4): random-effects pooling of the excess E of the consensus pairs over matched random pairs across the human public
## islet datasets, the pre-specified decision rule (margin +/-0.10 on the E scale), leave-one-out and estimator sensitivities (PRESPEC_WP4_choices.md section B), and pair-level pooled rho for consensus pairs
## evaluable in at least 3 human datasets. Estimators (fixed in advance): tau-squared by REML; Hartung-Knapp standard error with t(k-1) quantile and no floor = PRIMARY interval; sensitivity (never used for the rule):
## Knapp-Hartung factor floored at 1, DerSimonian-Laird tau-squared with the Hartung-Knapp interval, DerSimonian-Laird with the normal interval, fixed effect. Datasets qualify with >= 3 evaluable consensus pairs;
## with < 2 qualifying datasets the rule is not applied. Input: <WP4_RES>/<tag>__<variant>__summary.json and <tag>__<variant>__adjusted__cons.csv from run_datasets (default: the archived wp4_islets/results).
## Output (WP4_OUT, default results_R/wp4/meta): study_effects__*.csv, pair_level_pooled_rho.csv, meta_summary.json, meta_report.md.  Usage: Rscript wp4_meta_analysis.R
## Sums follow numpy (np_sum: sequential below 8 values), so that the iteration of the REML estimate gives the same digits as the Python run.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- HERE
source(file.path(R_PORT, "pycompat.R")); suppressMessages(library(jsonlite)); options(py_json_nan = TRUE)
RES <- Sys.getenv("WP4_RES", normalizePath(file.path(R_PORT, "..", "wp4_islets", "results"))); OUT <- Sys.getenv("WP4_OUT", file.path(R_PORT, "results_R", "wp4", "meta")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
HUMAN <- c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061")
LABEL <- c(gse81608 = "GSE81608 (Xin 2016)", gse84133_human = "GSE84133 human (Baron 2016)", gse86469 = "GSE86469 (Lawlor 2017)", gse81547 = "GSE81547 (Enge 2017)", emtab5061 = "E-MTAB-5061 (Segerstolpe 2016)"); MARGIN <- 0.10
source(file.path(R_PORT, "wp4_islets", "wp4_meta_lib.R"))      # ssum, reml_tau2, dl_tau2, pool, verdict (shared with the other WP4 scripts that pool)
study_table <- function(variant, mode) {
  rows <- list()
  for (tag in HUMAN) {
    p <- file.path(RES, sprintf("%s__%s__summary.json", tag, variant)); if (!file.exists(p)) next; s <- fromJSON(p, simplifyVector = FALSE)[[mode]]; g <- function(k) if (is.null(s[[k]])) NA_real_ else s[[k]]
    rows[[length(rows) + 1]] <- list(tag = tag, label = LABEL[[tag]], variant = variant, mode = mode, beta_cells = s$beta_cells, donors = s$donors, genes_evaluable = s$genes_evaluable, pairs = s$consensus_pairs_evaluable, E = g("E"), SE = g("SE_E"), median_rho = g("median_rho_consensus"), null_rho = g("null_median_rho_mean"))
  }
  T <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)); T$qualifies <- T$pairs >= 3; T
}
fs <- function(x, nd) paste0(if (x < 0) "-" else "+", py_fmt_f(abs(x), nd))       # '%+.3f'
run <- function() {
  variants <- list()
  for (variant in c("primary", "nondiabetic")) for (mode in c("adjusted", "unadjusted")) {
    T <- study_table(variant, mode)
    if (variant == "nondiabetic") { P <- study_table("primary", mode); extra <- P[!(P$tag %in% T$tag), , drop = FALSE]; extra$variant <- "nondiabetic"; T <- rbind(T, extra); T <- T[match(HUMAN[HUMAN %in% T$tag], T$tag), , drop = FALSE]; rownames(T) <- NULL }
    py_write_csv_repr(T, file.path(OUT, sprintf("study_effects__%s__%s.csv", variant, mode)), eol = "\r\n"); Q <- if (sum(T$qualifies) >= 2) T[T$qualifies, , drop = FALSE] else NULL
    rec <- list(variant = variant, mode = mode, qualifying = as.list(if (is.null(Q)) character(0) else Q$tag), n_qualifying = if (is.null(Q)) 0L else nrow(Q))
    if (!is.null(Q)) {
      y <- Q$E; se <- Q$SE; rec$estimators <- list()
      for (e in list(list("REML+HK (primary)", list(tau = "reml", ci = "hk", floor = FALSE)), list("REML+HK, factor floored at 1", list(tau = "reml", ci = "hk", floor = TRUE)), list("DL+HK", list(tau = "dl", ci = "hk")), list("DL, normal CI", list(tau = "dl", ci = "normal")), list("fixed effect", list(ci = "fixed")))) {
        r <- do.call(pool, c(list(y, se), e[[2]])); r$estimator <- e[[1]]; r$verdict <- verdict(r$lo, r$hi); rec$estimators[[length(rec$estimators) + 1]] <- r
      }
      if (nrow(Q) >= 3) { rec$leave_one_out <- list(); for (t in Q$tag) { R <- Q[Q$tag != t, , drop = FALSE]; r <- pool(R$E, R$SE); r$left_out <- t; r$verdict <- verdict(r$lo, r$hi); rec$leave_one_out[[length(rec$leave_one_out) + 1]] <- r } }
    }
    variants[[length(variants) + 1]] <- rec
  }
  out <- list(variants = variants)
  pairs <- list()      # pair-level pooled rho (adjusted, primary, human datasets)
  for (tag in HUMAN) {
    p <- file.path(RES, sprintf("%s__primary__adjusted__cons.csv", tag)); sp <- file.path(RES, sprintf("%s__primary__summary.json", tag)); if (!(file.exists(p) && file.exists(sp))) next
    df <- fromJSON(sp, simplifyVector = FALSE)$adjusted$df_for_p; C <- read.csv(p, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(C))) { k <- paste(C$gene1[i], C$gene2[i], sep = "\t"); pairs[[k]] <- c(pairs[[k]], list(list(tag = tag, rho = C$rho[i], df = df))) }
  }
  rows <- list()
  for (k in names(pairs)) {
    L <- pairs[[k]]; if (length(L) < 3) next; g <- strsplit(k, "\t", fixed = TRUE)[[1]]; z <- vapply(L, function(x) atanh(x$rho), 1); se <- vapply(L, function(x) sqrt(1.06 / (x$df - 1)), 1); r <- pool(z, se)
    rows[[length(rows) + 1]] <- list(gene1 = g[1], gene2 = g[2], datasets = length(L), tags = paste(vapply(L, function(x) x$tag, ""), collapse = ";"), rho_by_dataset = paste(vapply(L, function(x) py_fmt_f(x$rho, 3), ""), collapse = ";"), pooled_rho = tanh(r$mu), lo = tanh(r$lo), hi = tanh(r$hi), tau2 = r$tau2, I2 = r$I2)
  }
  py_write_csv_repr(if (length(rows)) do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)) else data.frame(), file.path(OUT, "pair_level_pooled_rho.csv"), eol = "\r\n"); out$pairs_evaluable_in_ge3 <- length(rows)
  write_json_typed(out, file.path(OUT, "meta_summary.json"))
  L <- c("# WP4 meta-analysis (protocol 4.4)", "", sprintf("Primary = depth-adjusted E, all frozen beta-cells, REML tau-squared + Hartung-Knapp 95%% CI (t, k - 1 df), margin +/-%s on the E (Fisher z) scale.", py_fmt_f(MARGIN, 2)), "")
  for (rec in variants) {
    L <- c(L, sprintf("## %s, %s mode: %d qualifying dataset(s): %s", rec$variant, rec$mode, rec$n_qualifying, if (length(rec$qualifying)) paste(unlist(rec$qualifying), collapse = ", ") else "-"))
    if (!is.null(rec$estimators)) {
      for (r in rec$estimators) L <- c(L, sprintf("* %-30s pooled E %s  95%% CI [%s, %s]  tau2 %s  I2 %s%%  -> %s", r$estimator, fs(r$mu, 3), fs(r$lo, 3), fs(r$hi, 3), py_fmt_f(r$tau2, 4), py_fmt_f(100 * r$I2, 0), r$verdict))
      for (r in rec$leave_one_out) L <- c(L, sprintf("  leave out %-15s pooled E %s  [%s, %s] -> %s", r$left_out, fs(r$mu, 3), fs(r$lo, 3), fs(r$hi, 3), r$verdict))
    } else L <- c(L, "* decision rule not applied (fewer than 2 qualifying datasets)")
    L <- c(L, "")
  }
  write_text_bytes(paste(L, collapse = .NL), file.path(OUT, "meta_report.md")); cat(paste(L, collapse = "\n"), "\n")
}
if (sys.nframe() == 0L) run()
