## perturbation_gse60803.R - R port of wp4_islets/perturbation_gse60803.py: protocol 4.6, step 3: log2 fold changes of the 17 genes in GSE60803 (adult human beta-cell-derived cells, beta-catenin shRNA; Affymetrix HuGene-1.1 ST, GPL11532)
## and their agreement in sign with the DKK-1 arm of the INS-1 experiment. Design and contrast fixed in PRESPEC_WP4_choices.md addendum G BEFORE any expression value was read: per donor: control = the untreated array (S108, S109, S112) or the
## nontarget array (S114); perturbed = mean of the donor's beta-catenin shRNA arrays (9 and 45; both replicates of S114); log2 fold change of a gene = mean over the 4 donors of (perturbed - control) on the log2 scale of the series matrix; the probe
## set of a gene (Entrez id match in the GPL11532 annotation) with the highest mean control expression is used. Secondary, descriptive: donor S114 only. INS-1 side: log2 fold change DKK-1 vs CTRL, 'as reported', RT-qPCR panel, from
## wp5_reporting/effect_sizes.csv (ES_CSV). Agreement: genes with both fold changes nonzero, exact two-sided binomial test of the number of equal signs against 0.5, Spearman correlation. DESCRIPTIVE.
## Output (WP4_PERT_OUT, default results_R/wp4/perturbation): gse60803_fold_changes.csv, gse60803_concordance.json.   Usage: Rscript perturbation_gse60803.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp4_core.R")); source(file.path(HERE, "wp4_stats.R")); options(py_json_nan = TRUE)
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); OUT <- Sys.getenv("WP4_PERT_OUT", file.path(R_PORT, "results_R", "wp4", "perturbation")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ES_CSV <- Sys.getenv("ES_CSV", normalizePath(file.path(WP4_DIR, "..", "wp5_reporting", "effect_sizes.csv")))
run <- function() {
  D <- file.path(DATA, "GSE60803"); sm <- read_series_matrix(file.path(D, "GSE60803_series_matrix.txt.gz")); X <- as.matrix(sm$table[, -1]); rownames(X) <- sm$table[[1]]; storage.mode(X) <- "double"
  stopifnot("sample order differs" = identical(colnames(X), sm$gsm)); S <- data.frame(gsm = sm$gsm, title = sm$title, stringsAsFactors = FALSE); S$donor <- sub("^.*donor (S[0-9]+).*$", "\\1", S$title); S$donor[!grepl("donor (S[0-9]+)", S$title)] <- NA
  S$arm <- ifelse(grepl("untreated", S$title, fixed = TRUE), "control", ifelse(grepl("nontarget", S$title, fixed = TRUE), "control", ifelse(grepl("β-catenin", S$title, fixed = TRUE), "bcat", "unknown"))); stopifnot(all(S$arm != "unknown"), length(unique(S$donor)) == 4)
  cat(sprintf("expression range of the series matrix: %.2f to %.2f (median %.2f), probe sets %d\n", min(X), max(X), median(X), nrow(X))); stopifnot("the matrix is not on a log2 scale" = max(X) < 30)
  con <- gzfile(file.path(D, "GPL11532.annot.gz"), "rt"); A <- readLines(con, warn = FALSE, encoding = "UTF-8"); close(con); A <- A[!(startsWith(A, "#") | startsWith(A, "!") | startsWith(A, "^"))]; hdr <- strsplit(A[1], "\t", fixed = TRUE)[[1]]; iG <- match("Gene ID", hdr)
  want <- as.character(ENTREZ[GENES]); P <- strsplit(A[-1], "\t", fixed = TRUE); ann <- setNames(vector("list", length(want)), want)      # Entrez id -> probe sets, only the ids of the 17 genes (file order)
  for (p in P) { if (length(p) < iG) next; ids <- trimws(strsplit(p[iG], "\\s*///\\s*")[[1]]); for (gid in intersect(ids, want)) ann[[gid]] <- c(ann[[gid]], p[1]) }
  ctrl <- S[S$arm == "control", ]; rownames(ctrl) <- ctrl$donor; pert <- S[S$arm == "bcat", ]; donors <- sort(unique(S$donor), method = "radix"); res <- list()
  for (g in GENES) {
    probes <- intersect(ann[[as.character(ENTREZ[[g]])]], rownames(X)); probes <- ann[[as.character(ENTREZ[[g]])]][ann[[as.character(ENTREZ[[g]])]] %in% rownames(X)]
    if (!length(probes)) { res[[length(res) + 1]] <- list(gene = g); next }
    cm <- rowMeans(X[probes, ctrl$gsm, drop = FALSE]); best <- probes[which.max(cm)]; d <- vapply(donors, function(dn) mean(X[best, pert$gsm[pert$donor == dn]]) - X[best, ctrl[dn, "gsm"]], 1); s114 <- mean(X[best, pert$gsm[pert$donor == "S114"]]) - X[best, ctrl["S114", "gsm"]]
    res[[length(res) + 1]] <- c(list(gene = g, probe = best, probes_for_gene = length(probes), mean_control_log2 = unname(cm[best])), setNames(as.list(unname(d)), paste0("delta_", donors)), list(log2fc_all_donors = mean(d), sd_across_donors = sd(d), log2fc_S114_only = s114))
  }
  cols <- c("gene", "probe", "probes_for_gene", "mean_control_log2", paste0("delta_", donors), "log2fc_all_donors", "sd_across_donors", "log2fc_S114_only")
  R <- do.call(rbind, lapply(res, function(r) { x <- lapply(cols, function(k) if (is.null(r[[k]])) NA else r[[k]]); names(x) <- cols; as.data.frame(x, stringsAsFactors = FALSE) }))
  E <- read.csv(ES_CSV, stringsAsFactors = FALSE, check.names = FALSE); E <- E[E$readout == "RT-qPCR" & E$contrast == "DKK-1 vs CTRL" & E$variant == "as reported", ]; ev <- setNames(E$log2fc_geo, E$feature)
  R$ins1_log2fc_DKK1 <- unname(ev[R$gene]); agree <- function(a, b) { x <- sign(a) == sign(b); x[is.na(x)] <- FALSE; x }; R$sign_agree_all <- agree(R$log2fc_all_donors, R$ins1_log2fc_DKK1); R$sign_agree_S114 <- agree(R$log2fc_S114_only, R$ins1_log2fc_DKK1)
  py_write_csv_repr(R, file.path(OUT, "gse60803_fold_changes.csv"), eol = "\r\n")
  out <- list(series = "GSE60803", platform = "GPL11532", donors = as.list(donors), arrays = nrow(S), genes_with_probe = sum(!is.na(R$probe)), genes_without_probe = as.list(R$gene[is.na(R$probe)]), ins1_arm = "DKK-1 vs CTRL, as reported",
              confounded = "controls are not expanded cells (untreated and nontarget arrays are described as intact islet cells; nontarget exists for one donor)")
  for (spec in list(list("all_donors", "log2fc_all_donors", "sign_agree_all"), list("S114_only", "log2fc_S114_only", "sign_agree_S114"))) {
    col <- spec[[2]]; ok <- !is.na(R[[col]]) & !is.na(R$ins1_log2fc_DKK1); ok <- ok & R[[col]] != 0 & R$ins1_log2fc_DKK1 != 0; Dd <- R[ok, ]; k <- sum(Dd[[spec[[3]]]]); n <- nrow(Dd); sp <- if (n > 2) py_spearmanr(Dd[[col]], Dd$ins1_log2fc_DKK1) else list(rho = NaN, p = NaN)
    out[[spec[[1]]]] <- list(n_genes = as.integer(n), n_same_sign = as.integer(k), binomial_p = if (n) py_binom_p(k, n) else NULL, spearman_rho = sp$rho, spearman_p = sp$p, genes = as.list(Dd$gene), same_sign_genes = as.list(Dd$gene[Dd[[spec[[3]]]]]), opposite_sign_genes = as.list(Dd$gene[!Dd[[spec[[3]]]]]))
  }
  write_json_typed(out, file.path(OUT, "gse60803_concordance.json")); print(R); cat(py_json_typed(out[c("all_donors", "S114_only")]), "\n")
}
if (sys.nframe() == 0L) run()
