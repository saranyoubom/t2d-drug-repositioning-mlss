## check_wp4b_numbers.R - R port of wp4_islets/check_wp4b_numbers.py: every number of README_WP4b.md and WORDING_islet_crosscheck_v3.md is rebuilt from positive_control/, wider_search/ and meta/ and searched for in the format used in the text
## (a typing slip is caught). Prints one line per check and the number of failures (exit status 1 if any). No analysis. The two texts are the archived ones (wp4_islets/); the numbers come from the R run (results_R/wp4/...) or, with
## CHECK_SOURCE=archive, from the archived Python results. The record of the live queries (wider_search/queries.json, counts of 6 Nov 2026) is always the archived one: the databases have moved on since.
## Usage: Rscript check_wp4b_numbers.R          CHECK_SOURCE = R (default) | archive
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); suppressMessages(library(jsonlite))
ARCH <- normalizePath(file.path(R_PORT, "..", "wp4_islets")); W4 <- file.path(R_PORT, "results_R", "wp4"); SRC <- Sys.getenv("CHECK_SOURCE", "R")
at <- function(rdir, adir, file) if (SRC == "archive") file.path(ARCH, adir, file) else file.path(W4, rdir, file)
read_text <- function(p) { b <- readBin(p, "raw", file.size(p)); s <- rawToChar(b); Encoding(s) <- "UTF-8"; s }
RD <- read_text(file.path(ARCH, "README_WP4b.md")); WD <- read_text(file.path(ARCH, "WORDING_islet_crosscheck_v3.md")); fails <- character(0)
has <- function(text, label, s) { ok <- grepl(enc2utf8(s), text, fixed = TRUE, useBytes = TRUE); cat(if (ok) "ok   " else "FAIL ", label, "   [", substr(s, 1, 80), "]\n", sep = ""); if (!ok) fails <<- c(fails, label) }
sg <- function(x, d = 3) { s <- py_fmt_f(x, d); ifelse(startsWith(s, "-"), s, paste0("+", s)) }
ff <- function(x, d) py_fmt_f(x, as.integer(d)); ri <- function(x) as.integer(round(x))
rcsv <- function(p) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
P <- rcsv(at("positive_control", "positive_control", "positive_control_pooled.csv")); A <- rcsv(at("positive_control", "positive_control", "positive_control_by_dataset.csv")); RP <- rcsv(at("positive_control", "positive_control", "rp_control.csv"))
pick <- function(D, mode, rho) D[D$mode == mode & D$rho == rho, ][1, ]
for (rho in c(0.05, 0.1, 0.2, 0.3, 0.5)) {
  d <- pick(P, "detected", rho); a <- pick(P, "all", rho)
  has(RD, sprintf("pooled row rho %s", format(rho)), sprintf("%s / %s / %s | %s", ff(d$share_no_meaningful_excess, 2), ff(d$share_inconclusive, 2), ff(d$share_co_expressed_more, 2), sg(d$pooled_E_mean)))
  has(RD, sprintf("pooled row rho %s (all)", format(rho)), sprintf("%s / %s / %s | %s", ff(a$share_no_meaningful_excess, 2), ff(a$share_inconclusive, 2), ff(a$share_co_expressed_more, 2), sg(a$pooled_E_mean)))
}
d3 <- pick(P, "detected", 0.3)
has(RD, "false reassurance at 0.3", sprintf("in only %s%% of 200 simulations (%d of 200", ff(100 * d3$share_no_meaningful_excess, 1), ri(200 * d3$share_no_meaningful_excess)))
has(RD, "inconclusive 97.5%", sprintf("\"inconclusive\" in %s%%", ff(100 * d3$share_inconclusive, 1)))
has(RD, "R1 value", sprintf("MET: %s%% (%d of 200)", ff(100 * d3$share_no_meaningful_excess, 1), ri(200 * d3$share_no_meaningful_excess)))
for (rho in c(0.2, 0.1)) has(RD, sprintf("summary share at %s", format(rho)), sprintf("at %s in %d%%", format(rho), ri(100 * pick(P, "detected", rho)$share_no_meaningful_excess)))
x <- pick(P, "detected", 0.05); has(RD, "summary share at 0.05", sprintf("at 0.05 in %s%% (%d of 200)", ff(100 * x$share_no_meaningful_excess, 1), ri(200 * x$share_no_meaningful_excess)))
ord <- c("emtab5061", "gse81547", "gse81608", "gse86469", "gse84133_mouse", "gse84133_human")
nm <- c(emtab5061 = "E-MTAB-5061", gse81547 = "GSE81547", gse81608 = "GSE81608", gse86469 = "GSE86469", gse84133_mouse = "GSE84133 mouse", gse84133_human = "GSE84133 human")
D3 <- A[A$mode == "detected" & A$rho == 0.3, ]; rownames(D3) <- D3$tag
for (t in ord) has(RD, paste("recovered E", t), sprintf("%s %s (%s)", nm[[t]], sg(D3[t, "E_recovered_mean"]), sg(D3[t, "E_primary"])))
has(RD, "attenuation list", "attenuation of the injected effect")
for (k in list(c("gse81608", "GSE81608"), c("gse86469", "GSE86469"), c("gse81547", "GSE81547"))) has(RD, paste("power", k[1]), sprintf("%s %s", k[2], ff(D3[k[1], "share_empirical_P_lt_0_05"], 3)))
has(RD, "power mouse", sprintf("mouse %s", ff(D3["gse84133_mouse", "share_empirical_P_lt_0_05"], 3)))
HUM5 <- c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061"); LAB <- c(gse81608 = "GSE81608", gse84133_human = "GSE84133 human", gse86469 = "GSE86469", gse81547 = "GSE81547", emtab5061 = "E-MTAB-5061", gse84133_mouse = "mouse")
for (i in seq_len(nrow(RP))) { r <- RP[i, ]; has(RD, paste("RP E", r$tag), if (r$tag %in% HUM5) sprintf("%s (%s)", sg(r$E_RP_mean), LAB[[r$tag]]) else sprintf("%s (mouse)", sg(r$E_RP_mean))) }
has(RD, "RP SE units", paste0("in units of the SE of the consensus E: ", paste(ff(RP$E_RP_in_SE_units, 1), collapse = ", ")))
H <- RP[RP$tag != "gse84133_mouse", ]; has(RD, "RP range", sprintf("give an excess of %s to %s in all five human datasets", sg(min(H$E_RP_mean)), sg(max(H$E_RP_mean)))); has(RD, "RP 52%", "52% of its sets")
ms <- fromJSON(at("meta", "meta", "meta_summary.json"), simplifyVector = FALSE); vv <- Filter(function(v) v$variant == "primary" && v$mode == "adjusted", ms$variants)[[1]]; e <- Filter(function(z) startsWith(z$estimator, "REML+HK (primary)"), vv$estimators)[[1]]
has(RD, "margin", sprintf("is %s", ff(max(abs(e$lo), abs(e$hi)), 3))); has(RD, "primary interval", sprintf("[%s, %s]", sg(e$lo), sg(e$hi))); has(RD, "primary pooled E", sprintf("pooled E %s", sg(e$mu)))
q <- fromJSON(file.path(ARCH, "wider_search", "queries.json"), simplifyVector = FALSE)
Fz <- rcsv(at("wider_search_freeze", "wider_search", "frozen_candidates.csv")); Cc <- rcsv(at("wider_search", "wider_search", "concordance_results.csv")); M <- rcsv(at("wider_search", "wider_search", "concordance_meta.csv"))
has(RD, "S1 count", sprintf("GEO DataSets %d records", as.integer(q$S1$count))); has(RD, "S3 counts", sprintf("PubMed %d PMIDs linked to %d GEO records", as.integer(q$S3$pmids), as.integer(q$S3$linked_gds_uids)))
has(RD, "S2 hits", sprintf("%d hits (9 outside GEO)", as.integer(q$S2$hits))); has(RD, "unique", sprintf("%d unique records", as.integer(q$unique_records)))
has(RD, "kept/dropped", sprintf("%d kept, %d dropped", sum(Fz$decision == "keep"), sum(Fz$decision == "drop")))
An <- Cc[!is.na(Cc$status) & Cc$status == "analysed", ]
for (i in seq_len(nrow(An))) { r <- An[i, ]; has(RD, paste("concordance", r$accession), sprintf("| %d of 17 | %s | %s (%s) | %d of %d |", as.integer(r$n_same_sign), ff(r$binomial_p, 3), sg(r$spearman_rho), ff(r$spearman_p, 3), as.integer(r$strong_same_sign), as.integer(r$strong_ins1_genes))) }
m <- M[1, ]; has(RD, "pooled rho", sprintf("rho %s, 95%% CI %s to %s", sg(m$pooled_rho), sg(m$lo), sg(m$hi)))
mn <- function(s) gsub("-", "−", s, fixed = TRUE)
has(WD, "wording pooled rho", sprintf("ρ %s, 95%% CI %s to %s", mn(sg(m$pooled_rho, 2)), mn(sg(m$lo, 2)), mn(sg(m$hi, 2))))
has(WD, "wording shares", "97.5% of simulations (70% at 0.2, 13% at 0.1)")
has(WD, "wording ribosomal", "+0.10 to +0.44 in all five human datasets"); has(WD, "wording design numbers", "0.78 to 0.96"); has(WD, "wording design numbers 2", "52% to 95%")
cat("FAILURES: ", length(fails), "\n", sep = ""); if (sys.nframe() == 0L) quit(status = if (length(fails)) 1L else 0L)
