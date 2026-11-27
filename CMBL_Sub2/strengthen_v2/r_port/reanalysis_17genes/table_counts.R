## table_counts.R - R port of reanalysis_17genes/table_counts.py: counts behind the manuscript's pairs table for the 20 and 17 genes: per scheme the pairs, pooled p < 0.05, pooled BH q < 0.05,
## treatment-adjusted p < 0.05 and q < 0.05; over the pairs common to all schemes the numbers significant in all four (p, pooled q, adjusted p, adjusted q) and the sign-consistent ones. Writes table_counts.json.
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R")); source(file.path(HERE, "pipeline.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
counts <- function(genes) {
  M <- P_load(genes); sch <- P_schemes(M); res <- P_sig_sets(sch); adj <- P_adjusted(sch); out <- list(); pooled_q <- list()
  for (k in names(res)) { up <- upper.tri(res[[k]]$p); q <- matrix(NA_real_, nrow(up), ncol(up)); q[up] <- p.adjust(res[[k]]$p[up], "BH"); q[lower.tri(q)] <- t(q)[lower.tri(q)]; pooled_q[[k]] <- list(genes = res[[k]]$genes, q = q) }   # symmetric, see P_adjusted
  for (k in names(res)) {
    up <- upper.tri(res[[k]]$p)
    out[[k]] <- list(pairs = sum(up), pooled_p = sum(res[[k]]$p[up] < 0.05), pooled_q = sum(pooled_q[[k]]$q[up] < 0.05), adj_p = sum(adj[[k]]$p[up] < 0.05), adj_q = sum(adj[[k]]$q[up] < 0.05))
  }
  cg <- .P_common(res); allf <- function(L, f) { m <- Reduce(`&`, lapply(L, f)); m[is.na(m)] <- FALSE; m }
  Pp <- .P_sub(res, "p", cg); Rr <- .P_sub(res, "r", cg); Qq <- .P_sub(pooled_q, "q", cg); AP <- .P_sub(adj, "p", cg); AQ <- .P_sub(adj, "q", cg)
  allp <- allf(Pp, function(x) x < 0.05); allq <- allf(Qq, function(x) x < 0.05); adjp <- allf(AP, function(x) x < 0.05); adjq <- allf(AQ, function(x) x < 0.05)
  s1 <- sign(Rr[[1]]); same <- allf(Rr[-1], function(r) sign(r) == s1); up <- upper.tri(allp)
  out$all <- list(pooled_p = sum(allp & up), pooled_q = sum(allq & up), adj_p = sum(adjp & up), adj_q = sum(adjq & up), signcons = sum(allp & same & up))
  out
}
o20 <- counts(P_GENES20); o17 <- counts(P_GENES17)
pr <- function(o) for (k in names(o)) cat("  ", k, paste0("{", paste(sprintf("'%s': %d", names(o[[k]]), unlist(o[[k]])), collapse = ", "), "}"), "\n")
cat("20 genes (compare with the existing table: raw 190/80/52/15/0; Hes1 171/169/169/43/12; Kcnj11 171/118/111/54/5; Ins1 171/98/89/76/49; all-four 27/18/5/0; sign-cons 21)\n"); pr(o20)
cat("17 genes\n"); pr(o17)
write_json_typed(list(`20` = o20, `17` = o17), file.path(OUT, "table_counts.json"))
