## netw.R - R port of reanalysis_17genes/netw.py: degree, betweenness (normalized) and closeness centrality of the consensus networks of pipeline_20genes.json and pipeline_17genes.json
## (any-sign and sign-consistent consensus pairs). networkx is not needed: Brandes' betweenness and the closeness of networkx (wf_improved=False) are written out below with the same iteration
## order as networkx (nodes in order of first appearance in the edge list, neighbours in order of edge insertion), so ties are broken as in Python. Prints the same table.
suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "rng", "numpy_rng.R")); source(file.path(R_PORT, "pycompat.R"))
OUT <- Sys.getenv("RA_OUT", file.path(R_PORT, "results_R", "reanalysis_17genes"))
nx_graph <- function(pairs) {                         # nx.Graph(); G.add_edges_from(pairs)
  nodes <- character(0); adj <- list()
  for (e in pairs) {
    u <- e[[1]]; v <- e[[2]]
    for (x in c(u, v)) if (!(x %in% nodes)) { nodes <- c(nodes, x); adj[[x]] <- character(0) }
    if (!(v %in% adj[[u]])) adj[[u]] <- c(adj[[u]], v)
    if (!(u %in% adj[[v]])) adj[[v]] <- c(adj[[v]], u)
  }
  list(nodes = nodes, adj = adj)
}
nx_betweenness <- function(G, normalized = TRUE) {    # networkx betweenness_centrality (undirected, k = None, endpoints = False)
  nodes <- G$nodes; n <- length(nodes); bt <- setNames(numeric(n), nodes)
  for (s in nodes) {
    S <- character(0); P <- setNames(replicate(n, character(0), simplify = FALSE), nodes); sigma <- setNames(numeric(n), nodes); D <- setNames(rep(NA_real_, n), nodes); sigma[s] <- 1; D[s] <- 0; Q <- s
    while (length(Q)) {
      v <- Q[1]; Q <- Q[-1]; S <- c(S, v); Dv <- D[[v]]; sigmav <- sigma[[v]]
      for (w in G$adj[[v]]) {
        if (is.na(D[[w]])) { Q <- c(Q, w); D[w] <- Dv + 1 }
        if (D[[w]] == Dv + 1) { sigma[w] <- sigma[[w]] + sigmav; P[[w]] <- c(P[[w]], v) }
      }
    }
    delta <- setNames(numeric(length(S)), S)
    while (length(S)) {
      w <- S[length(S)]; S <- S[-length(S)]; coeff <- (1 + delta[[w]]) / sigma[[w]]
      for (v in P[[w]]) delta[v] <- delta[[v]] + sigma[[v]] * coeff
      if (w != s) bt[w] <- bt[[w]] + delta[[w]]
    }
  }
  if (normalized && n > 2) bt <- bt * (1 / ((n - 1) * (n - 2)))
  bt
}
nx_distances <- function(G, s) {                      # single_source_shortest_path_length
  D <- setNames(rep(NA_real_, length(G$nodes)), G$nodes); D[s] <- 0; Q <- s
  while (length(Q)) { v <- Q[1]; Q <- Q[-1]; for (w in G$adj[[v]]) if (is.na(D[[w]])) { D[w] <- D[[v]] + 1; Q <- c(Q, w) } }
  D[!is.na(D)]
}
nx_closeness <- function(G) vapply(G$nodes, function(u) { sp <- nx_distances(G, u); tot <- sum(sp); if (tot > 0 && length(G$nodes) > 1) (length(sp) - 1) / tot else 0 }, 1)   # wf_improved = False
nx_components <- function(G) { seen <- character(0); k <- 0L; for (s in G$nodes) if (!(s %in% seen)) { k <- k + 1L; seen <- c(seen, names(nx_distances(G, s))) }; k }
metrics <- function(pairs) {
  G <- nx_graph(pairs); deg <- setNames(vapply(G$nodes, function(x) length(G$adj[[x]]), 1L), G$nodes); bt <- nx_betweenness(G); cl <- nx_closeness(G)
  rows <- G$nodes[order(-deg[G$nodes], -bt[G$nodes], method = "radix")]                                  # sorted(G.nodes, key=(-deg, -bt)), stable
  list(G = G, rows = data.frame(gene = rows, degree = unname(deg[rows]), betweenness = round(unname(bt[rows]), 4), closeness = round(unname(cl[rows]), 4), stringsAsFactors = FALSE), ncomp = nx_components(G))
}
all_metrics <- list()
for (n in c(20, 17)) {
  o <- fromJSON(file.path(OUT, sprintf("pipeline_%dgenes.json", n)), simplifyVector = FALSE)
  for (key in c("consensus", "consensus_signcons")) {
    if (is.null(o[[key]])) next
    pairs <- lapply(o[[key]], function(p) unlist(p)); m <- metrics(pairs); all_metrics[[sprintf("%d|%s", n, key)]] <- m$rows
    cat(sprintf("== %d genes | %s: %d pairs, %d genes, components %d\n", n, key, length(pairs), length(m$G$nodes), m$ncomp))
    for (i in seq_len(nrow(m$rows))) cat(sprintf("   %-8s degree %d  betweenness %.4f  closeness %.4f\n", m$rows$gene[i], m$rows$degree[i], m$rows$betweenness[i], m$rows$closeness[i]))
    cat("   degree>=6 hubs:", paste0("[", paste(sprintf("'%s'", m$rows$gene[m$rows$degree >= 6]), collapse = ", "), "]"), "| degree>=5:", paste0("[", paste(sprintf("'%s'", m$rows$gene[m$rows$degree >= 5]), collapse = ", "), "]"), "\n")
  }
}
saveRDS(all_metrics, file.path(OUT, "network_metrics.rds"))
