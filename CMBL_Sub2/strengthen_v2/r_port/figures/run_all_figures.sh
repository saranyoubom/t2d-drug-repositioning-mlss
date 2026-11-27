#!/bin/sh
# run_all_figures.sh - draws every figure of the study in R (PNG 300 dpi and PDF) from the saved results: Figs 7, 8, 9, 10, A5.1, S1 and the route-3 benchmark figure.
# Usage: sh run_all_figures.sh [Rscript]   (default: Rscript on the PATH). Inputs default to the archived result folders; set WP4_RES, WP4_META, WP4_WS, WP4_PC, WP2_RES, POSTHOC_RES, RA_RES, ROUTE3_RES to use the results of an R run
# (results_R/...); FIG_OUT is the output folder (default results_R/figures). The figures were compared with the archived matplotlib figures (see R_PORT_PLAN.md, step 4).
RS="${1:-Rscript}"; cd "$(dirname "$0")" || exit 1
for f in make_fig7.R make_fig8.R make_fig9.R make_fig10.R make_figA51.R make_figS1.R make_fig_benchmark.R; do
  echo "== $f"; "$RS" "$f" || { echo "FAILED: $f"; exit 1; }
done
echo "all figures done"
