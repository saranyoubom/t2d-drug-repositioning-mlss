#!/bin/bash
# 2b (annotate route for Kruskal-Wallis) and 3c-f (needs patchwork >= 1.3.2 for ggplot2 4.0; installed in ./Rlib only, the global R library is untouched)
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
export R_LIBS="$(pwd)/Rlib"
echo "=== 2b $(date +%H:%M:%S)"; "$RS" Fig_2b_glucose_uptake_revised.r > logs/Fig_2b_glucose_uptake_revised.r.log 2>&1; echo "    exit $?"
echo "=== 3ef $(date +%H:%M:%S)"; "$RS" Fig_3ef_Boxplot_script.r > logs/Fig_3ef_Boxplot_script.r.log 2>&1; echo "    exit $?"
echo "RERUN3 DONE $(date +%H:%M:%S)"
