#!/bin/bash
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
mkdir -p logs
for f in Fig_1d_Heatmap_Fold_change_mRNA_expression.R Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.R Fig_3f_PCA_Relative_mRNA_expression.R Fig_6b_Consensus_Hub_Network.r; do
  echo "=== $f $(date +%H:%M:%S)"
  "$RS" "$f" > "logs/${f%.*}.log" 2>&1; rc=$?; echo "    exit $rc"
  if [ $rc -ne 0 ]; then echo "FAILED: $f"; tail -12 "logs/${f%.*}.log"; fi
done
echo "ALL DONE $(date +%H:%M:%S)"
