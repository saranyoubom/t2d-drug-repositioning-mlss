#!/bin/bash
# Re-runs only the scripts whose figure text was changed for the italics fixes (gene symbols and statistical symbols). Data outputs must stay byte-identical.
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
mkdir -p logs
for f in Fig_3f_PCA_Relative_mRNA_expression.R Fig_4b_script.r Fig_4c_script.r Fig_4d_script.r Fig_4e_sig_consensus.r Fig_5abcd_Gene_Interaction_Network_Analysis.r; do
  echo "=== $f $(date +%H:%M:%S)"
  "$RS" "$f" > "logs/italics_${f%.*}.log" 2>&1
  rc=$?; echo "    exit $rc"
  if [ $rc -ne 0 ]; then echo "FAILED: $f"; tail -15 "logs/italics_${f%.*}.log"; fi
done
echo "ALL DONE $(date +%H:%M:%S)"
