#!/bin/bash
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
echo "=== 1c $(date +%H:%M:%S)"; "$RS" R_panel_label_fixes.R . out_1c > logs/R_panel_label_fixes.R.log 2>&1; echo "    exit $?"
for f in Fig_2b_glucose_uptake_revised.r Fig_2cdef_functional_genes_revised.r Fig_3ef_Boxplot_script.r; do
  echo "=== $f $(date +%H:%M:%S)"; "$RS" "$f" > "logs/${f}.log" 2>&1; echo "    exit $?"
done
echo "RERUN2 DONE $(date +%H:%M:%S)"
