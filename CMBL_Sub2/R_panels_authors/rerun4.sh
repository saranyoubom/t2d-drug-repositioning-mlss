#!/bin/bash
cd "$(dirname "$0")"; RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
echo "=== 1c $(date +%H:%M:%S)"; "$RS" R_panel_label_fixes.R . out_1c > logs/R_panel_label_fixes.R.log 2>&1; echo "    exit $?"
echo "=== 2b $(date +%H:%M:%S)"; "$RS" Fig_2b_glucose_uptake_revised.r > logs/Fig_2b_glucose_uptake_revised.r.log 2>&1; echo "    exit $?"
echo "=== 2cdef $(date +%H:%M:%S)"; "$RS" Fig_2cdef_functional_genes_revised.r > logs/Fig_2cdef_functional_genes_revised.r.log 2>&1; echo "    exit $?"
echo "RERUN4 DONE $(date +%H:%M:%S)"
