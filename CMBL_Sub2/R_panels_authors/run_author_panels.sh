#!/bin/bash
# Re-draws the author panels whose omnibus-test annotation carried a plain "p" (Fig. 1c, 2b, 2c-f, 3b, 3c-f). Same data and statistics; only the label is changed to an italic p.
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
run() { echo "=== $* $(date +%H:%M:%S)"; "$RS" "$@" > "logs/$(echo $1 | tr -c 'A-Za-z0-9._\n' '_').log" 2>&1; echo "    exit $?"; }
run Fig_2b_glucose_uptake_revised.r
run Fig_2cdef_functional_genes_revised.r
run Fig_3ef_Boxplot_script.r
run R_panel_label_fixes.R . out_1c
run R_Fig3ab_cpeptide_descriptive.R Cpeptide_raw_elisa_tidy.csv out_3ab
echo "ALL DONE $(date +%H:%M:%S)"
