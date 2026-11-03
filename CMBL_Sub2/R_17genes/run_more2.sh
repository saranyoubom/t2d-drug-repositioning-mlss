#!/bin/bash
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
"$RS" -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if(!requireNamespace("ggfortify",quietly=TRUE)) install.packages("ggfortify"); cat("ggfortify:",requireNamespace("ggfortify",quietly=TRUE),"\n")' > logs/install_ggfortify.log 2>&1
for f in Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.R Fig_3f_PCA_Relative_mRNA_expression.R; do
  echo "=== $f $(date +%H:%M:%S)"
  "$RS" "$f" > "logs/${f%.*}.log" 2>&1; rc=$?; echo "    exit $rc"
  if [ $rc -ne 0 ]; then echo "FAILED: $f"; tail -8 "logs/${f%.*}.log"; fi
done
echo "ALL DONE $(date +%H:%M:%S)"
