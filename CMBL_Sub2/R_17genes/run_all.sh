#!/bin/bash
# Runs the Fig. 4-6 R pipeline on the 17-gene matrix (R 4.4.0). Stops at the first failing script.
cd "$(dirname "$0")"
RS="/c/Program Files/R/R-4.4.0/bin/Rscript.exe"
mkdir -p logs
START="${1:-Fig_4b_script.r}"; started=0
for f in Fig_4b_script.r Fig_4c_script.r Fig_4d_script.r Fig_4e_script.r Fig_4e_sig_consensus.r \
         Fig_5abcd_Gene_Interaction_Network_Analysis.r Fig_6b_Consensus_Hub_Network.r Fig_6c_MLSS_v4.0.r \
         Fig_6d_S2_Validation_Top15.r Fig_6e_Drug_Combination_Network.R Fig_6f_Regulatory_Mechanistic_Sankey.r; do
  [ "$f" = "$START" ] && started=1
  [ $started -eq 0 ] && continue
  echo "=== $f $(date +%H:%M:%S)"
  "$RS" "$f" > "logs/${f%.*}.log" 2>&1
  rc=$?
  echo "    exit $rc"
  if [ $rc -ne 0 ]; then echo "FAILED: $f"; tail -15 "logs/${f%.*}.log"; exit $rc; fi
done
echo "ALL DONE $(date +%H:%M:%S)"
