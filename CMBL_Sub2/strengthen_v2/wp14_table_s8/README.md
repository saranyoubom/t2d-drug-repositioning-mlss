# WP14: Table S8 recomputed from the Cq values of the analyzed matrix (10 Nov 2026)

**Strength of this evidence.** A descriptive recomputation that replaces a table built on the wrong data (finding F09 of the independent review, accepted by the user on 10 Nov 2026). No hypothesis, no decision rule; the reported analysis is unchanged. Independent check: `check_wp14_numbers.py` (`check_wp14_last_run.txt`).

## 1. What was wrong with the old Table S8
`reanalysis_17genes/effsens.py` (the old Table S8) took the Cq values from the 2025 re-exports of every plate on which a gene occurs and averaged the per-plate fold changes. For the nine June genes that meant the two 3 June plates, which supplied none of the analyzed values (WP11), and for the five 26 August genes it used all four Gapdh wells of the DAPT group including the failed wells A10 and B10 (WP13). Its group-level fold changes therefore did not equal the reported ones (for example Cacna1d under DAPT 1.46 in the old table against 0.76 reported).

## 2. What was done
Same method and flag logic as `effsens.py`, but on the Cq values that the reported analysis used (`../wp1_ins1_structure/ins1_sheet_cq.csv`: the target and reference Cq of every replicate row of the 17 analyzed genes). For each gene and treatment (DAPT, DKK-1) against CTRL: dT = mean target Cq (treated) minus mean target Cq (CTRL); dG = the same for the reference (Gapdh) rows of that gene; FC(E = 100%) = 2^-(dT - dG); FC(E amplicon) = (1 + E_gene)^-dT / (1 + E_Gapdh)^-dG with the mean per-amplicon efficiencies of Additional file 1 Table S1 (Gapdh 70.1%); ratio = FC(E amplicon) / FC(E = 100%). Flag: "direction flips" if the two fold changes lie on different sides of 1; ">1.5x shift" if |log2 ratio| > 0.585; ">1.25x shift" if > 0.32. The reported fold change (mean of the per-replicate 2^-dCt of the treated group divided by that of CTRL) is listed for comparison. Output: `results/efficiency_sensitivity_analyzed_matrix.csv` (17 genes).

## 3. Result
* **The group-level fold change at E = 100% reproduces the reported fold change** (for example Tcf7l2 DAPT 0.776 against 0.778; Kcnj11 DAPT 0.393 against 0.396; Hes1 DAPT 0.137 against 0.136); 33 of the 34 gene-by-treatment pairs lie on the same side of 1 (the exception, Lef1 under DKK-1, is 0.993 against 1.007).
* **Efficiency correction (mean amplicon efficiencies 53.9 to 74.6%) moves the fold changes by up to 2.0-fold** (Hes1 DAPT from 0.137 to 0.275; ratios 1.33 for Wnt5a DAPT, 1.27 for Glp1r DAPT, 1.26 for Kcnj11 DAPT) and **changes the direction of exactly one fold change, Tcf7l2 under DKK-1, which is 1.002 without and 0.943 with the correction (no change from CTRL in either)**. Five of 17 genes carry a flag (four with a shift above 1.25-fold, one direction flip).
* The DAPT decreases of Hes1, Kcnj11, Ptbp1, Itpr1, Ryr2, Ryr3, Cacna1c, Cacna1d and Tcf7l2 stay below 1 under the correction, and the DAPT increase of Wnt9a stays above 1 (table in `results/`); the comparison with every change that the manuscript describes is left to the text pass.

## 4. Text for the manuscript (WP8 text pass)
* **Limitations sentence (replaces "A group-level sensitivity analysis changed no fold change in direction (Additional file 1: Table S8)"):** "A group-level sensitivity analysis with the mean amplicon efficiencies, computed from the Cq values used in the analysis, shifted the fold changes by up to 2.0-fold (Hes1 under DAPT, 0.14 to 0.28) and changed the direction of one fold change that was not different from control (Tcf7l2 under DKK-1, 1.00 to 0.94; Additional file 1: Table S8)."
* **Table S8** (Additional file 1): title "Table S8. Group-level sensitivity of the fold changes to the amplification efficiency, computed from the Cq values of the analyzed matrix"; columns: gene; mean amplicon efficiency (%); for DAPT and for DKK-1 each: reported fold change; fold change at E = 100%; fold change with the amplicon efficiencies; ratio; last column: flag. Rows: the 17 analyzed genes in the order of Table S1. Note under the table: the method of section 2 in two sentences, the Gapdh efficiency (70.1%), and that fold changes are relative to CTRL.
* Wherever the Methods or the Limitations mention that the sensitivity analysis used the 3 June plates or "per-plate" fold changes, the wording goes.

## 5. Limits
Group-level (Ct-space) means, not the per-replicate analysis; post hoc efficiencies estimated from the 2025 exports (mean per amplicon, not per well); the Gapdh efficiency (70.1%) is applied to all genes' references; one preparation.

## 6. Files
`wp14_effsens_analyzed_matrix.py` (reads Additional file 1 Table S1 and `ins1_sheet_cq.csv`), `check_wp14_numbers.py`, `results/efficiency_sensitivity_analyzed_matrix.csv`. Use `python -I`.
