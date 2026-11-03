# Gate B re-analysis, 2 Nov 2026

Python re-implementation of the pooled co-expression pipeline (there is no R on this machine). Input: `SBA/Figures/Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx` (12 observations x 21 genes; rows 0-3 CTRL, 4-7 DAPT, 8-11 DKK-1).

## Validation on the 20-gene set (reproduces every number in the manuscript)
| Quantity | Manuscript | This pipeline |
|---|---|---|
| Significant pairs: raw / Ins1 / Kcnj11 / Hes1 | 80 / 98 / 118 / 169 | 80 / 98 / 118 / 169 |
| Pooled consensus (significant in all four schemes, any sign) | 27 | 27 (21 sign-consistent, 6 Hey1 pairs flip sign) |
| Permutation null (2000), mean; P(null >= 27) | 29.8; 0.85 | 29.7; 0.849 |
| Treatment-adjusted q < 0.05: raw / Hes1 / Kcnj11 / Ins1 | 0 / 12 / 5 / 49 | 0 / 12 / 5 / 49 |
| Consensus pairs significant after adjustment in >= 1 scheme | 13 of 27 | 13 of 27 |
| Pairs significant in three adjusted schemes | Cacna1d-Wnt5a, Cacna1d-Ptbp1 | same two |
| Hub table (degree, betweenness, closeness) | `Fig_6b_Gene_Centrality_Summary.csv` | identical |

Note: the Methods sentence "same sign of r under all four conditions" does not describe the 27 (which includes the six sign-flipping Hey1 pairs). The 27 is "significant (p <= 0.05) in all four schemes".

## 17 genes (Hey1, Ryr1, Wnt2, Wnt2b excluded)
- Pairs tested: raw 136, each normalized scheme 120. Significant: raw 66, Ins1 80, Kcnj11 84, Hes1 120 of 120 (117 with |r| >= 0.8).
- Pooled consensus: 21 pairs, none sign-inconsistent. Permutation null mean 22.1 (SD 3.1), P(null >= 21) = 0.69: treatment structure alone still reproduces the count.
- Treatment-adjusted (df = 8, BH within scheme): q < 0.05 for raw 0, Ins1 51, Kcnj11 11, Hes1 15; none in all four schemes. 13 of 21 consensus pairs significant in at least one scheme; three in three schemes (Cacna1d-Ptbp1, Cacna1d-Wnt5a, Tcf7-Wnt9a). Smallest |r| for p < 0.05 at df = 8 is 0.632; the single best pair needs |r| >= 0.90 for q < 0.05.
- Pooled consensus network: 13 genes, 21 edges. Degree: Glp1r 6; Tcf7l2, Cacna1d, Ptbp1 5; Cacna1c, Wnt5a 4. **Only Glp1r reaches degree 6**, the hub rule used for the MLSS. With a degree >= 5 rule the hubs are Glp1r, Tcf7l2, Cacna1d, Ptbp1; the drug targets among them are still only CACNA1D and GLP1R, so MLSS values would be numerically unchanged (not re-run: R script).

## Findings that affect the text
1. **Replicate structure.** On the 3 Jun 2020 plates (Gr1, Gr2) each treatment has 2 wells per gene per plate, so the four values per group are 2 wells (Gr1) + 2 wells (Gr2). The RNA sheet lists two RNA samples per treatment (G1, G2). The 26 Aug and 9 Nov plates have 4 wells per group on one plate. Methods say "four wells from a single cell preparation".
2. **Gapdh dilution.** Gapdh Cq is about 15.6 on Gr1 and about 22.5 on Gr2 (a 6.7-cycle, ~100-fold difference), consistent with Gapdh at 1:100 on Gr2 only. Within-plate fold changes are not affected; cross-plate ΔCt values are.
3. **26 Aug plate, Gapdh in the DAPT group.** Gapdh Cq in wells A10 24.99, B10 24.13, C10 22.57, D10 22.80 (Control 23.6-23.8; DKK-1 22.5-22.9). The manuscript's DAPT fold changes for Kcnj11 (0.40), Cacna1d (0.76), Cacna1c (0.71), Ptbp1 (0.62) and Tcf7l2 (0.78) are reproduced only if Gapdh is averaged over C10 and D10 (0.39, 0.76, 0.64, 0.59, 0.73). Using all four wells gives 0.76, 1.46, 1.23, 1.13, 1.39. The handling of A10/B10 is not documented.
4. **Other genes** (Wnt5a/5b/9a, Itpr1, Ryr2/3 and the DKK-1 arm generally) agree between the manuscript matrix and the exports within ~10%.
5. **Efficiency sensitivity** (`efficiency_sensitivity_group_level.csv`): no direction flips. Hes1 DAPT 0.07 -> 0.17 and DKK-1 0.34 -> 0.49 with per-amplicon efficiency; Wnt5a DAPT 0.54 -> 0.68; all other changes < 1.25-fold.
6. **Primers** (`../primer_check/`): Hes1 and Hey1 as printed in Table S1 match canine accessions; the measured melt temperatures fit the rat `_NEW` Hes1 and Tcf7 amplicons and not the printed ones. See `primer_evidence.csv`.

## Not done here
- Figures 4-6, the MLSS tables and the STRING classification were not regenerated (they need R).
- The PAGE pre-test evidence for seven genes could not be added (no record found).

## Gapdh reference of the 26 August DAPT group: provenance, QC and exact sensitivity (3 Nov 2026; `gapdh_exact.py`, outputs `gapdh_exact.json`, `gapdh_exact_pairs.csv`)
Source: the authors' original calculation sheet `Undergrad_Glucose uptake signaling/RT-PCR/2020-08-26_... Quantification Summary.xlsx` (formulas: HK = Gapdh row, dCt = Cq - HK, relative expression = 2^-dCt). This replaces the constant-factor approximation of the first version of this section (deleted; it overstated the effect on Glp1r and the hubs).
- **What was done originally.** Raw Gapdh Cq of the DAPT group: A10 24.99, B10 24.13, C10 22.57, D10 22.80. In the sheet the four DAPT Gapdh rows hold 22.5685, 22.7983, 22.5685, 22.7983 (C10 and D10 entered twice; the A10/B10 values sit in a scratch column V). DAPT replicate k of every target gene is paired with row k. The matrix `Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx` equals the sheet's relative-expression column exactly (max difference 0.0) and the reported DAPT fold changes are reproduced exactly (Kcnj11 0.396, Cacna1d 0.763, Cacna1c 0.712, Ptbp1 0.616, Tcf7l2 0.778). The earlier "recalculated" columns of Table S1c used group means in Ct space (a different basis) and were replaced.
- **QC of A10/B10** (`gapdh_qc_amplification.py`, `gapdh_qc_replicates.py`; raw CFX exports): end-point fluorescence 513 and 725 RFU vs 1,279-1,507 RFU in the other ten Gapdh wells (no target well on the plate is below 700 RFU); melt-peak height 169 and 245 vs 453-544; melt Tm unchanged (82.5 C). Lower template would delay the Cq but not lower the plateau (CTRL, 1 cycle later, has a normal plateau), so these look like failed wells. Gapdh replicate SD: DAPT 1.14 cycles vs 0.10 (CTRL) and 0.15 (DKK-1) in the same run, and 0.03-0.43 in all groups of the other runs. No reason is recorded anywhere in the folder (calculation sheet has no comments; senior-project documents do not mention it).
- **Exact recalculation** (CTRL and DKK-1 untouched; reference applied to all four DAPT replicates): C0 = mean of the two valid wells (22.68); C1 = mean of all four wells (23.62).
  - Per gene (DAPT fold change, reported / C0 / C1): Tcf7l2 0.78 / 0.78 / 1.49; Kcnj11 0.40 / 0.40 / 0.76; Cacna1d 0.76 / 0.76 / 1.46; Cacna1c 0.71 / 0.72 / 1.37; Ptbp1 0.62 / 0.62 / 1.18. C0 leaves every significance category unchanged; C1 reverses four genes (Cacna1c p 0.054, Ptbp1 0.026) and removes the DAPT-vs-DKK-1 difference of Kcnj11.
  - Pooled network, reported / C0 / C1: significant pairs raw 66 / 63 / 62, Ins1 80 / 80 / 49, Kcnj11 84 / 86 / 44, Hes1 120 in all; consensus pairs 21 / 24 / 21 (C0 keeps 20, loses 1, adds 4; C1 keeps 11, loses 10, adds 10); adjusted q < 0.05 Kcnj11-normalized 11 / 36 / 12; hubs: Glp1r 6, Tcf7l2/Cacna1d/Ptbp1 5 (reported); C0 adds Wnt5a (5); C1 leaves Glp1r 5 and Cacna1d 5, Ptbp1 and Tcf7l2 4.
  - Cacna1d-Glp1r: the MLSS Synergy Bonus uses its Hes1-normalized r (0.960 reported, 0.976 C0, 0.974 C1), so the bonus would barely change; what changes is consensus status (C1: raw, Ins1- and Kcnj11-normalized r become negative and non-significant). The MLSS was not recomputed. (An earlier version of this line said the edge was 'lost' for the MLSS; that was wrong.)
- Not done: permutation null for C0/C1.
