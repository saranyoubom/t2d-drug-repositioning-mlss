# Wording proposals for the plate provenance (WP11, 8 Nov 2026; for the WP8 fix pass, not yet applied to any manuscript)

Rules: American spelling, no thousands commas, no figure keys in images, every number below comes from `README.md` of this folder (checked by `check_wp11_numbers.py`), items the files cannot settle are marked **[AUTHOR]** and go in yellow in the manuscript.

## 1. Table S1c: plate and run of each gene of the matrix
| run | file (run start, UTC) | genes of the 17 analyzed | other matrix columns |
|---|---|---|---|
| 24 June 2020, plate 1 | `2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1.pcrd` (07:07) | Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7 | Wnt2 (Wnt2b has its own DAPT block) |
| 24 June 2020, plate 2 | `2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2.pcrd` (08:26) | Hes1, Glut2 (Slc2a2), Ins1, Glp1r | Hey1 |
| 26 August 2020 | `2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1.pcrd` (16:49) | Tcf7l2, Cacna1c, Cacna1d, Kcnj11, Ptbp1 | none |
| 9 November 2020 | `2020-11-09_ITPR, RYR1-3.pcrd` | Itpr1, Ryr2, Ryr3 | Ryr1 |
| 3 June 2020, Gr1 and Gr2 | `2020-Jun-03_Gr1_...pcrd` (05:56), `..._Gr2_...pcrd` (07:27) | none | none |
The present Table S1c lists Glp1r, Glut2, Hes1 and Ins1 under the 3 June plates; that is wrong (fix list A12).

## 2. Methods, RT-qPCR data handling (replaces the sentences about the plates and replicates)
"Reactions of the 24 June 2020 experiment were run on two plates; each gene was normalized to the Gapdh wells of its own plate. Cq values were read from the baseline-subtracted amplification curves with one manual threshold per plate (138.06 RFU for the second plate and about 210 RFU for the first, recovered by recomputing the Cq values from the exported curves; the sheet values reproduce to 0.04 cycles on average). **[AUTHOR: confirm the thresholds and the baseline setting of the first plate from the original file.]** The two plates of 3 June 2020 contributed no value to the analysis. For the 26 August 2020 run, the four Cq values of each gene and group were sorted in ascending order in the calculation workbook and the k-th smallest target value was paired with the k-th smallest Gapdh value, so that replicate k of this run denotes a rank and not a well; Gapdh wells A10 and B10 of the DAPT group were replaced by copies of wells C10 and D10 **[AUTHOR: reason]**. The 24 June and 9 November values are in the order of the calculation sheet, which is not sorted."

## 3. Limitations (one sentence, replaces any wording about paired replicates)
"Because target and Gapdh Cq values of the 26 August run were paired by rank, the replicate spread of the five genes of that run is smaller than a well-based pairing would give, and the correlations among those genes (mean +0.54, WP1) follow from this pairing and from their shared Gapdh vector (see the row-order permutation in the Methods); they are not read as co-regulation."
(The number +0.54 is the within-run correlation already reported in the manuscript from WP1; no new number is introduced.)

## 4. Optional, not recommended without a new asset
The comparison of the three pairings (README section 3) could become a short table in Additional file 5; it would need a WP6a round, because Additional file 5 has no part for it. The manuscript does not need it.

## 5. What the authors confirm before submission
1. The threshold and baseline setting of 24 June plate 1 (about 210 RFU, residual 0.04 cycles) and that one manual threshold per plate was the working practice.
2. The reason for setting aside Gapdh wells A10 and B10 (26 August, DAPT).
3. That the first plate shows an automatic threshold (45.70 RFU) in the 2025 re-export while the sheet corresponds to about 210 RFU.
