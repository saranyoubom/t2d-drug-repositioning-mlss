# Fig. 5 corrected for the Glut2 / Slc2a2 naming error (8 Nov 2026)

**Decision.** The user chose "correct it" on 8 Nov 2026 (fix list C8). This folder holds the corrected figure for approval. Nothing in `R_17genes/`, `Figures/`, `Submission_package_JTM/` or the manuscripts was changed; everything is in scratch copies here.

## What was done
1. `as_reported/` and `alias/`: copies of `R_17genes/Fig_5abcd_Gene_Interaction_Network_Analysis.r` with the hard-coded `setwd()` removed (the original would overwrite `R_17genes/figures_publication`), run with R 4.4.0 (`Rscript`) in their own folders on copies of the three input files. `alias/` differs from `as_reported/` by one line (STRING node names Slc2a2 read as Glut2, marked `ALIAS FIX`). `patch_r_copies.py` adds an environment patch only: ragg as the implicit device (Rscript prints plots to its pdf device, which does not know the font "Arial" and stops with "invalid font type"), `windowsFonts(Arial = ...)`, and svg saves skipped (the composite uses the png panels).
2. **Gate 1 (data):** the six data exports of `as_reported/` (`Fig5B_Coverage_Summary`, `Fig5C_All_Pairs_Classified`, `Fig5D_Novel_Discoveries`, `Fig5_Network_Nodes`, `Fig5_Network_Edges`, `Fig5_Summary_Statistics`) are byte-identical to those of the published run in `R_17genes/figures_publication/data_export/` (5 / 16 / 29 / 70).
3. `compose_fig5.py` (the project's own `compose_figures_17genes.py` with the four panels replaced; the label patch of `patch_fig_labels.py` is NOT applied because the R script already sets the statistical symbol in italic and the raster patch would draw a second r at a position fitted to the old rendering). **Gate 2 (figure):** the `as_reported` composite differs from the published Figure 5 in 3.3% of the pixels, all inside the four panel slots (anti-aliasing of the new rendering) and in 0 pixels outside them.
4. Result: `composite_alias/Figure_5_final.png` (2000 x 2000 px, 300 dpi, same size as the published figure) and `Figure_5_final_600dpi.png`; `Fig5_published_vs_corrected.png` for review. Between the two composites: panel (a) 0.6% of the pixels (four edges change from grey to cyan), (b) 1.4% (pie), (c) 11.1% (histograms), (d) 0.0% (unchanged).

## What changes in the numbers (STRING 12.0 as used; `alias/figures_publication/data_export/Fig5_Summary_Statistics.csv`)
| | published | corrected |
|---|---|---|
| STRING Only | 29 (24.2%), mean |r| 0.951, median 0.958 | **33 (27.5%)**, mean 0.955, median 0.966 |
| Not analytically consistent | 70 (58.3%), mean 0.934, median 0.944 | **66 (55.0%)**, mean 0.932, median 0.943 |
| pairs in STRING | 34 of 120 (28.3%) | **38 of 120 (31.7%)** |
| Known Validated / Novel Discovery | 5 (4.2%) / 16 (13.3%), mean 0.968 | unchanged |
With STRING 12.5 (WP10 part A, alias): STRING Only 34, Not Validated 65, 39 of 120 (32.5%).

## Status 9 Nov 2026: APPROVED by the user at 00:26 +0700 (`APPROVED.txt`, with checksums). WP6b was released to apply item 15 of its brief; the items below are the list of what follows from the approval (manuscript by WP6b; Additional file 4, JTM package and `Figures/Figure_5.png` in the WP8 fix pass; WP10 README and Additional file 5 part 5.6 in WP6a round 5).
## To do after the user's approval of the figure (list; status in the line above)
* Manuscript v2 (WP6b brief item 15): replace the embedded Fig. 5 image by `Figure_5_final.png` (same pixel size, `embed_composites_in_manuscript.py` style replacement of `word/media/image5.png` in the v2 build); text numbers in the Results (the four-category sentence, the mean and median |r| sentence, the STRING coverage sentence) and in the legend of Fig. 5b; replace the Glut2 disclosure sentence by one Methods sentence on the name mapping; the sensitivity sentence uses 39 of 120 (32.5%) against 38 (31.7%).
* Additional file 4 (`Submission_package_JTM/Additional file 4.xlsx`, built by `Analysis_outputs/build_additional_file4_17genes.py`): the four sheets `Fig_5_All_Classified_Pairs`, `Fig_5_Summary_Statistics`, `Fig_5b_Coverage_Summary`, `Fig_5d_Top_unregistered_pairs` from the corrected exports (the fourth is unchanged in content).
* JTM package (`Submission_package_JTM/Figures/`) and `Figures/Figure_5.png`: the corrected figure (WP8, in a copy first).
* WP10 README and Additional file 5 part 5.6: state the decision and give the corrected numbers as the reported ones (WP6a round, after the figure is approved).
