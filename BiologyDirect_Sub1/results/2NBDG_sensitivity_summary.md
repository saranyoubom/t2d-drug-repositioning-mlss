# 2-NBDG sensitivity analysis (for the Supplement)

Source: `2021-04-07_Gluclose uptake assay.xlsx` and `_Original.xlsx` (assay date 7 Apr 2021; both files last saved 21 Apr 2026). Script: `R_2NBDG_sensitivity.R`; table: `2NBDG_sensitivity_table.csv`; console output: `2NBDG_sensitivity_output.txt`; input values with cell addresses: `input_data/nbdg_blocks.csv`.

## What the source workbook contains
- Two read-outs ("1st", "2nd") of the same wells at two seeding densities (25K and 50K cells per well), each with 4 wells per group and 4 negative-control (blank) wells.
- The published figure uses the 25K, 2nd read, with CTRL = 195, 331, 331, 3382 (analysis block, rows 20-23), the highest CTRL well (3382) excluded by a Grubbs test on CTRL only, and fold change = intensity / mean of the remaining three CTRL values. Blank wells were not subtracted.
- CTRL replicate 2 reads 184 in the reading block (row 12) and 331 in every analysis block (rows 21-22), where it duplicates replicate 3. No comment, highlight or note explains the difference; the workbook was created 7 Apr 2021 and last saved 21 Apr 2026.

## Result
- Scenario 1 reproduces the published values exactly (CTRL 1.000 ± 0.159, DAPT 0.946 ± 0.198, DKK-1 1.288 ± 0.388; Kruskal-Wallis p = 0.839; all Bonferroni-adjusted p = 1.0).
- Across 12 scenarios (CTRL replicate 2 = 184 or 331; high CTRL well excluded or kept; blank not subtracted, subtracted by row or by mean; 1st or 2nd read; 25K or 50K cells), no comparison was significant: the smallest Bonferroni-adjusted p was 0.733 (log10-scale ANOVA with t-tests) and omnibus p ranged from 0.48 to 0.87.
- The direction and size of the differences are not stable. Mean fold change relative to CTRL ranged from 0.19 to 1.25 for DAPT and from 0.29 to 1.92 for DKK-1 (medians 0.76 to 1.30 and 0.96 to 1.73). With the high CTRL well kept, both treated groups appear 65-80% lower than CTRL; with it excluded, they appear equal or higher. Using 184 instead of 331 raises the published estimates from 0.95 to 1.14 (DAPT) and from 1.29 to 1.56 (DKK-1).
- Blank wells (74-122 AU) were of the same order as treated-well signals (157-681 AU in the 25K block).

## Reading
The assay detected no significant difference under any analysis choice, but it cannot support a statement about the direction of any effect on glucose uptake. The earlier statement that DAPT and DKK-1 impair glucose uptake (paired percent-of-CTRL analysis in the 2022 senior-project report and the earlier manuscript) rested on one of these analysis variants and is not supported.
