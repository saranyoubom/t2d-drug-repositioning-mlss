# O6 disclosure: draft wording (27 Nov 2026)

Decision of the user, 27 Nov 2026: option 1, keep the as-run (pre-specified) calibration result as primary (Tr2 type I error 7.5%, FAIL) and disclose the exact-rule result (4.8%, PASS) as a post hoc correction found after the results were known.

Status: DRAFT, not applied. Nothing in the package or the manuscript source has been edited. The wording below waits for the user's approval; it will be applied together with the Python-to-R software wording and the package rebuild.

Numbers are read from `results_R/o6_anova_exact/wp2_exact/calibration_verdict.txt`, `o6_calibration_summary.md` and `o6_summary.md`. "Tr2", "Tp" keep the manuscript's subscript formatting; "10-29" is written with a Word superscript, not with Unicode superscript characters.

## 1. Manuscript, Results, body paragraph [137] (replaces the text from "The excess came from" to the end)

Current:

> The excess came from one series (GSE20677) in which 13 of the 17 genes of gene set A had identical values in all 12 samples; the analysis of variance used to select genes flags such genes as differentially expressed because their within-group sum of squares is exactly zero. Without the two series with constant genes the type I error was 5.0% (175 of 3500) and the power 91% (post hoc). Tp (type I error 5.8%, power 91%) met the requirement and is the secondary statistic named by the fixed rule.

Proposed (the C7 correction of the audit plus the disclosure):

> The excess came from one series (GSE20677) in which 13 of the 17 genes of gene set A had identical values in all 12 samples; the analysis of variance used to select genes flagged such genes as differentially expressed, an artifact of floating-point arithmetic: their within-group sum of squares is exactly zero, and rounding of the grand mean leaves a between-group sum of squares of about 4 × 10^-29, so F is infinite. Without the two series with constant genes the type I error was 5.0% (175 of 3500) and the power 91% (post hoc). When the analysis was reimplemented in R, after all results were known, we found that the same artifact also decided which series were analyzable. With genes that are constant in all 12 samples given P = 1 (post hoc), two development series dropped out, one entered and gene set A changed in three, and Tr2 met the requirement (type I error 4.8%, 172 of 3600; power 91%, 656 of 720; Additional file 5, part 4.1); the Tr2 results in the validation sets did not change. The pre-specified result remains the primary one. Tp (type I error 5.8%, power 91%) met the requirement and is the secondary statistic named by the fixed rule.

## 2. Manuscript, Methods, body paragraph [46] (one sentence added after "...without within-group variation (Results).")

> A rounding artifact in this gene selection, found later, was examined in a second post hoc analysis in which genes with identical values in all 12 samples were given P = 1 (Results; Additional file 5, part 4.1).

## 3. Manuscript, Discussion, body paragraph [166] (parenthesis extended)

Current: "(type I error 7.5%, above the 6% limit set in advance)"

Proposed: "(type I error 7.5%, above the 6% limit set in advance; 4.8% after a post hoc correction of a rounding artifact in the gene selection)"

## 4. Cover letter, paragraph [10] (parenthesis added)

Current: "A within-group coupling test, whose type I error of 7.5% is above the 6% limit set in advance in our protocol (its SHA-256 fingerprints ...)"

Proposed: "A within-group coupling test, whose type I error of 7.5% is above the 6% limit set in advance in our protocol (4.8% after a post hoc correction of a rounding artifact in the gene selection; its SHA-256 fingerprints ...)"

## 5. Fig. 9 legend [141] and Results [138]: unchanged

Fig. 9 and paragraph [138] show the pre-specified analysis, and the legend already points to the Results for the type I error.

## 6. Additional file 5, part 4.1: new block after the post hoc sensitivity block (after document paragraph [541])

> POST HOC CORRECTION (found on 24 Nov 2026 when the analysis was reimplemented in R, after all results were known; not the pre-specified analysis): genes with identical values in all 12 samples given P = 1 in the analysis of variance of part 1, everything else as pre-specified; development set of 36 series (16 on GPL570, 20 on GPL1261): GSE20677 and GSE49079 have no differentially expressed gene and are no longer analyzable, GSE30164 enters, gene set A is redrawn in GSE15892, GSE28477 and GSE53183
>
>   T_r2 type I error: 172/3600 = 0.0478 (Wilson 95% CI 0.0413 to 0.0552) -> PASS
>
>   T_r2 power at r = 0.7, m = 6: 656/720 = 0.9111 (Wilson 95% CI 0.8881 to 0.9298) -> PASS; per series: min 0.75, median 0.90, 35 of 36 series at or above 0.80
>
>   OVERALL: PASS
>
>   All pre-listed statistics (type I error, power at r = 0.7, m = 6): T_r2 0.048, 0.911; T_p 0.031, 0.918; T_eig 0.045, 0.908; count 0.043, 0.275

## 7. Additional file 5, part 8 Deviations: new entry after the entry of 20 Nov 2026 (document paragraph [72])

> 27 Nov 2026 (POST HOC, found when the analysis was reimplemented in R, after all results were known; no pre-specified rule, threshold or definition was changed and the pre-specified result stays primary). (1) The rounding artifact diagnosed for GSE20677 (entry of 5 Nov 2026, item 2) also sets the gene selection of other series. A probe with identical values in all 12 samples gets P = 0 because the computed grand mean differs from the exact value by about one unit in the last place, and these P values make the Benjamini-Hochberg adjustment lenient: GSE49079 had 262 differentially expressed genes as run and none in exact arithmetic. (2) Corrected analysis: such probes are given P = 1, everything else is unchanged. Development set: 36 series (GSE20677 and GSE49079 out, GSE30164 in; gene set A redrawn in GSE15892, GSE28477 and GSE53183; gene set B unchanged). Calibration: part 4.1. Development application: T_r2 flagged 22 of 36 series (61%) for gene set A (25 of 37 as run). Validation: only the gene set A of GSE68424 changes; T_r2 flagged 47 of 59 sets for gene set A and 48 of 59 for gene set B, as run; the adjusted count of the third ratio scheme was above its null in 15 of 59 sets (14 as run); the medians of Fig. 9d would be 1.03, 1.11 and 1.10 over 95 data sets (1.02, 1.12 and 1.11 over 96 as run). (3) The two post hoc checks of 20 Nov 2026 (removal of the first principal component and the expected share of the group means) were not repeated with the corrected selection.

## Open points for the user

- Item 7 (3): the post hoc PC1 removal and expected share were not recomputed under the exact rule. Option 1 does not need them, so the entry says so.
- Additional file 6: no new sheet proposed; the exact-rule tables stay in the repository (`results_R/o6_anova_exact/`).
- The sentence "When the analysis was reimplemented in R" presumes the software wording of the manuscript changes from Python to R (pending).
