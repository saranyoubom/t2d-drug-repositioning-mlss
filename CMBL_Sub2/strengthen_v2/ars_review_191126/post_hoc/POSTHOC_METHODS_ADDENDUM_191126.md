# Addendum to the post hoc method of 19 Nov 2026: the INS-1 matrix (written after the development results were seen)

**Status: post hoc, and added after the results of the 37 development series were known** (their summary is in `results_dev/summary.md`). The main method file (`POSTHOC_METHODS_191126.md`) did not list the INS-1 matrix. It is added because the manuscript quotes 97% for INS-1 (Abstract and Results) and a reader needs the value that independence would give for the same matrix. Nothing in the method changes; only the data set and the permutation count differ.

- Data: the 12 x 17 INS-1 matrix of the manuscript (3 groups x 4 wells; linear scale, relative expression 2^-dCt as in the manuscript), read from the file that `test_coexp_checks.py` reads (`SBA/Figures/Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx`, read only), gene order of the manuscript.
- Same code and the same two checks, 10000 permutations and seed 2026 (the setting of the manuscript's INS-1 coupling test: T_r2 0.133, null mean 0.114, P = 0.061).
- Gates before the result is read: observed share 0.971 (manuscript: 97.1%), observed T_r2 0.133 and its P within 0.004 of 0.061 (manuscript).
- Reading rule as in the main file: reported whichever direction; no variant is added; the number goes into the text only through `results_ins1/summary.md`.
