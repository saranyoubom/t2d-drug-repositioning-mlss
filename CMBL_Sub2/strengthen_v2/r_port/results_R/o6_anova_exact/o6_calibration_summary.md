# Impact check O6, part 2: WP2 numbers with the exact ANOVA rule

Written 2026-11-27 00:58 by o6_anova_exact/compare_o6_calib.R. Exact rule: development set of 36 series (GSE20677 and GSE49079 out, GSE30164 in; set A redrawn in GSE15892, GSE28477, GSE53183); held-out GSE68424 set A redrawn; everything else reused bit for bit (prep_exact_calib.log, merge_heldout.log).

## 1. Calibration verdict (protocol 3.2)

### R run, numpy order (37 series; equals the archived Python verdict except T_p and T_eig, BLAS noise)

```
PROTOCOL 3.2 REQUIREMENT (T_r2 type I error <= 0.06 and power >= 80% at r = 0.7, m = 6; pooled over 37 series)

T_r2 type I error: 276/3700 = 0.0746 (Wilson 95% CI 0.0666 to 0.0835) -> FAIL
T_r2 power at r = 0.7, m = 6: 671/740 = 0.9068 (Wilson 95% CI 0.8837 to 0.9257) -> PASS; per series: min 0.75, median 0.90, 35 of 37 series at or above 0.80
OVERALL: FAIL

All pre-listed statistics (type I error, power at r = 0.7, m = 6): T_r2 0.075, 0.907; T_p 0.058, 0.914; T_eig 0.074, 0.909; count 0.041, 0.282
Secondary statistic named by the fixed rule (highest power at r = 0.7, m = 6 among statistics with type I error <= 0.06): T_p

POST HOC SENSITIVITY (not the pre-specified analysis): without the 2 series whose gene set A has genes with no within-group variation (GSE20677, GSE28477); 35 series remain
  T_r2   type I error 175/3500 = 0.0500 (0.0433 to 0.0577); power at r = 0.7, m = 6: 636/700 = 0.9086 (0.8849 to 0.9278)
  T_p    type I error 115/3500 = 0.0329 (0.0274 to 0.0393); power at r = 0.7, m = 6: 638/700 = 0.9114 (0.8881 to 0.9303)
  T_eig  type I error 171/3500 = 0.0489 (0.0422 to 0.0565); power at r = 0.7, m = 6: 635/700 = 0.9071 (0.8834 to 0.9265)
  count  type I error 30/700 = 0.0429 (0.0302 to 0.0605); power at r = 0.7, m = 6: 205/700 = 0.2929 (0.2603 to 0.3276)
```

### Exact rule (36 series)

```
PROTOCOL 3.2 REQUIREMENT (T_r2 type I error <= 0.06 and power >= 80% at r = 0.7, m = 6; pooled over 36 series)

T_r2 type I error: 172/3600 = 0.0478 (Wilson 95% CI 0.0413 to 0.0552) -> PASS
T_r2 power at r = 0.7, m = 6: 656/720 = 0.9111 (Wilson 95% CI 0.8881 to 0.9298) -> PASS; per series: min 0.75, median 0.90, 35 of 36 series at or above 0.80
OVERALL: PASS

All pre-listed statistics (type I error, power at r = 0.7, m = 6): T_r2 0.048, 0.911; T_p 0.031, 0.918; T_eig 0.045, 0.908; count 0.043, 0.275
Pre-listed statistics meeting both requirements: T_r2, T_p, T_eig

POST HOC SENSITIVITY (not the pre-specified analysis): without the 1 series whose gene set A has genes with no within-group variation (GSE28477); 35 series remain
  T_r2   type I error 171/3500 = 0.0489 (0.0422 to 0.0565); power at r = 0.7, m = 6: 637/700 = 0.9100 (0.8865 to 0.9290)
  T_p    type I error 112/3500 = 0.0320 (0.0267 to 0.0384); power at r = 0.7, m = 6: 643/700 = 0.9186 (0.8960 to 0.9366)
  T_eig  type I error 161/3500 = 0.0460 (0.0395 to 0.0535); power at r = 0.7, m = 6: 637/700 = 0.9100 (0.8865 to 0.9290)
  count  type I error 31/700 = 0.0443 (0.0314 to 0.0622); power at r = 0.7, m = 6: 197/700 = 0.2814 (0.2494 to 0.3159)
```

## 2. T_r2 per series: the changed series and every series with a type I error above 0.06

| series | platform | type I error, numpy | type I error, exact | power r = 0.7, m = 6, numpy | exact |
|---|---|---|---|---|---|
| GSE20677 | GPL570 | 100/100 | - | 20/20 | - |
| GSE49079 | GPL1261 | 4/100 | - | 18/20 | - |
| GSE15892 | GPL1261 | 7/100 | 5/100 | 17/20 | 19/20 |
| GSE28477 | GPL1261 | 1/100 | 1/100 | 15/20 | 19/20 |
| GSE53183 | GPL570 | 4/100 | 5/100 | 19/20 | 18/20 |
| GSE30164 | GPL1261 | - | 1/100 | - | 18/20 |
| GSE13274 | GPL570 | 7/100 | 7/100 | 20/20 | 20/20 |
| GSE26656 | GPL570 | 8/100 | 8/100 | 15/20 | 15/20 |
| GSE36176 | GPL570 | 8/100 | 8/100 | 18/20 | 18/20 |
| GSE39820 | GPL1261 | 7/100 | 7/100 | 19/20 | 19/20 |
| GSE53480 | GPL1261 | 7/100 | 7/100 | 17/20 | 17/20 |
| GSE58093 | GPL570 | 9/100 | 9/100 | 18/20 | 18/20 |
| GSE7890 | GPL570 | 7/100 | 7/100 | 16/20 | 16/20 |

## 3. Pooled calibration table (all conditions)

| r | m | statistic | numpy | exact |
|---|---|---|---|---|
| 0.0 | 0 | count | 30/740 = 0.041 | 31/720 = 0.043 |
| 0.0 | 0 | T_eig | 272/3700 = 0.074 | 162/3600 = 0.045 |
| 0.0 | 0 | T_p | 216/3700 = 0.058 | 113/3600 = 0.031 |
| 0.0 | 0 | T_r2 | 276/3700 = 0.075 | 172/3600 = 0.048 |
| 0.3 | 3 | T_eig | 63/740 = 0.085 | 39/720 = 0.054 |
| 0.3 | 3 | T_p | 51/740 = 0.069 | 33/720 = 0.046 |
| 0.3 | 3 | T_r2 | 67/740 = 0.091 | 47/720 = 0.065 |
| 0.5 | 3 | T_eig | 81/740 = 0.109 | 60/720 = 0.083 |
| 0.5 | 3 | T_p | 68/740 = 0.092 | 46/720 = 0.064 |
| 0.5 | 3 | T_r2 | 89/740 = 0.120 | 68/720 = 0.094 |
| 0.7 | 3 | T_eig | 140/740 = 0.189 | 117/720 = 0.163 |
| 0.7 | 3 | T_p | 134/740 = 0.181 | 108/720 = 0.150 |
| 0.7 | 3 | T_r2 | 171/740 = 0.231 | 140/720 = 0.194 |
| 0.9 | 3 | T_eig | 222/740 = 0.300 | 204/720 = 0.283 |
| 0.9 | 3 | T_p | 157/740 = 0.212 | 140/720 = 0.194 |
| 0.9 | 3 | T_r2 | 257/740 = 0.347 | 244/720 = 0.339 |
| 0.3 | 6 | count | 51/740 = 0.069 | 48/720 = 0.067 |
| 0.3 | 6 | T_eig | 98/740 = 0.132 | 84/720 = 0.117 |
| 0.3 | 6 | T_p | 75/740 = 0.101 | 56/720 = 0.078 |
| 0.3 | 6 | T_r2 | 98/740 = 0.132 | 83/720 = 0.115 |
| 0.5 | 6 | count | 99/740 = 0.134 | 89/720 = 0.124 |
| 0.5 | 6 | T_eig | 326/740 = 0.441 | 300/720 = 0.417 |
| 0.5 | 6 | T_p | 223/740 = 0.301 | 197/720 = 0.274 |
| 0.5 | 6 | T_r2 | 316/740 = 0.427 | 292/720 = 0.406 |
| 0.7 | 6 | count | 209/740 = 0.282 | 198/720 = 0.275 |
| 0.7 | 6 | T_eig | 673/740 = 0.909 | 654/720 = 0.908 |
| 0.7 | 6 | T_p | 676/740 = 0.914 | 661/720 = 0.918 |
| 0.7 | 6 | T_r2 | 671/740 = 0.907 | 656/720 = 0.911 |
| 0.9 | 6 | count | 500/740 = 0.676 | 478/720 = 0.664 |
| 0.9 | 6 | T_eig | 739/740 = 0.999 | 719/720 = 0.999 |
| 0.9 | 6 | T_p | 739/740 = 0.999 | 718/720 = 0.997 |
| 0.9 | 6 | T_r2 | 740/740 = 1.000 | 720/720 = 1.000 |
| 0.3 | 10 | T_eig | 226/740 = 0.305 | 200/720 = 0.278 |
| 0.3 | 10 | T_p | 147/740 = 0.199 | 129/720 = 0.179 |
| 0.3 | 10 | T_r2 | 223/740 = 0.301 | 198/720 = 0.275 |
| 0.5 | 10 | T_eig | 711/740 = 0.961 | 690/720 = 0.958 |
| 0.5 | 10 | T_p | 595/740 = 0.804 | 570/720 = 0.792 |
| 0.5 | 10 | T_r2 | 703/740 = 0.950 | 683/720 = 0.949 |
| 0.7 | 10 | T_eig | 740/740 = 1.000 | 720/720 = 1.000 |
| 0.7 | 10 | T_p | 740/740 = 1.000 | 720/720 = 1.000 |
| 0.7 | 10 | T_r2 | 740/740 = 1.000 | 720/720 = 1.000 |
| 0.9 | 10 | T_eig | 740/740 = 1.000 | 720/720 = 1.000 |
| 0.9 | 10 | T_p | 740/740 = 1.000 | 720/720 = 1.000 |
| 0.9 | 10 | T_r2 | 740/740 = 1.000 | 720/720 = 1.000 |

## 4. Development application (dev_application_summary.md)

33 of 59 lines differ (left: numpy order, R run inputs; right: exact rule).

- line 7
  - numpy: | Within-group coupling, T_r2 (primary), P < 0.05 | 25/37 = 68% (51-80%) |
  - exact: | Within-group coupling, T_r2 (primary), P < 0.05 | 22/36 = 61% (45-75%) |
- line 8
  - numpy: | T_p, P < 0.05 | 21/37 = 57% (41-71%) |
  - exact: | T_p, P < 0.05 | 20/36 = 56% (40-70%) |
- line 9
  - numpy: | T_eig, P < 0.05 | 21/37 = 57% (41-71%) |
  - exact: | T_eig, P < 0.05 | 19/36 = 53% (37-68%) |
- line 10
  - numpy: | Consensus count not above its permutation null (E1; P >= 0.05), this run | 32/37 = 86% (72-94%) |
  - exact: | Consensus count not above its permutation null (E1; P >= 0.05), this run | 30/36 = 83% (68-92%) |
- line 12
  - numpy: | E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 37) | median 0% [0%, 6%] (n=37) |
  - exact: | E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 36) | median 0% [0%, 3%] (n=36) |
- line 13
  - numpy: | E2 any scheme (part 1 definition), this run / part 1 | 60% [33%, 90%] (n=37) / 60% [34%, 90%] (n=37) |
  - exact: | E2 any scheme (part 1 definition), this run / part 1 | 72% [32%, 93%] (n=36) / 60% [34%, 90%] (n=37) |
- line 14
  - numpy: | E2 all four schemes, this run / part 1 | 0% [0%, 2%] (n=37) / 0% [0%, 2%] (n=37) |
  - exact: | E2 all four schemes, this run / part 1 | 0% [0%, 0%] (n=36) / 0% [0%, 2%] (n=37) |
- line 15
  - numpy: | Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.0] (n=37) / 0.0 [0.0, 1.0] (n=37) |
  - exact: | Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.0] (n=36) / 0.0 [0.0, 1.0] (n=37) |
- line 16
  - numpy: | E3 share of pooled covariance carried by group means | 98% [96%, 99%] (n=37) |
  - exact: | E3 share of pooled covariance carried by group means | 98% [96%, 99%] (n=36) |
- line 17
  - numpy: | Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.81 [0.62, 0.91] (n=35); above 0.61 in 27/35 = 77% (61-88%) |
  - exact: | Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.81 [0.65, 0.91] (n=36); above 0.61 in 30/36 = 83% (68-92%) |
- line 23
  - numpy: | raw | 9/37 = 24% (13-40%) | 0.0 [0.0, 2.0] (n=37) | 0.1 [0.1, 0.2] (n=37) | 0.00 [0.00, 13.16] (n=37) (defined in 37) | -0.0 [-0.1, 1.9] (n=37) |
  - exact: | raw | 9/36 = 25% (14-41%) | 0.0 [0.0, 2.0] (n=36) | 0.1 [0.1, 0.2] (n=36) | 0.00 [0.00, 13.20] (n=36) (defined in 36) | -0.1 [-0.1, 1.7] (n=36) |
- line 24
  - numpy: | ratio to gene 1 | 8/37 = 22% (11-37%) | 6.0 [2.0, 49.0] (n=37) | 5.0 [0.7, 27.3] (n=37) | 1.01 [0.80, 2.51] (n=37) (defined in 37) | 0.6 [-0.7, 3.5] (n=37) |
  - exact: | ratio to gene 1 | 8/36 = 22% (12-38%) | 12.0 [1.8, 52.0] (n=36) | 4.3 [0.7, 46.8] (n=36) | 1.02 [0.71, 2.99] (n=36) (defined in 36) | 0.7 [-0.7, 4.0] (n=36) |
- line 25
  - numpy: | ratio to gene 2 | 11/37 = 30% (17-46%) | 12.0 [3.0, 38.0] (n=37) | 8.2 [2.5, 23.9] (n=37) | 0.99 [0.39, 2.49] (n=37) (defined in 37) | -0.2 [-1.9, 7.2] (n=37) |
  - exact: | ratio to gene 2 | 9/36 = 25% (14-41%) | 13.5 [3.0, 39.2] (n=36) | 9.6 [1.4, 26.2] (n=36) | 0.99 [0.50, 2.33] (n=36) (defined in 36) | -0.2 [-1.8, 5.6] (n=36) |
- line 26
  - numpy: | ratio to gene 3 | 4/37 = 11% (4-25%) | 14.0 [4.0, 51.0] (n=37) | 13.9 [2.2, 49.7] (n=37) | 1.03 [0.45, 1.88] (n=37) (defined in 37) | 0.6 [-6.8, 6.8] (n=37) |
  - exact: | ratio to gene 3 | 3/36 = 8% (3-22%) | 16.5 [4.8, 57.5] (n=36) | 25.5 [3.0, 54.5] (n=36) | 1.01 [0.43, 1.54] (n=36) (defined in 36) | 0.4 [-5.7, 5.4] (n=36) |
- line 28
  - numpy: By platform, human GPL570: T_r2 P < 0.05 in 12/17 = 71% (47-87%); consensus count not above null in 14/17 = 82% (59-94%).
  - exact: By platform, human GPL570: T_r2 P < 0.05 in 11/16 = 69% (44-86%); consensus count not above null in 13/16 = 81% (57-93%).
- line 29
  - numpy: By platform, mouse GPL1261: T_r2 P < 0.05 in 13/20 = 65% (43-82%); consensus count not above null in 18/20 = 90% (70-97%).
  - exact: By platform, mouse GPL1261: T_r2 P < 0.05 in 11/20 = 55% (34-74%); consensus count not above null in 17/20 = 85% (64-95%).
- line 30
  - numpy: Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 23/35 = 66% (49-79%); T_p in 20/35 = 57% (41-72%); consensus count not above null in 30/35 = 86% (71-94%); E2 raw scheme alone 0% [0%, 5%] (n=35); covariance share carried by group means 98% [96%, 99%] (n=35).
  - exact: Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 22/35 = 63% (46-77%); T_p in 20/35 = 57% (41-72%); consensus count not above null in 29/35 = 83% (67-92%); E2 raw scheme alone 0% [0%, 3%] (n=35); covariance share carried by group means 98% [96%, 99%] (n=35).
- line 36
  - numpy: | Within-group coupling, T_r2 (primary), P < 0.05 | 29/37 = 78% (63-89%) |
  - exact: | Within-group coupling, T_r2 (primary), P < 0.05 | 28/36 = 78% (62-88%) |
- line 37
  - numpy: | T_p, P < 0.05 | 28/37 = 76% (60-87%) |
  - exact: | T_p, P < 0.05 | 27/36 = 75% (59-86%) |
- line 38
  - numpy: | T_eig, P < 0.05 | 27/37 = 73% (57-85%) |
  - exact: | T_eig, P < 0.05 | 26/36 = 72% (56-84%) |
- line 39
  - numpy: | Consensus count not above its permutation null (E1; P >= 0.05), this run | 24/37 = 65% (49-78%) |
  - exact: | Consensus count not above its permutation null (E1; P >= 0.05), this run | 24/36 = 67% (50-80%) |
- line 41
  - numpy: | E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 37) | median 0% [0%, 20%] (n=37) |
  - exact: | E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 36) | median 0% [0%, 21%] (n=36) |
- line 42
  - numpy: | E2 any scheme (part 1 definition), this run / part 1 | 100% [67%, 100%] (n=37) / 100% [67%, 100%] (n=37) |
  - exact: | E2 any scheme (part 1 definition), this run / part 1 | 98% [65%, 100%] (n=36) / 100% [67%, 100%] (n=37) |
- line 43
  - numpy: | E2 all four schemes, this run / part 1 | 0% [0%, 0%] (n=37) / 0% [0%, 0%] (n=37) |
  - exact: | E2 all four schemes, this run / part 1 | 0% [0%, 2%] (n=36) / 0% [0%, 0%] (n=37) |
- line 44
  - numpy: | Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.0] (n=37) / 0.0 [0.0, 1.0] (n=37) |
  - exact: | Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.5] (n=36) / 0.0 [0.0, 1.0] (n=37) |
- line 45
  - numpy: | E3 share of pooled covariance carried by group means | 69% [42%, 86%] (n=37) |
  - exact: | E3 share of pooled covariance carried by group means | 69% [48%, 87%] (n=36) |
- line 46
  - numpy: | Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.70 [0.40, 0.94] (n=27); above 0.61 in 15/27 = 56% (37-72%) |
  - exact: | Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.69 [0.37, 0.91] (n=28); above 0.61 in 15/28 = 54% (36-70%) |
- line 52
  - numpy: | raw | 14/37 = 38% (24-54%) | 1.0 [0.0, 6.0] (n=37) | 0.1 [0.1, 0.2] (n=37) | 4.64 [0.00, 31.25] (n=37) (defined in 37) | 0.8 [-0.1, 5.9] (n=37) |
  - exact: | raw | 13/36 = 36% (22-52%) | 1.0 [0.0, 3.8] (n=36) | 0.1 [0.1, 0.2] (n=36) | 4.24 [0.00, 29.87] (n=36) (defined in 36) | 0.8 [-0.1, 3.7] (n=36) |
- line 53
  - numpy: | ratio to gene 1 | 12/37 = 32% (20-49%) | 19.0 [2.0, 68.0] (n=37) | 7.8 [1.2, 50.8] (n=37) | 1.61 [0.96, 5.60] (n=37) (defined in 37) | 2.5 [-0.1, 12.0] (n=37) |
  - exact: | ratio to gene 1 | 11/36 = 31% (18-47%) | 17.5 [2.0, 67.2] (n=36) | 5.9 [1.1, 45.2] (n=36) | 1.53 [0.95, 3.98] (n=36) (defined in 36) | 2.1 [-0.3, 10.2] (n=36) |
- line 54
  - numpy: | ratio to gene 2 | 8/37 = 22% (11-37%) | 29.0 [4.0, 67.0] (n=37) | 10.0 [2.5, 52.7] (n=37) | 1.08 [0.98, 1.81] (n=37) (defined in 37) | 1.3 [-0.1, 7.3] (n=37) |
  - exact: | ratio to gene 2 | 7/36 = 19% (10-35%) | 14.0 [3.8, 67.0] (n=36) | 10.3 [2.4, 53.3] (n=36) | 1.05 [0.98, 1.76] (n=36) (defined in 36) | 1.3 [-0.1, 7.1] (n=36) |
- line 55
  - numpy: | ratio to gene 3 | 7/37 = 19% (9-34%) | 20.0 [4.0, 49.0] (n=37) | 7.0 [1.8, 22.7] (n=37) | 1.23 [0.84, 2.83] (n=37) (defined in 37) | 2.3 [-0.3, 12.5] (n=37) |
  - exact: | ratio to gene 3 | 6/36 = 17% (8-32%) | 20.0 [5.5, 49.5] (n=36) | 8.6 [2.1, 33.2] (n=36) | 1.17 [0.87, 2.80] (n=36) (defined in 36) | 2.1 [-1.5, 12.0] (n=36) |
- line 57
  - numpy: By platform, human GPL570: T_r2 P < 0.05 in 13/17 = 76% (53-90%); consensus count not above null in 11/17 = 65% (41-83%).
  - exact: By platform, human GPL570: T_r2 P < 0.05 in 12/16 = 75% (51-90%); consensus count not above null in 11/16 = 69% (44-86%).
- line 59
  - numpy: Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 28/35 = 80% (64-90%); T_p in 27/35 = 77% (61-88%); consensus count not above null in 23/35 = 66% (49-79%); E2 raw scheme alone 0% [0%, 17%] (n=35); covariance share carried by group means 69% [45%, 85%] (n=35).
  - exact: Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 28/35 = 80% (64-90%); T_p in 27/35 = 77% (61-88%); consensus count not above null in 23/35 = 66% (49-79%); E2 raw scheme alone 0% [0%, 22%] (n=35); covariance share carried by group means 69% [48%, 85%] (n=35).

## 5. Validation summary (validation_summary.md)

44 of 55 lines differ (left: numpy order, R run inputs; right: exact rule).

- line 5
  - numpy: | Endpoint | Development arrays (37) (n = 37) | Held-out GPL1261 (n = 20) | Held-out GPL570 (n = 20) | Held-out arrays pooled (n = 40) | RNA-seq human (n = 19) | RNA-seq pooled (n = 19) | All validation sets pooled (n = 59) | Development arrays without degenerate series (post hoc) (n = 35) | All validation sets pooled without degenerate series (post hoc) (n = 59) |
  - exact: | Endpoint | Development arrays (37) (n = 36) | Held-out GPL1261 (n = 20) | Held-out GPL570 (n = 20) | Held-out arrays pooled (n = 40) | RNA-seq human (n = 19) | RNA-seq pooled (n = 19) | All validation sets pooled (n = 59) | Development arrays without degenerate series (post hoc) (n = 35) | All validation sets pooled without degenerate series (post hoc) (n = 59) |
- line 7
  - numpy: | Within-group coupling: T_r2 P < 0.05 (primary) | 25/37 = 68% (51-80%) | 15/20 = 75% (53-89%) | 17/20 = 85% (64-95%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 23/35 = 66% (49-79%) | 47/59 = 80% (68-88%) |
  - exact: | Within-group coupling: T_r2 P < 0.05 (primary) | 22/36 = 61% (45-75%) | 15/20 = 75% (53-89%) | 17/20 = 85% (64-95%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 22/35 = 63% (46-77%) | 47/59 = 80% (68-88%) |
- line 8
  - numpy: | T_p P < 0.05 | 21/37 = 57% (41-71%) | 8/20 = 40% (22-61%) | 15/20 = 75% (53-89%) | 23/40 = 57% (42-71%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 37/59 = 63% (50-74%) | 20/35 = 57% (41-72%) | 37/59 = 63% (50-74%) |
  - exact: | T_p P < 0.05 | 20/36 = 56% (40-70%) | 8/20 = 40% (22-61%) | 15/20 = 75% (53-89%) | 23/40 = 57% (42-71%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 37/59 = 63% (50-74%) | 20/35 = 57% (41-72%) | 37/59 = 63% (50-74%) |
- line 9
  - numpy: | T_eig P < 0.05 | 21/37 = 57% (41-71%) | 10/20 = 50% (30-70%) | 13/20 = 65% (43-82%) | 23/40 = 57% (42-71%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 38/59 = 64% (52-75%) | 20/35 = 57% (41-72%) | 38/59 = 64% (52-75%) |
  - exact: | T_eig P < 0.05 | 19/36 = 53% (37-68%) | 10/20 = 50% (30-70%) | 13/20 = 65% (43-82%) | 23/40 = 57% (42-71%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 38/59 = 64% (52-75%) | 19/35 = 54% (38-70%) | 38/59 = 64% (52-75%) |
- line 10
  - numpy: | E1: consensus count not above null (P >= 0.05) | 32/37 = 86% (72-94%) | 16/20 = 80% (58-92%) | 18/20 = 90% (70-97%) | 34/40 = 85% (71-93%) | 11/19 = 58% (36-77%) | 11/19 = 58% (36-77%) | 45/59 = 76% (64-85%) | 30/35 = 86% (71-94%) | 45/59 = 76% (64-85%) |
  - exact: | E1: consensus count not above null (P >= 0.05) | 30/36 = 83% (68-92%) | 16/20 = 80% (58-92%) | 18/20 = 90% (70-97%) | 34/40 = 85% (71-93%) | 11/19 = 58% (36-77%) | 11/19 = 58% (36-77%) | 45/59 = 76% (64-85%) | 29/35 = 83% (67-92%) | 45/59 = 76% (64-85%) |
- line 11
  - numpy: | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 6%] (n=37) | 0% [0%, 1%] (n=20) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=40) | 0% [0%, 19%] (n=19) | 0% [0%, 19%] (n=19) | 0% [0%, 4%] (n=59) | 0% [0%, 5%] (n=35) | 0% [0%, 4%] (n=59) |
  - exact: | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 3%] (n=36) | 0% [0%, 1%] (n=20) | 0% [0%, 1%] (n=20) | 0% [0%, 1%] (n=40) | 0% [0%, 19%] (n=19) | 0% [0%, 19%] (n=19) | 0% [0%, 4%] (n=59) | 0% [0%, 3%] (n=35) | 0% [0%, 4%] (n=59) |
- line 12
  - numpy: | E2 any scheme (part 1 definition) | 60% [33%, 90%] (n=37) | 69% [58%, 85%] (n=20) | 73% [53%, 100%] (n=20) | 72% [58%, 95%] (n=40) | 84% [59%, 100%] (n=19) | 84% [59%, 100%] (n=19) | 77% [57%, 98%] (n=59) | 60% [29%, 89%] (n=35) | 77% [57%, 98%] (n=59) |
  - exact: | E2 any scheme (part 1 definition) | 72% [32%, 93%] (n=36) | 69% [58%, 85%] (n=20) | 73% [53%, 100%] (n=20) | 72% [58%, 93%] (n=40) | 84% [59%, 100%] (n=19) | 84% [59%, 100%] (n=19) | 77% [57%, 97%] (n=59) | 69% [29%, 91%] (n=35) | 77% [57%, 97%] (n=59) |
- line 13
  - numpy: | E2 all four schemes | 0% [0%, 2%] (n=37) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=40) | 0% [0%, 6%] (n=19) | 0% [0%, 6%] (n=19) | 0% [0%, 0%] (n=59) | 0% [0%, 0%] (n=35) | 0% [0%, 0%] (n=59) |
  - exact: | E2 all four schemes | 0% [0%, 0%] (n=36) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=40) | 0% [0%, 6%] (n=19) | 0% [0%, 6%] (n=19) | 0% [0%, 0%] (n=59) | 0% [0%, 0%] (n=35) | 0% [0%, 0%] (n=59) |
- line 14
  - numpy: | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.0] (n=37) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=40) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 1.0] (n=59) | 0.0 [0.0, 1.0] (n=35) | 0.0 [0.0, 1.0] (n=59) |
  - exact: | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.0] (n=36) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=40) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 1.0] (n=59) | 0.0 [0.0, 1.0] (n=35) | 0.0 [0.0, 1.0] (n=59) |
- line 15
  - numpy: | E3 covariance share carried by group means | 98% [96%, 99%] (n=37) | 97% [95%, 98%] (n=20) | 98% [96%, 100%] (n=20) | 97% [95%, 99%] (n=40) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=59) | 98% [96%, 99%] (n=35) | 98% [95%, 99%] (n=59) |
  - exact: | E3 covariance share carried by group means | 98% [96%, 99%] (n=36) | 97% [95%, 98%] (n=20) | 98% [96%, 100%] (n=20) | 97% [95%, 99%] (n=40) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=59) | 98% [96%, 99%] (n=35) | 98% [95%, 99%] (n=59) |
- line 16
  - numpy: | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.81 [0.62, 0.91] (n=35) | 0.73 [0.68, 0.88] (n=20) | 0.73 [0.60, 0.88] (n=20) | 0.73 [0.65, 0.88] (n=40) | 0.80 [0.62, 0.88] (n=19) | 0.80 [0.62, 0.88] (n=19) | 0.73 [0.63, 0.88] (n=59) | 0.81 [0.62, 0.91] (n=35) | 0.73 [0.63, 0.88] (n=59) |
  - exact: | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.81 [0.65, 0.91] (n=36) | 0.73 [0.68, 0.88] (n=20) | 0.73 [0.60, 0.88] (n=20) | 0.73 [0.63, 0.88] (n=40) | 0.80 [0.62, 0.88] (n=19) | 0.80 [0.62, 0.88] (n=19) | 0.73 [0.62, 0.88] (n=59) | 0.81 [0.64, 0.91] (n=35) | 0.73 [0.62, 0.88] (n=59) |
- line 17
  - numpy: | Series with share above 0.61 | 27/35 = 77% (61-88%) | 18/20 = 90% (70-97%) | 14/20 = 70% (48-85%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 27/35 = 77% (61-88%) | 47/59 = 80% (68-88%) |
  - exact: | Series with share above 0.61 | 30/36 = 83% (68-92%) | 18/20 = 90% (70-97%) | 14/20 = 70% (48-85%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 29/35 = 83% (67-92%) | 47/59 = 80% (68-88%) |
- line 18
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 9/37 = 24% (13-40%) | 3/20 = 15% (5-36%) | 7/20 = 35% (18-57%) | 10/40 = 25% (14-40%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 20/59 = 34% (23-47%) | 9/35 = 26% (14-42%) | 20/59 = 34% (23-47%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 9/36 = 25% (14-41%) | 3/20 = 15% (5-36%) | 8/20 = 40% (22-61%) | 11/40 = 28% (16-43%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 21/59 = 36% (25-48%) | 9/35 = 26% (14-42%) | 21/59 = 36% (25-48%) |
- line 19
  - numpy: | Median observed / null mean adjusted count, raw scheme | 0.00 [0.00, 13.16] (n=37) | 2.51 [0.00, 8.56] (n=20) | 0.16 [0.00, 9.71] (n=20) | 0.86 [0.00, 9.71] (n=40) | 4.79 [0.00, 30.93] (n=19) | 4.79 [0.00, 30.93] (n=19) | 1.39 [0.00, 14.60] (n=59) | 0.00 [0.00, 13.25] (n=35) | 1.39 [0.00, 14.60] (n=59) |
  - exact: | Median observed / null mean adjusted count, raw scheme | 0.00 [0.00, 13.20] (n=36) | 2.51 [0.00, 8.56] (n=20) | 0.99 [0.00, 10.76] (n=20) | 1.52 [0.00, 10.76] (n=40) | 4.79 [0.00, 30.93] (n=19) | 4.79 [0.00, 30.93] (n=19) | 1.64 [0.00, 14.60] (n=59) | 0.00 [0.00, 13.25] (n=35) | 1.64 [0.00, 14.60] (n=59) |
- line 20
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 8/37 = 22% (11-37%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 5/19 = 26% (12-49%) | 5/19 = 26% (12-49%) | 13/59 = 22% (13-34%) | 8/35 = 23% (12-39%) | 13/59 = 22% (13-34%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 8/36 = 22% (12-38%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 5/19 = 26% (12-49%) | 5/19 = 26% (12-49%) | 13/59 = 22% (13-34%) | 8/35 = 23% (12-39%) | 13/59 = 22% (13-34%) |
- line 21
  - numpy: | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.01 [0.80, 2.51] (n=37) | 0.98 [0.68, 2.35] (n=20) | 0.93 [0.38, 1.26] (n=20) | 0.96 [0.57, 1.74] (n=40) | 1.42 [0.97, 3.50] (n=19) | 1.42 [0.97, 3.50] (n=19) | 1.05 [0.73, 2.38] (n=59) | 1.01 [0.63, 3.46] (n=35) | 1.05 [0.73, 2.38] (n=59) |
  - exact: | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.02 [0.71, 2.99] (n=36) | 0.98 [0.68, 2.35] (n=20) | 0.93 [0.38, 1.94] (n=20) | 0.96 [0.57, 2.35] (n=40) | 1.42 [0.97, 3.50] (n=19) | 1.42 [0.97, 3.50] (n=19) | 1.05 [0.73, 2.78] (n=59) | 1.03 [0.63, 3.46] (n=35) | 1.05 [0.73, 2.78] (n=59) |
- line 22
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 11/37 = 30% (17-46%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 15/59 = 25% (16-38%) | 10/35 = 29% (16-45%) | 15/59 = 25% (16-38%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 9/36 = 25% (14-41%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 15/59 = 25% (16-38%) | 9/35 = 26% (14-42%) | 15/59 = 25% (16-38%) |
- line 23
  - numpy: | Median observed / null mean adjusted count, ratio to gene 2 scheme | 0.99 [0.39, 2.49] (n=37) | 1.30 [0.72, 1.96] (n=20) | 0.96 [0.34, 1.35] (n=20) | 1.04 [0.45, 1.81] (n=40) | 2.03 [1.03, 7.07] (n=19) | 2.03 [1.03, 7.07] (n=19) | 1.15 [0.72, 2.70] (n=59) | 0.96 [0.36, 2.52] (n=35) | 1.15 [0.72, 2.70] (n=59) |
  - exact: | Median observed / null mean adjusted count, ratio to gene 2 scheme | 0.99 [0.50, 2.33] (n=36) | 1.30 [0.72, 1.96] (n=20) | 0.96 [0.14, 1.35] (n=20) | 1.04 [0.37, 1.81] (n=40) | 2.03 [1.03, 7.07] (n=19) | 2.03 [1.03, 7.07] (n=19) | 1.15 [0.66, 2.70] (n=59) | 0.99 [0.46, 2.38] (n=35) | 1.15 [0.66, 2.70] (n=59) |
- line 24
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 4/37 = 11% (4-25%) | 3/20 = 15% (5-36%) | 3/20 = 15% (5-36%) | 6/40 = 15% (7-29%) | 8/19 = 42% (23-64%) | 8/19 = 42% (23-64%) | 14/59 = 24% (15-36%) | 4/35 = 11% (5-26%) | 14/59 = 24% (15-36%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 3/36 = 8% (3-22%) | 3/20 = 15% (5-36%) | 4/20 = 20% (8-42%) | 7/40 = 18% (9-32%) | 8/19 = 42% (23-64%) | 8/19 = 42% (23-64%) | 15/59 = 25% (16-38%) | 3/35 = 9% (3-22%) | 15/59 = 25% (16-38%) |
- line 25
  - numpy: | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.03 [0.45, 1.88] (n=37) | 1.26 [0.95, 2.14] (n=20) | 1.00 [0.87, 1.85] (n=20) | 1.01 [0.87, 2.14] (n=40) | 2.07 [1.10, 11.21] (n=19) | 2.07 [1.10, 11.21] (n=19) | 1.19 [0.94, 2.88] (n=59) | 1.00 [0.42, 2.02] (n=35) | 1.19 [0.94, 2.88] (n=59) |
  - exact: | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.01 [0.43, 1.54] (n=36) | 1.26 [0.95, 2.14] (n=20) | 1.00 [0.87, 1.85] (n=20) | 1.07 [0.87, 2.14] (n=40) | 2.07 [1.10, 11.21] (n=19) | 2.07 [1.10, 11.21] (n=19) | 1.20 [0.94, 2.88] (n=59) | 1.03 [0.42, 1.65] (n=35) | 1.20 [0.94, 2.88] (n=59) |
- line 27
  - numpy: T_r2 prevalence against the 5% expected under no coupling (exact binomial, one-sided): Development arrays (37) P = 3.1e-24; Held-out GPL1261 P = 3.7e-16; Held-out GPL570 P = 7.5e-20; Held-out arrays pooled P = 1.2e-34; RNA-seq human P = 9.8e-17; RNA-seq pooled P = 9.8e-17; All validation sets pooled P = 4.4e-50; Development arrays without degenerate series (post hoc) P = 5.5e-22; All validation sets pooled without degenerate series (post hoc) P = 4.4e-50
  - exact: T_r2 prevalence against the 5% expected under no coupling (exact binomial, one-sided): Development arrays (37) P = 4.6e-20; Held-out GPL1261 P = 3.7e-16; Held-out GPL570 P = 7.5e-20; Held-out arrays pooled P = 1.2e-34; RNA-seq human P = 9.8e-17; RNA-seq pooled P = 9.8e-17; All validation sets pooled P = 4.4e-50; Development arrays without degenerate series (post hoc) P = 1.9e-20; All validation sets pooled without degenerate series (post hoc) P = 4.4e-50
- line 28
  - numpy: Held-out arrays versus development (Fisher exact, two-sided): T_r2 prevalence P = 0.3; consensus count not above null P = 1.
  - exact: Held-out arrays versus development (Fisher exact, two-sided): T_r2 prevalence P = 0.082; consensus count not above null P = 1.
- line 32
  - numpy: | Endpoint | Development arrays (37) (n = 37) | Held-out GPL1261 (n = 20) | Held-out GPL570 (n = 20) | Held-out arrays pooled (n = 40) | RNA-seq human (n = 19) | RNA-seq pooled (n = 19) | All validation sets pooled (n = 59) | Development arrays without degenerate series (post hoc) (n = 35) | All validation sets pooled without degenerate series (post hoc) (n = 59) |
  - exact: | Endpoint | Development arrays (37) (n = 36) | Held-out GPL1261 (n = 20) | Held-out GPL570 (n = 20) | Held-out arrays pooled (n = 40) | RNA-seq human (n = 19) | RNA-seq pooled (n = 19) | All validation sets pooled (n = 59) | Development arrays without degenerate series (post hoc) (n = 35) | All validation sets pooled without degenerate series (post hoc) (n = 59) |
- line 34
  - numpy: | Within-group coupling: T_r2 P < 0.05 (primary) | 29/37 = 78% (63-89%) | 15/20 = 75% (53-89%) | 16/20 = 80% (58-92%) | 31/40 = 78% (62-88%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 48/59 = 81% (70-89%) | 28/35 = 80% (64-90%) | 48/59 = 81% (70-89%) |
  - exact: | Within-group coupling: T_r2 P < 0.05 (primary) | 28/36 = 78% (62-88%) | 15/20 = 75% (53-89%) | 16/20 = 80% (58-92%) | 31/40 = 78% (62-88%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 48/59 = 81% (70-89%) | 28/35 = 80% (64-90%) | 48/59 = 81% (70-89%) |
- line 35
  - numpy: | T_p P < 0.05 | 28/37 = 76% (60-87%) | 10/20 = 50% (30-70%) | 16/20 = 80% (58-92%) | 26/40 = 65% (50-78%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 43/59 = 73% (60-83%) | 27/35 = 77% (61-88%) | 43/59 = 73% (60-83%) |
  - exact: | T_p P < 0.05 | 27/36 = 75% (59-86%) | 10/20 = 50% (30-70%) | 16/20 = 80% (58-92%) | 26/40 = 65% (50-78%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 43/59 = 73% (60-83%) | 27/35 = 77% (61-88%) | 43/59 = 73% (60-83%) |
- line 36
  - numpy: | T_eig P < 0.05 | 27/37 = 73% (57-85%) | 14/20 = 70% (48-85%) | 16/20 = 80% (58-92%) | 30/40 = 75% (60-86%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 47/59 = 80% (68-88%) | 26/35 = 74% (58-86%) | 47/59 = 80% (68-88%) |
  - exact: | T_eig P < 0.05 | 26/36 = 72% (56-84%) | 14/20 = 70% (48-85%) | 16/20 = 80% (58-92%) | 30/40 = 75% (60-86%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 47/59 = 80% (68-88%) | 26/35 = 74% (58-86%) | 47/59 = 80% (68-88%) |
- line 37
  - numpy: | E1: consensus count not above null (P >= 0.05) | 24/37 = 65% (49-78%) | 12/20 = 60% (39-78%) | 13/20 = 65% (43-82%) | 25/40 = 62% (47-76%) | 9/19 = 47% (27-68%) | 9/19 = 47% (27-68%) | 34/59 = 58% (45-69%) | 23/35 = 66% (49-79%) | 34/59 = 58% (45-69%) |
  - exact: | E1: consensus count not above null (P >= 0.05) | 24/36 = 67% (50-80%) | 12/20 = 60% (39-78%) | 13/20 = 65% (43-82%) | 25/40 = 62% (47-76%) | 9/19 = 47% (27-68%) | 9/19 = 47% (27-68%) | 34/59 = 58% (45-69%) | 23/35 = 66% (49-79%) | 34/59 = 58% (45-69%) |
- line 38
  - numpy: | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 20%] (n=37) | 0% [0%, 15%] (n=20) | 0% [0%, 41%] (n=20) | 0% [0%, 29%] (n=40) | 21% [4%, 42%] (n=19) | 21% [4%, 42%] (n=19) | 3% [0%, 35%] (n=59) | 0% [0%, 17%] (n=35) | 3% [0%, 35%] (n=59) |
  - exact: | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 21%] (n=36) | 0% [0%, 15%] (n=20) | 0% [0%, 41%] (n=20) | 0% [0%, 29%] (n=40) | 21% [4%, 42%] (n=19) | 21% [4%, 42%] (n=19) | 3% [0%, 35%] (n=59) | 0% [0%, 22%] (n=35) | 3% [0%, 35%] (n=59) |
- line 39
  - numpy: | E2 any scheme (part 1 definition) | 100% [67%, 100%] (n=37) | 95% [69%, 100%] (n=20) | 96% [74%, 100%] (n=20) | 96% [73%, 100%] (n=40) | 100% [84%, 100%] (n=19) | 100% [84%, 100%] (n=19) | 100% [75%, 100%] (n=59) | 95% [64%, 100%] (n=35) | 100% [75%, 100%] (n=59) |
  - exact: | E2 any scheme (part 1 definition) | 98% [65%, 100%] (n=36) | 95% [69%, 100%] (n=20) | 96% [74%, 100%] (n=20) | 96% [73%, 100%] (n=40) | 100% [84%, 100%] (n=19) | 100% [84%, 100%] (n=19) | 100% [75%, 100%] (n=59) | 95% [64%, 100%] (n=35) | 100% [75%, 100%] (n=59) |
- line 40
  - numpy: | E2 all four schemes | 0% [0%, 0%] (n=37) | 0% [0%, 13%] (n=20) | 0% [0%, 22%] (n=20) | 0% [0%, 17%] (n=40) | 8% [0%, 17%] (n=19) | 8% [0%, 17%] (n=19) | 0% [0%, 17%] (n=59) | 0% [0%, 0%] (n=35) | 0% [0%, 17%] (n=59) |
  - exact: | E2 all four schemes | 0% [0%, 2%] (n=36) | 0% [0%, 13%] (n=20) | 0% [0%, 22%] (n=20) | 0% [0%, 17%] (n=40) | 8% [0%, 17%] (n=19) | 8% [0%, 17%] (n=19) | 0% [0%, 17%] (n=59) | 0% [0%, 5%] (n=35) | 0% [0%, 17%] (n=59) |
- line 41
  - numpy: | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.0] (n=37) | 0.0 [0.0, 2.5] (n=20) | 0.5 [0.0, 4.0] (n=20) | 0.0 [0.0, 4.0] (n=40) | 3.0 [1.0, 8.5] (n=19) | 3.0 [1.0, 8.5] (n=19) | 1.0 [0.0, 5.0] (n=59) | 0.0 [0.0, 1.0] (n=35) | 1.0 [0.0, 5.0] (n=59) |
  - exact: | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.5] (n=36) | 0.0 [0.0, 2.5] (n=20) | 0.5 [0.0, 4.0] (n=20) | 0.0 [0.0, 4.0] (n=40) | 3.0 [1.0, 8.5] (n=19) | 3.0 [1.0, 8.5] (n=19) | 1.0 [0.0, 5.0] (n=59) | 0.0 [0.0, 2.0] (n=35) | 1.0 [0.0, 5.0] (n=59) |
- line 42
  - numpy: | E3 covariance share carried by group means | 69% [42%, 86%] (n=37) | 90% [64%, 95%] (n=20) | 82% [33%, 93%] (n=20) | 85% [41%, 95%] (n=40) | 80% [31%, 97%] (n=19) | 80% [31%, 97%] (n=19) | 85% [40%, 95%] (n=59) | 69% [45%, 85%] (n=35) | 85% [40%, 95%] (n=59) |
  - exact: | E3 covariance share carried by group means | 69% [48%, 87%] (n=36) | 90% [64%, 95%] (n=20) | 82% [33%, 93%] (n=20) | 85% [41%, 95%] (n=40) | 80% [31%, 97%] (n=19) | 80% [31%, 97%] (n=19) | 85% [40%, 95%] (n=59) | 69% [48%, 85%] (n=35) | 85% [40%, 95%] (n=59) |
- line 43
  - numpy: | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.70 [0.40, 0.94] (n=27) | 0.60 [0.55, 0.80] (n=17) | 0.70 [0.51, 0.94] (n=15) | 0.66 [0.53, 0.87] (n=32) | 0.71 [0.46, 0.79] (n=15) | 0.71 [0.46, 0.79] (n=15) | 0.67 [0.50, 0.83] (n=47) | 0.70 [0.39, 0.96] (n=26) | 0.67 [0.50, 0.83] (n=47) |
  - exact: | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.69 [0.37, 0.91] (n=28) | 0.60 [0.55, 0.80] (n=17) | 0.70 [0.51, 0.94] (n=15) | 0.66 [0.53, 0.87] (n=32) | 0.71 [0.46, 0.79] (n=15) | 0.71 [0.46, 0.79] (n=15) | 0.67 [0.50, 0.83] (n=47) | 0.67 [0.36, 0.94] (n=27) | 0.67 [0.50, 0.83] (n=47) |
- line 44
  - numpy: | Series with share above 0.61 | 15/27 = 56% (37-72%) | 8/17 = 47% (26-69%) | 10/15 = 67% (42-85%) | 18/32 = 56% (39-72%) | 10/15 = 67% (42-85%) | 10/15 = 67% (42-85%) | 28/47 = 60% (45-72%) | 14/26 = 54% (35-71%) | 28/47 = 60% (45-72%) |
  - exact: | Series with share above 0.61 | 15/28 = 54% (36-70%) | 8/17 = 47% (26-69%) | 10/15 = 67% (42-85%) | 18/32 = 56% (39-72%) | 10/15 = 67% (42-85%) | 10/15 = 67% (42-85%) | 28/47 = 60% (45-72%) | 14/27 = 52% (34-69%) | 28/47 = 60% (45-72%) |
- line 45
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 14/37 = 38% (24-54%) | 8/20 = 40% (22-61%) | 9/20 = 45% (26-66%) | 17/40 = 42% (29-58%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 31/59 = 53% (40-65%) | 13/35 = 37% (23-54%) | 31/59 = 53% (40-65%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 13/36 = 36% (22-52%) | 8/20 = 40% (22-61%) | 9/20 = 45% (26-66%) | 17/40 = 42% (29-58%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 31/59 = 53% (40-65%) | 13/35 = 37% (23-54%) | 31/59 = 53% (40-65%) |
- line 46
  - numpy: | Median observed / null mean adjusted count, raw scheme | 4.64 [0.00, 31.25] (n=37) | 0.00 [0.00, 36.49] (n=20) | 5.38 [0.00, 45.86] (n=20) | 3.37 [0.00, 40.89] (n=40) | 18.87 [6.52, 79.70] (n=19) | 18.87 [6.52, 79.70] (n=19) | 8.77 [0.00, 50.96] (n=59) | 4.64 [0.00, 30.33] (n=35) | 8.77 [0.00, 50.96] (n=59) |
  - exact: | Median observed / null mean adjusted count, raw scheme | 4.24 [0.00, 29.87] (n=36) | 0.00 [0.00, 36.49] (n=20) | 5.38 [0.00, 45.86] (n=20) | 3.37 [0.00, 40.89] (n=40) | 18.87 [6.52, 79.70] (n=19) | 18.87 [6.52, 79.70] (n=19) | 8.77 [0.00, 50.96] (n=59) | 4.64 [0.00, 30.33] (n=35) | 8.77 [0.00, 50.96] (n=59) |
- line 47
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 12/37 = 32% (20-49%) | 5/20 = 25% (11-47%) | 9/20 = 45% (26-66%) | 14/40 = 35% (22-50%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 21/59 = 36% (25-48%) | 11/35 = 31% (19-48%) | 21/59 = 36% (25-48%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 11/36 = 31% (18-47%) | 5/20 = 25% (11-47%) | 9/20 = 45% (26-66%) | 14/40 = 35% (22-50%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 21/59 = 36% (25-48%) | 11/35 = 31% (19-48%) | 21/59 = 36% (25-48%) |
- line 48
  - numpy: | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.61 [0.96, 5.60] (n=37) | 1.13 [0.43, 2.93] (n=20) | 2.50 [1.18, 4.38] (n=20) | 1.56 [0.78, 4.04] (n=40) | 1.24 [1.00, 2.53] (n=19) | 1.24 [1.00, 2.53] (n=19) | 1.51 [0.96, 3.77] (n=59) | 1.61 [0.94, 4.52] (n=35) | 1.51 [0.96, 3.77] (n=59) |
  - exact: | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.53 [0.95, 3.98] (n=36) | 1.13 [0.43, 2.93] (n=20) | 2.50 [1.18, 4.38] (n=20) | 1.56 [0.78, 4.04] (n=40) | 1.24 [1.00, 2.53] (n=19) | 1.24 [1.00, 2.53] (n=19) | 1.51 [0.96, 3.77] (n=59) | 1.61 [0.94, 4.52] (n=35) | 1.51 [0.96, 3.77] (n=59) |
- line 49
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 8/37 = 22% (11-37%) | 5/20 = 25% (11-47%) | 11/20 = 55% (34-74%) | 16/40 = 40% (26-55%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 26/59 = 44% (32-57%) | 7/35 = 20% (10-36%) | 26/59 = 44% (32-57%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 7/36 = 19% (10-35%) | 5/20 = 25% (11-47%) | 11/20 = 55% (34-74%) | 16/40 = 40% (26-55%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 26/59 = 44% (32-57%) | 7/35 = 20% (10-36%) | 26/59 = 44% (32-57%) |
- line 50
  - numpy: | Median observed / null mean adjusted count, ratio to gene 2 scheme | 1.08 [0.98, 1.81] (n=37) | 1.10 [0.98, 3.68] (n=20) | 2.54 [1.08, 20.11] (n=20) | 1.38 [1.01, 7.25] (n=40) | 1.92 [1.02, 9.19] (n=19) | 1.92 [1.02, 9.19] (n=19) | 1.40 [1.01, 7.40] (n=59) | 1.08 [0.99, 1.77] (n=35) | 1.40 [1.01, 7.40] (n=59) |
  - exact: | Median observed / null mean adjusted count, ratio to gene 2 scheme | 1.05 [0.98, 1.76] (n=36) | 1.10 [0.98, 3.68] (n=20) | 2.54 [1.08, 20.11] (n=20) | 1.38 [1.01, 7.25] (n=40) | 1.92 [1.02, 9.19] (n=19) | 1.92 [1.02, 9.19] (n=19) | 1.40 [1.01, 7.40] (n=59) | 1.08 [0.99, 1.77] (n=35) | 1.40 [1.01, 7.40] (n=59) |
- line 51
  - numpy: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 7/37 = 19% (9-34%) | 5/20 = 25% (11-47%) | 8/20 = 40% (22-61%) | 13/40 = 32% (20-48%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 20/59 = 34% (23-47%) | 6/35 = 17% (8-33%) | 20/59 = 34% (23-47%) |
  - exact: | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 6/36 = 17% (8-32%) | 5/20 = 25% (11-47%) | 8/20 = 40% (22-61%) | 13/40 = 32% (20-48%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 20/59 = 34% (23-47%) | 6/35 = 17% (8-33%) | 20/59 = 34% (23-47%) |
- line 52
  - numpy: | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.23 [0.84, 2.83] (n=37) | 0.99 [0.69, 2.30] (n=20) | 1.35 [0.99, 3.12] (n=20) | 1.15 [0.92, 2.83] (n=40) | 1.35 [0.96, 3.95] (n=19) | 1.35 [0.96, 3.95] (n=19) | 1.20 [0.94, 2.89] (n=59) | 1.23 [0.79, 2.81] (n=35) | 1.20 [0.94, 2.89] (n=59) |
  - exact: | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.17 [0.87, 2.80] (n=36) | 0.99 [0.69, 2.30] (n=20) | 1.35 [0.99, 3.12] (n=20) | 1.15 [0.92, 2.83] (n=40) | 1.35 [0.96, 3.95] (n=19) | 1.35 [0.96, 3.95] (n=19) | 1.20 [0.94, 2.89] (n=59) | 1.23 [0.86, 2.81] (n=35) | 1.20 [0.94, 2.89] (n=59) |
- line 54
  - numpy: T_r2 prevalence against the 5% expected under no coupling (exact binomial, one-sided): Development arrays (37) P = 4.8e-31; Held-out GPL1261 P = 3.7e-16; Held-out GPL570 P = 6.1e-18; Held-out arrays pooled P = 8.1e-33; RNA-seq human P = 1.2e-20; RNA-seq pooled P = 1.2e-20; All validation sets pooled P = 5.7e-52; Development arrays without degenerate series (post hoc) P = 1.8e-30; All validation sets pooled without degenerate series (post hoc) P = 5.7e-52
  - exact: T_r2 prevalence against the 5% expected under no coupling (exact binomial, one-sided): Development arrays (37) P = 7.6e-30; Held-out GPL1261 P = 3.7e-16; Held-out GPL570 P = 6.1e-18; Held-out arrays pooled P = 8.1e-33; RNA-seq human P = 1.2e-20; RNA-seq pooled P = 1.2e-20; All validation sets pooled P = 5.7e-52; Development arrays without degenerate series (post hoc) P = 1.8e-30; All validation sets pooled without degenerate series (post hoc) P = 5.7e-52
- line 55
  - numpy: Held-out arrays versus development (Fisher exact, two-sided): T_r2 prevalence P = 1; consensus count not above null P = 1.
  - exact: Held-out arrays versus development (Fisher exact, two-sided): T_r2 prevalence P = 1; consensus count not above null P = 0.81.

## 6. WP2 summary table (WP2_summary_table.md)

38 of 40 lines differ (left: numpy order, R run inputs; right: exact rule).

- line 3
  - numpy: | A | Within-group coupling: T_r2 P < 0.05 (primary) | 25/37 = 68% (51-80%) | 15/20 = 75% (53-89%) | 17/20 = 85% (64-95%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 23/35 = 66% (49-79%) | 47/59 = 80% (68-88%) |
  - exact: | A | Within-group coupling: T_r2 P < 0.05 (primary) | 22/36 = 61% (45-75%) | 15/20 = 75% (53-89%) | 17/20 = 85% (64-95%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 22/35 = 63% (46-77%) | 47/59 = 80% (68-88%) |
- line 4
  - numpy: | A | T_p P < 0.05 | 21/37 = 57% (41-71%) | 8/20 = 40% (22-61%) | 15/20 = 75% (53-89%) | 23/40 = 57% (42-71%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 37/59 = 63% (50-74%) | 20/35 = 57% (41-72%) | 37/59 = 63% (50-74%) |
  - exact: | A | T_p P < 0.05 | 20/36 = 56% (40-70%) | 8/20 = 40% (22-61%) | 15/20 = 75% (53-89%) | 23/40 = 57% (42-71%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 37/59 = 63% (50-74%) | 20/35 = 57% (41-72%) | 37/59 = 63% (50-74%) |
- line 5
  - numpy: | A | T_eig P < 0.05 | 21/37 = 57% (41-71%) | 10/20 = 50% (30-70%) | 13/20 = 65% (43-82%) | 23/40 = 57% (42-71%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 38/59 = 64% (52-75%) | 20/35 = 57% (41-72%) | 38/59 = 64% (52-75%) |
  - exact: | A | T_eig P < 0.05 | 19/36 = 53% (37-68%) | 10/20 = 50% (30-70%) | 13/20 = 65% (43-82%) | 23/40 = 57% (42-71%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 38/59 = 64% (52-75%) | 19/35 = 54% (38-70%) | 38/59 = 64% (52-75%) |
- line 6
  - numpy: | A | E1: consensus count not above null (P >= 0.05) | 32/37 = 86% (72-94%) | 16/20 = 80% (58-92%) | 18/20 = 90% (70-97%) | 34/40 = 85% (71-93%) | 11/19 = 58% (36-77%) | 11/19 = 58% (36-77%) | 45/59 = 76% (64-85%) | 30/35 = 86% (71-94%) | 45/59 = 76% (64-85%) |
  - exact: | A | E1: consensus count not above null (P >= 0.05) | 30/36 = 83% (68-92%) | 16/20 = 80% (58-92%) | 18/20 = 90% (70-97%) | 34/40 = 85% (71-93%) | 11/19 = 58% (36-77%) | 11/19 = 58% (36-77%) | 45/59 = 76% (64-85%) | 29/35 = 83% (67-92%) | 45/59 = 76% (64-85%) |
- line 7
  - numpy: | A | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 6%] (n=37) | 0% [0%, 1%] (n=20) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=40) | 0% [0%, 19%] (n=19) | 0% [0%, 19%] (n=19) | 0% [0%, 4%] (n=59) | 0% [0%, 5%] (n=35) | 0% [0%, 4%] (n=59) |
  - exact: | A | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 3%] (n=36) | 0% [0%, 1%] (n=20) | 0% [0%, 1%] (n=20) | 0% [0%, 1%] (n=40) | 0% [0%, 19%] (n=19) | 0% [0%, 19%] (n=19) | 0% [0%, 4%] (n=59) | 0% [0%, 3%] (n=35) | 0% [0%, 4%] (n=59) |
- line 8
  - numpy: | A | E2 any scheme (part 1 definition) | 60% [33%, 90%] (n=37) | 69% [58%, 85%] (n=20) | 73% [53%, 100%] (n=20) | 72% [58%, 95%] (n=40) | 84% [59%, 100%] (n=19) | 84% [59%, 100%] (n=19) | 77% [57%, 98%] (n=59) | 60% [29%, 89%] (n=35) | 77% [57%, 98%] (n=59) |
  - exact: | A | E2 any scheme (part 1 definition) | 72% [32%, 93%] (n=36) | 69% [58%, 85%] (n=20) | 73% [53%, 100%] (n=20) | 72% [58%, 93%] (n=40) | 84% [59%, 100%] (n=19) | 84% [59%, 100%] (n=19) | 77% [57%, 97%] (n=59) | 69% [29%, 91%] (n=35) | 77% [57%, 97%] (n=59) |
- line 9
  - numpy: | A | E2 all four schemes | 0% [0%, 2%] (n=37) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=40) | 0% [0%, 6%] (n=19) | 0% [0%, 6%] (n=19) | 0% [0%, 0%] (n=59) | 0% [0%, 0%] (n=35) | 0% [0%, 0%] (n=59) |
  - exact: | A | E2 all four schemes | 0% [0%, 0%] (n=36) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=20) | 0% [0%, 0%] (n=40) | 0% [0%, 6%] (n=19) | 0% [0%, 6%] (n=19) | 0% [0%, 0%] (n=59) | 0% [0%, 0%] (n=35) | 0% [0%, 0%] (n=59) |
- line 10
  - numpy: | A | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.0] (n=37) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=40) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 1.0] (n=59) | 0.0 [0.0, 1.0] (n=35) | 0.0 [0.0, 1.0] (n=59) |
  - exact: | A | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.0] (n=36) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=20) | 0.0 [0.0, 1.0] (n=40) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 5.5] (n=19) | 0.0 [0.0, 1.0] (n=59) | 0.0 [0.0, 1.0] (n=35) | 0.0 [0.0, 1.0] (n=59) |
- line 11
  - numpy: | A | E3 covariance share carried by group means | 98% [96%, 99%] (n=37) | 97% [95%, 98%] (n=20) | 98% [96%, 100%] (n=20) | 97% [95%, 99%] (n=40) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=59) | 98% [96%, 99%] (n=35) | 98% [95%, 99%] (n=59) |
  - exact: | A | E3 covariance share carried by group means | 98% [96%, 99%] (n=36) | 97% [95%, 98%] (n=20) | 98% [96%, 100%] (n=20) | 97% [95%, 99%] (n=40) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=19) | 98% [95%, 99%] (n=59) | 98% [96%, 99%] (n=35) | 98% [95%, 99%] (n=59) |
- line 12
  - numpy: | A | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.81 [0.62, 0.91] (n=35) | 0.73 [0.68, 0.88] (n=20) | 0.73 [0.60, 0.88] (n=20) | 0.73 [0.65, 0.88] (n=40) | 0.80 [0.62, 0.88] (n=19) | 0.80 [0.62, 0.88] (n=19) | 0.73 [0.63, 0.88] (n=59) | 0.81 [0.62, 0.91] (n=35) | 0.73 [0.63, 0.88] (n=59) |
  - exact: | A | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.81 [0.65, 0.91] (n=36) | 0.73 [0.68, 0.88] (n=20) | 0.73 [0.60, 0.88] (n=20) | 0.73 [0.63, 0.88] (n=40) | 0.80 [0.62, 0.88] (n=19) | 0.80 [0.62, 0.88] (n=19) | 0.73 [0.62, 0.88] (n=59) | 0.81 [0.64, 0.91] (n=35) | 0.73 [0.62, 0.88] (n=59) |
- line 13
  - numpy: | A | Series with share above 0.61 | 27/35 = 77% (61-88%) | 18/20 = 90% (70-97%) | 14/20 = 70% (48-85%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 27/35 = 77% (61-88%) | 47/59 = 80% (68-88%) |
  - exact: | A | Series with share above 0.61 | 30/36 = 83% (68-92%) | 18/20 = 90% (70-97%) | 14/20 = 70% (48-85%) | 32/40 = 80% (65-90%) | 15/19 = 79% (57-91%) | 15/19 = 79% (57-91%) | 47/59 = 80% (68-88%) | 29/35 = 83% (67-92%) | 47/59 = 80% (68-88%) |
- line 14
  - numpy: | A | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 9/37 = 24% (13-40%) | 3/20 = 15% (5-36%) | 7/20 = 35% (18-57%) | 10/40 = 25% (14-40%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 20/59 = 34% (23-47%) | 9/35 = 26% (14-42%) | 20/59 = 34% (23-47%) |
  - exact: | A | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 9/36 = 25% (14-41%) | 3/20 = 15% (5-36%) | 8/20 = 40% (22-61%) | 11/40 = 28% (16-43%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 21/59 = 36% (25-48%) | 9/35 = 26% (14-42%) | 21/59 = 36% (25-48%) |
- line 15
  - numpy: | A | Median observed / null mean adjusted count, raw scheme | 0.00 [0.00, 13.16] (n=37) | 2.51 [0.00, 8.56] (n=20) | 0.16 [0.00, 9.71] (n=20) | 0.86 [0.00, 9.71] (n=40) | 4.79 [0.00, 30.93] (n=19) | 4.79 [0.00, 30.93] (n=19) | 1.39 [0.00, 14.60] (n=59) | 0.00 [0.00, 13.25] (n=35) | 1.39 [0.00, 14.60] (n=59) |
  - exact: | A | Median observed / null mean adjusted count, raw scheme | 0.00 [0.00, 13.20] (n=36) | 2.51 [0.00, 8.56] (n=20) | 0.99 [0.00, 10.76] (n=20) | 1.52 [0.00, 10.76] (n=40) | 4.79 [0.00, 30.93] (n=19) | 4.79 [0.00, 30.93] (n=19) | 1.64 [0.00, 14.60] (n=59) | 0.00 [0.00, 13.25] (n=35) | 1.64 [0.00, 14.60] (n=59) |
- line 16
  - numpy: | A | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 8/37 = 22% (11-37%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 5/19 = 26% (12-49%) | 5/19 = 26% (12-49%) | 13/59 = 22% (13-34%) | 8/35 = 23% (12-39%) | 13/59 = 22% (13-34%) |
  - exact: | A | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 8/36 = 22% (12-38%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 5/19 = 26% (12-49%) | 5/19 = 26% (12-49%) | 13/59 = 22% (13-34%) | 8/35 = 23% (12-39%) | 13/59 = 22% (13-34%) |
- line 17
  - numpy: | A | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.01 [0.80, 2.51] (n=37) | 0.98 [0.68, 2.35] (n=20) | 0.93 [0.38, 1.26] (n=20) | 0.96 [0.57, 1.74] (n=40) | 1.42 [0.97, 3.50] (n=19) | 1.42 [0.97, 3.50] (n=19) | 1.05 [0.73, 2.38] (n=59) | 1.01 [0.63, 3.46] (n=35) | 1.05 [0.73, 2.38] (n=59) |
  - exact: | A | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.02 [0.71, 2.99] (n=36) | 0.98 [0.68, 2.35] (n=20) | 0.93 [0.38, 1.94] (n=20) | 0.96 [0.57, 2.35] (n=40) | 1.42 [0.97, 3.50] (n=19) | 1.42 [0.97, 3.50] (n=19) | 1.05 [0.73, 2.78] (n=59) | 1.03 [0.63, 3.46] (n=35) | 1.05 [0.73, 2.78] (n=59) |
- line 18
  - numpy: | A | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 11/37 = 30% (17-46%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 15/59 = 25% (16-38%) | 10/35 = 29% (16-45%) | 15/59 = 25% (16-38%) |
  - exact: | A | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 9/36 = 25% (14-41%) | 4/20 = 20% (8-42%) | 4/20 = 20% (8-42%) | 8/40 = 20% (10-35%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 15/59 = 25% (16-38%) | 9/35 = 26% (14-42%) | 15/59 = 25% (16-38%) |
- line 19
  - numpy: | A | Median observed / null mean adjusted count, ratio to gene 2 scheme | 0.99 [0.39, 2.49] (n=37) | 1.30 [0.72, 1.96] (n=20) | 0.96 [0.34, 1.35] (n=20) | 1.04 [0.45, 1.81] (n=40) | 2.03 [1.03, 7.07] (n=19) | 2.03 [1.03, 7.07] (n=19) | 1.15 [0.72, 2.70] (n=59) | 0.96 [0.36, 2.52] (n=35) | 1.15 [0.72, 2.70] (n=59) |
  - exact: | A | Median observed / null mean adjusted count, ratio to gene 2 scheme | 0.99 [0.50, 2.33] (n=36) | 1.30 [0.72, 1.96] (n=20) | 0.96 [0.14, 1.35] (n=20) | 1.04 [0.37, 1.81] (n=40) | 2.03 [1.03, 7.07] (n=19) | 2.03 [1.03, 7.07] (n=19) | 1.15 [0.66, 2.70] (n=59) | 0.99 [0.46, 2.38] (n=35) | 1.15 [0.66, 2.70] (n=59) |
- line 20
  - numpy: | A | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 4/37 = 11% (4-25%) | 3/20 = 15% (5-36%) | 3/20 = 15% (5-36%) | 6/40 = 15% (7-29%) | 8/19 = 42% (23-64%) | 8/19 = 42% (23-64%) | 14/59 = 24% (15-36%) | 4/35 = 11% (5-26%) | 14/59 = 24% (15-36%) |
  - exact: | A | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 3/36 = 8% (3-22%) | 3/20 = 15% (5-36%) | 4/20 = 20% (8-42%) | 7/40 = 18% (9-32%) | 8/19 = 42% (23-64%) | 8/19 = 42% (23-64%) | 15/59 = 25% (16-38%) | 3/35 = 9% (3-22%) | 15/59 = 25% (16-38%) |
- line 21
  - numpy: | A | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.03 [0.45, 1.88] (n=37) | 1.26 [0.95, 2.14] (n=20) | 1.00 [0.87, 1.85] (n=20) | 1.01 [0.87, 2.14] (n=40) | 2.07 [1.10, 11.21] (n=19) | 2.07 [1.10, 11.21] (n=19) | 1.19 [0.94, 2.88] (n=59) | 1.00 [0.42, 2.02] (n=35) | 1.19 [0.94, 2.88] (n=59) |
  - exact: | A | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.01 [0.43, 1.54] (n=36) | 1.26 [0.95, 2.14] (n=20) | 1.00 [0.87, 1.85] (n=20) | 1.07 [0.87, 2.14] (n=40) | 2.07 [1.10, 11.21] (n=19) | 2.07 [1.10, 11.21] (n=19) | 1.20 [0.94, 2.88] (n=59) | 1.03 [0.42, 1.65] (n=35) | 1.20 [0.94, 2.88] (n=59) |
- line 22
  - numpy: | B | Within-group coupling: T_r2 P < 0.05 (primary) | 29/37 = 78% (63-89%) | 15/20 = 75% (53-89%) | 16/20 = 80% (58-92%) | 31/40 = 78% (62-88%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 48/59 = 81% (70-89%) | 28/35 = 80% (64-90%) | 48/59 = 81% (70-89%) |
  - exact: | B | Within-group coupling: T_r2 P < 0.05 (primary) | 28/36 = 78% (62-88%) | 15/20 = 75% (53-89%) | 16/20 = 80% (58-92%) | 31/40 = 78% (62-88%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 48/59 = 81% (70-89%) | 28/35 = 80% (64-90%) | 48/59 = 81% (70-89%) |
- line 23
  - numpy: | B | T_p P < 0.05 | 28/37 = 76% (60-87%) | 10/20 = 50% (30-70%) | 16/20 = 80% (58-92%) | 26/40 = 65% (50-78%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 43/59 = 73% (60-83%) | 27/35 = 77% (61-88%) | 43/59 = 73% (60-83%) |
  - exact: | B | T_p P < 0.05 | 27/36 = 75% (59-86%) | 10/20 = 50% (30-70%) | 16/20 = 80% (58-92%) | 26/40 = 65% (50-78%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 43/59 = 73% (60-83%) | 27/35 = 77% (61-88%) | 43/59 = 73% (60-83%) |
- line 24
  - numpy: | B | T_eig P < 0.05 | 27/37 = 73% (57-85%) | 14/20 = 70% (48-85%) | 16/20 = 80% (58-92%) | 30/40 = 75% (60-86%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 47/59 = 80% (68-88%) | 26/35 = 74% (58-86%) | 47/59 = 80% (68-88%) |
  - exact: | B | T_eig P < 0.05 | 26/36 = 72% (56-84%) | 14/20 = 70% (48-85%) | 16/20 = 80% (58-92%) | 30/40 = 75% (60-86%) | 17/19 = 89% (69-97%) | 17/19 = 89% (69-97%) | 47/59 = 80% (68-88%) | 26/35 = 74% (58-86%) | 47/59 = 80% (68-88%) |
- line 25
  - numpy: | B | E1: consensus count not above null (P >= 0.05) | 24/37 = 65% (49-78%) | 12/20 = 60% (39-78%) | 13/20 = 65% (43-82%) | 25/40 = 62% (47-76%) | 9/19 = 47% (27-68%) | 9/19 = 47% (27-68%) | 34/59 = 58% (45-69%) | 23/35 = 66% (49-79%) | 34/59 = 58% (45-69%) |
  - exact: | B | E1: consensus count not above null (P >= 0.05) | 24/36 = 67% (50-80%) | 12/20 = 60% (39-78%) | 13/20 = 65% (43-82%) | 25/40 = 62% (47-76%) | 9/19 = 47% (27-68%) | 9/19 = 47% (27-68%) | 34/59 = 58% (45-69%) | 23/35 = 66% (49-79%) | 34/59 = 58% (45-69%) |
- line 26
  - numpy: | B | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 20%] (n=37) | 0% [0%, 15%] (n=20) | 0% [0%, 41%] (n=20) | 0% [0%, 29%] (n=40) | 21% [4%, 42%] (n=19) | 21% [4%, 42%] (n=19) | 3% [0%, 35%] (n=59) | 0% [0%, 17%] (n=35) | 3% [0%, 35%] (n=59) |
  - exact: | B | E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3] | 0% [0%, 21%] (n=36) | 0% [0%, 15%] (n=20) | 0% [0%, 41%] (n=20) | 0% [0%, 29%] (n=40) | 21% [4%, 42%] (n=19) | 21% [4%, 42%] (n=19) | 3% [0%, 35%] (n=59) | 0% [0%, 22%] (n=35) | 3% [0%, 35%] (n=59) |
- line 27
  - numpy: | B | E2 any scheme (part 1 definition) | 100% [67%, 100%] (n=37) | 95% [69%, 100%] (n=20) | 96% [74%, 100%] (n=20) | 96% [73%, 100%] (n=40) | 100% [84%, 100%] (n=19) | 100% [84%, 100%] (n=19) | 100% [75%, 100%] (n=59) | 95% [64%, 100%] (n=35) | 100% [75%, 100%] (n=59) |
  - exact: | B | E2 any scheme (part 1 definition) | 98% [65%, 100%] (n=36) | 95% [69%, 100%] (n=20) | 96% [74%, 100%] (n=20) | 96% [73%, 100%] (n=40) | 100% [84%, 100%] (n=19) | 100% [84%, 100%] (n=19) | 100% [75%, 100%] (n=59) | 95% [64%, 100%] (n=35) | 100% [75%, 100%] (n=59) |
- line 28
  - numpy: | B | E2 all four schemes | 0% [0%, 0%] (n=37) | 0% [0%, 13%] (n=20) | 0% [0%, 22%] (n=20) | 0% [0%, 17%] (n=40) | 8% [0%, 17%] (n=19) | 8% [0%, 17%] (n=19) | 0% [0%, 17%] (n=59) | 0% [0%, 0%] (n=35) | 0% [0%, 17%] (n=59) |
  - exact: | B | E2 all four schemes | 0% [0%, 2%] (n=36) | 0% [0%, 13%] (n=20) | 0% [0%, 22%] (n=20) | 0% [0%, 17%] (n=40) | 8% [0%, 17%] (n=19) | 8% [0%, 17%] (n=19) | 0% [0%, 17%] (n=59) | 0% [0%, 5%] (n=35) | 0% [0%, 17%] (n=59) |
- line 29
  - numpy: | B | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.0] (n=37) | 0.0 [0.0, 2.5] (n=20) | 0.5 [0.0, 4.0] (n=20) | 0.0 [0.0, 4.0] (n=40) | 3.0 [1.0, 8.5] (n=19) | 3.0 [1.0, 8.5] (n=19) | 1.0 [0.0, 5.0] (n=59) | 0.0 [0.0, 1.0] (n=35) | 1.0 [0.0, 5.0] (n=59) |
  - exact: | B | Pairs significant in all schemes after adjustment | 0.0 [0.0, 1.5] (n=36) | 0.0 [0.0, 2.5] (n=20) | 0.5 [0.0, 4.0] (n=20) | 0.0 [0.0, 4.0] (n=40) | 3.0 [1.0, 8.5] (n=19) | 3.0 [1.0, 8.5] (n=19) | 1.0 [0.0, 5.0] (n=59) | 0.0 [0.0, 2.0] (n=35) | 1.0 [0.0, 5.0] (n=59) |
- line 30
  - numpy: | B | E3 covariance share carried by group means | 69% [42%, 86%] (n=37) | 90% [64%, 95%] (n=20) | 82% [33%, 93%] (n=20) | 85% [41%, 95%] (n=40) | 80% [31%, 97%] (n=19) | 80% [31%, 97%] (n=19) | 85% [40%, 95%] (n=59) | 69% [45%, 85%] (n=35) | 85% [40%, 95%] (n=59) |
  - exact: | B | E3 covariance share carried by group means | 69% [48%, 87%] (n=36) | 90% [64%, 95%] (n=20) | 82% [33%, 93%] (n=20) | 85% [41%, 95%] (n=40) | 80% [31%, 97%] (n=19) | 80% [31%, 97%] (n=19) | 85% [40%, 95%] (n=59) | 69% [48%, 85%] (n=35) | 85% [40%, 95%] (n=59) |
- line 31
  - numpy: | B | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.70 [0.40, 0.94] (n=27) | 0.60 [0.55, 0.80] (n=17) | 0.70 [0.51, 0.94] (n=15) | 0.66 [0.53, 0.87] (n=32) | 0.71 [0.46, 0.79] (n=15) | 0.71 [0.46, 0.79] (n=15) | 0.67 [0.50, 0.83] (n=47) | 0.70 [0.39, 0.96] (n=26) | 0.67 [0.50, 0.83] (n=47) |
  - exact: | B | Raw pooled-significant share among group-effect genes (expectation 0.61) | 0.69 [0.37, 0.91] (n=28) | 0.60 [0.55, 0.80] (n=17) | 0.70 [0.51, 0.94] (n=15) | 0.66 [0.53, 0.87] (n=32) | 0.71 [0.46, 0.79] (n=15) | 0.71 [0.46, 0.79] (n=15) | 0.67 [0.50, 0.83] (n=47) | 0.67 [0.36, 0.94] (n=27) | 0.67 [0.50, 0.83] (n=47) |
- line 32
  - numpy: | B | Series with share above 0.61 | 15/27 = 56% (37-72%) | 8/17 = 47% (26-69%) | 10/15 = 67% (42-85%) | 18/32 = 56% (39-72%) | 10/15 = 67% (42-85%) | 10/15 = 67% (42-85%) | 28/47 = 60% (45-72%) | 14/26 = 54% (35-71%) | 28/47 = 60% (45-72%) |
  - exact: | B | Series with share above 0.61 | 15/28 = 54% (36-70%) | 8/17 = 47% (26-69%) | 10/15 = 67% (42-85%) | 18/32 = 56% (39-72%) | 10/15 = 67% (42-85%) | 10/15 = 67% (42-85%) | 28/47 = 60% (45-72%) | 14/27 = 52% (34-69%) | 28/47 = 60% (45-72%) |
- line 33
  - numpy: | B | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 14/37 = 38% (24-54%) | 8/20 = 40% (22-61%) | 9/20 = 45% (26-66%) | 17/40 = 42% (29-58%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 31/59 = 53% (40-65%) | 13/35 = 37% (23-54%) | 31/59 = 53% (40-65%) |
  - exact: | B | Adjusted q < 0.05 count above its null (P < 0.05), raw scheme | 13/36 = 36% (22-52%) | 8/20 = 40% (22-61%) | 9/20 = 45% (26-66%) | 17/40 = 42% (29-58%) | 14/19 = 74% (51-88%) | 14/19 = 74% (51-88%) | 31/59 = 53% (40-65%) | 13/35 = 37% (23-54%) | 31/59 = 53% (40-65%) |
- line 34
  - numpy: | B | Median observed / null mean adjusted count, raw scheme | 4.64 [0.00, 31.25] (n=37) | 0.00 [0.00, 36.49] (n=20) | 5.38 [0.00, 45.86] (n=20) | 3.37 [0.00, 40.89] (n=40) | 18.87 [6.52, 79.70] (n=19) | 18.87 [6.52, 79.70] (n=19) | 8.77 [0.00, 50.96] (n=59) | 4.64 [0.00, 30.33] (n=35) | 8.77 [0.00, 50.96] (n=59) |
  - exact: | B | Median observed / null mean adjusted count, raw scheme | 4.24 [0.00, 29.87] (n=36) | 0.00 [0.00, 36.49] (n=20) | 5.38 [0.00, 45.86] (n=20) | 3.37 [0.00, 40.89] (n=40) | 18.87 [6.52, 79.70] (n=19) | 18.87 [6.52, 79.70] (n=19) | 8.77 [0.00, 50.96] (n=59) | 4.64 [0.00, 30.33] (n=35) | 8.77 [0.00, 50.96] (n=59) |
- line 35
  - numpy: | B | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 12/37 = 32% (20-49%) | 5/20 = 25% (11-47%) | 9/20 = 45% (26-66%) | 14/40 = 35% (22-50%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 21/59 = 36% (25-48%) | 11/35 = 31% (19-48%) | 21/59 = 36% (25-48%) |
  - exact: | B | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 1 scheme | 11/36 = 31% (18-47%) | 5/20 = 25% (11-47%) | 9/20 = 45% (26-66%) | 14/40 = 35% (22-50%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 21/59 = 36% (25-48%) | 11/35 = 31% (19-48%) | 21/59 = 36% (25-48%) |
- line 36
  - numpy: | B | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.61 [0.96, 5.60] (n=37) | 1.13 [0.43, 2.93] (n=20) | 2.50 [1.18, 4.38] (n=20) | 1.56 [0.78, 4.04] (n=40) | 1.24 [1.00, 2.53] (n=19) | 1.24 [1.00, 2.53] (n=19) | 1.51 [0.96, 3.77] (n=59) | 1.61 [0.94, 4.52] (n=35) | 1.51 [0.96, 3.77] (n=59) |
  - exact: | B | Median observed / null mean adjusted count, ratio to gene 1 scheme | 1.53 [0.95, 3.98] (n=36) | 1.13 [0.43, 2.93] (n=20) | 2.50 [1.18, 4.38] (n=20) | 1.56 [0.78, 4.04] (n=40) | 1.24 [1.00, 2.53] (n=19) | 1.24 [1.00, 2.53] (n=19) | 1.51 [0.96, 3.77] (n=59) | 1.61 [0.94, 4.52] (n=35) | 1.51 [0.96, 3.77] (n=59) |
- line 37
  - numpy: | B | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 8/37 = 22% (11-37%) | 5/20 = 25% (11-47%) | 11/20 = 55% (34-74%) | 16/40 = 40% (26-55%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 26/59 = 44% (32-57%) | 7/35 = 20% (10-36%) | 26/59 = 44% (32-57%) |
  - exact: | B | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 2 scheme | 7/36 = 19% (10-35%) | 5/20 = 25% (11-47%) | 11/20 = 55% (34-74%) | 16/40 = 40% (26-55%) | 10/19 = 53% (32-73%) | 10/19 = 53% (32-73%) | 26/59 = 44% (32-57%) | 7/35 = 20% (10-36%) | 26/59 = 44% (32-57%) |
- line 38
  - numpy: | B | Median observed / null mean adjusted count, ratio to gene 2 scheme | 1.08 [0.98, 1.81] (n=37) | 1.10 [0.98, 3.68] (n=20) | 2.54 [1.08, 20.11] (n=20) | 1.38 [1.01, 7.25] (n=40) | 1.92 [1.02, 9.19] (n=19) | 1.92 [1.02, 9.19] (n=19) | 1.40 [1.01, 7.40] (n=59) | 1.08 [0.99, 1.77] (n=35) | 1.40 [1.01, 7.40] (n=59) |
  - exact: | B | Median observed / null mean adjusted count, ratio to gene 2 scheme | 1.05 [0.98, 1.76] (n=36) | 1.10 [0.98, 3.68] (n=20) | 2.54 [1.08, 20.11] (n=20) | 1.38 [1.01, 7.25] (n=40) | 1.92 [1.02, 9.19] (n=19) | 1.92 [1.02, 9.19] (n=19) | 1.40 [1.01, 7.40] (n=59) | 1.08 [0.99, 1.77] (n=35) | 1.40 [1.01, 7.40] (n=59) |
- line 39
  - numpy: | B | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 7/37 = 19% (9-34%) | 5/20 = 25% (11-47%) | 8/20 = 40% (22-61%) | 13/40 = 32% (20-48%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 20/59 = 34% (23-47%) | 6/35 = 17% (8-33%) | 20/59 = 34% (23-47%) |
  - exact: | B | Adjusted q < 0.05 count above its null (P < 0.05), ratio to gene 3 scheme | 6/36 = 17% (8-32%) | 5/20 = 25% (11-47%) | 8/20 = 40% (22-61%) | 13/40 = 32% (20-48%) | 7/19 = 37% (19-59%) | 7/19 = 37% (19-59%) | 20/59 = 34% (23-47%) | 6/35 = 17% (8-33%) | 20/59 = 34% (23-47%) |
- line 40
  - numpy: | B | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.23 [0.84, 2.83] (n=37) | 0.99 [0.69, 2.30] (n=20) | 1.35 [0.99, 3.12] (n=20) | 1.15 [0.92, 2.83] (n=40) | 1.35 [0.96, 3.95] (n=19) | 1.35 [0.96, 3.95] (n=19) | 1.20 [0.94, 2.89] (n=59) | 1.23 [0.79, 2.81] (n=35) | 1.20 [0.94, 2.89] (n=59) |
  - exact: | B | Median observed / null mean adjusted count, ratio to gene 3 scheme | 1.17 [0.87, 2.80] (n=36) | 0.99 [0.69, 2.30] (n=20) | 1.35 [0.99, 3.12] (n=20) | 1.15 [0.92, 2.83] (n=40) | 1.35 [0.96, 3.95] (n=19) | 1.35 [0.96, 3.95] (n=19) | 1.20 [0.94, 2.89] (n=59) | 1.23 [0.86, 2.81] (n=35) | 1.20 [0.94, 2.89] (n=59) |

Note: summarize_validation.R labels the development arm "Development arrays (37)" whatever the number of series; under the exact rule it holds 36 series.

