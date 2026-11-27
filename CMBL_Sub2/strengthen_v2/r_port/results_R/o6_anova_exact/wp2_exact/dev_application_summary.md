# Development series: protocol 3.3 and 3.4 (WP2 deliverable C)

## Gene set A (17 random differentially expressed genes, FDR < 0.05)

| Endpoint | Result (37 series) |
|---|---|
| Within-group coupling, T_r2 (primary), P < 0.05 | 22/36 = 61% (45-75%) |
| T_p, P < 0.05 | 20/36 = 56% (40-70%) |
| T_eig, P < 0.05 | 19/36 = 53% (37-68%) |
| Consensus count not above its permutation null (E1; P >= 0.05), this run | 30/36 = 83% (68-92%) |
| same, part 1 (loop null) | 30/37 = 81% (66-91%) |
| E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 36) | median 0% [0%, 3%] (n=36) |
| E2 any scheme (part 1 definition), this run / part 1 | 72% [32%, 93%] (n=36) / 60% [34%, 90%] (n=37) |
| E2 all four schemes, this run / part 1 | 0% [0%, 0%] (n=36) / 0% [0%, 2%] (n=37) |
| Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.0] (n=36) / 0.0 [0.0, 1.0] (n=37) |
| E3 share of pooled covariance carried by group means | 98% [96%, 99%] (n=36) |
| Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.81 [0.65, 0.91] (n=36); above 0.61 in 30/36 = 83% (68-92%) |

Ratio-scheme diagnostics (observed adjusted q < 0.05 count against its null from 500 within-group permutations):

| Scheme | Series with P < 0.05 | Median observed | Median null mean | Median observed / null mean | Median excess |
|---|---|---|---|---|---|
| raw | 9/36 = 25% (14-41%) | 0.0 [0.0, 2.0] (n=36) | 0.1 [0.1, 0.2] (n=36) | 0.00 [0.00, 13.20] (n=36) (defined in 36) | -0.1 [-0.1, 1.7] (n=36) |
| ratio to gene 1 | 8/36 = 22% (12-38%) | 12.0 [1.8, 52.0] (n=36) | 4.3 [0.7, 46.8] (n=36) | 1.02 [0.71, 2.99] (n=36) (defined in 36) | 0.7 [-0.7, 4.0] (n=36) |
| ratio to gene 2 | 9/36 = 25% (14-41%) | 13.5 [3.0, 39.2] (n=36) | 9.6 [1.4, 26.2] (n=36) | 0.99 [0.50, 2.33] (n=36) (defined in 36) | -0.2 [-1.8, 5.6] (n=36) |
| ratio to gene 3 | 3/36 = 8% (3-22%) | 16.5 [4.8, 57.5] (n=36) | 25.5 [3.0, 54.5] (n=36) | 1.01 [0.43, 1.54] (n=36) (defined in 36) | 0.4 [-5.7, 5.4] (n=36) |

By platform, human GPL570: T_r2 P < 0.05 in 11/16 = 69% (44-86%); consensus count not above null in 13/16 = 81% (57-93%).
By platform, mouse GPL1261: T_r2 P < 0.05 in 11/20 = 55% (34-74%); consensus count not above null in 17/20 = 85% (64-95%).
Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 22/35 = 63% (46-77%); T_p in 20/35 = 57% (41-72%); consensus count not above null in 29/35 = 83% (67-92%); E2 raw scheme alone 0% [0%, 3%] (n=35); covariance share carried by group means 98% [96%, 99%] (n=35).

## Gene set B (17 random expressed genes)

| Endpoint | Result (37 series) |
|---|---|
| Within-group coupling, T_r2 (primary), P < 0.05 | 28/36 = 78% (62-88%) |
| T_p, P < 0.05 | 27/36 = 75% (59-86%) |
| T_eig, P < 0.05 | 26/36 = 72% (56-84%) |
| Consensus count not above its permutation null (E1; P >= 0.05), this run | 24/36 = 67% (50-80%) |
| same, part 1 (loop null) | 25/37 = 68% (51-80%) |
| E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 36) | median 0% [0%, 21%] (n=36) |
| E2 any scheme (part 1 definition), this run / part 1 | 98% [65%, 100%] (n=36) / 100% [67%, 100%] (n=37) |
| E2 all four schemes, this run / part 1 | 0% [0%, 2%] (n=36) / 0% [0%, 0%] (n=37) |
| Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.5] (n=36) / 0.0 [0.0, 1.0] (n=37) |
| E3 share of pooled covariance carried by group means | 69% [48%, 87%] (n=36) |
| Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.69 [0.37, 0.91] (n=28); above 0.61 in 15/28 = 54% (36-70%) |

Ratio-scheme diagnostics (observed adjusted q < 0.05 count against its null from 500 within-group permutations):

| Scheme | Series with P < 0.05 | Median observed | Median null mean | Median observed / null mean | Median excess |
|---|---|---|---|---|---|
| raw | 13/36 = 36% (22-52%) | 1.0 [0.0, 3.8] (n=36) | 0.1 [0.1, 0.2] (n=36) | 4.24 [0.00, 29.87] (n=36) (defined in 36) | 0.8 [-0.1, 3.7] (n=36) |
| ratio to gene 1 | 11/36 = 31% (18-47%) | 17.5 [2.0, 67.2] (n=36) | 5.9 [1.1, 45.2] (n=36) | 1.53 [0.95, 3.98] (n=36) (defined in 36) | 2.1 [-0.3, 10.2] (n=36) |
| ratio to gene 2 | 7/36 = 19% (10-35%) | 14.0 [3.8, 67.0] (n=36) | 10.3 [2.4, 53.3] (n=36) | 1.05 [0.98, 1.76] (n=36) (defined in 36) | 1.3 [-0.1, 7.1] (n=36) |
| ratio to gene 3 | 6/36 = 17% (8-32%) | 20.0 [5.5, 49.5] (n=36) | 8.6 [2.1, 33.2] (n=36) | 1.17 [0.87, 2.80] (n=36) (defined in 36) | 2.1 [-1.5, 12.0] (n=36) |

By platform, human GPL570: T_r2 P < 0.05 in 12/16 = 75% (51-90%); consensus count not above null in 11/16 = 69% (44-86%).
By platform, mouse GPL1261: T_r2 P < 0.05 in 16/20 = 80% (58-92%); consensus count not above null in 13/20 = 65% (43-82%).
Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 28/35 = 80% (64-90%); T_p in 27/35 = 77% (61-88%); consensus count not above null in 23/35 = 66% (49-79%); E2 raw scheme alone 0% [0%, 22%] (n=35); covariance share carried by group means 69% [48%, 85%] (n=35).
