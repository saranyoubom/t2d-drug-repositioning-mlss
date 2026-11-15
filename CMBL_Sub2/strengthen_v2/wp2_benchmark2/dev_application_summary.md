# Development series: protocol 3.3 and 3.4 (WP2 deliverable C)

## Gene set A (17 random differentially expressed genes, FDR < 0.05)

| Endpoint | Result (37 series) |
|---|---|
| Within-group coupling, T_r2 (primary), P < 0.05 | 25/37 = 68% (51-80%) |
| T_p, P < 0.05 | 21/37 = 57% (41-71%) |
| T_eig, P < 0.05 | 21/37 = 57% (41-71%) |
| Consensus count not above its permutation null (E1; P >= 0.05), this run | 32/37 = 86% (72-94%) |
| same, part 1 (loop null) | 30/37 = 81% (66-91%) |
| E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 37) | median 0% [0%, 6%] (n=37) |
| E2 any scheme (part 1 definition), this run / part 1 | 60% [24%, 88%] (n=37) / 60% [34%, 90%] (n=37) |
| E2 all four schemes, this run / part 1 | 0% [0%, 0%] (n=37) / 0% [0%, 2%] (n=37) |
| Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.0] (n=37) / 0.0 [0.0, 1.0] (n=37) |
| E3 share of pooled covariance carried by group means | 98% [96%, 99%] (n=37) |
| Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.81 [0.62, 0.91] (n=35); above 0.61 in 27/35 = 77% (61-88%) |

Ratio-scheme diagnostics (observed adjusted q < 0.05 count against its null from 500 within-group permutations):

| Scheme | Series with P < 0.05 | Median observed | Median null mean | Median observed / null mean | Median excess |
|---|---|---|---|---|---|
| raw | 10/37 = 27% (15-43%) | 0.0 [0.0, 2.0] (n=37) | 0.1 [0.1, 0.2] (n=37) | 0.00 [0.00, 13.16] (n=37) (defined in 37) | -0.0 [-0.1, 1.9] (n=37) |
| ratio to gene 1 | 8/37 = 22% (11-37%) | 6.0 [2.0, 49.0] (n=37) | 5.0 [0.7, 27.3] (n=37) | 1.01 [0.80, 2.51] (n=37) (defined in 37) | 0.6 [-0.7, 3.5] (n=37) |
| ratio to gene 2 | 11/37 = 30% (17-46%) | 12.0 [3.0, 38.0] (n=37) | 8.2 [2.5, 23.9] (n=37) | 0.99 [0.39, 2.49] (n=37) (defined in 37) | -0.2 [-1.9, 7.2] (n=37) |
| ratio to gene 3 | 5/37 = 14% (6-28%) | 14.0 [4.0, 51.0] (n=37) | 13.9 [2.2, 49.7] (n=37) | 1.03 [0.45, 1.88] (n=37) (defined in 37) | 0.6 [-6.8, 6.8] (n=37) |

By platform, human GPL570: T_r2 P < 0.05 in 12/17 = 71% (47-87%); consensus count not above null in 14/17 = 82% (59-94%).
By platform, mouse GPL1261: T_r2 P < 0.05 in 13/20 = 65% (43-82%); consensus count not above null in 18/20 = 90% (70-97%).
Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 23/35 = 66% (49-79%); T_p in 20/35 = 57% (41-72%); consensus count not above null in 30/35 = 86% (71-94%); E2 raw scheme alone 0% [0%, 5%] (n=35); covariance share carried by group means 98% [96%, 99%] (n=35).

## Gene set B (17 random expressed genes)

| Endpoint | Result (37 series) |
|---|---|
| Within-group coupling, T_r2 (primary), P < 0.05 | 29/37 = 78% (63-89%) |
| T_p, P < 0.05 | 28/37 = 76% (60-87%) |
| T_eig, P < 0.05 | 27/37 = 73% (57-85%) |
| Consensus count not above its permutation null (E1; P >= 0.05), this run | 24/37 = 65% (49-78%) |
| same, part 1 (loop null) | 25/37 = 68% (51-80%) |
| E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = 37) | median 0% [0%, 20%] (n=37) |
| E2 any scheme (part 1 definition), this run / part 1 | 100% [67%, 100%] (n=37) / 100% [67%, 100%] (n=37) |
| E2 all four schemes, this run / part 1 | 0% [0%, 0%] (n=37) / 0% [0%, 0%] (n=37) |
| Pairs significant in all schemes after adjustment, this run / part 1 | 0.0 [0.0, 1.0] (n=37) / 0.0 [0.0, 1.0] (n=37) |
| E3 share of pooled covariance carried by group means | 69% [42%, 86%] (n=37) |
| Raw pooled-significant share among genes with a group effect (expectation 0.61) | 0.70 [0.40, 0.94] (n=27); above 0.61 in 15/27 = 56% (37-72%) |

Ratio-scheme diagnostics (observed adjusted q < 0.05 count against its null from 500 within-group permutations):

| Scheme | Series with P < 0.05 | Median observed | Median null mean | Median observed / null mean | Median excess |
|---|---|---|---|---|---|
| raw | 14/37 = 38% (24-54%) | 1.0 [0.0, 6.0] (n=37) | 0.1 [0.1, 0.2] (n=37) | 4.64 [0.00, 31.25] (n=37) (defined in 37) | 0.8 [-0.1, 5.9] (n=37) |
| ratio to gene 1 | 12/37 = 32% (20-49%) | 19.0 [2.0, 68.0] (n=37) | 7.8 [1.2, 50.8] (n=37) | 1.61 [0.96, 5.60] (n=37) (defined in 37) | 2.5 [-0.1, 12.0] (n=37) |
| ratio to gene 2 | 8/37 = 22% (11-37%) | 29.0 [4.0, 67.0] (n=37) | 10.0 [2.5, 52.7] (n=37) | 1.08 [0.98, 1.81] (n=37) (defined in 37) | 1.3 [-0.1, 7.3] (n=37) |
| ratio to gene 3 | 7/37 = 19% (9-34%) | 20.0 [4.0, 49.0] (n=37) | 7.0 [1.8, 22.7] (n=37) | 1.23 [0.84, 2.83] (n=37) (defined in 37) | 2.3 [-0.3, 12.5] (n=37) |

By platform, human GPL570: T_r2 P < 0.05 in 13/17 = 76% (53-90%); consensus count not above null in 11/17 = 65% (41-83%).
By platform, mouse GPL1261: T_r2 P < 0.05 in 16/20 = 80% (58-92%); consensus count not above null in 13/20 = 65% (43-82%).
Post hoc sensitivity, without series with degenerate genes (post hoc sensitivity) (n = 35): T_r2 P < 0.05 in 28/35 = 80% (64-90%); T_p in 27/35 = 77% (61-88%); consensus count not above null in 23/35 = 66% (49-79%); E2 raw scheme alone 0% [0%, 17%] (n=35); covariance share carried by group means 69% [45%, 85%] (n=35).
