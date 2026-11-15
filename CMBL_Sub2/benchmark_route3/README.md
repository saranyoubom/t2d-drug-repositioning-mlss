# Route 3: benchmark of the pooled co-expression pipeline on public perturbation datasets (4 Nov 2026)

Question: is the manuscript's finding (pooled consensus pairs are what the 3-group x 4-replicate design produces, and they vanish after treatment adjustment) specific to the INS-1 data or general?
Everything was pre-specified before any benchmark expression data were read: `PROTOCOL_prespecified.md` (fingerprints in `PROTOCOL_SHA256.txt`; frozen dataset list in `frozen_datasets.csv`, `FROZEN_SHA256.txt`). One deviation (D1: duplicate series excluded, list not extended) is recorded in the protocol before any correlation was computed.

## What was run
* Code: `bench_pipeline.py` = the manuscript pipeline for any 12 x G matrix; `validate_pipeline.py` reproduces every INS-1 number first (21 consensus pairs; 66/80/84/120 pooled; 0/51/11/15 adjusted; null 22.2 vs 22.1; group-mean share 0.97).
* Data: GEO microarray series on GPL570 (human) and GPL1261 (mouse), 12-36 samples, three or more groups with at least 4 samples (`select_datasets.py`, metadata only); `run_benchmark.py` (first a dry run that only counts DE genes). 150 series screened per platform; 64 eligible (25 human, 39 mouse); 58 processed in the frozen order, 37 analysable: **17 human, 20 mouse** (20 had fewer than 34 DE genes, 1 was a duplicate; the human list was exhausted at 17, the mouse list stopped at the target of 20). Per series: 4 replicates of 3 groups, 17 genes on the linear scale, 3 of them denominators.
* Gene set A (primary): 17 random differentially expressed genes (FDR < 0.05). Gene set B: 17 random expressed genes. Positive control: shared within-group latent factor injected into 6 genes of set A (within-group r about 0.7).

## Results (`summary_endpoints.json`, `per_dataset_results.csv`, `Benchmark_figure.png`)
| Endpoint (gene set A, 37 series) | Result |
|---|---|
| E1: pooled consensus count not above the permutation null (P >= 0.05) | **30 of 37 = 81% (95% CI 66-91%)**; 7 series (19%, 9-34%) above the null |
| Pre-specified band | point estimate in the ">= 80%: generally explained by the design" band; the lower CI bound (66%) lies in the "often explained" band |
| Positive control: injected within-group coupling flagged (P < 0.05) | **10 of 37 = 27% (15-43%)** vs 19% without injection; below the pre-specified 80%, so "not above the null" is **weak evidence** |
| E2: consensus pairs with q < 0.05 after treatment adjustment | in at least one scheme: median 60% (IQR 34-90%); **in all four schemes: median 0% (IQR 0-1.6%); pairs significant in all four: median 0 (IQR 0-1)** |
| E3: share of pooled covariance carried by group means (raw, significant pairs) | **median 97.6% (IQR 95.5-99.2%)** (INS-1: 97%) |
| Gene set B (random expressed genes) | 25 of 37 (68%) not above null; 12 above; E3 median 69% (IQR 42-86%) |
| By platform (set A, P >= 0.05) | human 13 of 17 (76%); mouse 17 of 20 (85%) |
INS-1 for comparison: 21 vs null 22.1 (P 0.69), 0 pairs in all four schemes after adjustment, share 0.97: typical of the benchmark.

## Reading
* Supported: in differentially expressed gene sets the pooled correlation is almost entirely group-mean covariance (about 98%), and recurrence across the four normalization schemes does not survive treatment adjustment (median 0 pairs). The INS-1 result is typical, not special.
* Not supported: reading "the consensus count is within the null" as evidence that there is no co-regulation. The count test is insensitive to injected within-group coupling (27% detection), and in 19% of series the count is above the null, so real within-group structure exists in some data. The adjusted (within-group) analysis, not the count, is the informative step.
* Limits: 37 microarray series from two platforms; designs chosen mechanically (some compare cell types, tissues or genotypes, e.g. GSE30164, GSE73265); DE genes are selected on the same 12 samples (as in a typical study); the primary estimate sits at the edge of the pre-specified band; only the first 4 replicates per group were used.
