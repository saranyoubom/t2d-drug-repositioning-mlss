# Benchmark kit: 96 public designs of 3 groups x 4 replicates and 17 genes

`series_manifest.csv` lists the 96 public series on which the study benchmarked its checks: 37 development arrays, 40 held-out arrays and 19 held-out RNA-seq series. For each series it gives the accession, the platform (GPL570 human arrays, GPL1261 mouse arrays) or organism (RNA-seq), the labels of the groups in GEO, the 12 samples of the design in order (3 groups x the first 4 samples of each group), and the 17 gene symbols of each of two gene sets. With it anyone can rebuild the same 12 x 17 matrices from GEO and run their own method on them.

**Not included:** expression data. GEO data are not redistributed here; cite GEO and the original series.

## How the designs were chosen
- Series were selected from the Gene Expression Omnibus on sample metadata alone (titles and characteristics, never expression values): single-platform series on GPL570 (human) or GPL1261 (mouse) with 12 to 36 samples and a perturbation key that defines at least three groups of at least four samples; the first four samples of each group were used. Some designs compare cell types, tissues or genotypes rather than treatments, because the selection was mechanical.
- The 37 development series were used to calibrate the checks. The held-out lists (40 arrays, 19 RNA-seq series with NCBI-generated raw counts) were frozen before any expression value was read; the mouse RNA-seq arm could not be run because no counts were available.
- Gene set A: 17 genes drawn at random from the expressed genes with a significant group effect (one-way ANOVA, Benjamini-Hochberg q < 0.05; at least 34 such genes were required for a series to be analysed). Gene set B: 17 random expressed genes. Expressed = mean log expression above the median. One probe per gene symbol (highest mean expression). Seeds: 2026 + the GSE number (gene set B: seed sequence [2026 + GSE number, 1]).
- Linear scale: arrays 2^(log2 expression); RNA-seq counts per million of the genes with a count of at least 1 in all 12 samples.

## Columns of the manifest
`arm`, `accession`, `platform_or_organism`, `frozen_list_order` (position in the frozen candidate list), `group_labels_in_geo`, `samples_in_design_order` (12 GSM ids, semicolon separated; groups are consecutive blocks of 4), `gene_set_A_symbols`, `gene_set_B_symbols` (17 symbols each, semicolon separated), `n_de_genes` (genes with a significant group effect in the series), `scale_in_file`.

## Rebuilding the matrices
The loaders of the study are in the study repository (`CMBL_Sub2/strengthen_v2/wp2_benchmark2/`: `bench_data_io.py` for arrays, `run_rnaseq.py` for RNA-seq). They download the series matrix or count file from NCBI, pick the 12 samples, map probes to gene symbols and draw the gene sets with the seeds above. `build_manifest_provenance.py` shows how this manifest was assembled from the study's frozen lists and stored results; it refers to the study's folders and is not meant to be run elsewhere.

## Limits
One design (3 groups x 4 replicates), random gene sets and 17 genes; human and mouse arrays and human RNA-seq only (no rat and no mouse RNA-seq); results on this kit say how often a check fires in such designs and with what power in simulations, not whether a particular network is right.
