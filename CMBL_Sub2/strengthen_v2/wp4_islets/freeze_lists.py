"""freeze_lists.py - build the frozen dataset table and the frozen beta-cell lists of protocol 4.3 from METADATA ONLY (GEO series matrices, ArrayExpress SDRF), before any expression value of
these datasets is read. Output (beside this script): frozen_islet_datasets.csv (order, accession, decision, label source, beta cells, donors, reason) and frozen_beta_cells.csv (one row per author-labelled
beta cell: dataset, cell_key = join key to the expression file, gsm, donor, disease, label, label_source). The caller fingerprints both files (sha256) into ../PROTOCOL_v2_SHA256.txt.
Label rules (fixed here): the label is the cell-type annotation written by the data authors in the public record; no label is derived from expression by us; cells without that label are not beta-cells.
Usage: python -I freeze_lists.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd
from read_geo_meta import sample_table

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('V2_DATA_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_cache')) + os.sep + 'v2_data' + os.sep + 'wp4'
rows = []; summary = []

# 1. GSE86469 (Lawlor 2017): 'cell type' = Beta in the GEO sample characteristics; join key = sample title (column name of the expression csv)
T = sample_table(os.path.join(ROOT, 'GSE86469', 'GSE86469_series_matrix.txt.gz'))
B = T[T['ch:cell type'] == 'Beta']
for _, r in B.iterrows():
    rows.append(dict(dataset='GSE86469', cell_key=r['title'], gsm=r['geo_accession'], donor=r['ch:islet unos id'], disease='T2D' if r['ch:disease'].startswith('Type 2') else 'non-diabetic', label='Beta', label_source='GEO sample characteristic "cell type" (authors)'))
summary.append(dict(order=1, dataset='GSE86469', source='GEO ftp (series matrix + RSEM expected counts csv.gz)', decision='attempt', label_source='GEO sample characteristic "cell type" (authors)', beta_cells=len(B), donors=B['ch:islet unos id'].nunique(), reason=''))

# 2. GSE81547 (Enge 2017): 'inferred_cell_type' = beta in the GEO sample characteristics; join key = GSM accession (file name prefix inside GSE81547_RAW.tar)
T = sample_table(os.path.join(ROOT, 'GSE81547', 'GSE81547_series_matrix.txt.gz'))
B = T[T['ch:inferred_cell_type'] == 'beta']
for _, r in B.iterrows():
    rows.append(dict(dataset='GSE81547', cell_key=r['geo_accession'], gsm=r['geo_accession'], donor='age%s_%s' % (r['ch:donor_age'], r['ch:gender']), disease='not stated (series summary: previously healthy donors)', label='beta',
                     label_source='GEO sample characteristic "inferred_cell_type" (authors)'))
summary.append(dict(order=2, dataset='GSE81547', source='GEO ftp (GSE81547_RAW.tar, 176 MB, per-cell csv.gz)', decision='attempt', label_source='GEO sample characteristic "inferred_cell_type" (authors)', beta_cells=len(B),
                    donors=len(set('age%s_%s' % (a, g) for a, g in zip(B['ch:donor_age'], B['ch:gender']))), reason=''))

# 3. E-MTAB-5061 (Segerstolpe 2016): 'inferred cell type' = beta cell in the SDRF; join key = Source Name; only cells with single cell quality OK carry an inferred type
S = pd.read_csv(os.path.join(ROOT, 'E-MTAB-5061', 'E-MTAB-5061.sdrf.txt'), sep='\t', dtype=str, low_memory=False)
B = S[S['Characteristics [inferred cell type]'] == 'beta cell']
for _, r in B.iterrows():
    rows.append(dict(dataset='E-MTAB-5061', cell_key=r['Source Name'], gsm='', donor=r['Characteristics [individual]'], disease='T2D' if r['Characteristics [disease]'].startswith('type II') else 'non-diabetic', label='beta cell',
                     label_source='ArrayExpress SDRF "inferred cell type" (authors)'))
summary.append(dict(order=3, dataset='E-MTAB-5061', source='ArrayExpress / BioStudies (processed.1.zip 159 MB holds the 723 MB rpkm+counts text; the zip is streamed, the 723 MB text is not downloaded)', decision='attempt',
                    label_source='ArrayExpress SDRF "inferred cell type" (authors)', beta_cells=len(B), donors=B['Characteristics [individual]'].nunique(), reason=''))

# 4. GSE85241 (Muraro 2016): the GEO record has 32 library-level samples and no per-cell annotation; clusters were found by the authors' StemID and are not in the GEO files
summary.append(dict(order=4, dataset='GSE85241', source='GEO ftp (counts csv.gz 17 MB, barcode list, demultiplexing readme)', decision='skipped before expression download',
                    label_source='none in the public record (32 library-level GEO samples; supplementary files: counts, barcodes, readme)', beta_cells=0, donors=0,
                    reason='no author-provided beta-cell labels in the GEO or repository files; protocol 4.3 requires author-provided labels and plain public download, and we do not assign cell types ourselves'))
F = pd.DataFrame(rows); F.to_csv(os.path.join(HERE, 'frozen_beta_cells.csv'), index=False)
P = pd.DataFrame(summary); P.to_csv(os.path.join(HERE, 'frozen_islet_datasets.csv'), index=False)
print(P.to_string(index=False, max_colwidth=60))
print(F.groupby(['dataset', 'donor', 'disease']).size().to_string())
