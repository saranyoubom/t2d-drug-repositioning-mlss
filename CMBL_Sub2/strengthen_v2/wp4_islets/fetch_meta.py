"""fetch_meta.py - download ONLY the small metadata files of the four additional human islet datasets (protocol 4.3) into one scratch folder per source, and record every attempt in attempts.csv.
Files: GSE86469 series matrix; GSE81547 series matrix + filelist.txt; GSE85241 series matrix + barcode file; E-MTAB-5061 sdrf + idf. No expression data are fetched here (expression is downloaded only after the label files
have been read and the dataset / beta-cell lists are frozen). Folders: <data root>/wp4/<source>/ (untrusted data, never executed).
Usage: python -I fetch_meta.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import wp4_http as H
from attempts_log import record

ROOT = os.environ.get('V2_DATA_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_cache')) + os.sep + 'v2_data' + os.sep + 'wp4'
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
EBI = 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/061/E-MTAB-5061/Files/'
ITEMS = [
    ('GSE86469', 'series matrix', GEO + 'GSE86nnn/GSE86469/matrix/GSE86469_series_matrix.txt.gz', 'GSE86469/GSE86469_series_matrix.txt.gz', 'sample characteristics in the GEO series matrix'),
    ('GSE81547', 'series matrix', GEO + 'GSE81nnn/GSE81547/matrix/GSE81547_series_matrix.txt.gz', 'GSE81547/GSE81547_series_matrix.txt.gz', 'sample characteristics in the GEO series matrix'),
    ('GSE81547', 'suppl filelist', GEO + 'GSE81nnn/GSE81547/suppl/filelist.txt', 'GSE81547/filelist.txt', 'file names only'),
    ('GSE85241', 'series matrix', GEO + 'GSE85nnn/GSE85241/matrix/GSE85241_series_matrix.txt.gz', 'GSE85241/GSE85241_series_matrix.txt.gz', 'sample characteristics in the GEO series matrix'),
    ('GSE85241', 'barcode file', GEO + 'GSE85nnn/GSE85241/suppl/GSE85241_cel-seq_barcodes.csv.gz', 'GSE85241/GSE85241_cel-seq_barcodes.csv.gz', 'none (barcode list)'),
    ('E-MTAB-5061', 'sdrf', EBI + 'E-MTAB-5061.sdrf.txt', 'E-MTAB-5061/E-MTAB-5061.sdrf.txt', 'SDRF characteristics (authors, ArrayExpress)'),
    ('E-MTAB-5061', 'idf', EBI + 'E-MTAB-5061.idf.txt', 'E-MTAB-5061/E-MTAB-5061.idf.txt', 'IDF (description)'),
]
for ds, item, url, rel, lab in ITEMS:
    path = os.path.join(ROOT, rel.replace('/', os.sep))
    try:
        n = H.download(url, path, max_mb=50)
        record(ds, item, url, n, lab, 'downloaded (metadata only)', ''); print('ok', ds, item, n)
    except Exception as e:
        record(ds, item, url, None, lab, 'failed', repr(e)[:300]); print('FAILED', ds, item, repr(e)[:200])
