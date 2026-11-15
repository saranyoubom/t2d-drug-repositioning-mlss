"""fetch_expression.py - download the EXPRESSION files of the additional human islet datasets of protocol 4.3, one source per scratch folder, after the dataset table and beta-cell lists were frozen
(frozen_islet_datasets.csv, sha256 in ../PROTOCOL_v2_SHA256.txt). Every attempt, size and outcome is appended to attempts.csv. Sizes are checked with HEAD first; nothing above 300 MB is downloaded
(the 723 MB E-MTAB-5061 text file is NOT downloaded: its 159 MB zip is, and is streamed later). Files are untrusted data (never executed).
Usage: python -I fetch_expression.py GSE86469 | GSE81547 | E-MTAB-5061"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import wp4_http as H
from attempts_log import record

ROOT = os.environ.get('V2_DATA_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_cache')) + os.sep + 'v2_data' + os.sep + 'wp4'
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
FILES = {
    'GSE86469': [('RSEM expected counts csv.gz', GEO + 'GSE86nnn/GSE86469/suppl/GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz', 'GSE86469/GSE86469_RSEM_expected_counts.csv.gz', 'GEO sample characteristic "cell type" (authors)')],
    'GSE81547': [('RAW.tar (2544 per-cell csv.gz)', GEO + 'GSE81nnn/GSE81547/suppl/GSE81547_RAW.tar', 'GSE81547/GSE81547_RAW.tar', 'GEO sample characteristic "inferred_cell_type" (authors)')],
    'E-MTAB-5061': [('processed.1.zip (rpkm + counts text, 723 MB unzipped)', 'https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-5061/E-MTAB-5061.processed.1.zip', 'E-MTAB-5061/E-MTAB-5061.processed.1.zip',
                     'ArrayExpress SDRF "inferred cell type" (authors)')],
}
ds = sys.argv[1]
for item, url, rel, lab in FILES[ds]:
    path = os.path.join(ROOT, rel.replace('/', os.sep)); h = H.head(url)
    print(ds, item, 'HEAD', h)
    try:
        n = H.download(url, path, max_mb=300)
        record(ds, item, url, n, lab, 'downloaded', 'size checked by HEAD before download (limit 300 MB)'); print('ok', ds, item, n)
    except Exception as e:
        record(ds, item, url, h.get('size'), lab, 'failed or refused', repr(e)[:300]); print('FAILED', ds, item, repr(e)[:300])
