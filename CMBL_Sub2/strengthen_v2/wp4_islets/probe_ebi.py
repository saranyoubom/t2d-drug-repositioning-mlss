"""probe_ebi.py - read-only reconnaissance of the public download location of E-MTAB-5061 (Segerstolpe 2016) at EBI: tries the ArrayExpress/BioStudies ftp index over https and HEADs the files.
Nothing is downloaded. Usage: python -I probe_ebi.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import wp4_http as H

CANDIDATE_INDEXES = [
    'https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/MTAB/E-MTAB-5061/',
    'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/061/E-MTAB-5061/Files/',
    'https://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/061/E-MTAB-5061/Files/',
]
for url in CANDIDATE_INDEXES:
    print('INDEX', url)
    try:
        items = H.list_index(url)
    except Exception as e:
        print('   listing failed:', repr(e)); continue
    for name, href in items:
        if name.lower().startswith('parent'): continue
        h = H.head(url + href)
        print(f'   {href}\t{h["status"]}\t{h["size"]}\t{h["type"]}')
for u in ['https://www.ebi.ac.uk/biostudies/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt',
          'https://www.ebi.ac.uk/biostudies/files/E-MTAB-5061/E-MTAB-5061.idf.txt',
          'https://www.ebi.ac.uk/biostudies/files/E-MTAB-5061/pancreas_refseq_rpkms_counts_3514sc.txt']:
    print('HEAD', u, H.head(u))
