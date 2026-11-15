"""parse_gse81608.py - stream the GSE81608 RPKM table (Xin 2016; 39851 genes x 1600 cells) and keep ONLY the 503 beta-cells (type 'beta' in the GEO sample records) as a sparse float64 matrix, with the number of
detected genes per cell over ALL genes (the depth covariate for RPKM data, protocol 4.2). The same pandas C parser as the original script reads the file in chunks of rows, so the stored values are the same numbers.
Output: <PARSED>/gse81608/beta.npz (CSR, cells x genes, gene labels = Entrez ids) and meta.json (genes, donor, disease, tot = sum of RPKM, ndet, cell ids).
Usage: python -I parse_gse81608.py <dir with series_matrix.txt.gz and human_islets_rpkm.txt.gz> <output dir>"""
import sys, os, gzip, json
import numpy as np, pandas as pd
from scipy import sparse
sys.stdout.reconfigure(encoding='utf8')
D, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
rows = {}
with gzip.open(os.path.join(D, 'series_matrix.txt.gz'), 'rt', encoding='utf8', errors='replace') as f:
    for l in f:
        if l.startswith(('!Sample_title', '!Sample_characteristics_ch1')): rows.setdefault(l.split('\t')[0], []).append([x.strip('"') for x in l.rstrip('\n').split('\t')[1:]])
ch = rows['!Sample_characteristics_ch1']; title = rows['!Sample_title'][0]
col = lambda i: [x.split(': ', 1)[1] for x in ch[i]]
meta = pd.DataFrame({'sample': ['Sample_' + t.split()[-1] for t in title], 'donor': col(1), 'condition': col(2), 'type': col(6)})
b = meta[meta.type == 'beta'].reset_index(drop=True)
keep = []
for k, part in enumerate(pd.read_csv(os.path.join(D, 'human_islets_rpkm.txt.gz'), sep='\t', index_col=0, chunksize=2500)):
    if k == 0: assert list(part.columns) == list(meta['sample']), 'column order differs from the metadata'
    keep.append(part[list(b['sample'])].astype('float64'))
Xb = pd.concat(keep); del keep
A = Xb.to_numpy()                                                                      # genes x beta-cells
ndet = (A > 0).sum(axis=0); tot = A.sum(axis=0)
B = sparse.csr_matrix(A.T); del A
sparse.save_npz(os.path.join(OUT, 'beta.npz'), B)
json.dump(dict(name='GSE81608', species='human', kind='rpkm', id_type='entrez', genes=[int(g) for g in Xb.index], donor=b['donor'].tolist(), disease=b['condition'].tolist(), tot=tot.tolist(), ndet=ndet.tolist(),
               cell_ids=b['sample'].tolist()), open(os.path.join(OUT, 'meta.json'), 'w'))
print('GSE81608 beta cells', B.shape[0], 'genes', B.shape[1], 'nonzero', B.nnz, '| donors', b.donor.nunique(), '| T2D cells', int((b.condition == 'T2D').sum()), '| median detected genes %d' % np.median(ndet))
