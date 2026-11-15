"""parse_new_datasets.py - stream the expression files of the additional human islet datasets of protocol 4.3 and keep ONLY the author-labelled beta-cells of the frozen list (frozen_beta_cells.csv), written to the
unified sparse format of wp4_core.load_parsed: <PARSED>/<tag>/beta.npz (CSR, cells x genes) + meta.json (genes, donor, disease, tot, ndet, cell_ids, note).
Per dataset (fixed in PRESPEC_WP4_choices.md before any expression value was read):
  GSE86469    RSEM expected counts csv.gz (rows = Ensembl gene ids, columns = cell names)              kind 'count', id_type 'ensembl'; tot = column sum over all rows, ndet = rows > 0
  GSE81547    RAW.tar of per-cell htseq count tables (gene symbol, count); ERCC and the 5 htseq special rows (no_feature, ambiguous, too_low_aQual, not_aligned, alignment_not_unique) are dropped before tot / ndet
              kind 'count', id_type 'symbol'
  E-MTAB-5061 see parse_emtab5061() (the count block of the authors' rpkm+counts text file inside processed.1.zip; the 723 MB text is streamed from the zip, never written out)
No expression-derived label is created; cells not on the frozen list are discarded while streaming. Data are untrusted: only read, never executed. Memory: one dense block of (genes x beta-cells) at a time.
Usage: python -I parse_new_datasets.py GSE86469 | GSE81547 | E-MTAB-5061"""
import sys, os, gzip, json, tarfile, zipfile, io, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
from scipy import sparse
import wp4_core as C

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(C.PARSED)
FROZEN = pd.read_csv(os.path.join(HERE, 'frozen_beta_cells.csv'), dtype=str)
SPECIAL = {'no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned', 'alignment_not_unique'}


def save(tag, name, kind, id_type, genes, A, F, note):
    """A: dense float64 (genes x cells) in the order of the frozen rows F; tot and ndet are taken over all rows of A"""
    out = os.path.join(C.PARSED, tag); os.makedirs(out, exist_ok=True)
    tot, ndet = A.sum(axis=0), (A > 0).sum(axis=0); B = sparse.csr_matrix(A.T)
    sparse.save_npz(os.path.join(out, 'beta.npz'), B)
    json.dump(dict(name=name, species='human', kind=kind, id_type=id_type, genes=[str(g) for g in genes], donor=F.donor.tolist(), disease=F.disease.tolist(), tot=tot.tolist(), ndet=ndet.tolist(),
                   cell_ids=F.cell_key.tolist(), note=note), open(os.path.join(out, 'meta.json'), 'w'))
    print(name, 'beta cells', B.shape[0], 'rows kept', B.shape[1], 'nonzero', B.nnz, '| donors', F.donor.nunique(), '| median detected genes %d' % np.median(ndet), '| median library size %.0f' % np.median(tot), flush=True)


def parse_gse86469():
    F = FROZEN[FROZEN.dataset == 'GSE86469'].reset_index(drop=True); keys = list(F.cell_key); parts = []
    for k, part in enumerate(pd.read_csv(os.path.join(ROOT, 'GSE86469', 'GSE86469_RSEM_expected_counts.csv.gz'), index_col=0, chunksize=4000)):
        if k == 0: assert set(keys) <= set(part.columns), 'frozen cell names missing from the file: %s' % sorted(set(keys) - set(part.columns))[:5]
        parts.append(part[keys].astype('float64'))
    X = pd.concat(parts); genes = [str(g) for g in X.index]; nonens = [g for g in genes if not g.startswith('ENSG')]
    assert not nonens, 'rows that are not Ensembl genes: %s' % nonens[:5]
    assert len(set(genes)) == len(genes)
    save('gse86469', 'GSE86469', 'count', 'ensembl', genes, X.to_numpy(), F, 'RSEM expected counts; library size = sum over all %d Ensembl rows' % len(genes))


def parse_gse81547():
    F = FROZEN[FROZEN.dataset == 'GSE81547'].reset_index(drop=True); want = set(F.gsm); cols = {}; genes0 = None
    with tarfile.open(os.path.join(ROOT, 'GSE81547', 'GSE81547_RAW.tar'), 'r|') as t:
        for m in t:
            gsm = m.name.split('_')[0]
            if gsm not in want: continue
            lines = gzip.decompress(t.extractfile(m).read()).decode('utf8', 'replace').splitlines(); names = [l.split('\t')[0] for l in lines]
            if genes0 is None: genes0 = names
            assert names == genes0, 'row order differs in ' + m.name
            cols[gsm] = np.array([float(l.split('\t')[1]) for l in lines])
    assert set(cols) == want, 'frozen GSM missing from the tar: %s' % sorted(want - set(cols))[:5]
    keep = [i for i, n in enumerate(genes0) if not n.startswith('ERCC') and n not in SPECIAL]
    A = np.column_stack([cols[g] for g in F.gsm])[keep]; genes = [genes0[i] for i in keep]
    assert len(set(genes)) == len(genes)
    save('gse81547', 'GSE81547', 'count', 'symbol', genes, A, F, 'htseq read counts; ERCC (%d rows) and the 5 htseq special rows dropped; library size = sum over %d gene rows' % (len(genes0) - len(keep) - sum(n in SPECIAL for n in genes0), len(genes)))


def parse_emtab5061():
    """layout documented in the IDF (E-MTAB-5061.idf.txt): header '#samples' + 3514 sample ids; rows = gene symbol, RefSeq id, 3514 RPKM, 3514 read counts (same order); rows 1-26179 RefSeq genes, then 92 ERCC, then eGFP.
    The count block of the 26179 RefSeq rows of the frozen beta-cells is kept. Identity check on the full sample set (INS in beta cells, GCG in alpha cells, against all other author-labelled types) before saving."""
    F = FROZEN[FROZEN.dataset == 'E-MTAB-5061'].reset_index(drop=True); NS, NG = 3514, 26179
    S = pd.read_csv(os.path.join(ROOT, 'E-MTAB-5061', 'E-MTAB-5061.sdrf.txt'), sep='\t', dtype=str, low_memory=False).set_index('Source Name')['Characteristics [inferred cell type]']
    zp = zipfile.ZipFile(os.path.join(ROOT, 'E-MTAB-5061', 'E-MTAB-5061.processed.1.zip')); name = zp.infolist()[0].filename; assert name == 'pancreas_refseq_rpkms_counts_3514sc.txt'
    with zp.open(name) as f:
        t = io.TextIOWrapper(f, encoding='utf8', errors='replace', newline=''); ids = t.readline().rstrip('\r\n').split('\t')[1:]; assert len(ids) == NS, len(ids)
        pos = {k: i for i, k in enumerate(ids)}; cols = [2 + NS + pos[k] for k in F.cell_key]; genes, rows, ins, gcg = [], [], None, None
        for i, line in enumerate(t):
            p = line.rstrip('\r\n').split('\t'); assert len(p) == 2 + 2 * NS, 'row %d has %d fields' % (i, len(p))
            if i == NG: assert p[0].startswith('ERCC'), 'row %d is not the first ERCC row: %s' % (i, p[0]); break
            genes.append(p[0]); rows.append(np.array([p[c] for c in cols], dtype=np.float64))
            if p[0] == 'INS': ins = np.array(p[2 + NS:], dtype=float)
            if p[0] == 'GCG': gcg = np.array(p[2 + NS:], dtype=float)
    assert ins is not None and gcg is not None and len(genes) == NG
    # identity check (not a selection step)
    lab = np.array([S.get(k, 'unlabelled') for k in ids]); chk = {}
    for gene, vec, own in (('INS', ins, 'beta cell'), ('GCG', gcg, 'alpha cell')):
        # the comparison types exclude the authors' categories that carry several hormones by definition ('co-expression cell', 'unclassified endocrine cell') and cells without a type ('not applicable');
        # the first version of this check (all types) stopped on them and is logged as a deviation in the protocol, section 8 (6 Nov 2026)
        own_mean = float(vec[lab == own].mean()); others = {l: float(vec[lab == l].mean()) for l in sorted(set(lab)) if l not in (own, 'unlabelled', 'nan', 'co-expression cell', 'unclassified endocrine cell', 'not applicable') and (lab == l).sum() >= 10}
        chk[gene] = dict(cell_type=own, mean_count_own=own_mean, max_mean_count_other=max(others.values()), other_means=others, passed=bool(own_mean > 2 * max(others.values())))
    json.dump(chk, open(os.path.join(HERE, 'logs', 'identity_check_E-MTAB-5061.json'), 'w'), indent=1); print('identity check', {g: (round(v['mean_count_own'], 1), round(v['max_mean_count_other'], 1), v['passed']) for g, v in chk.items()}, flush=True)
    assert all(v['passed'] for v in chk.values()), 'identity check failed: the column alignment is wrong; dataset reported as failed, not repaired'
    X = pd.DataFrame(np.vstack(rows), index=genes); ndup = int(X.index.duplicated().sum()); X = X.groupby(level=0, sort=False).sum()
    save('emtab5061', 'E-MTAB-5061', 'count', 'symbol', list(X.index), X.to_numpy(), F, 'read counts of the 26179 RefSeq rows (ERCC and eGFP excluded); %d rows with a repeated gene symbol were summed; library size = sum over the %d remaining rows' % (ndup, len(X)))


if __name__ == '__main__':
    {'GSE86469': parse_gse86469, 'GSE81547': parse_gse81547, 'E-MTAB-5061': parse_emtab5061}[sys.argv[1]]()
