"""run_rnaseq.py - protocol v2 section 3.5, held-out RNA-seq arm (WP2 deliverable E): counts, dry run and full run on the frozen list (frozen_heldout_rnaseq.csv, fingerprinted before any count was read).
Written by the planner in the main session after the WP2 agent stopped; same structure and code path as run_heldout.py (bench_data_io.prepare_rnaseq -> prepare_from_log -> bench2_core.analyse_matrix).
Usage: python -I run_rnaseq.py dry | full
  For every frozen series in order (per organism): download the NCBI-generated raw count matrix <GSE>_raw_counts_<build>_NCBI.tsv.gz (human GRCh38.p13, mouse GRCm39; cached under
  <data root>/v2_data/rnaseq/counts, untrusted data read as numbers only). A 404 means the series has no NCBI-generated counts and is recorded as such. A reply that is a web page instead of a count
  file is recorded as blocked and is NOT worked around: after 5 such replies in a row the organism is stopped ('not done'). The 12 selected samples are taken by GSM accession; CPM on all genes of
  the file, genes need count >= 1 in all 12 samples, expressed and DE as in part 1 on log2(CPM + 1); a series needs >= 34 DE genes. Gene ids (Entrez GeneID) are the gene names: the draws of the
  17 genes are random, so no annotation is needed. Stop at 20 analysable series per organism or when the list is exhausted. Duplicate rule as in part 1: a series sharing a selected sample with an
  earlier analysable series is excluded.
  dry   decides analysability only (counts available, samples present, DE genes >= 34); writes heldout_rnaseq_dry.jsonl; no correlation is computed
  full  same loop with the analysis of gene sets A and B; one JSON line per series in heldout_rnaseq_results.jsonl; every decision is checked against the dry run (a mismatch is logged)"""
import sys, os, csv, json, time, gzip, urllib.error
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
MAXP = 20; MAX_BLOCKED_IN_A_ROW = 5


def counts_path(root, gse, build): return os.path.join(root, 'v2_data', 'rnaseq', 'counts', '%s_%s.tsv.gz' % (gse, build))


def fetch_counts(S, H, root, gse, org):
    """(path, None) or (None, reason); the count file is cached"""
    path = counts_path(root, gse, S.BUILD[org]); os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and os.path.getsize(path) > 0: return path, None
    try: data = H.http_get(S.COUNT_URL.format(acc=gse, build=S.BUILD[org]), binary=True, timeout=180)
    except urllib.error.HTTPError as e: return None, 'no NCBI-generated counts (HTTP %d)' % e.code
    if data[:2] != b'\x1f\x8b': return None, 'blocked: the reply was a web page, not a count file'
    open(path + '.part', 'wb').write(data); os.replace(path + '.part', path); return path, None


def read_counts(path, gsms):
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f: T = pd.read_csv(f, sep='\t', index_col=0, low_memory=False)
    miss = [g for g in gsms if g not in T.columns]
    if miss: raise ValueError('%d of 12 selected samples are missing in the count file' % len(miss))
    C = T[gsms].apply(pd.to_numeric, errors='coerce')
    if C.isna().any().any(): raise ValueError('non-numeric entries in the count file')
    return C.values.astype(float), [str(i) for i in C.index]


def main(dry):
    import bench_data_io as D, bench2_core as B, select_rnaseq as S, ncbi_http as H
    root = D.data_root(); rows = list(csv.DictReader(open(os.path.join(HERE, 'frozen_heldout_rnaseq.csv'), encoding='utf8'))); orgs = sorted({r['organism'] for r in rows})
    dry_path = os.path.join(HERE, 'heldout_rnaseq_dry.jsonl'); ref = {}
    if not dry and os.path.exists(dry_path): ref = {j['gse']: j for j in (json.loads(l) for l in open(dry_path, encoding='utf8'))}
    out = open(dry_path if dry else os.path.join(HERE, 'heldout_rnaseq_results.jsonl'), 'w', encoding='utf8')
    done = {o: 0 for o in orgs}; used = {}; blocked = {o: 0 for o in orgs}; stopped = {o: False for o in orgs}; mism = 0; tally = {o: {} for o in orgs}
    for row in rows:
        org = row['organism']
        if done[org] >= MAXP or stopped[org]: continue
        t0 = time.time(); gs = row['gsms'].split(';'); clash = sorted({used[g] for g in gs if g in used}); sets = None
        rec = dict(gse=row['gse'], platform=org, organism=org, order=int(row['order']), how=row['how'], groups=row['groups'])
        if clash: rec['status'] = 'not analysable: duplicate of ' + clash[0]
        else:
            path, why = fetch_counts(S, H, root, row['gse'], org)
            if path is None:
                rec['status'] = 'not analysable: ' + why; blocked[org] = blocked[org] + 1 if why.startswith('blocked') else 0
                if blocked[org] >= MAX_BLOCKED_IN_A_ROW: stopped[org] = True; rec['status'] += ' (organism stopped after %d in a row: not done, no work-around)' % MAX_BLOCKED_IN_A_ROW
            else:
                blocked[org] = 0
                try: C, ids = read_counts(path, gs); rec, sets = D.prepare_rnaseq(rec, C, ids, [None] * len(ids), int(row['gse'][3:]), dry)
                except Exception as e: rec['status'] = 'not analysable: %s %s' % (type(e).__name__, str(e)[:120])
        if sets is not None and not dry:
            try:
                for tag in ('A', 'B'):
                    flat, detail = B.analyse_matrix(sets[tag]['X'], sets[tag]['genes'], int(row['gse'][3:]))
                    rec[tag] = dict(flat=flat, detail=detail, genes=sets[tag]['genes'], n_de_in_set=sets[tag]['n_de_in_set'])
                rec['status'] = 'analysed'
            except Exception as e: rec['status'] = 'not analysable: analysis error %s %s' % (type(e).__name__, str(e)[:120]); rec.pop('A', None); rec.pop('B', None)
        if rec['status'] in ('analysed', 'analysable (dry run)'): done[org] += 1; used.update({g: row['gse'] for g in gs})
        tally[org][rec['status'][:60]] = tally[org].get(rec['status'][:60], 0) + 1
        if ref:
            r0 = ref.get(row['gse']); a0 = None if r0 is None else (r0['status'] == 'analysable (dry run)'); a1 = rec['status'] == 'analysed'
            if r0 is None or a0 != a1 or r0.get('n_de') != rec.get('n_de'): mism += 1; rec['dry_run_mismatch'] = True
        out.write(json.dumps(rec, default=lambda x: x.item() if hasattr(x, 'item') else (x.tolist() if hasattr(x, 'tolist') else str(x))) + '\n'); out.flush()
        print(f"{org} {row['gse']} (#{row['order']}): {rec['status'][:110]} | DE {rec.get('n_de', '-')} | done {done} | {time.time() - t0:.0f}s", flush=True)
    print('finished', done, 'dry-run mismatches:', mism, flush=True)
    for o in orgs: print('status counts', o, json.dumps(tally[o]), flush=True)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    if mode in ('dry', 'full'): main(mode == 'dry')
    else: raise SystemExit('usage: run_rnaseq.py dry|full')
