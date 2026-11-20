"""post_hoc_validation_191126.py - main session, 19 Nov 2026. The two post hoc checks of post_hoc_checks_191126.py on the 59 validation sets (40 held-out arrays, 19 RNA-seq series).
The stored results hold scalars only, so the 12 x 17 matrices are rebuilt from the public GEO files with the benchmark's own loaders (wp2_benchmark2/bench_data_io.py and run_rnaseq.py, read only), on the series that
were 'analysed' in heldout_arrays_results.jsonl and heldout_rnaseq_results.jsonl. Gates (before any summary is written): the rebuilt gene sets equal the stored ones, and the observed share and the original T_r2 P
equal the stored share_between and T_r2_p of every matrix.
NO FILE IS DOWNLOADED unless --allow-download is given, and that flag is given only after the user has permitted the download. Downloaded files are untrusted data: read as text or numbers only,
kept in their own folders under DATA_DIR (v2_data/arrays_heldout/matrices, bench_data/platforms, v2_data/rnaseq/counts), never executed. NCBI politeness as in ncbi_http.py (<= 3 requests per second, no e-mail address).
Usage (python -I):  python -I post_hoc_validation_191126.py DATA_DIR --list            the files that would be needed (no network)
                    python -I post_hoc_validation_191126.py DATA_DIR --allow-download   rebuild, run the two checks, write results_validation/"""
import os, sys, csv, json, argparse
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import post_hoc_checks_191126 as P
W2 = P.W2
import bench_data_io as D, select_rnaseq as S, ncbi_http as H, run_rnaseq as RR


def load_analysed():
    arr = [json.loads(l) for l in open(os.path.join(W2, 'heldout_arrays_results.jsonl'), encoding='utf8')]; rna = [json.loads(l) for l in open(os.path.join(W2, 'heldout_rnaseq_results.jsonl'), encoding='utf8')]
    return [r for r in arr if r.get('status') == 'analysed'], [r for r in rna if r.get('status') == 'analysed']


def needed_files(root, arr, rna):
    mat = os.path.join(root, 'v2_data', 'arrays_heldout', 'matrices'); plat = os.path.join(root, 'bench_data', 'platforms'); out = []
    for r in arr:
        n = int(r['gse'][3:]); folder = 'GSE%dnnn' % (n // 1000) if n >= 1000 else 'GSEnnn'
        out.append((os.path.join(mat, '%s_%s.txt.gz' % (r['gse'], r['platform'])), '%sseries/%s/%s/matrix/%s_series_matrix.txt.gz (the file name carries the platform for multi-platform series)' % (D.FTP, folder, r['gse'], r['gse'])))
    for gpl in sorted({r['platform'] for r in arr}):
        n = int(gpl[3:]); folder = 'GPL%dnnn' % (n // 1000) if n >= 1000 else 'GPLnnn'; out.append((os.path.join(plat, gpl + '.annot.gz'), '%splatforms/%s/%s/annot/%s.annot.gz' % (D.FTP, folder, gpl, gpl)))
    for r in rna:
        org = r.get('organism', 'human'); out.append((RR.counts_path(root, r['gse'], S.BUILD[org]), S.COUNT_URL.format(acc=r['gse'], build=S.BUILD[org])))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('data_dir'); ap.add_argument('--list', action='store_true'); ap.add_argument('--allow-download', action='store_true'); a = ap.parse_args()
    root = a.data_dir; arr, rna = load_analysed(); need = needed_files(root, arr, rna)
    if a.list:
        for p, u in need: print('present' if os.path.exists(p) and os.path.getsize(p) > 0 else 'missing', '|', p, '|', u)
        print(len(need), 'files;', sum(1 for p, _ in need if os.path.exists(p) and os.path.getsize(p) > 0), 'present'); return
    D.configure(os.path.join(root, 'v2_data', 'arrays_heldout', 'matrices'), os.path.join(root, 'bench_data', 'platforms')); orig_fetch = D.fetch

    def guarded(url, path):
        if os.path.exists(path) and os.path.getsize(path) > 0: return path
        if not a.allow_download: raise SystemExit('download not allowed: %s is missing (%s); rerun with --allow-download only after the user has permitted it' % (path, url))
        return orig_fetch(url, path)
    D.fetch = guarded
    frozen_a = {r['gse']: r for r in csv.DictReader(open(os.path.join(W2, 'frozen_heldout_arrays.csv'), encoding='utf8'))}
    frozen_r = {r['gse']: r for r in csv.DictReader(open(os.path.join(W2, 'frozen_heldout_rnaseq.csv'), encoding='utf8'))}
    mats, stored, arm_of, gene_mismatch = [], {}, {}, []
    for rec in arr:
        gse = rec['gse']; rec2, sets = D.prepare_array(frozen_a[gse])
        if sets is None: raise SystemExit('%s: not analysable when rebuilt: %s' % (gse, rec2.get('status')))
        for tag in ('A', 'B'):
            if sets[tag]['genes'] != rec[tag]['genes']: gene_mismatch.append((gse, tag))
            mats.append((gse, rec['platform'], tag, sets[tag]['X'])); stored[(gse, tag)] = rec[tag]['flat']; arm_of[gse] = 'arrays ' + rec['platform']
    for rec in rna:
        gse = rec['gse']; row = frozen_r[gse]; org = row['organism']; path = RR.counts_path(root, gse, S.BUILD[org])
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            if not a.allow_download: raise SystemExit('download not allowed: %s is missing; rerun with --allow-download only after the user has permitted it' % path)
            path, why = RR.fetch_counts(S, H, root, gse, org)
            if path is None: raise SystemExit('%s: %s' % (gse, why))
        C, ids = RR.read_counts(path, row['gsms'].split(';')); r0 = dict(gse=gse, platform=org, organism=org, order=int(row['order']), how=row['how'], groups=row['groups'])
        rec2, sets = D.prepare_rnaseq(r0, C, ids, [None] * len(ids), int(gse[3:]))
        if sets is None: raise SystemExit('%s: not analysable when rebuilt: %s' % (gse, rec2.get('status')))
        for tag in ('A', 'B'):
            if sets[tag]['genes'] != rec[tag]['genes']: gene_mismatch.append((gse, tag))
            mats.append((gse, 'rnaseq ' + org, tag, sets[tag]['X'])); stored[(gse, tag)] = rec[tag]['flat']; arm_of[gse] = 'rnaseq ' + org
    print('rebuilt %d matrices; gene sets different from the stored ones: %d' % (len(mats), len(gene_mismatch)), flush=True)
    out = os.path.join(HERE, 'results_validation'); os.makedirs(out, exist_ok=True)
    if gene_mismatch: json.dump(gene_mismatch, open(os.path.join(out, 'GATE_FAILED_genes.json'), 'w')); raise SystemExit('GATE FAILED: gene sets differ from the stored ones')
    rows, gate = P.analyse(mats, 'validation sets (59), rebuilt matrices', stored)
    if gate['share_max_abs_diff'] > 1e-9 or gate['tr2_p_mismatch'] != 0 or gate['n_checked'] != len(mats):
        json.dump(dict(gate=gate, n=len(mats)), open(os.path.join(out, 'GATE_FAILED.json'), 'w')); print('GATE FAILED', gate); sys.exit(2)
    for r in rows: r['arm'] = arm_of[r['gse']]
    with open(os.path.join(out, 'per_matrix.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    S_all = P.summarize(rows, 'validation sets (59), rebuilt matrices', gate); json.dump(S_all, open(os.path.join(out, 'summary.json'), 'w'), indent=1, default=float)
    txt = [P.write_md(S_all, os.path.join(out, 'summary.md'))]
    for arm in sorted(set(arm_of.values())):
        sub = [r for r in rows if r['arm'] == arm]; Sa = P.summarize(sub, 'validation sets, ' + arm, gate); txt.append(P.write_md(Sa, os.path.join(out, 'summary_%s.md' % arm.replace(' ', '_')))); json.dump(Sa, open(os.path.join(out, 'summary_%s.json' % arm.replace(' ', '_')), 'w'), indent=1, default=float)
    print('\n\n'.join(txt))


if __name__ == '__main__': main()
