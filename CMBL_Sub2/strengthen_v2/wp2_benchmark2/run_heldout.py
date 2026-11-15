"""run_heldout.py - dry run and full run on the frozen held-out list (protocol v2 section 3.5, WP2 deliverable D), same code path as the development series
(bench_data_io.prepare_array -> bench2_core.analyse_matrix; the loop, the duplicate rule and the stopping rule are those of run_benchmark.main).
Usage: python -I run_heldout.py arrays dry | arrays full
  dry   decides analysability only (file integrity, number of DE genes >= 34), writes heldout_arrays_dry.jsonl; no correlation is computed
  full  processes the frozen list in order until 20 series per platform are analysable or the list is exhausted; one JSON line per series (status or reason) in heldout_arrays_results.jsonl with the
        flat endpoints and details of gene sets A and B; every decision is checked against the dry run (a mismatch is logged, never silently accepted)
Duplicate rule (as in part 1, deviation D1): a series sharing at least one of its 12 selected samples with an earlier analysable held-out series is excluded. Downloads are throttled (ncbi_http) and kept
under <data root>/v2_data/arrays_heldout/matrices (untrusted data, read as numbers only)."""
import sys, os, csv, json, time
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
MAXP = 20


def arrays(dry):
    import bench_data_io as D, bench2_core as C
    root = D.data_root(); D.configure(os.path.join(root, 'v2_data', 'arrays_heldout', 'matrices'), os.path.join(root, 'bench_data', 'platforms'))
    rows = list(csv.DictReader(open(os.path.join(HERE, 'frozen_heldout_arrays.csv'), encoding='utf8'))); plats = sorted({r['platform'] for r in rows})
    dry_path = os.path.join(HERE, 'heldout_arrays_dry.jsonl'); ref = {}
    if not dry and os.path.exists(dry_path): ref = {j['gse']: j for j in (json.loads(l) for l in open(dry_path, encoding='utf8'))}
    out = open(dry_path if dry else os.path.join(HERE, 'heldout_arrays_results.jsonl'), 'w', encoding='utf8'); done = {p: 0 for p in plats}; used = {}; mism = 0
    for row in rows:
        gpl = row['platform']
        if done[gpl] >= MAXP: continue
        t0 = time.time(); gs = set(row['gsms'].split(';')); clash = sorted({used[g] for g in gs if g in used}); sets = None
        if clash: rec = dict(gse=row['gse'], platform=gpl, order=int(row['order']), status='not analysable: duplicate of ' + clash[0])
        else:
            try: rec, sets = D.prepare_array(row, dry)
            except Exception as e: rec = dict(gse=row['gse'], platform=gpl, order=int(row['order']), status='not analysable: %s %s' % (type(e).__name__, str(e)[:120]))
        if sets is not None and not dry:
            try:
                for tag in ('A', 'B'):
                    flat, detail = C.analyse_matrix(sets[tag]['X'], sets[tag]['genes'], int(row['gse'][3:]))
                    rec[tag] = dict(flat=flat, detail=detail, genes=sets[tag]['genes'], n_de_in_set=sets[tag]['n_de_in_set'])
                rec['status'] = 'analysed'
            except Exception as e: rec['status'] = 'not analysable: analysis error %s %s' % (type(e).__name__, str(e)[:120]); rec.pop('A', None); rec.pop('B', None)
        if rec['status'] in ('analysed', 'analysable (dry run)'): done[gpl] += 1; used.update({g: row['gse'] for g in gs})
        if ref:
            r0 = ref.get(row['gse']); a0 = None if r0 is None else (r0['status'] == 'analysable (dry run)'); a1 = rec['status'] == 'analysed'
            if r0 is None or a0 != a1 or r0.get('n_de') != rec.get('n_de'): mism += 1; rec['dry_run_mismatch'] = True
        out.write(json.dumps(rec, default=lambda x: x.item() if hasattr(x, 'item') else (x.tolist() if hasattr(x, 'tolist') else str(x))) + '\n'); out.flush()
        print(f"{gpl} {row['gse']} (#{row['order']}): {rec['status']} | DE {rec.get('n_de', '-')} | done {done} | {time.time() - t0:.0f}s", flush=True)
    print('finished', done, 'dry-run mismatches:', mism, flush=True)


if __name__ == '__main__':
    arm = sys.argv[1] if len(sys.argv) > 1 else ''; mode = sys.argv[2] if len(sys.argv) > 2 else ''
    if arm == 'arrays' and mode in ('dry', 'full'): arrays(mode == 'dry')
    else: raise SystemExit('usage: run_heldout.py arrays dry|full')
