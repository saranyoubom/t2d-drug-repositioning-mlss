"""gate_part1.py - WP2 deliverable A gate: reload the development series from the cached part 1 data and reproduce the part 1 per-dataset numbers before anything new is done.
Steps (python -I gate_part1.py <step>):
  extract      re-run the part 1 loop (frozen order, duplicate rule, 20 analysable per platform) with bench_data_io on bench_data/{matrices,platforms}; compare status, DE counts, scale and the two
               17-gene draws with results_per_dataset.jsonl; save the 12 x 17 matrices (compact/dev_sets.npz, compact/dev_prepared.json)
  compare      bench_pipeline.analyse (the part 1 function, 2000 permutations, seed 2026 + GSE number: the part 1 call) on gene sets A and B; compare observed consensus count, null mean, null SD, P
               with per_dataset_results.csv (exact match required); workers = argv[2] (default 2); resumes from compact/gate_bench_pipeline.jsonl
  consistency  coexp_checks (bench2_core.analyse_matrix) on the same 74 matrices: observed values must be identical to part 1, null means and P agree within Monte Carlo error (not a gate);
               the results are cached in compact/dev_core_results.jsonl and reused by apply_dev.py
  report       writes gate_part1_reproduction.txt from the three result files
Downloaded data are untrusted: read as numbers only."""
import sys, os, csv, json, time
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); AO = os.path.dirname(os.path.dirname(HERE)); P1 = os.path.join(AO, 'benchmark_route3')
for p in (HERE, P1): sys.path.insert(0, p)
import numpy as np
COMPACT = os.path.join(HERE, 'compact'); os.makedirs(COMPACT, exist_ok=True)
F_EXTRACT = os.path.join(COMPACT, 'extract_check.json'); F_NPZ = os.path.join(COMPACT, 'dev_sets.npz'); F_PREP = os.path.join(COMPACT, 'dev_prepared.json')
F_BENCH = os.path.join(COMPACT, 'gate_bench_pipeline.jsonl'); F_CORE = os.path.join(COMPACT, 'dev_core_results.jsonl'); F_REPORT = os.path.join(HERE, 'gate_part1_reproduction.txt')


def part1_records():
    return [json.loads(l) for l in open(os.path.join(P1, 'results_per_dataset.jsonl'), encoding='utf8')]


def extract():
    import bench_data_io as D
    root = D.data_root(); D.configure(os.path.join(root, 'bench_data', 'matrices'), os.path.join(root, 'bench_data', 'platforms'))
    rows = list(csv.DictReader(open(os.path.join(P1, 'frozen_datasets.csv'), encoding='utf8'))); ref = {r['gse']: r for r in part1_records()}
    done = {'GPL570': 0, 'GPL1261': 0}; used = {}; checks = []; arrays = {}; prep = {}
    for row in rows:
        gpl = row['platform']
        if done[gpl] >= 20: continue
        t0 = time.time(); gs = set(row['gsms'].split(';')); clash = sorted({used[g] for g in gs if g in used}); sets = None
        if clash: rec = dict(gse=row['gse'], platform=gpl, order=int(row['order']), status='not analysable: duplicate of ' + clash[0])
        else:
            try: rec, sets = D.prepare_array(row)
            except Exception as e: rec = dict(gse=row['gse'], platform=gpl, order=int(row['order']), status='not analysable: %s %s' % (type(e).__name__, str(e)[:120]))
        if rec['status'] == 'prepared': rec['status'] = 'analysed'
        if rec['status'] in ('analysed', 'analysable (dry run)'): done[gpl] += 1; used.update({g: row['gse'] for g in gs})
        r0 = ref.get(row['gse']); ok = {}
        if r0 is None: ok['in_part1_log'] = False
        else:
            ok['status'] = rec['status'] == r0['status']; ok['n_de'] = rec.get('n_de') == r0.get('n_de'); ok['n_expressed'] = rec.get('n_expressed') == r0.get('n_expressed')
            ok['n_genes_mapped'] = rec.get('n_genes_mapped') == r0.get('n_genes_mapped'); ok['scale'] = rec.get('scale') == r0.get('scale')
            ok['q99'] = abs(rec.get('q99', 0) - r0.get('q99', 0)) < 1e-9 if 'q99' in r0 else ('q99' not in rec)
            if sets is not None:
                for tag in ('A', 'B'): ok['genes_' + tag] = sets[tag]['genes'] == r0[tag]['genes']; ok['n_de_in_set_' + tag] = sets[tag]['n_de_in_set'] == r0[tag]['n_de_in_set']
        checks.append(dict(gse=row['gse'], platform=gpl, order=int(row['order']), status=rec['status'], n_de=rec.get('n_de'), ok=ok, all_ok=all(ok.values()), seconds=round(time.time() - t0, 1)))
        if sets is not None:
            for tag in ('A', 'B'): arrays['%s__%s' % (row['gse'], tag)] = sets[tag]['X']
            prep[row['gse']] = dict(platform=gpl, order=int(row['order']), n_de=rec['n_de'], genes_A=sets['A']['genes'], genes_B=sets['B']['genes'], n_de_in_set_A=sets['A']['n_de_in_set'], n_de_in_set_B=sets['B']['n_de_in_set'],
                                    scale=rec['scale'], groups=rec['groups'])
        print('%s %s (#%s): %s | DE %s | all matches %s | %ss' % (gpl, row['gse'], row['order'], rec['status'], rec.get('n_de', '-'), checks[-1]['all_ok'], checks[-1]['seconds']), flush=True)
    n_ref = len(ref); n_hit = sum(1 for c in checks if c['ok'].get('in_part1_log', True))
    summary = dict(n_rows_processed=len(checks), n_in_part1_log=n_hit, n_part1_records=n_ref, n_all_ok=sum(c['all_ok'] for c in checks), n_analysed=sum(c['status'] == 'analysed' for c in checks),
                   analysed_per_platform=done, checks=checks)
    json.dump(summary, open(F_EXTRACT, 'w'), indent=1); np.savez_compressed(F_NPZ, **arrays); json.dump(prep, open(F_PREP, 'w'), indent=1)
    print('extract finished:', {k: v for k, v in summary.items() if k != 'checks'}, flush=True)


def _work(task):
    import warnings; warnings.filterwarnings('ignore')
    import bench_pipeline as B
    gse, tag, X, n = task; t0 = time.time(); r = B.analyse(np.array(X), [0, 1, 2], nperm=2000, seed=2026 + n)
    return dict(gse=gse, tag=tag, res=r, sec=round(time.time() - t0, 1))


def compare(nw=2):
    import multiprocessing as mp
    Z = np.load(F_NPZ); prep = json.load(open(F_PREP)); tasks = []
    have = set()
    if os.path.exists(F_BENCH): have = {(j['gse'], j['tag']) for j in (json.loads(l) for l in open(F_BENCH, encoding='utf8'))}
    for gse in prep:
        for tag in ('A', 'B'):
            if (gse, tag) not in have: tasks.append((gse, tag, Z['%s__%s' % (gse, tag)].tolist(), int(gse[3:])))
    print('tasks to run:', len(tasks), 'already done:', len(have), 'workers:', nw, flush=True)
    with open(F_BENCH, 'a', encoding='utf8') as out, mp.get_context('spawn').Pool(nw) as pool:
        for j in pool.imap_unordered(_work, tasks):
            out.write(json.dumps(j) + '\n'); out.flush(); print('done', j['gse'], j['tag'], 'obs', j['res']['observed'], 'P', j['res']['perm_p_plain'], '%ss' % j['sec'], flush=True)
    print('compare finished', flush=True)


def consistency():
    import bench2_core as C
    Z = np.load(F_NPZ); prep = json.load(open(F_PREP)); have = set()
    if os.path.exists(F_CORE): have = {(j['gse'], j['tag']) for j in (json.loads(l) for l in open(F_CORE, encoding='utf8'))}
    with open(F_CORE, 'a', encoding='utf8') as out:
        for gse, p in prep.items():
            for tag in ('A', 'B'):
                if (gse, tag) in have: continue
                t0 = time.time(); X = Z['%s__%s' % (gse, tag)]; flat, detail = C.analyse_matrix(X, p['genes_' + tag], int(gse[3:]))
                out.write(json.dumps(dict(gse=gse, platform=p['platform'], tag=tag, flat=flat, detail=detail), default=float) + '\n'); out.flush()
                print('core', gse, tag, 'obs', flush=True, end=' '); print(flat['observed'], 'T_r2 P %.3f' % flat['T_r2_p'], '%.1fs' % (time.time() - t0), flush=True)
    print('consistency finished', flush=True)


def report():
    import pandas as pd
    ex = json.load(open(F_EXTRACT)); T = pd.read_csv(os.path.join(P1, 'per_dataset_results.csv')); ref = {(r['gse'], t): r[t] for r in part1_records() if r['status'] == 'analysed' for t in ('A', 'B')}
    new = {(j['gse'], j['tag']): j['res'] for j in (json.loads(l) for l in open(F_BENCH, encoding='utf8'))}
    L = ['GATE A: reproduction of the part 1 per-dataset numbers on the development series', '', 'Run: %s; numpy %s' % (time.strftime('%Y-%m-%d %H:%M'), np.__version__), '',
         '1. Reload of the cached series (bench_data), same loop as run_benchmark.py (frozen order, duplicate rule, at most 20 analysable per platform)',
         '   rows processed %d (part 1 log: %d records); status, DE counts, expressed counts, scale, 99th percentile, gene lists of sets A and B and DE-in-set counts identical for %d of %d rows; analysed per platform %s'
         % (ex['n_rows_processed'], ex['n_part1_records'], ex['n_all_ok'], ex['n_rows_processed'], ex['analysed_per_platform'])]
    bad = [c for c in ex['checks'] if not c['all_ok']]
    L.append('   mismatching rows: ' + ('none' if not bad else '; '.join('%s %s' % (c['gse'], [k for k, v in c['ok'].items() if not v]) for c in bad)))
    L += ['', '2. bench_pipeline.analyse (part 1 function; 2000 permutations; seed 2026 + GSE number) against per_dataset_results.csv, gene sets A and B',
          '   criterion: observed count identical; null mean within 1e-9; permutation P (plain fraction) identical', '',
          'series        set  observed(p1/now)  null_mean(p1/now)         P(p1/now)       null_sd(p1/now)       match']
    n_ok = 0; n_tot = 0; maxd = 0.0; rows = []
    for (gse, tag) in sorted(ref, key=lambda k: (k[0], k[1])):
        n_tot += 1; t = T[(T.gse == gse) & (T.set == tag)]
        if (gse, tag) not in new or len(t) != 1: L.append('%-12s %s  MISSING (new %s, table rows %d)' % (gse, tag, (gse, tag) in new, len(t))); continue
        r = new[(gse, tag)]; t = t.iloc[0]; ok = (int(t.observed) == r['observed']) and abs(t.null_mean - r['null_mean']) <= 1e-9 and abs(t.perm_p - r['perm_p_plain']) <= 1e-12 and abs(t.null_sd - r['null_sd']) <= 1e-9
        maxd = max(maxd, abs(t.null_mean - r['null_mean'])); n_ok += ok
        L.append('%-12s %s  %3d / %3d        %9.4f / %9.4f   %6.3f / %6.3f    %8.4f / %8.4f   %s' % (gse, tag, int(t.observed), r['observed'], t.null_mean, r['null_mean'], t.perm_p, r['perm_p_plain'], t.null_sd, r['null_sd'], 'OK' if ok else 'MISMATCH'))
    L += ['', '   matches: %d of %d (series x set); largest absolute difference in the null mean %.2e' % (n_ok, n_tot, maxd)]
    gate = (n_ok == n_tot and n_tot == 2 * ex['n_analysed'] and ex['n_all_ok'] == ex['n_rows_processed'] and ex['n_analysed'] == 37)
    L += ['', 'GATE RESULT: %s' % ('PASS (all %d analysed series reload, the draws are identical and the part 1 numbers are reproduced exactly)' % ex['n_analysed'] if gate else 'FAIL'), '']
    if os.path.exists(F_CORE):
        core = {(j['gse'], j['tag']): j for j in (json.loads(l) for l in open(F_CORE, encoding='utf8'))}; rows = []
        for (gse, tag), r0 in ref.items():
            if (gse, tag) not in core: continue
            f = core[(gse, tag)]['flat']; rows.append(dict(gse=gse, tag=tag, obs_eq=f['observed'] == r0['observed'], same_eq=f['observed_same_sign'] == r0['observed_same_sign'],
                                                      pooled_eq=[f['pooled_raw'], f['pooled_den1'], f['pooled_den2'], f['pooled_den3']] == r0['pooled_sig_per_scheme'],
                                                      adj_eq=[f['adjq_raw_obs'], f['adjq_den1_obs'], f['adjq_den2_obs'], f['adjq_den3_obs']] == r0['adj_q05_per_scheme'], any_eq=f['adj_any'] == r0['cons_adj_any'],
                                                      all4_eq=f['adj_all4'] == r0['cons_adj_all4'], pairs4_eq=f['adj_pairs_all4'] == r0['adj_pairs_all4'], share_eq=abs(f['share_between'] - r0['share_between_raw']) < 1e-9 if r0['share_between_raw'] == r0['share_between_raw'] else True,
                                                      dmean=f['null_mean'] - r0['null_mean'], dsd=f['null_sd'] - r0['null_sd'], dp=f['perm_p'] - r0['perm_p_plain'], seflag=abs(f['null_mean'] - r0['null_mean']) / max(r0['null_sd'], 1e-9) * np.sqrt(2000)))
        R = pd.DataFrame(rows)
        L += ['3. Consistency of the new module (coexp_checks via bench2_core, 2000 permutations) with part 1 on the same %d matrices (not a gate: random streams differ by design)' % len(R),
              '   identical observed values: consensus count %d/%d, same-sign count %d/%d, pooled counts per scheme %d/%d, adjusted q<0.05 counts per scheme %d/%d, consensus pairs retained (any scheme) %d/%d, (all schemes) %d/%d,' %
              (R.obs_eq.sum(), len(R), R.same_eq.sum(), len(R), R.pooled_eq.sum(), len(R), R.adj_eq.sum(), len(R), R.any_eq.sum(), len(R), R.all4_eq.sum(), len(R)),
              '   pairs significant in all schemes after adjustment %d/%d, covariance share %d/%d' % (R.pairs4_eq.sum(), len(R), R.share_eq.sum(), len(R)),
              '   null mean difference (new - part 1): mean %.3f, SD %.3f, largest absolute %.3f; null SD difference: mean %.3f; P difference: mean %.4f, largest absolute %.3f;' % (R.dmean.mean(), R.dmean.std(), R.dmean.abs().max(), R.dsd.mean(), R.dp.mean(), R.dp.abs().max()),
              '   |difference of null means| in Monte Carlo SEs of the difference (sqrt(2/2000) x SD): median %.2f, 95th percentile %.2f, maximum %.2f' % (np.nanmedian(R.seflag / np.sqrt(2)), np.nanpercentile(R.seflag / np.sqrt(2), 95), np.nanmax(R.seflag / np.sqrt(2))), '']
    open(F_REPORT, 'w', encoding='utf8').write('\n'.join(L) + '\n'); print('\n'.join(L[:12])); print('...'); print('\n'.join(L[-14:]))


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'report'
    if step == 'extract': extract()
    elif step == 'compare': compare(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    elif step == 'consistency': consistency()
    elif step == 'report': report()
    else: raise SystemExit('step: extract | compare [workers] | consistency | report')
