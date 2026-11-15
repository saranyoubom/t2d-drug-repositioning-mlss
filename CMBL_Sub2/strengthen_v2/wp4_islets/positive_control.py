"""positive_control.py - WP4 extension A1-A3 of PRESPEC_WP4b_positive_control_and_wider_search.md (fixed before this script ran): what could the islet cross-check have detected?
A1  synthetic spike-in: one latent factor per replicate is added (donor-centred log scale, before the depth adjustment) to the distinct genes of the evaluable consensus pairs of a dataset with loading
    lambda_g = sqrt(rho / (1 - rho)) x sd_g, nominal rho in {0.05, 0.1, 0.2, 0.3, 0.5}; mode 'detected' (primary: only in cells where the gene is detected, zeros stay zeros) and 'all' (upper bound); the unchanged
    pipeline (depth-adjusted partial Spearman, excess E over the 2000 matched random pair sets of the primary run) is applied; 200 replicates; the five human datasets are pooled per replicate with the primary
    estimator (REML + Hartung-Knapp, margin +/-0.10). Gate: without injection the code must reproduce E and SE of results/<tag>__primary__summary.json (adjusted) to 1e-9.
A2  real-gene control: 100 random sets of 11 pairs of cytosolic ribosomal protein genes (rp_list.py); E_RP per set against matched random pairs (deciles of detection; pool without any ribosomal gene and without
    the 17 genes; 2000 draws per decile combination).
A3  margin sensitivity of the primary pooled interval (descriptive).
Output: positive_control/inj_<tag>.npz, rp_<tag>.json, positive_control_by_dataset.csv, positive_control_pooled.csv, rp_control.csv, margin_sensitivity.csv, positive_control_report.md.
Usage: python -I positive_control.py [tag ...]   run the datasets named (default all six) and then the report;   python -I positive_control.py report   only pool and write the report from the saved files."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
import wp4_core as C
import meta_analysis as MA
from metacell import fast_draws                    # same draws as wp4_core.draws_main (tested identical)
from run_datasets import load
from rp_list import RP, mouse_symbol

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'positive_control'); os.makedirs(OUT, exist_ok=True)
TAGS = ['gse81608', 'gse84133_human', 'gse86469', 'gse81547', 'emtab5061', 'gse84133_mouse']; HUMAN = TAGS[:5]
RHOS = [0.05, 0.1, 0.2, 0.3, 0.5]; MODES = ['detected', 'all']; NREP = 200; CLIP = 0.999999
NSETS, NPAIRS = 100, 11
if os.environ.get('PC_TEST'): NREP = 5; NSETS = 5; OUT = os.path.join(HERE, 'positive_control_test'); os.makedirs(OUT, exist_ok=True)      # code test only (5 replicates, separate folder); the real run never sets it


def prepare(cells):
    cons, consset, redge, cen = C.consensus_tables(); sym = cells.label
    det = np.asarray((cells.B > 0).mean(axis=0)).ravel(); elig = np.nonzero(det >= 0.10)[0]
    labels_e = [cells.genes[j] for j in elig]; idx = {g: r for r, g in enumerate(labels_e)}; det_e = det[elig]; own = {sym(g) for g in C.GENES}
    bins, pool = C.make_pool(det_e, labels_e, own)
    Xc = C.center_by_donor(C.lognorm(cells, elig), cells.donor); Z = C.depth_design(cells.depth, cells.donor); Rk = C.adjusted_ranks(Xc, Z)
    pairs = [(a, c) for i, a in enumerate(C.GENES) for c in C.GENES[i + 1:] if frozenset((a, c)) in consset and sym(a) in idx and sym(c) in idx]
    pidx = [(idx[sym(a)], idx[sym(c)]) for a, c in pairs]
    IA, IC = fast_draws(pool, bins, pidx, C.NN, C.SEED_NULL); V = C.null_values(Rk, IA, IC); zm = np.arctanh(np.clip(V, -CLIP, CLIP)).mean(axis=1)
    genes = [g for g in C.GENES if any(g in p for p in pairs)]; cols = [idx[sym(g)] for g in genes]; gi = {g: k for k, g in enumerate(genes)}
    G = np.asarray(Xc[:, cols], np.float64); M = np.asarray(cells.B[:, elig[cols]].todense() > 0)
    rho0 = np.array([float(Rk[a] @ Rk[c]) for a, c in pidx]); z_obs = float(np.arctanh(np.clip(rho0, -CLIP, CLIP)).mean())
    return dict(n=cells.n, donor=cells.donor, Z=Z, G=G, M=M, sd=G.std(axis=0), pair_k=[(gi[a], gi[c]) for a, c in pairs], zm=zm, zm_mean=float(zm.mean()), SE=float(zm.std(ddof=1)), E_primary=z_obs - float(zm.mean()),
                npairs=len(pairs), ngenes=len(genes), own=own, idx=idx, labels_e=labels_e, det_e=det_e, Rk=Rk)


def one_rep(P, rho, mode, rng):
    z = rng.standard_normal(P['n']); z = C.center_by_donor(z[:, None], P['donor'])[:, 0]; z = z / z.std()
    add = z[:, None] * (np.sqrt(rho / (1 - rho)) * P['sd'])[None, :]
    if mode == 'detected': add = add * P['M']
    R = C.adjusted_ranks(C.center_by_donor(P['G'] + add, P['donor']), P['Z'])
    zobs = float(np.arctanh(np.clip(np.array([R[a] @ R[c] for a, c in P['pair_k']]), -CLIP, CLIP)).mean())
    return zobs - P['zm_mean'], (1 + int((P['zm'] >= zobs).sum())) / (C.NN + 1)


def null_cell_mean(Rk, pool, da, dc, seed, n=C.NN):
    rng = np.random.default_rng(np.random.SeedSequence(seed)); pa, pc = pool[da], pool[dc]
    ia = pa[rng.integers(0, len(pa), n)]; ic = pc[rng.integers(0, len(pc), n)]; bad = ia == ic; k = 0
    while bad.any() and k < 1000: ic[bad] = pc[rng.integers(0, len(pc), int(bad.sum()))]; bad = ia == ic; k += 1
    return float(np.arctanh(np.clip(np.einsum('ij,ij->i', Rk[ia], Rk[ic]), -CLIP, CLIP)).mean())


def rp_control(cells, P, di):
    ids = pd.read_csv(os.path.join(HERE, 'rp_ids_human.csv'))
    lab = {s: int(e) for s, e in zip(ids.symbol, ids.entrez)} if cells.id_type == 'entrez' else {s: e for s, e in zip(ids.symbol, ids.ensembl)} if cells.id_type == 'ensembl' else \
          {s: mouse_symbol(s) for s in RP} if cells.species == 'mouse' else {s: s for s in RP}
    present = [s for s in RP if lab[s] in cells.gi]; eligible = [lab[s] for s in present if lab[s] in P['idx']]
    bins, pool = C.make_pool(P['det_e'], P['labels_e'], P['own'] | {lab[s] for s in RP})
    rng = np.random.default_rng(2030 + di); table = {}; Es = []; rhos = []
    for s in range(NSETS):
        zs, ns = [], []
        for _ in range(NPAIRS):
            a, c = rng.choice(len(eligible), 2, replace=False); ia, ic = P['idx'][eligible[a]], P['idx'][eligible[c]]; r = float(P['Rk'][ia] @ P['Rk'][ic]); zs.append(np.arctanh(np.clip(r, -CLIP, CLIP))); rhos.append(r)
            key = (int(bins[ia]), int(bins[ic]))
            if key not in table: table[key] = null_cell_mean(P['Rk'], pool, key[0], key[1], [2030, di, key[0], key[1]])
            ns.append(table[key])
        Es.append(float(np.mean(zs) - np.mean(ns)))
    Es = np.array(Es); return dict(rp_genes_in_data=len(present), rp_genes_eligible=len(eligible), E_RP_mean=float(Es.mean()), E_RP_sd=float(Es.std(ddof=1)), share_E_RP_ge_0_10=float((Es >= 0.10).mean()), mean_rho_RP_pairs=float(np.mean(rhos)))


def run_dataset(tag):
    di = TAGS.index(tag); t0 = time.time(); cells = load(tag); P = prepare(cells); s = json.load(open(os.path.join(HERE, 'results', f'{tag}__primary__summary.json')))['adjusted']
    assert abs(P['E_primary'] - s['E']) < 1e-9 and abs(P['SE'] - s['SE_E']) < 1e-9, 'GATE FAILED %s: E %.12f vs %.12f, SE %.12f vs %.12f' % (tag, P['E_primary'], s['E'], P['SE'], s['SE_E'])
    print('%-15s gate ok (E %+.4f, SE %.4f; %d pairs, %d genes) | prepared in %.0f s' % (tag, P['E_primary'], P['SE'], P['npairs'], P['ngenes'], time.time() - t0), flush=True)
    Erep = np.zeros((len(RHOS), len(MODES), NREP)); Pemp = np.zeros_like(Erep)
    for ri, rho in enumerate(RHOS):
        for mi, mode in enumerate(MODES):
            for rep in range(NREP): Erep[ri, mi, rep], Pemp[ri, mi, rep] = one_rep(P, rho, mode, np.random.default_rng(np.random.SeedSequence([2029, di, ri, mi, rep])))
    np.savez(os.path.join(OUT, f'inj_{tag}.npz'), Erep=Erep, Pemp=Pemp, E_primary=P['E_primary'], SE=P['SE'], rhos=np.array(RHOS), nrep=NREP)
    print('%-15s injection done (%.0f s): recovered E at nominal rho 0.3: detected %+.3f, all %+.3f' % (tag, time.time() - t0, Erep[RHOS.index(0.3), 0].mean(), Erep[RHOS.index(0.3), 1].mean()), flush=True)
    rp = rp_control(cells, P, di); json.dump(rp, open(os.path.join(OUT, f'rp_{tag}.json'), 'w'), indent=1)
    print('%-15s ribosomal control: %d genes eligible | E_RP mean %+.3f (SD %.3f), share >= 0.10: %.2f | %.0f s' % (tag, rp['rp_genes_eligible'], rp['E_RP_mean'], rp['E_RP_sd'], rp['share_E_RP_ge_0_10'], time.time() - t0), flush=True)


def report():
    D = {t: np.load(os.path.join(OUT, f'inj_{t}.npz')) for t in TAGS if os.path.exists(os.path.join(OUT, f'inj_{t}.npz'))}; rows = []
    for t, d in D.items():
        for mi, mode in enumerate(MODES):
            for ri, rho in enumerate(RHOS):
                E = d['Erep'][ri, mi]; zn = float(np.arctanh(rho))
                rows.append(dict(tag=t, mode=mode, rho=rho, nominal_z=zn, E_primary=float(d['E_primary']), SE_E=float(d['SE']), E_recovered_mean=float(E.mean()), E_recovered_sd=float(E.std(ddof=1)),
                                 attenuation=float((E.mean() - float(d['E_primary'])) / zn), share_empirical_P_lt_0_05=float((d['Pemp'][ri, mi] < 0.05).mean())))
    A = pd.DataFrame(rows); A.to_csv(os.path.join(OUT, 'positive_control_by_dataset.csv'), index=False)
    prow = []
    if all(t in D for t in HUMAN):
        se = np.array([float(D[t]['SE']) for t in HUMAN])
        for mi, mode in enumerate(MODES):
            for ri, rho in enumerate(RHOS):
                V, mu, lo, hi = [], [], [], []
                for rep in range(NREP):
                    r = MA.pool(np.array([D[t]['Erep'][ri, mi, rep] for t in HUMAN]), se); V.append(MA.verdict(r['lo'], r['hi'])); mu.append(r['mu']); lo.append(r['lo']); hi.append(r['hi'])
                V = np.array(V); prow.append(dict(mode=mode, rho=rho, replicates=NREP, share_no_meaningful_excess=float((V == 'no meaningful excess co-expression').mean()), share_inconclusive=float((V == 'inconclusive').mean()),
                                                  share_co_expressed_more=float((V == 'consensus pairs are co-expressed more than random pairs').mean()), pooled_E_mean=float(np.mean(mu)), pooled_E_sd=float(np.std(mu, ddof=1)),
                                                  pooled_lo_mean=float(np.mean(lo)), pooled_hi_mean=float(np.mean(hi))))
    Pp = pd.DataFrame(prow); Pp.to_csv(os.path.join(OUT, 'positive_control_pooled.csv'), index=False)
    rp = [dict(tag=t, **json.load(open(os.path.join(OUT, f'rp_{t}.json')))) for t in TAGS if os.path.exists(os.path.join(OUT, f'rp_{t}.json'))]; R = pd.DataFrame(rp)
    R['SE_E_consensus'] = [float(D[t]['SE']) if t in D else np.nan for t in R.tag]; R['E_RP_in_SE_units'] = R.E_RP_mean / R.SE_E_consensus; R.to_csv(os.path.join(OUT, 'rp_control.csv'), index=False)
    ms = json.load(open(os.path.join(HERE, 'meta', 'meta_summary.json'))); e = [x for x in [v for v in ms['variants'] if v['variant'] == 'primary' and v['mode'] == 'adjusted'][0]['estimators'] if x['estimator'].startswith('REML+HK (primary)')][0]
    M = pd.DataFrame([dict(margin=m, verdict=MA.verdict(e['lo'], e['hi'], m)) for m in (0.05, 0.075, 0.10, 0.15, 0.20)]); M.to_csv(os.path.join(OUT, 'margin_sensitivity.csv'), index=False); mmin = max(abs(e['lo']), abs(e['hi']))
    L = ['# WP4 extension A1-A3: positive control of the islet cross-check (reading rules fixed in PRESPEC_WP4b_positive_control_and_wider_search.md before this ran)', '']
    if len(Pp):
        L += ['## A1. Pooled verdict over the five human datasets, per nominal rho (200 replicates; no injection = the primary result, pooled E %+.3f [%+.3f, %+.3f])' % (e['mu'], e['lo'], e['hi']), '',
              '| mode | nominal rho | no meaningful excess | inconclusive | co-expressed more | pooled E mean (SD) | pooled CI mean |', '|---|---|---|---|---|---|---|']
        for r in Pp.itertuples(): L.append('| %s | %.2f | %.2f | %.2f | %.2f | %+.3f (%.3f) | [%+.3f, %+.3f] |' % (r.mode, r.rho, r.share_no_meaningful_excess, r.share_inconclusive, r.share_co_expressed_more, r.pooled_E_mean, r.pooled_E_sd, r.pooled_lo_mean, r.pooled_hi_mean))
        d = Pp[(Pp['mode'] == 'detected')]; fr = {r.rho: r.share_no_meaningful_excess for r in d.itertuples()}; r1 = fr[0.3] <= 0.20
        ok = [r for r in sorted(fr) if fr[r] <= 0.20]
        L += ['', '**Reading rule R1 (mode "detected", nominal rho 0.3, false-reassurance rate <= 20%%): %s** (false-reassurance rate at 0.3: %.2f).' % ('MET: the cross-check is SENSITIVE' if r1 else 'NOT MET', fr[0.3]),
              ('R2 not needed.' if r1 else 'R2: the smallest nominal rho of the grid with a false-reassurance rate <= 20%% is %s.' % (('%.2f' % ok[0]) if ok else 'not reached up to 0.5')), '']
    L += ['## A1. Per dataset (recovered excess and attenuation; mode "detected" unless stated)', '', '| dataset | mode | nominal rho | E primary | E recovered mean (SD) | attenuation | share with empirical P < 0.05 |', '|---|---|---|---|---|---|---|']
    for r in A.itertuples(): L.append('| %s | %s | %.2f | %+.3f | %+.3f (%.3f) | %.2f | %.2f |' % (r.tag, r.mode, r.rho, r.E_primary, r.E_recovered_mean, r.E_recovered_sd, r.attenuation, r.share_empirical_P_lt_0_05))
    if len(R):
        n_ok = int(sum(1 for r in R.itertuples() if r.tag in HUMAN and r.E_RP_mean > 0.10))
        L += ['', '## A2. Ribosomal protein pairs (100 random sets of 11 pairs per dataset; excess E_RP over matched random pairs)', '', '| dataset | RP genes eligible | mean E_RP (SD) | E_RP in SE units of the consensus E | share of sets with E_RP >= 0.10 | mean rho of the RP pairs |', '|---|---|---|---|---|---|']
        for r in R.itertuples(): L.append('| %s | %d | %+.3f (%.3f) | %.1f | %.2f | %+.3f |' % (r.tag, r.rp_genes_eligible, r.E_RP_mean, r.E_RP_sd, r.E_RP_in_SE_units, r.share_E_RP_ge_0_10, r.mean_rho_RP_pairs))
        L += ['', '**Reading rule (mean E_RP > +0.10 in at least 4 of the 5 human datasets): %s** (%d of 5).' % ('MET' if n_ok >= 4 else 'NOT MET', n_ok)]
    L += ['', '## A3. Margin sensitivity of the primary pooled interval [%+.3f, %+.3f] (the pre-set margin stays +/-0.10)' % (e['lo'], e['hi']), '', '| margin | verdict of the primary interval |', '|---|---|'] + ['| %.3f | %s |' % (r.margin, r.verdict) for r in M.itertuples()]
    L += ['', 'The smallest margin at which the verdict is "no meaningful excess" is %.3f: the data exclude an excess in Fisher z above that value.' % mmin]
    open(os.path.join(OUT, 'positive_control_report.md'), 'w', encoding='utf8').write('\n'.join(L) + '\n'); print('\n'.join(L))


if __name__ == '__main__':
    a = sys.argv[1:]
    if a != ['report']:
        for t in (a or TAGS): run_dataset(t)
    report()
