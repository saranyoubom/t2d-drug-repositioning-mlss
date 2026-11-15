"""drugcentral_sensitivity.py - WP10 part B: the MLSS of the reported run (DrugCentral file of the reported analysis) versus the same MLSS computed with the current DrugCentral
drug-target file (database dump of 25 Sep 2026, downloaded by the user into ../../../Database/). The Python port ../wp3_mlss/mlss_port.py is used unchanged (rules: exact GENE match to the reported hubs,
human target, pIC50 = largest human IC50/Ki/EC50 with relation '=', same Hes1-normalized correlation lookup, same weights); only the interaction file differs.
Gate first: the port on the OLD file reproduces the reported scores (its own gate). Pre-specification: PRESPEC_WP10_database_versions.md (section 4).
Usage: python -I drugcentral_sensitivity.py     (reads only; writes results/mlss_new_all_combinations.csv, mlss_comparison.csv, pic50_changes.csv, drugcentral_summary.json)"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'wp3_mlss'))
import numpy as np, pandas as pd
from scipy import stats
import mlss_port as M

PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
NEW = os.path.join(PROJECT, 'Database', 'drug.target.interaction.tsv.gz')
RES = os.path.join(HERE, 'results')
os.makedirs(RES, exist_ok=True)


def read_new(path):
    dc = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, na_values=[''], low_memory=False, encoding='utf-8-sig', compression='gzip')
    dc['ACT_VALUE'] = pd.to_numeric(dc['ACT_VALUE'], errors='coerce')
    dc['G'] = dc['GENE'].str.upper()
    return dc


def key(a, b):
    return frozenset((a, b))


def main():
    lines = []
    ok, P0, _ = M.gate(lines, None)
    print('GATE on the old file:', 'PASS' if ok else 'FAIL')
    if not ok:
        sys.exit(1)
    hubs, hub_df = M.hubs_from_centrality()
    corr = M.corr_from_matrix_csv()
    old = M.read_drugcentral()
    new = read_new(NEW)
    P1 = M.run_reported(new, hubs, corr)
    P0 = M.run_reported(old, hubs, corr)
    out = {'hubs': hubs, 'rows_old': int(len(old)), 'rows_new': int(len(new)), 'distinct_drugs_old': int(old['DRUG_NAME'].nunique()), 'distinct_drugs_new': int(new['DRUG_NAME'].nunique())}
    t0, t1 = M.targets_from(old, hubs), M.targets_from(new, hubs)
    out['mlss_drugs_old'], out['mlss_drugs_new'] = len(t0), len(t1)
    out['drugs_added'], out['drugs_removed'] = sorted(set(t1) - set(t0)), sorted(set(t0) - set(t1))
    byhub = lambda t, g: sorted(d for d, gs in t.items() if g in gs)
    out['drugs_per_hub_old'] = {g: len(byhub(t0, g)) for g in hubs}
    out['drugs_per_hub_new'] = {g: len(byhub(t1, g)) for g in hubs}
    out['combinations_old'], out['combinations_new'] = int(len(P0)), int(len(P1))
    out['cross_target_old'], out['cross_target_new'] = int((P0['Complementarity_C'] == 1).sum()), int((P1['Complementarity_C'] == 1).sum())
    P0['key'] = [key(a, b) for a, b in zip(P0['Drug_A'], P0['Drug_B'])]
    P1['key'] = [key(a, b) for a, b in zip(P1['Drug_A'], P1['Drug_B'])]
    r0 = dict(zip(P0['key'], P0.to_dict('records')))
    r1 = dict(zip(P1['key'], P1.to_dict('records')))
    shared = [k for k in r0 if k in r1]
    out['combinations_shared'] = len(shared)
    s0 = np.array([r0[k]['MLSS_v4.0'] for k in shared])
    s1 = np.array([r1[k]['MLSS_v4.0'] for k in shared])
    k0 = np.array([r0[k]['Rank'] for k in shared])
    k1 = np.array([r1[k]['Rank'] for k in shared])
    out['spearman_scores_shared'] = round(float(stats.spearmanr(s0, s1)[0]), 4)
    out['spearman_ranks_shared'] = round(float(stats.spearmanr(k0, k1)[0]), 4)
    out['score_changed_shared'] = int((np.abs(s0 - s1) > 1e-9).sum())
    out['max_abs_score_change_shared'] = round(float(np.abs(s0 - s1).max()), 4)
    top0 = [k for k in P0.sort_values('Rank')['key'][:15]]
    top1 = [k for k in P1.sort_values('Rank')['key'][:15]]
    out['top15_old_still_in_top15_new'] = int(sum(1 for k in top0 if k in set(top1)))
    out['top15_old_ranks_in_new'] = [[r0[k]['Combination'], int(r0[k]['Rank']), int(r1[k]['Rank']) if k in r1 else None] for k in top0]
    out['top15_new'] = [[int(r1[k]['Rank']), r1[k]['Combination'], round(float(r1[k]['MLSS_v4.0']), 4)] for k in top1]
    # tie group straddling rank 15 in the new ranking
    sc = P1.sort_values('Rank')['MLSS_v4.0'].round(9).tolist()
    out['tie_straddles_rank15_new'] = bool(sc[14] == sc[15])
    out['tie_straddles_rank15_old'] = bool(P0.sort_values('Rank')['MLSS_v4.0'].round(9).tolist()[14] == P0.sort_values('Rank')['MLSS_v4.0'].round(9).tolist()[15])
    # new drugs: ranks of their combinations
    nd = out['drugs_added']
    mask = np.array([d in nd for d in P1['Drug_A']]) | np.array([d in nd for d in P1['Drug_B']])
    rows_new = P1[mask] if nd else P1.iloc[0:0]
    out['combinations_with_new_drugs'] = int(len(rows_new))
    out['new_drug_combination_ranks'] = [[int(r['Rank']), r['Combination'], round(float(r['MLSS_v4.0']), 4)] for _, r in rows_new.sort_values('Rank').iterrows()][:25]
    out['new_drug_in_top15'] = int(sum(1 for k in top1 if k not in r0))
    # potency changes for drugs and hub targets of the reported run
    pot0, pot1 = M.potency_from(old), M.potency_from(new)
    hubU = {h.upper() for h in hubs}
    chg = []
    for (d, g) in sorted(set(pot0) | set(pot1)):
        if g in hubU and (d in t0 or d in t1):
            a, b = pot0.get((d, g)), pot1.get((d, g))
            if a is None or b is None or abs(a - b) > 1e-9:
                chg.append([d, g, a, b])
    out['pic50_changes'] = len(chg)
    pd.DataFrame(chg, columns=['drug', 'gene', 'pIC50_old', 'pIC50_new']).to_csv(os.path.join(RES, 'pic50_changes.csv'), index=False)
    # full tables
    cols = [c for c in ['Rank', 'Combination', 'MLSS_v4.0', 'Complementarity_C', 'Balance_B', 'Potency_P', 'Mean_pIC50', 'Targets_A', 'Targets_B', 'Drug_A', 'Drug_B'] if c in P1.columns]
    P1[cols].to_csv(os.path.join(RES, 'mlss_new_all_combinations.csv'), index=False)
    union = sorted(set(r0) | set(r1), key=lambda k: (r1[k]['Rank'] if k in r1 else 10 ** 6, r0[k]['Rank'] if k in r0 else 10 ** 6))
    cmp_rows = []
    for k in union:
        a = r0.get(k)
        b = r1.get(k)
        cmp_rows.append([(b if b is not None else a)['Combination'], int(a['Rank']) if a is not None else '', '%.6f' % a['MLSS_v4.0'] if a is not None else '', int(b['Rank']) if b is not None else '', '%.6f' % b['MLSS_v4.0'] if b is not None else '',
                         'both' if (a is not None and b is not None) else ('new only' if b is not None else 'old only')])
    pd.DataFrame(cmp_rows, columns=['combination', 'rank_old', 'MLSS_old', 'rank_new', 'MLSS_new', 'status']).to_csv(os.path.join(RES, 'mlss_comparison.csv'), index=False)
    with open(os.path.join(RES, 'drugcentral_summary.json'), 'w', encoding='utf8') as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:6000])


if __name__ == '__main__':
    main()
