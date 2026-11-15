"""run_datasets.py - protocol 4.2 / 4.3: cell-level analysis of every public islet dataset, unadjusted and depth-adjusted, with the expression-matched null (wp4_core.analyse, gate-verified).
Datasets (all frozen before analysis): gse81608, gse86469, gse81547, emtab5061 (parsed by parse_*.py, human) and GSE84133 human and mouse (earlier loader). Variants: 'primary' (all frozen beta-cells) and, for every
dataset with donors that have type 2 diabetes, 'nondiabetic' (those donors removed, same rules). Output in results/: <tag>__<variant>__<mode>__{pairs,cons,detection}.csv, <tag>__<variant>__summary.json (both modes),
and results/index.csv (one row per dataset, variant and mode). One dataset at a time is held in memory.
Usage: python -I run_datasets.py [tag ...]      tags: gse81608 gse86469 gse81547 emtab5061 gse84133_human gse84133_mouse (default: all available)"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
import wp4_core as C

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'results'); os.makedirs(OUT, exist_ok=True)
BARON = os.path.join(C.DATA_ROOT, 'gse84133', 'parsed')
TAGS = ['gse81608', 'gse84133_human', 'gse86469', 'gse81547', 'emtab5061', 'gse84133_mouse']


def load(tag):
    if tag == 'gse84133_human': return C.load_baron('human', BARON)
    if tag == 'gse84133_mouse': return C.load_baron('mouse', BARON)
    return C.load_parsed(tag)


def is_t2d(d): return np.array([str(x).strip().lower().startswith(('t2d', 'type 2', 'type ii')) for x in d])


def run(tag):
    t0 = time.time(); cells = load(tag); t2d = is_t2d(cells.disease); rows = []
    variants = [('primary', cells)] + ([('nondiabetic', cells.subset(~t2d))] if t2d.any() and (~t2d).sum() > 0 else [])
    for vname, cv in variants:
        res = C.analyse(cv, 0.10, ('unadjusted', 'adjusted'))
        for mode, r in res.items():
            for part in ('pairs', 'cons', 'detection'): r[part].to_csv(os.path.join(OUT, f'{tag}__{vname}__{mode}__{part}.csv'), index=False)
            s = dict(r['summary']); s.update(tag=tag, variant=vname, t2d_beta_cells=int(is_t2d(cv.disease).sum()), kind=cv.kind, depth_covariate='log10 library size' if cv.kind in ('umi', 'count') else 'log10 detected genes')
            rows.append(s)
        json.dump({m: r['summary'] for m, r in res.items()}, open(os.path.join(OUT, f'{tag}__{vname}__summary.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
        for s in rows[-2:]:
            print('%-16s %-11s %-10s cells %5d donors %2d | evaluable genes %2d pairs %2d | median rho %s | null %s | E %s SE %s' % (tag, vname, s['mode'], s['beta_cells'], s['donors'], s['genes_evaluable'], s['consensus_pairs_evaluable'],
                  '%.3f' % s['median_rho_consensus'] if s['median_rho_consensus'] is not None else '-', '%.3f' % s['null_median_rho_mean'] if 'null_median_rho_mean' in s else '-', '%.3f' % s['E'] if 'E' in s else '-', '%.3f' % s['SE_E'] if 'SE_E' in s else '-'), flush=True)
    keep = ['tag', 'variant', 'mode', 'kind', 'depth_covariate', 'beta_cells', 'donors', 't2d_beta_cells', 'df_for_p', 'genes_evaluable', 'consensus_pairs_evaluable', 'other_pairs_evaluable', 'median_rho_consensus', 'median_rho_other', 'null_median_rho_mean',
            'null_median_rho_95pct', 'empirical_p_median_ge_observed', 'empirical_p_mean_ge_observed', 'mean_z_observed', 'mean_z_null', 'E', 'SE_E', 'empirical_p_E', 'consensus_positive', 'consensus_sig_positive_q05', 'consensus_sig_negative_q05', 'edges_q05_rho02']
    I = pd.DataFrame(rows).reindex(columns=keep); print('%s done in %.0f s' % (tag, time.time() - t0), flush=True); return I


if __name__ == '__main__':
    tags = sys.argv[1:] or TAGS
    for tag in tags:
        I = run(tag); p = os.path.join(OUT, f'index_{tag}.csv'); I.to_csv(p, index=False)
    parts = [pd.read_csv(os.path.join(OUT, f'index_{t}.csv')) for t in TAGS if os.path.exists(os.path.join(OUT, f'index_{t}.csv'))]
    pd.concat(parts).to_csv(os.path.join(OUT, 'index.csv'), index=False); print('index.csv written with', len(parts), 'datasets')
