"""WP1 / protocol v2 section 2.1: gene-to-run map of the INS-1 matrix, established by matching VALUES (not labels) against every block of the authors' calculation sheets (read only).
Every column of Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx (21 genes x 12 observations: rows 0-3 CTRL, 4-7 DAPT, 8-11 DKK-1) is compared with every candidate series of
the three workbooks, in sheet row order within group (no well-to-replicate mapping is reconstructed). Writes next to this script: gene_run_map.csv (the map), gene_run_evidence.csv
(every gene x candidate comparison), ins1_sheet_cq.csv (target and reference Cq of the matched sheet rows, used by the shared-reference analysis 2.4).
Usage: python -I build_gene_run_map.py      (env RTPCR_DIR and INS1_MATRIX override the default paths; takes about 10 s)"""
import sys, os, re, itertools
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd, openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
RT = os.environ.get('RTPCR_DIR', 'RT-PCR')
MATRIX = os.environ.get('INS1_MATRIX', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))
WB = {'jun': '2020-Jun-03+24_RT-PCR_Manual Cal_RESULT_WR.xlsx', 'aug': '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx',
      'nov': '2020-Nov-09_RT-PCR_Manual Cal_RESULT (version 1).xlsb.xlsx'}
GENES17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']; DEN = ('Ins1', 'Kcnj11', 'Hes1')
GROUPS = ('CTRL', 'DAPT', 'DKK-1')
RUN_OF = {'jun03': 'Jun03', 'jun24': 'Jun24', 'mix': 'Jun24', 'aug': 'Aug26', 'nov': 'Nov09', 'novbox': 'Nov09'}


def norm(s): return re.sub(r'[^a-z0-9]', '', str(s).lower())


def group_of(label):
    s = str(label).strip().lower()
    return 0 if s.startswith(('ctrl', 'control')) else 1 if s.startswith('dapt') else 2 if s.startswith('dkk') else None


def order_by_group(labels):
    """positions in matrix order: CTRL replicates, DAPT replicates, DKK-1 replicates, each in sheet row order"""
    g = [group_of(l) for l in labels]
    if None in g: return None
    return [i for k in range(3) for i, x in enumerate(g) if x == k]


def num(v): return isinstance(v, (int, float)) and not isinstance(v, bool)


def cand(source, sheet, gene, rows, labels, tcq, ref, cached, formulas):
    o = order_by_group(labels)
    if o is None or len(o) != 12 or not all(num(tcq[i]) and num(ref[i]) for i in range(12)): return None
    t = np.array([tcq[i] for i in o], float); rf = np.array([ref[i] for i in o], float)
    return dict(source=source, sheet=sheet, gene=gene, rows=[rows[i] for i in o], labels=[labels[i] for i in o], tcq=t, ref=rf, rel_re=2.0 ** -(t - rf),
                rel_cached=np.array([cached[i] for i in o], float) if cached is not None and all(num(cached[i]) for i in range(12)) else None, formulas=formulas)


def blocks_by_header(ws, wf, hcol, lab_col, t_col, r_col, src, sheet, cached_col=None, gene_cols=('Q',), formula_cols=()):
    out = []
    for r in range(1, ws.max_row + 1):
        if ws[f'{hcol}{r}'].value != 'Target': continue
        gene = next((ws[f'{c}{r}'].value for c in gene_cols if ws[f'{c}{r}'].value not in (None, '')), None)
        rows = list(range(r + 1, r + 13)); labels = [ws[f'{lab_col}{i}'].value for i in rows]
        if gene is None or any(l is None for l in labels): continue
        f = {c: wf[f'{c}{rows[0]}'].value for c in formula_cols}; f['target_cq_cell'] = wf[f'{t_col}{rows[0]}'].value; f['ref_cq_cell'] = wf[f'{r_col}{rows[0]}'].value
        c = cand(src, sheet, str(gene), rows, labels, [ws[f'{t_col}{i}'].value for i in rows], [ws[f'{r_col}{i}'].value for i in rows],
                 None if cached_col is None else [ws[f'{cached_col}{i}'].value for i in rows], f)
        if c: out.append(c)
    return out


def parse_all():
    wj_v = openpyxl.load_workbook(os.path.join(RT, WB['jun']), data_only=True); wj_f = openpyxl.load_workbook(os.path.join(RT, WB['jun']), data_only=False)
    wa_v = openpyxl.load_workbook(os.path.join(RT, WB['aug']), data_only=True); wa_f = openpyxl.load_workbook(os.path.join(RT, WB['aug']), data_only=False)
    wn_v = openpyxl.load_workbook(os.path.join(RT, WB['nov']), data_only=True); wn_f = openpyxl.load_workbook(os.path.join(RT, WB['nov']), data_only=False)
    C = []
    C += blocks_by_header(wj_v['03 June 2020'], wj_f['03 June 2020'], 'E', 'D', 'E', 'F', 'jun03', '03 June 2020', None, ('Q',), ('G',))
    C += blocks_by_header(wj_v['24 June 2020'], wj_f['24 June 2020'], 'D', 'C', 'D', 'E', 'jun24', '24 June 2020', None, ('Q',), ('F',))
    mix = blocks_by_header(wj_v['Mix_3+24 June 2020'], wj_f['Mix_3+24 June 2020'], 'D', 'C', 'D', 'E', 'mix', 'Mix_3+24 June 2020', 'M', ('Q',), ('F', 'M'))
    C += mix
    ws, wf = wn_v['Box Plot'], wn_f['Box Plot']
    for r in range(1, ws.max_row + 1):                                                       # 9 Nov sheet: gene name in column A (header row or the row below), cached 2^-dCq in column N
        if ws[f'E{r}'].value != 'Target': continue
        gene = ws[f'A{r}'].value or ws[f'A{r + 1}'].value or ws[f'Q{r}'].value; rows = list(range(r + 1, r + 13)); labels = [ws[f'D{i}'].value for i in rows]
        if any(l is None for l in labels): continue
        cached = [ws[f'N{i}'].value for i in rows] if ws[f'N{rows[0]}'].value is not None else None
        f = dict(G=wf[f'G{rows[0]}'].value, N=wf[f'N{rows[0]}'].value, target_cq_cell=wf[f'E{rows[0]}'].value, ref_cq_cell=wf[f'F{rows[0]}'].value)
        c = cand('nov' if cached else 'novbox', "Box Plot", str(gene), rows, labels, [ws[f'E{i}'].value for i in rows], [ws[f'F{i}'].value for i in rows], cached, f)
        if c: C.append(c)
    ws, wf = wa_v['0'], wa_f['0']; blocks = {}
    for r in range(2, ws.max_row + 1):                                                       # 26 Aug summary: one row per well; target in D, sample in F, Cq in G, reference (HK) in J, 2^-dCq in M
        if ws[f'D{r}'].value in (None, 'Gapdh'): continue
        blocks.setdefault(ws[f'D{r}'].value, []).append(r)
    for gene, rows in blocks.items():
        f = dict(J=wf[f'J{rows[0]}'].value, K=wf[f'K{rows[0]}'].value, M=wf[f'M{rows[0]}'].value, target_cq_cell=wf[f'G{rows[0]}'].value, ref_cq_cell=wf[f'J{rows[0]}'].value)
        c = cand('aug', '0', gene, rows, [ws[f'F{i}'].value for i in rows], [ws[f'G{i}'].value for i in rows], [ws[f'J{i}'].value for i in rows], [ws[f'M{i}'].value for i in rows], f)
        if c: C.append(c)
    return C


def ascending(v): return bool(np.all(np.diff(v) >= 0))


def main():
    X = pd.read_excel(MATRIX); C = parse_all(); print('candidate blocks parsed: %d (%s)' % (len(C), ', '.join('%s %d' % (s, sum(c['source'] == s for c in C)) for s in RUN_OF)))
    ev, mp, long = [], [], []
    for gene in X.columns:
        m = X[gene].astype(float).values; rows = []
        for c in C:
            for kind in (['cached'] if c['rel_cached'] is not None else []) + ['recomputed']:
                v = c['rel_cached'] if kind == 'cached' else c['rel_re']
                rows.append(dict(matrix_gene=gene, source=c['source'], run=RUN_OF[c['source']], sheet=c['sheet'], block_gene=c['gene'], first_row=c['rows'][0], value_kind=kind,
                                 max_abs_diff=float(np.abs(v - m).max()), max_rel_diff=float((np.abs(v - m) / np.abs(m)).max())))
        e = pd.DataFrame(rows); ev.append(e); ex = e[e.max_rel_diff <= 1e-9]
        ok = ex.copy(); lab_ok = ok[ok.block_gene.map(norm) == norm(gene)]
        non = e[e.max_rel_diff > 1e-9].sort_values('max_rel_diff'); best = non.iloc[0]
        runs_exact = sorted(set(ok.run)); status = 'resolved' if len(runs_exact) == 1 and len(lab_ok) > 0 else ('unresolved: no exact match' if len(ok) == 0 else 'ambiguous: exact match in %s' % '+'.join(runs_exact))
        # the sheet series used downstream: the cached cell if one exists, else the recomputed 2^-(target - reference) of the same block
        pick = None
        for c in C:
            if RUN_OF[c['source']] in runs_exact and norm(c['gene']) == norm(gene) and (np.abs(c['rel_re'] - m) / m).max() <= 1e-9 and c['source'] in ('mix', 'aug', 'nov'): pick = c; break
        if pick is None:
            for c in C:
                if RUN_OF[c['source']] in runs_exact and norm(c['gene']) == norm(gene) and (np.abs(c['rel_re'] - m) / m).max() <= 1e-9: pick = c; break
        if pick is None: mp.append(dict(gene=gene, status=status)); continue
        used = pick['rel_cached'] if pick['rel_cached'] is not None else pick['rel_re']
        tc, rf = pick['tcq'], pick['ref']; srt_t = sum(ascending(tc[4 * k:4 * k + 4]) for k in range(3)); srt_r = sum(ascending(rf[4 * k:4 * k + 4]) for k in range(3))
        mp.append(dict(gene=gene, analysed_in_17=gene in GENES17, role='denominator' if gene in DEN else ('analysed' if gene in GENES17 else 'excluded from the 17-gene analysis'), run=runs_exact[0] if len(runs_exact) == 1 else '', status=status,
                       workbook=WB['jun'] if pick['source'] in ('jun03', 'jun24', 'mix') else WB[pick['source'][:3]], sheet=pick['sheet'], sheet_rows='%d-%d' % (min(pick['rows']), max(pick['rows'])),
                       value_used='cached column (formula 2^-dCq)' if pick['rel_cached'] is not None else '2^-(target Cq - reference Cq) recomputed from the typed Cq',
                       exact_sources='; '.join(sorted(set('%s %s (%s)' % (r.sheet, r.block_gene, r.value_kind) for r in ok.itertuples()))), n_exact_sources=int(len(ok)),
                       max_abs_diff_used=float(np.abs(used - m).max()), max_rel_diff_used=float((np.abs(used - m) / m).max()), max_rel_diff_recomputed=float((np.abs(pick['rel_re'] - m) / m).max()),
                       best_nonmatching_source='%s %s' % (best.sheet, best.block_gene), best_nonmatching_max_rel_diff=float(best.max_rel_diff),
                       best_3June_plate_block='%s' % e[e.run == 'Jun03'].sort_values('max_rel_diff').iloc[0].block_gene, best_3June_plate_max_rel_diff=float(e[e.run == 'Jun03'].max_rel_diff.min()),
                       target_cq_min=float(tc.min()), target_cq_max=float(tc.max()), ref_cq_min=float(rf.min()), ref_cq_max=float(rf.max()),
                       target_cq_ascending_groups=int(srt_t), ref_cq_ascending_groups=int(srt_r), formulas=str({k: v for k, v in pick['formulas'].items()})))
        for i in range(12):
            long.append(dict(gene=gene, replicate=i, group=GROUPS[i // 4], replicate_in_group=i % 4 + 1, sheet=pick['sheet'], sheet_row=pick['rows'][i], sheet_label=pick['labels'][i], target_cq=tc[i], reference_cq=rf[i],
                             matrix_value=float(m[i]), sheet_value=float(used[i]), recomputed_2_pow_minus_dCq=float(pick['rel_re'][i])))
    mp = pd.DataFrame(mp); long = pd.DataFrame(long); ev = pd.concat(ev, ignore_index=True)
    # reference sets: genes of one run whose 12 per-replicate reference Cq values are identical share one reference vector (one Gapdh value per replicate)
    refset = {}; sets = []; runof = dict(zip(mp.gene, mp.get('run', pd.Series(dtype=str))))
    for g in mp.gene[mp.status == 'resolved']:
        v = long[long.gene == g].reference_cq.values; hit = next((i for i, (r, vec) in enumerate(sets) if r == runof[g] and np.allclose(vec, v, atol=1e-9)), None)
        if hit is None: sets.append((runof[g], v)); hit = len(sets) - 1
        refset[g] = hit
    cnt = {}; name = {}
    for i, (r, _) in enumerate(sets): cnt[r] = cnt.get(r, 0) + 1; name[i] = '%s-ref%s' % (r, 'ABCDEFG'[cnt[r] - 1])
    mp['reference_set'] = [name[refset[g]] if g in refset else '' for g in mp.gene]
    mp.to_csv(os.path.join(HERE, 'gene_run_map.csv'), index=False); ev.to_csv(os.path.join(HERE, 'gene_run_evidence.csv'), index=False); long.to_csv(os.path.join(HERE, 'ins1_sheet_cq.csv'), index=False)
    pd.set_option('display.width', 250, 'display.max_columns', 30, 'display.max_colwidth', 60)
    print(mp[['gene', 'role', 'run', 'status', 'sheet', 'sheet_rows', 'n_exact_sources', 'max_rel_diff_used', 'max_abs_diff_used', 'best_nonmatching_source', 'best_nonmatching_max_rel_diff', 'reference_set', 'target_cq_ascending_groups', 'ref_cq_ascending_groups']].to_string())
    ok = mp.status == 'resolved'; print('\nresolved %d of %d matrix genes; largest relative difference over all resolved genes %.3g, largest absolute difference %.3g' % (ok.sum(), len(mp), mp.loc[ok, 'max_rel_diff_used'].max(), mp.loc[ok, 'max_abs_diff_used'].max()))
    print('runs: ' + ', '.join('%s %d genes' % (r, n) for r, n in mp[ok].run.value_counts().items()) + ' | reference sets: ' + ', '.join('%s %s' % (s, list(mp[mp.reference_set == s].gene)) for s in sorted(set(mp.reference_set))))
    return 0 if ok.sum() >= 11 else 1


if __name__ == '__main__':
    sys.exit(main())
