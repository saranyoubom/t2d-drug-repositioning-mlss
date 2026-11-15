"""wp4_core.py - shared loaders and statistics of WP4 (public islet single-cell cross-check, version 2; protocol v2 section 4). No side effects at import.
Same statistics as ../../human_islet_GSE81608/human_islet_check.py and ../../islet_GSE84133/islet_check_baron.py (log-normalised expression of beta-cells, donor-centring, rank correlation, genes evaluable when
detected in >= 10% of beta-cells, expression-matched null of random gene pairs from the same detection deciles, 2000 draws seed 2026, pair-level seed 2027), with three additions of protocol 4.2 / 4.4:
  * depth adjustment: per gene, the donor-centred log expression is regressed on the donor-centred depth covariate and its square (common slopes, i.e. donor fixed effects + depth + depth^2), the residuals are
    ranked and correlated (partial Spearman); depth = log10 library size for count data (UMI or read counts) and log10 number of detected genes for RPKM data. df of the p values = cells - donors - 4.
  * excess E = mean Fisher z of the evaluable consensus pairs minus the mean over the 2000 matched random pair sets of their mean Fisher z; SE = SD over the 2000 draws.
  * the matched null is rebuilt on the adjusted ranks with the SAME random draws as the unadjusted one (the draws depend only on the gene pools).
The random draws replicate the loop of the earlier scripts call for call (so that the earlier numbers are reproduced exactly); legacy numeric types are kept (float32 for UMI data).
Data container Cells: raw beta-cell values (sparse, cells x genes), kind 'umi' | 'count' | 'rpkm', donor, disease, library size, number of detected genes."""
import os, json
import numpy as np, pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
R17 = os.environ.get('R17_DIR', os.path.normpath(os.path.join(HERE, '..', '..', 'R_17genes')))
DATA_ROOT = os.environ.get('V2_DATA_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_cache'))
PARSED = os.path.join(DATA_ROOT, 'v2_data', 'wp4', 'parsed')
NN, SEED_NULL, SEED_PAIR = 2000, 2026, 2027

GENES = ['Hes1', 'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']     # order as in the earlier scripts
ENTREZ = {'Hes1': 3280, 'Wnt5a': 7474, 'Wnt5b': 81029, 'Wnt9a': 7483, 'Lef1': 51176, 'Tcf7': 6932, 'Tcf7l2': 6934, 'Itpr1': 3708, 'Ryr2': 6262, 'Ryr3': 6263,
          'Cacna1c': 775, 'Cacna1d': 776, 'Kcnj11': 3767, 'Ptbp1': 5725, 'Glut2': 6514, 'Ins1': 3630, 'Glp1r': 2740}      # human orthologs (Glut2 = SLC2A2, Ins1 = INS)
SYMBOL = {'human': lambda g: {'Glut2': 'SLC2A2', 'Ins1': 'INS'}.get(g, g.upper()), 'mouse': lambda g: {'Glut2': 'Slc2a2'}.get(g, g)}


_ENS = {}


def ensembl_map():
    """human Ensembl gene ids of the 17 genes (ensembl_ids_17genes.csv, from the NCBI Datasets gene report; used for datasets whose gene labels are Ensembl ids, id_type 'ensembl')"""
    if not _ENS:
        d = pd.read_csv(os.path.join(HERE, 'ensembl_ids_17genes.csv')); _ENS.update(dict(zip(d.gene, d.ensembl)))
    return _ENS


def consensus_tables():
    """the 21 consensus pairs (in the order of Fig_4d_consensus_pairs.csv), their set, and the in vitro edge r (Hes1-normalised) and centrality degree"""
    cons = [tuple(x.split('_')) for x in pd.read_csv(os.path.join(R17, 'Fig_4d_consensus_pairs.csv')).iloc[:, 0]]
    edge = pd.read_csv(os.path.join(R17, 'Fig_6b_Edge_Correlation_Statistics.csv')); redge = {frozenset((r['Gene1'], r['Gene2'])): r['Correlation (r)'] for _, r in edge.iterrows()}
    cen = pd.read_csv(os.path.join(R17, 'Fig_6b_Gene_Centrality_Summary.csv')).set_index('Gene')['Degree']
    return cons, {frozenset(c) for c in cons}, redge, cen


class Cells:
    """beta-cells of one dataset: B sparse (cells x genes) of raw counts (kind umi / count) or RPKM; id_type 'entrez' (gene labels are Entrez ids), 'ensembl' (Ensembl gene ids) or 'symbol'"""
    def __init__(self, name, species, kind, id_type, genes, B, donor, disease, tot, ndet):
        self.name, self.species, self.kind, self.id_type = name, species, kind, id_type
        self.genes = list(genes); self.gi = {g: i for i, g in enumerate(self.genes)}; self.B = sparse.csr_matrix(B)
        self.donor = np.asarray(donor); self.disease = np.asarray(disease); self.tot = np.asarray(tot, float); self.ndet = np.asarray(ndet, float)
        assert self.B.shape == (len(self.donor), len(self.genes)) and len(self.tot) == len(self.donor) == len(self.ndet) == len(self.disease)

    @property
    def n(self): return self.B.shape[0]

    @property
    def n_donors(self): return len(set(self.donor.tolist()))

    @property
    def depth(self):
        """protocol 4.2: log10 library size for UMI / read-count data, log10 number of detected genes for RPKM data"""
        return np.log10(self.tot) if self.kind in ('umi', 'count') else np.log10(self.ndet)

    def label(self, g):
        if self.id_type == 'entrez': return ENTREZ[g]
        if self.id_type == 'ensembl': return ensembl_map()[g]
        return SYMBOL[self.species](g)

    def subset(self, mask, name=None):
        m = np.asarray(mask, bool)
        return Cells(name or self.name, self.species, self.kind, self.id_type, self.genes, self.B[m], self.donor[m], self.disease[m], self.tot[m], self.ndet[m])


def load_parsed(tag):
    """unified parsed format written by the parse_*.py scripts: <PARSED>/<tag>/beta.npz (CSR, cells x genes) + meta.json"""
    d = os.path.join(PARSED, tag); m = json.load(open(os.path.join(d, 'meta.json')))
    return Cells(m['name'], m['species'], m['kind'], m['id_type'], m['genes'], sparse.load_npz(os.path.join(d, 'beta.npz')), m['donor'], m['disease'], m['tot'], m['ndet'])


def load_baron(species, base):
    """GSE84133 from the files written by ../../islet_GSE84133/load_beta_counts.py (<donor>_beta_counts.npz + <donor>_meta.json): same donor order and identity check as islet_check_baron.run;
    human4 has type 2 diabetes according to the GEO sample records"""
    donors = {'human': ['human1', 'human2', 'human3', 'human4'], 'mouse': ['mouse1', 'mouse2']}[species]; mats, genes0, tots, donor = [], None, [], []
    for d in donors:
        m = json.load(open(os.path.join(base, f'{d}_meta.json'))); B = sparse.load_npz(os.path.join(base, f'{d}_beta_counts.npz')); cl = m['clusters']
        ins_b, ins_o = cl['beta']['INS_cp10k'], max(v['INS_cp10k'] for c, v in cl.items() if c != 'beta'); gcg_a, gcg_o = cl['alpha']['GCG_cp10k'], max(v['GCG_cp10k'] for c, v in cl.items() if c != 'alpha')
        assert ins_b > 2 * ins_o and gcg_a > 2 * gcg_o, f'identity check failed in {d}'
        if genes0 is not None: assert m['genes'] == genes0, 'gene lists differ between donors'
        genes0 = m['genes']; mats.append(B); tots += m['beta_total_umi']; donor += [d] * B.shape[0]
    B = sparse.vstack(mats).tocsr(); donor = np.array(donor); disease = np.where(donor == 'human4', 'T2D', 'non-diabetic')
    return Cells('GSE84133_' + species, species, 'umi', 'symbol', genes0, B, donor, disease, np.array(tots, float), np.asarray((B > 0).sum(axis=1)).ravel())


# ------------------------------------------------------------------------------------------------------------------------------------------ expression frame
def lognorm(cells, elig):
    """dense log2-normalised expression of the eligible genes (cells x genes) in the numeric type of the earlier scripts: UMI -> log2(CP10K + 1) in float32; read counts -> log2(CPM + 1); RPKM -> log2(RPKM + 1)"""
    D = cells.B[:, elig].toarray()
    if cells.kind == 'umi': return np.log2(D.astype(np.float32) / cells.tot[:, None].astype(np.float32) * 1e4 + 1.0)
    if cells.kind == 'count': return np.log2(D.astype(np.float64) / cells.tot[:, None] * 1e6 + 1.0)
    return np.log2(D.astype(np.float64) + 1.0)


def center_by_donor(X, donor):
    Xc = X.copy()
    for d in dict.fromkeys(donor.tolist()):
        m = donor == d; Xc[m] -= Xc[m].mean(axis=0, keepdims=True)
    return Xc


def rank_std(Xc_part):
    """cells x genes chunk -> genes x cells matrix of mean-centred, unit-norm ranks (average ranks for ties); the dot product of two rows is the Spearman rho (the 1e-12 of the earlier scripts is kept)"""
    R = stats.rankdata(np.asarray(Xc_part, np.float64), axis=0).T
    R = R - R.mean(1, keepdims=True)
    return R / (np.sqrt((R ** 2).sum(1, keepdims=True)) + 1e-12)


def ranks_of(Xc, chunk=1500):
    out = np.empty((Xc.shape[1], Xc.shape[0]), np.float64)
    for s in range(0, Xc.shape[1], chunk): out[s:s + chunk] = rank_std(Xc[:, s:s + chunk])
    return out


def depth_design(depth, donor):
    """Z (cells x 2): donor-centred depth and donor-centred squared depth (unit variance); with the donor-centred expression this is the donor-fixed-effects model with common depth slopes"""
    d = np.asarray(depth, float).copy()
    for g in dict.fromkeys(donor.tolist()):
        m = donor == g; d[m] -= d[m].mean()
    d2 = d ** 2
    for g in dict.fromkeys(donor.tolist()):
        m = donor == g; d2[m] -= d2[m].mean()
    return np.column_stack([d / d.std(), d2 / d2.std()])


def adjusted_ranks(Xc, Z, chunk=1500):
    """ranks of the residuals of Xc (cells x genes, donor-centred) after regression on Z, genes x cells, standardised as rank_std"""
    Zp = np.linalg.pinv(Z); out = np.empty((Xc.shape[1], Xc.shape[0]), np.float64)
    for s in range(0, Xc.shape[1], chunk):
        part = np.asarray(Xc[:, s:s + chunk], np.float64); out[s:s + chunk] = rank_std(part - Z @ (Zp @ part))
    return out


# ------------------------------------------------------------------------------------------------------------------------------------------ matched null
def make_pool(det_e, labels_e, own):
    """detection deciles of the eligible genes (rank-based qcut as before) and, per decile, the positions of the eligible genes that are not among the 17 analysed genes (ascending positions)"""
    bins = pd.qcut(pd.Series(det_e).rank(method='first'), 10, labels=False).to_numpy()
    pool = {k: np.array([i for i in range(len(labels_e)) if bins[i] == k and labels_e[i] not in own]) for k in range(10)}
    return bins, pool


def draws_main(pool, bins, pair_idx, nn=NN, seed=SEED_NULL):
    """index draws of the main null, in the order of the earlier loop: for each draw, for each evaluable consensus pair: ia, ic (ic redrawn while equal to ia); returns two (nn x pairs) arrays"""
    rng = np.random.default_rng(seed); IA = np.empty((nn, len(pair_idx)), np.int64); IC = np.empty_like(IA)
    for s in range(nn):
        for j, (a, c) in enumerate(pair_idx):
            ia = rng.choice(pool[bins[a]]); ic = rng.choice(pool[bins[c]])
            while ic == ia: ic = rng.choice(pool[bins[c]])
            IA[s, j], IC[s, j] = ia, ic
    return IA, IC


def draws_pair(pool, bins, pair_idx, nn=NN, seed=SEED_PAIR):
    """index draws of the pair-level null (pair-major order of the earlier loop, separate generator)"""
    rng = np.random.default_rng(seed); IA = np.empty((nn, len(pair_idx)), np.int64); IC = np.empty_like(IA)
    for j, (a, c) in enumerate(pair_idx):
        for k in range(nn):
            ia = rng.choice(pool[bins[a]]); ic = rng.choice(pool[bins[c]])
            while ic == ia: ic = rng.choice(pool[bins[c]])
            IA[k, j], IC[k, j] = ia, ic
    return IA, IC


def null_values(Rk, IA, IC):
    """Spearman rho of the drawn gene pairs (nn x pairs)"""
    out = np.empty(IA.shape)
    for j in range(IA.shape[1]): out[:, j] = np.einsum('ij,ij->i', Rk[IA[:, j]], Rk[IC[:, j]])
    return out


# ------------------------------------------------------------------------------------------------------------------------------------------ analysis of one dataset
def analyse(cells, thr=0.10, modes=('unadjusted', 'adjusted'), nn=NN):
    """Full cell-level analysis of one Cells object. Returns {mode: dict(pairs=DataFrame 136 rows, cons=DataFrame of evaluable consensus pairs with pair-level null, detection=DataFrame 17 genes, summary=dict)}.
    Adjusted mode needs cells.depth. The draws of the nulls are made once and used for both modes."""
    cons, consset, redge, cen = consensus_tables(); sym = cells.label; n = cells.n; nd = cells.n_donors
    det = np.asarray((cells.B > 0).mean(axis=0)).ravel(); elig = np.nonzero(det >= thr)[0]
    labels_e = [cells.genes[j] for j in elig]; idx = {g: r for r, g in enumerate(labels_e)}; det_e = det[elig]
    own = {sym(g) for g in GENES}; bins, pool = make_pool(det_e, labels_e, own)
    Xc = center_by_donor(lognorm(cells, elig), cells.donor)
    # detection of the 17 genes (all genes of the data, evaluable or not)
    D = pd.DataFrame([dict(gene_invitro=g, gene=sym(g), in_data=sym(g) in cells.gi, beta_cells=n, detected_pct=100 * det[cells.gi[sym(g)]] if sym(g) in cells.gi else np.nan, evaluable=sym(g) in idx) for g in GENES])
    genes_pairs = [(a, c) for i, a in enumerate(GENES) for c in GENES[i + 1:]]
    out = {}; draws = None
    for mode in modes:
        if mode == 'unadjusted': Rk = ranks_of(Xc); dfree = n - nd - 2
        elif mode == 'adjusted': Rk = adjusted_ranks(Xc, depth_design(cells.depth, cells.donor)); dfree = n - nd - 4
        else: raise ValueError(mode)
        rows = []
        for a, c in genes_pairs:
            ea, ec = sym(a), sym(c); evb = (ea in idx) and (ec in idx); r_ = float(Rk[idx[ea]] @ Rk[idx[ec]]) if evb else np.nan; t_ = r_ * np.sqrt(dfree / (1 - r_ ** 2)) if evb else np.nan
            rows.append(dict(gene1=a, gene2=c, consensus_invitro=frozenset((a, c)) in consset, r_invitro_Hes1norm=redge.get(frozenset((a, c)), np.nan), evaluable=evb,
                             detected_pct_1=100 * det[cells.gi[ea]] if ea in cells.gi else np.nan, detected_pct_2=100 * det[cells.gi[ec]] if ec in cells.gi else np.nan, rho=r_, p=2 * stats.t.sf(abs(t_), dfree) if evb else np.nan))
        P = pd.DataFrame(rows); ev = P[P.evaluable].copy()
        if len(ev): ev['q'] = multipletests(ev.p, method='fdr_bh')[1]; P = P.merge(ev[['gene1', 'gene2', 'q']], on=['gene1', 'gene2'], how='left')
        else: P['q'] = np.nan
        c_ev, o_ev = ev[ev.consensus_invitro], ev[~ev.consensus_invitro]; rc = c_ev.rho; ro = o_ev.rho
        s = dict(dataset=cells.name, mode=mode, detection_threshold=thr, beta_cells=int(n), donors=int(nd), df_for_p=int(dfree), genes_in_data=int(D.in_data.sum()), genes_evaluable=int(D.evaluable.sum()),
                 genes_not_evaluable=[g for g, e in zip(D.gene_invitro, D.evaluable) if not e], eligible_genes_in_pool=int(len(elig)), consensus_pairs_total=len(cons), consensus_pairs_evaluable=int(len(c_ev)),
                 other_pairs_evaluable=int(len(o_ev)), median_rho_consensus=float(rc.median()) if len(rc) else None, median_rho_other=float(ro.median()) if len(ro) else None,
                 mannwhitney_p_consensus_vs_other=float(stats.mannwhitneyu(rc, ro).pvalue) if len(rc) > 1 and len(ro) > 1 else None, consensus_positive=int((rc > 0).sum()),
                 consensus_sig_positive_q05=int(((rc > 0) & (c_ev.q < 0.05)).sum()) if len(rc) else 0, consensus_sig_negative_q05=int(((rc < 0) & (c_ev.q < 0.05)).sum()) if len(rc) else 0, consensus_rho_ge_0_2=int((rc >= 0.2).sum()))
        cn = c_ev[['gene1', 'gene2', 'r_invitro_Hes1norm', 'rho', 'p', 'q']].copy().reset_index(drop=True)
        if len(c_ev) >= 1:
            ev_pairs = [(sym(r.gene1), sym(r.gene2)) for r in c_ev.itertuples()]; pidx = [(idx[a], idx[c]) for a, c in ev_pairs]
            if draws is None: draws = (draws_main(pool, bins, pidx, nn, SEED_NULL), draws_pair(pool, bins, pidx, nn, SEED_PAIR))
            (IA, IC), (IA2, IC2) = draws
            V = null_values(Rk, IA, IC); meds, means = np.median(V, axis=1), V.mean(axis=1); zm = np.arctanh(V).mean(axis=1)
            obs_med, obs_mean = float(rc.median()), float(rc.mean()); z_obs = float(np.arctanh(rc.to_numpy()).mean())
            s.update(null_median_rho_mean=float(meds.mean()), null_median_rho_95pct=float(np.percentile(meds, 95)), observed_median_rho=obs_med,
                     empirical_p_median_ge_observed=float((1 + (meds >= obs_med).sum()) / (nn + 1)), empirical_p_mean_ge_observed=float((1 + (means >= obs_mean).sum()) / (nn + 1)),
                     mean_z_observed=z_obs, mean_z_null=float(zm.mean()), E=z_obs - float(zm.mean()), SE_E=float(zm.std(ddof=1)), empirical_p_E=float((1 + (zm >= z_obs).sum()) / (nn + 1)))
            V2 = null_values(Rk, IA2, IC2)
            cn['null_mean'] = V2.mean(axis=0); cn['null_95pct'] = np.percentile(V2, 95, axis=0); cn['empirical_p_ge_observed'] = (1 + (V2 >= cn.rho.to_numpy()[None, :]).sum(axis=0)) / (nn + 1)
        # hubs: degree in the beta-cell network (q < 0.05 and |rho| >= 0.2) vs in vitro degree
        E_ = ev[(ev.q < 0.05) & (ev.rho.abs() >= 0.2)] if len(ev) else ev; hd = {g: int(((E_.gene1 == g) | (E_.gene2 == g)).sum()) for g in GENES}
        s.update(edges_q05_rho02=int(len(E_)), invitro_hubs_degree={g: hd[g] for g in ('Glp1r', 'Tcf7l2', 'Cacna1d', 'Ptbp1')})
        pr = P[((P.gene1 == 'Glp1r') & (P.gene2 == 'Cacna1d')) | ((P.gene1 == 'Cacna1d') & (P.gene2 == 'Glp1r'))].iloc[0]
        s['Cacna1d_Glp1r'] = dict(evaluable=bool(pr.evaluable), rho=float(pr.rho) if pr.evaluable else None, p=float(pr.p) if pr.evaluable else None,
                                  detected_pct_Glp1r=float(pr.detected_pct_1 if pr.gene1 == 'Glp1r' else pr.detected_pct_2), detected_pct_Cacna1d=float(pr.detected_pct_1 if pr.gene1 == 'Cacna1d' else pr.detected_pct_2))
        out[mode] = dict(pairs=P, cons=cn, detection=D, summary=s); del Rk
    return out
