"""string_web_export_check.py - WP10 addendum (8 Nov 2026): the user's own STRING 12.5 web export of 8 Nov 2026 compared with the two STRING tables already used.
Inputs (read only):
  Database/STRING12.5/string_interactions_21_nodes.tsv   user's web export, 8 Nov 2026 (file times 13:16-13:23); network URL in STRING12.5.txt (version-12-5.string-db.org, networkId beSoaw2hyqcw)
  R_17genes/string_interactions.csv                      STRING 12.0 web export that Fig. 5 used (original download date not recorded; copy date 2 Nov 2026)
  wp10_versions/results/string_v12_5_network.tsv         12.5 network retrieved through the API by string_sensitivity.py (rat 10116, 17 + 4 genes, minimum score 0.4)
Questions: (1) do the web export and the API network of 12.5 agree (edges, scores)? (2) what changed from 12.0 to 12.5 in the web exports? (3) does the Fig. 5 classification computed from the web export equal the API-based result of string_sensitivity.py?
Usage: python -I string_web_export_check.py   (writes results/string_web_export_check.json)"""
import os, sys, csv, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import string_sensitivity as ss
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
WEB = os.path.join(P, 'Database', 'STRING12.5', 'string_interactions_21_nodes.tsv')
URLF = os.path.join(P, 'Database', 'STRING12.5', 'STRING12.5.txt')
API = os.path.join(HERE, 'results', 'string_v12_5_network.tsv')
OUT = os.path.join(HERE, 'results', 'string_web_export_check.json')
G17 = {'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r'}


def read_web(path):
    rows = []
    with open(path, encoding='utf-8-sig', newline='') as f:
        rd = csv.reader(f, delimiter='\t')
        hdr = [h.lstrip('#') for h in next(rd)]
        for r in rd:
            if r: rows.append(dict(zip(hdr, r)))
    return rows


web = read_web(WEB)
old = ss.load_inputs()[0]
api = [dict(zip(l.split('\t')[0:0] or [], [])) for l in []]
lines = open(API, encoding='utf8').read().strip().splitlines()
ah = lines[0].split('\t')
api = [dict(zip(ah, l.split('\t'))) for l in lines[1:]]
alias = lambda x: ss.ALIAS.get(x, x)
out = {'network_url_file': open(URLF, encoding='utf8').read().strip()}


def summarize(rows, a='node1', b='node2', sc='combined_score'):
    nodes = sorted({r[a] for r in rows} | {r[b] for r in rows})
    sc_by = {}
    for r in rows:
        k = ss.pname(r[a], r[b]); v = float(r[sc])
        sc_by.setdefault(k, set()).add(v)
    return nodes, sc_by


wn, ws = summarize(web)
on, os_ = summarize(old)
an, as_ = summarize(api, 'preferredName_A', 'preferredName_B', 'score')
out['rows'] = {'web_12_5': len(web), 'old_12_0': len(old), 'api_12_5': len(api)}
out['nodes'] = {'web_12_5': wn, 'old_12_0': on, 'api_12_5': an}
out['distinct_pairs'] = {'web_12_5': len(ws), 'old_12_0': len(os_), 'api_12_5': len(as_)}
out['pairs_with_two_different_scores_in_one_file'] = {'web_12_5': sorted(k for k, v in ws.items() if len(v) > 1), 'old_12_0': sorted(k for k, v in os_.items() if len(v) > 1), 'api_12_5': sorted(k for k, v in as_.items() if len(v) > 1)}
sc = lambda d: {k: sorted(v)[0] for k, v in d.items()}
W, O, A = sc(ws), sc(os_), sc(as_)
out['min_score_in_file'] = {'web_12_5': min(W.values()), 'old_12_0': min(O.values()), 'api_12_5': min(A.values())}
# (1) web 12.5 against API 12.5
out['web_vs_api_12_5'] = {'pairs_only_in_web': sorted(set(W) - set(A)), 'pairs_only_in_api': sorted(set(A) - set(W)), 'pairs_in_both': len(set(W) & set(A)),
                          'max_abs_score_difference_in_common_pairs': round(max([abs(W[k] - A[k]) for k in set(W) & set(A)] or [0]), 4)}
# (2) web 12.5 against old 12.0
common = set(W) & set(O)
d = sorted(((W[k] - O[k]), k) for k in common)
out['web_12_5_vs_old_12_0'] = {'pairs_gained_in_12_5': sorted(set(W) - set(O)), 'pairs_lost_in_12_5': sorted(set(O) - set(W)), 'pairs_in_both': len(common),
                              'score_change_in_common_pairs': {'n_changed': sum(1 for x, _ in d if abs(x) > 0.0005), 'min': round(d[0][0], 3), 'max': round(d[-1][0], 3), 'min_pair': d[0][1], 'max_pair': d[-1][1]}}
# (3) classification from the web export, as the R code does it (literal names) and with the alias
old_, cons, pairs = ss.load_inputs()
edges_web = [(r['node1'], r['node2']) for r in web]
res = {}
for name, al in [('as_reported', False), ('alias', True)]:
    cl = ss.classify(pairs, cons, ss.string_set(edges_web, al))
    c = ss.counts(cl); c['in_STRING'] = c['Known Validated'] + c['STRING Only']; c['coverage_percent'] = round(100.0 * c['in_STRING'] / len(pairs), 1)
    res[name] = c
out['classification_from_web_export_12_5'] = res
sj = json.load(open(os.path.join(HERE, 'results', 'string_summary.json'), encoding='utf8'))
out['classification_from_api_12_5_saved'] = {k: sj[k]['v12_5'] for k in ('as_reported', 'alias')}
out['classification_web_equals_api'] = all(res[k] == sj[k]['v12_5'] for k in ('as_reported', 'alias'))
# edges among the 17 analyzed genes (alias)
def e17(edges): return {ss.pname(alias(a), alias(b)) for a, b in edges if alias(a) in G17 and alias(b) in G17}
ew, ea, eo = e17(edges_web), e17([(r['preferredName_A'], r['preferredName_B']) for r in api]), e17([(r['node1'], r['node2']) for r in old])
out['edges_among_17_genes'] = {'web_12_5': len(ew), 'api_12_5': len(ea), 'old_12_0': len(eo), 'web_equals_api': ew == ea, 'gained_vs_12_0': sorted(ew - eo), 'lost_vs_12_0': sorted(eo - ew)}
json.dump(out, open(OUT, 'w', encoding='utf8'), indent=1)
print(json.dumps(out, indent=1)[:6000])
