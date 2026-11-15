"""string_sensitivity.py - WP10 part A: does the STRING classification of Fig. 5 change with the current STRING release (12.5) instead of 12.0?
Reproduces the classification of R_17genes/Fig_5abcd_Gene_Interaction_Network_Analysis.r in Python from the reported inputs (read only):
  string_interactions.csv (STRING 12.0 web export, rat, minimum score 0.4), Fig_4d_consensus_pairs.csv, Fig_4e_Hes1_significant_correlations_only.csv.
Categories (R code): Known Validated = in STRING and four-method consensus; Novel Discovery = not in STRING but consensus; STRING Only = in STRING, not consensus; Not Validated = neither,
over the pairs with |r| >= 0.7 in the Hes1-normalized significant correlations (120 pairs). Pair names are the two gene names sorted and joined by a hyphen.
NOTE (found while porting): the correlation matrix names the gene Glut2, the STRING export names it Slc2a2, and the R code matches names literally, so no pair with Glut2 can ever be 'in STRING' in the reported analysis.
Both versions are therefore reported: 'as reported' (literal names, the R behavior) and 'alias' (Slc2a2 read as Glut2).
Usage (python -I): python -I string_sensitivity.py gate | retrieve | analyze     (gate and analyze read only saved files; retrieve needs the network)"""
import os, sys, csv, json, time, urllib.request, urllib.parse
HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(os.path.dirname(HERE))
R = os.path.join(ANALYSIS, 'R_17genes')
RES = os.path.join(HERE, 'results')
os.makedirs(RES, exist_ok=True)
STABLE = 'https://version-12-5.string-db.org'
ALIAS = {'Slc2a2': 'Glut2'}


def pname(a, b):
    return '-'.join(sorted([a, b]))


def load_inputs():
    with open(os.path.join(R, 'string_interactions.csv'), encoding='utf-8-sig', newline='') as f:
        old = list(csv.DictReader(f))
    with open(os.path.join(R, 'Fig_4d_consensus_pairs.csv'), encoding='utf-8-sig', newline='') as f:
        cons = {pname(*r['Consensus_Pairs'].split('_')[:2]) for r in csv.DictReader(f)}
    with open(os.path.join(R, 'Fig_4e_Hes1_significant_correlations_only.csv'), encoding='utf-8-sig', newline='') as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
        mat = {}
        for row in rd:
            for g2, v in zip(hdr, row[1:]):
                if v not in ('', 'NA'):
                    mat[(row[0], g2)] = float(v)
    pairs = {}
    for (a, b), v in mat.items():
        if a != b and abs(v) >= 0.7:
            pairs[pname(a, b)] = v
    return old, cons, pairs


def string_set(edges, alias):
    s = set()
    for a, b in edges:
        if alias:
            a, b = ALIAS.get(a, a), ALIAS.get(b, b)
        s.add(pname(a, b))
    return s


def classify(pairs, cons, sset):
    out = {}
    for p in pairs:
        ins, val = p in sset, p in cons
        out[p] = 'Known Validated' if ins and val else ('Novel Discovery' if (not ins and val) else ('STRING Only' if ins else 'Not Validated'))
    return out


def counts(cl):
    c = {}
    for v in cl.values():
        c[v] = c.get(v, 0) + 1
    return {k: c.get(k, 0) for k in ['Known Validated', 'Novel Discovery', 'STRING Only', 'Not Validated']}


def gate():
    old, cons, pairs = load_inputs()
    edges = [(r['node1'], r['node2']) for r in old]
    out = {'pairs_abs_r_ge_0.7': len(pairs), 'consensus_pairs': len(cons), 'string_export_edges': len(edges), 'string_export_distinct_pairs': len(string_set(edges, False))}
    for name, alias in [('as_reported', False), ('alias', True)]:
        cl = classify(pairs, cons, string_set(edges, alias))
        out[name] = counts(cl)
        out[name]['in_STRING'] = out[name]['Known Validated'] + out[name]['STRING Only']
        out[name]['coverage_percent'] = round(100.0 * out[name]['in_STRING'] / len(pairs), 1)
    out['pairs_with_Glut2'] = sum(1 for p in pairs if 'Glut2' in p.split('-'))
    out['export_edges_with_Slc2a2_among_matrix_genes'] = sorted(pname(a, b) for a, b in edges if 'Slc2a2' in (a, b))
    with open(os.path.join(RES, 'string_gate.json'), 'w', encoding='utf8') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    rep = out['as_reported']
    ok = (len(pairs) == 120 and rep['Known Validated'] == 5 and rep['Novel Discovery'] == 16 and rep['STRING Only'] == 29 and rep['Not Validated'] == 70 and rep['in_STRING'] == 34)
    print('GATE (reported numbers 120 pairs; 5 / 16 / 29 / 70; 34 in STRING = 28.3%%): %s' % ('PASS' if ok else 'FAIL'))
    return ok


def get(url, tries=6):
    delay = 3
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'research-version-check/1.0'})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode('utf8')
        except Exception as e:
            print('  retry %d after %s' % (k + 1, str(e)[:100]), flush=True)
            time.sleep(delay)
            delay = min(delay * 2, 60)
    raise RuntimeError('giving up on ' + url)


GENES = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Slc2a2', 'Ins1', 'Glp1r', 'Wnt2', 'Wnt2b', 'Hey1', 'Ryr1']


def retrieve():
    v = get(STABLE + '/api/json/version')
    open(os.path.join(RES, 'string_version.json'), 'w', encoding='utf8').write(v)
    print('version endpoint:', v.strip())
    ids = '%0d'.join(urllib.parse.quote(g) for g in GENES)
    url = STABLE + '/api/tsv/network?identifiers=' + ids + '&species=10116&required_score=400&caller_identity=wp10_version_check'
    t = get(url)
    open(os.path.join(RES, 'string_v12_5_network.tsv'), 'w', encoding='utf8').write(t)
    print('rows (with header):', len(t.strip().splitlines()))


def analyze():
    old, cons, pairs = load_inputs()
    oedges = [(r['node1'], r['node2']) for r in old]
    lines = open(os.path.join(RES, 'string_v12_5_network.tsv'), encoding='utf8').read().strip().splitlines()
    hdr = lines[0].split('\t')
    ia, ib, isc = hdr.index('preferredName_A'), hdr.index('preferredName_B'), hdr.index('score')
    nedges = []
    nscore = {}
    for ln in lines[1:]:
        c = ln.split('\t')
        nedges.append((c[ia], c[ib]))
        nscore[pname(c[ia], c[ib])] = float(c[isc])
    oscore = {pname(r['node1'], r['node2']): float(r['combined_score']) for r in old}
    out = {'version_response': json.load(open(os.path.join(RES, 'string_version.json'), encoding='utf8')), 'edges_v12_0': len(set(oscore)), 'edges_v12_5': len(set(nscore))}
    rows = []
    cls = {}
    for name, alias in [('as_reported', False), ('alias', True)]:
        so, sn = string_set(oedges, alias), string_set(nedges, alias)
        c0, c1 = classify(pairs, cons, so), classify(pairs, cons, sn)
        cls[name] = (c0, c1)
        out[name] = {'v12_0': counts(c0), 'v12_5': counts(c1)}
        for v in ('v12_0', 'v12_5'):
            d = out[name][v]
            d['in_STRING'] = d['Known Validated'] + d['STRING Only']
            d['coverage_percent'] = round(100.0 * d['in_STRING'] / len(pairs), 1)
        out[name]['pairs_gained'] = sorted(p for p in pairs if p in sn and p not in so)
        out[name]['pairs_lost'] = sorted(p for p in pairs if p in so and p not in sn)
        out[name]['category_changes'] = {p: [c0[p], c1[p]] for p in pairs if c0[p] != c1[p]}
    # all STRING edges among the 17 analyzed genes (alias applied), both versions
    g17 = {'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r'}
    def e17(edges):
        return {pname(ALIAS.get(a, a), ALIAS.get(b, b)) for a, b in edges if ALIAS.get(a, a) in g17 and ALIAS.get(b, b) in g17}
    a0, a1 = e17(oedges), e17(nedges)
    out['edges_among_17_genes'] = {'v12_0': len(a0), 'v12_5': len(a1), 'gained': sorted(a1 - a0), 'lost': sorted(a0 - a1)}
    out['pairs'] = len(pairs)
    for p in sorted(pairs):
        rows.append([p, pairs[p], int(p in cons), '%.3f' % oscore.get(p, float('nan')) if p in oscore else '', '%.3f' % nscore.get(p, float('nan')) if p in nscore else '', cls['as_reported'][0][p], cls['as_reported'][1][p], cls['alias'][0][p], cls['alias'][1][p]])
    with open(os.path.join(RES, 'string_classification.csv'), 'w', encoding='utf8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['pair', 'r_Hes1_normalized', 'four_method_consensus', 'score_v12_0_literal_names', 'score_v12_5_literal_names', 'category_as_reported_v12_0', 'category_as_reported_v12_5', 'category_alias_v12_0', 'category_alias_v12_5'])
        w.writerows(rows)
    with open(os.path.join(RES, 'string_summary.json'), 'w', encoding='utf8') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ('version_response',)}, indent=1)[:3500])


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'gate':
        sys.exit(0 if gate() else 1)
    elif cmd == 'retrieve':
        retrieve()
    elif cmd == 'analyze':
        analyze()
    else:
        print(__doc__)
