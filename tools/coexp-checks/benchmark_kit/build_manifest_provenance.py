"""Builds benchmark_kit/series_manifest.csv from the frozen lists and the stored results of the study (read only). One row per analysed series (37 development + 59 validation = 96)."""
import os, sys, json, csv
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
A = os.environ.get('STUDY_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'CMBL_Sub2'))   # the folder CMBL_Sub2 of the study repository
V2 = os.path.join(A, 'strengthen_v2'); W2 = os.path.join(V2, 'wp2_benchmark2')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'series_manifest.csv')
rows = []
prep = json.load(open(os.path.join(W2, 'compact', 'before_lean', 'dev_prepared.json'), encoding='utf8'))
frozen = {r['gse']: r for r in csv.DictReader(open(os.path.join(A, 'benchmark_route3', 'frozen_datasets.csv'), encoding='utf8'))}
for gse, p in prep.items():
    f = frozen[gse]
    rows.append(dict(arm='development_arrays', accession=gse, platform_or_organism=p['platform'], frozen_list_order=f['order'], group_labels_in_geo=f['groups'], samples_in_design_order=f['gsms'],
                     gene_set_A_symbols=';'.join(p['genes_A']), gene_set_B_symbols=';'.join(p['genes_B']), n_de_genes=p['n_de'], scale_in_file=p['scale']))
fa = {r['gse']: r for r in csv.DictReader(open(os.path.join(W2, 'frozen_heldout_arrays.csv'), encoding='utf8'))}
for l in open(os.path.join(W2, 'heldout_arrays_results.jsonl'), encoding='utf8'):
    j = json.loads(l)
    if j.get('status') != 'analysed': continue
    f = fa[j['gse']]
    rows.append(dict(arm='heldout_arrays', accession=j['gse'], platform_or_organism=j['platform'], frozen_list_order=f['order'], group_labels_in_geo=f['groups'], samples_in_design_order=f['gsms'],
                     gene_set_A_symbols=';'.join(j['A']['genes']), gene_set_B_symbols=';'.join(j['B']['genes']), n_de_genes=j['n_de'], scale_in_file=j['scale']))
fr = {r['gse']: r for r in csv.DictReader(open(os.path.join(W2, 'frozen_heldout_rnaseq.csv'), encoding='utf8'))}
for l in open(os.path.join(W2, 'heldout_rnaseq_results.jsonl'), encoding='utf8'):
    j = json.loads(l)
    if j.get('status') != 'analysed': continue
    f = fr[j['gse']]
    rows.append(dict(arm='heldout_rnaseq', accession=j['gse'], platform_or_organism=j.get('organism', 'human'), frozen_list_order=f['order'], group_labels_in_geo=f['groups'], samples_in_design_order=f['gsms'],
                     gene_set_A_symbols=';'.join(j['A']['genes']), gene_set_B_symbols=';'.join(j['B']['genes']), n_de_genes=j.get('n_de', ''), scale_in_file=j.get('scale', 'counts->CPM')))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
from collections import Counter
print(len(rows), dict(Counter(r['arm'] for r in rows)))
print('samples per series (should be 12):', Counter(len(r['samples_in_design_order'].split(';')) for r in rows))
print('genes per set (should be 17):', Counter((len(r['gene_set_A_symbols'].split(';')), len(r['gene_set_B_symbols'].split(';'))) for r in rows))
