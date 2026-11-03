import sys,json
sys.path.insert(0,'.')
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pipeline as P
from statsmodels.stats.multitest import multipletests
def counts(genes):
    df=P.load(genes); sch=P.schemes(df); res=P.sig_sets(sch); adj=P.adjusted(sch); out={}
    pooled_q={}
    for k,d in res.items():
        keys=list(d); q=multipletests([d[x][1] for x in keys],method='fdr_bh')[1]; pooled_q[k]=dict(zip(keys,q))
    for k in res:
        out[k]=dict(pairs=len(res[k]),pooled_p=sum(v[1]<0.05 for v in res[k].values()),pooled_q=sum(q<0.05 for q in pooled_q[k].values()),
                    adj_p=sum(v[1]<0.05 for v in adj[k].values()),adj_q=sum(v[2]<0.05 for v in adj[k].values()))
    common=set.intersection(*[set(d) for d in res.values()])
    allp=[p for p in common if all(res[k][p][1]<0.05 for k in res)]
    allq=[p for p in common if all(pooled_q[k][p]<0.05 for k in res)]
    adjp=[p for p in common if all(adj[k][p][1]<0.05 for k in res)]
    adjq=[p for p in common if all(adj[k][p][2]<0.05 for k in res)]
    sc=[p for p in allp if len({np.sign(res[k][p][0]) for k in res})==1]
    out['all']=dict(pooled_p=len(allp),pooled_q=len(allq),adj_p=len(adjp),adj_q=len(adjq),signcons=len(sc))
    return out
o20=counts(P.GENES20); o17=counts(P.GENES17)
print('20 genes (compare with the existing table: raw 190/80/52/15/0; Hes1 171/169/169/43/12; Kcnj11 171/118/111/54/5; Ins1 171/98/89/76/49; all-four 27/18/5/0; sign-cons 21)'); 
for k,v in o20.items(): print('  ',k,v)
print('17 genes'); 
for k,v in o17.items(): print('  ',k,v)
json.dump({'20':o20,'17':o17},open('table_counts.json','w'),default=int,indent=1)
