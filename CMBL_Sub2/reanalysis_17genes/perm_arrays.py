import sys,numpy as np,pandas as pd
sys.path.insert(0,'.')
import pipeline as P
sys.stdout.reconfigure(encoding='utf8')
genes=P.GENES17; n=2000
df=P.load(genes); L=np.log(df.values); gm=np.vstack([L[P.GROUP==g].mean(0) for g in range(3)]); resid=L-gm[P.GROUP]
rng=np.random.default_rng(2026); rec={k:[] for k in ('Raw','Hes1-normalized','Kcnj11-normalized','Ins1-normalized','cons','cons_sign')}
for _ in range(n):
    R=resid.copy()
    for g in range(3):
        idx=np.where(P.GROUP==g)[0]
        for c in range(L.shape[1]): R[idx,c]=rng.permutation(R[idx,c])
    Pm=pd.DataFrame(np.exp(gm[P.GROUP]+R),columns=df.columns); res=P.sig_sets(P.schemes(Pm))
    for k,d in res.items(): rec[k].append(sum(v[1]<0.05 for v in d.values()))
    rec['cons'].append(len(P.consensus(res,use_sign=False))); rec['cons_sign'].append(len(P.consensus(res,use_sign=True)))
np.savez('perm_arrays_17.npz',**{k:np.array(v) for k,v in rec.items()})
print('saved',{k:float(np.mean(v)) for k,v in rec.items()})
