import sys,json,numpy as np,pandas as pd
sys.path.insert(0,'.')
import pipeline as P
sys.stdout.reconfigure(encoding='utf8')
def run(genes,n=2000,seed=2026):
    df=P.load(genes); L=np.log(df.values); gm=np.vstack([L[P.GROUP==g].mean(0) for g in range(3)]); resid=L-gm[P.GROUP]
    rng=np.random.default_rng(seed); rec={'Raw':[],'Ins1-normalized':[],'Kcnj11-normalized':[],'Hes1-normalized':[],'cons':[],'cons_sign':[]}
    for _ in range(n):
        R=resid.copy()
        for g in range(3):
            idx=np.where(P.GROUP==g)[0]
            for c in range(L.shape[1]): R[idx,c]=rng.permutation(R[idx,c])
        Pm=pd.DataFrame(np.exp(gm[P.GROUP]+R),columns=df.columns); res=P.sig_sets(P.schemes(Pm))
        for k,d in res.items(): rec[k].append(sum(v[1]<0.05 for v in d.values()))
        rec['cons'].append(len(P.consensus(res,use_sign=False))); rec['cons_sign'].append(len(P.consensus(res,use_sign=True)))
    obs=P.sig_sets(P.schemes(df)); out={}
    for k,v in rec.items():
        a=np.array(v); out[k]=dict(mean=float(a.mean()),lo=float(np.percentile(a,2.5)),hi=float(np.percentile(a,97.5)))
    out['observed']={k:int(sum(x[1]<0.05 for x in d.values())) for k,d in obs.items()}
    out['observed']['cons']=int(len(P.consensus(obs,use_sign=False))); out['observed']['cons_sign']=int(len(P.consensus(obs,use_sign=True)))
    out['P_ge_obs']={'cons':float((np.array(rec['cons'])>=out['observed']['cons']).mean()),'cons_sign':float((np.array(rec['cons_sign'])>=out['observed']['cons_sign']).mean())}
    return out
res={'20':run(P.GENES20),'17':run(P.GENES17)}
for k,v in res.items():
    print('==',k,'genes')
    for kk,vv in v.items(): print('  ',kk,vv)

json.dump(res,open('perm_detail.json','w'),indent=1,default=float)
