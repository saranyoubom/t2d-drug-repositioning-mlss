"""Pooled co-expression pipeline, re-implemented in Python from the manuscript Methods and the Fig_4c R script.
Input: Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx (12 observations x 21 genes, 2^-dCt, rows 0-3 CTRL, 4-7 DAPT, 8-11 DKK-1)."""
import numpy as np, pandas as pd, itertools, json
from scipy import stats
from statsmodels.stats.multitest import multipletests
SRC=r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\SBA\Figures\Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'
GROUP=np.array([0]*4+[1]*4+[2]*4); DEN=['Ins1','Kcnj11','Hes1']
GENES20=['Wnt2','Wnt2b','Wnt5a','Wnt5b','Wnt9a','Lef1','Tcf7','Tcf7l2','Hes1','Hey1','Itpr1','Ryr2','Ryr3','Cacna1c','Cacna1d','Kcnj11','Ptbp1','Glut2','Ins1','Glp1r']
EXCLUDED=['Hey1','Wnt2','Wnt2b','Ryr1']
GENES17=[g for g in GENES20 if g not in EXCLUDED]
def load(genes):
    df=pd.read_excel(SRC); return df[genes].astype(float).reset_index(drop=True)
def schemes(df):
    out={'Raw':df}
    for d in ['Ins1','Kcnj11','Hes1']:
        if d in df: out[d+'-normalized']=df.drop(columns=[d]).div(df[d],axis=0)
    return out
def corr_p(M,df_resid=None):
    X=np.asarray(M,float); n=X.shape[0]; r=np.corrcoef(X,rowvar=False)
    dfree=(n-2) if df_resid is None else df_resid
    with np.errstate(divide='ignore',invalid='ignore'):
        t=r*np.sqrt(dfree/(1-r**2)); p=2*stats.t.sf(np.abs(t),dfree)
    np.fill_diagonal(p,np.nan); return r,p
def pairs(genes): return list(itertools.combinations(range(len(genes)),2))
def sig_sets(sch,alpha=0.05):
    """per scheme: dict pair(gene names sorted)->(r,p)"""
    res={}
    for k,m in sch.items():
        r,p=corr_p(m.values); g=list(m.columns); d={}
        for i,j in pairs(g): d[tuple(sorted((g[i],g[j])))]=(r[i,j],p[i,j])
        res[k]=d
    return res
def consensus(res,alpha=0.05,use_sign=True):
    keys=list(res); common=set.intersection(*[set(res[k]) for k in keys]); out=[]
    for pr in common:
        ok=all(res[k][pr][1]<=alpha for k in keys); sgn={np.sign(res[k][pr][0]) for k in keys}
        if ok and (len(sgn)==1 or not use_sign): out.append(pr)
    return sorted(out)
def counts(res):
    c={}
    for k,d in res.items():
        sig=[v for v in d.values() if v[1]<0.05]; c[k]=dict(pairs=len(d),sig=len(sig),strong=sum(abs(v[0])>=0.7 for v in sig),very=sum(abs(v[0])>=0.8 for v in sig))
    return c
def adjusted(sch):
    X=np.column_stack([np.ones(12),(GROUP==1).astype(float),(GROUP==2).astype(float)]); out={}
    for k,m in sch.items():
        R=m.values-X@np.linalg.lstsq(X,m.values,rcond=None)[0]
        r,p=corr_p(R,df_resid=12-2-2); g=list(m.columns); items=[(tuple(sorted((g[i],g[j]))),r[i,j],p[i,j]) for i,j in pairs(g)]
        q=multipletests([it[2] for it in items],method='fdr_bh')[1]
        out[k]={it[0]:(it[1],it[2],qq) for it,qq in zip(items,q)}
    return out
def permutation_null(df,n=2000,seed=2026,observed=None):
    rng=np.random.default_rng(seed); L=np.log(df.values); gm=np.vstack([L[GROUP==g].mean(0) for g in range(3)])
    resid=L-gm[GROUP]; counts_=[]; counts_s=[]
    for _ in range(n):
        R=resid.copy()
        for g in range(3):
            idx=np.where(GROUP==g)[0]
            for c in range(L.shape[1]): R[idx,c]=rng.permutation(R[idx,c])
        P=pd.DataFrame(np.exp(gm[GROUP]+R),columns=df.columns)
        rs=sig_sets(schemes(P)); counts_.append(len(consensus(rs,use_sign=False))); counts_s.append(len(consensus(rs,use_sign=True)))
    return np.array(counts_),np.array(counts_s)
def rcrit(df_,alpha):
    t=stats.t.ppf(1-alpha/2,df_); return t/np.sqrt(df_+t**2)
def analyze(genes,nperm=2000,seed=2026,label=''):
    df=load(genes); sch=schemes(df); res=sig_sets(sch); cons=consensus(res,use_sign=False); cons_s=consensus(res,use_sign=True); cn=counts(res)
    adj=adjusted(sch)
    out=dict(label=label,genes=genes,n_genes=len(genes),counts=cn,consensus=cons,consensus_signcons=cons_s)
    out['adj_sig_q05']={k:sum(v[2]<0.05 for v in d.values()) for k,d in adj.items()}
    out['adj_sig_p05']={k:sum(v[1]<0.05 for v in d.values()) for k,d in adj.items()}
    out['adj_all_four']=[pr for pr in set.intersection(*[set(d) for d in adj.values()]) if all(adj[k][pr][2]<0.05 for k in adj)]
    out['cons_adj_any']=[pr for pr in cons if any(adj[k][pr][2]<0.05 for k in adj if pr in adj[k])]
    out['cons_adj_counts']={'%s|%s'%pr:sum(adj[k][pr][2]<0.05 for k in adj if pr in adj[k]) for pr in cons}
    m_norm=max(c['pairs'] for k,c in cn.items() if k!='Raw'); m_raw=cn['Raw']['pairs']
    out['rcrit_df8_q05_best']=dict(raw=rcrit(8,0.05/m_raw),norm=rcrit(8,0.05/m_norm)); out['rcrit_df8_p05']=rcrit(8,0.05)
    c,cs=permutation_null(df,nperm,seed)
    out['perm']=dict(n=nperm,mean=float(c.mean()),sd=float(c.std()),p_ge_obs=float((c>=len(cons)).mean()),observed=len(cons),q05=float(np.percentile(c,5)),q95=float(np.percentile(c,95)))
    out['perm_signcons']=dict(n=nperm,mean=float(cs.mean()),sd=float(cs.std()),p_ge_obs=float((cs>=len(cons_s)).mean()),observed=len(cons_s),q05=float(np.percentile(cs,5)),q95=float(np.percentile(cs,95)))
    # sign consistency of consensus pairs across schemes
    out['sign_inconsistent']=[pr for pr in set.intersection(*[set(d) for d in res.values()]) if all(res[k][pr][1]<=0.05 for k in res) and len({np.sign(res[k][pr][0]) for k in res})>1]
    out['cons_table']={ '%s|%s'%pr:{k:dict(r=float(res[k][pr][0]),p=float(res[k][pr][1])) for k in res} for pr in cons}
    out['adj_table']={ '%s|%s'%pr:{k:dict(r=float(adj[k][pr][0]),q=float(adj[k][pr][2])) for k in adj if pr in adj[k]} for pr in cons}
    return out,res,adj
if __name__=='__main__':
    import sys; sys.stdout.reconfigure(encoding='utf8')
    for lab,genes in (('20 genes (current manuscript)',GENES20),('17 genes (after exclusion)',GENES17)):
        o,res,adj=analyze(genes,label=lab)
        print('=====',lab,'| genes',o['n_genes'])
        for k,c in o['counts'].items(): print('  ',k,c)
        print('   pooled consensus (p<=0.05 in all 4 schemes, any sign):',len(o['consensus']),'| sign-consistent:',len(o['consensus_signcons']))
        for nm in ('perm','perm_signcons'): print('   permutation null (%s): mean %.1f SD %.1f, P(null>=%d)=%.3f, 5-95%%: %.0f-%.0f'%(nm,o[nm]['mean'],o[nm]['sd'],o[nm]['observed'],o[nm]['p_ge_obs'],o[nm]['q05'],o[nm]['q95']))
        print('   treatment-adjusted q<0.05 per scheme:',o['adj_sig_q05'],'| all four:',len(o['adj_all_four']),'| consensus pairs sig in >=1 scheme:',len(o['cons_adj_any']),'of',len(o['consensus']))
        print('   consensus pairs sig in >=3 schemes:',[k for k,v in o['cons_adj_counts'].items() if v>=3])
        print('   sign-inconsistent pairs (sig all 4 but sign differs):',len(o['sign_inconsistent']))
        print('   |r| needed at df=8: p<0.05 %.3f ; best pair q<0.05 raw %.3f / norm %.3f'%(o['rcrit_df8_p05'],o['rcrit_df8_q05_best']['raw'],o['rcrit_df8_q05_best']['norm']))
        json.dump(o,open('pipeline_%dgenes.json'%o['n_genes'],'w'),default=lambda x:list(x) if isinstance(x,tuple) else str(x))
