import sys,json
sys.path.insert(0,'.')
sys.stdout.reconfigure(encoding='utf8')
import pipeline as P
def run(genes,alts,orig_den='Hes1'):
    df=P.load(genes); base=P.sig_sets(P.schemes(df)); cons=P.consensus(base,use_sign=False)
    out={'n_consensus':len(cons)}
    for D in alts:
        d3=df.drop(columns=[D]) if False else df
        sch={'Raw':df}
        for dd in ('Ins1','Kcnj11'): sch[dd+'-normalized']=df.drop(columns=[dd]).div(df[dd],axis=0)
        sch[D+'-normalized']=df.drop(columns=[D]).div(df[D],axis=0)
        res=P.sig_sets(sch); newc=set(P.consensus(res,use_sign=False))
        evaluable=[p for p in cons if D not in p]
        kept=[p for p in evaluable if p in newc]
        out[D]=dict(evaluable=len(evaluable),retained=len(kept),new_consensus=len(newc),gained=len(newc-set(cons)))
    return out
r17=run(P.GENES17,['Tcf7','Cacna1c']); r20=run(P.GENES20,['Hey1','Tcf7','Cacna1c'])
print('20 genes (check vs manuscript: Hey1/Tcf7 retained all evaluable, Cacna1c 12 of 22):',json.dumps(r20))
print('17 genes:',json.dumps(r17))
json.dump({'20':r20,'17':r17},open('denom_sens.json','w'),indent=1)
