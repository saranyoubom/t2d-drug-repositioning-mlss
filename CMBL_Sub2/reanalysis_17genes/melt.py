import sys,re,json,statistics as st,math
sys.path.insert(0,'.')
sys.stdout.reconfigure(encoding='utf8')
import cfx
rows=[]   # run, well, gene, sample, cq, melt list
for run in cfx.runs():
    cq={}; mt={}
    for t,r in cfx.sheet_rows(run,'Quantification Cq Results').items():
        if t.startswith('Run'): continue
        h=r[0]; ix={k:i for i,k in enumerate(h) if k}
        for x in r[1:]:
            if x[ix['Well']]: cq[x[ix['Well']]]=(x[ix['Target']],x[ix['Sample']],x[ix['Cq']],x[ix['Content']])
    for t,r in cfx.sheet_rows(run,'Melt Curve Peak Results').items():
        if t.startswith('Run'): continue
        h=r[0]; ix={k:i for i,k in enumerate(h) if k}
        for x in r[1:]:
            if x[ix['Well']]: mt.setdefault(x[ix['Well']],[]).append((x[ix['Melt Temperature']],x[ix['Peak Height']])) if isinstance(x[ix['Melt Temperature']],(int,float)) else mt.setdefault(x[ix['Well']],[])
    for w,(tg,sm,c,ct) in cq.items():
        g=re.sub(r'^\d+\.\s*','',str(tg)).strip(); suf=re.search(r'_(Ctrl|DAPT|DKK1)$',g); g=re.sub(r'_(Ctrl|DAPT|DKK1)$','',g); sm=(sm if sm else (suf.group(1) if suf else sm))
        g={'GAPDH':'Gapdh','GLUT-2':'Glut2','Ins-1':'Ins1','ITPR1':'Itpr1','RYR1':'Ryr1','RYR2':'Ryr2','RYR3':'Ryr3'}.get(g,g)
        rows.append(dict(run=run,well=w,gene=g,sample=sm,cq=c,content=ct,melt=mt.get(w,[])))
json.dump(rows,open('cfx_long.json','w',encoding='utf8'),ensure_ascii=False)
print(len(rows),'wells;',len({r['gene'] for r in rows}),'targets:',sorted({r['gene'] for r in rows}))
from collections import defaultdict
G=defaultdict(list)
for r in rows:
    if r['content'].startswith('Unkn') or True: G[r['gene']].append(r)
print(f"{'gene':10s} {'wells':>5s} {'Cq>=35/NA':>9s} {'nomelt':>6s} {'multi':>5s} {'Tm mean':>8s} {'Tm SD':>6s}  runs")
summ={}
for g,rs in sorted(G.items()):
    n=len(rs); late=sum(1 for r in rs if not isinstance(r['cq'],(int,float)) or r['cq']>=35); nomelt=sum(1 for r in rs if not r['melt'])
    multi=sum(1 for r in rs if len(r['melt'])>1)
    tms=[r['melt'][0][0] for r in rs if r['melt']]
    summ[g]=dict(n=n,late=late,nomelt=nomelt,multi=multi,tm=(st.mean(tms) if tms else None),sd=(st.pstdev(tms) if len(tms)>1 else None),runs=sorted({r['run'][:16] for r in rs}))
    print(f"{g:10s} {n:5d} {late:9d} {nomelt:6d} {multi:5d} {summ[g]['tm'] or float('nan'):8.2f} {summ[g]['sd'] or float('nan'):6.2f}  {','.join(summ[g]['runs'])}")
json.dump(summ,open('cfx_summary.json','w',encoding='utf8'),ensure_ascii=False)
