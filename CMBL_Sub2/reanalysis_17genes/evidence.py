import sys,re,json,csv,os
sys.path.insert(0,'.')
sys.stdout.reconfigure(encoding='utf8')
exec(open('insilico.py',encoding='utf8').read().split("d=docx.Document")[0])
import numpy as np, docx
from tmnn import tm
summ=json.load(open('cfx_summary.json',encoding='utf8'))
d=docx.Document(r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Supplementary\Additional file 1.docx')
rows=[]
for r in d.tables[0].rows[2:]:
    c=[x.text.strip() for x in r.cells]
    if c[0]==c[1]: continue
    s=c[3].split()
    if len(s)>=2: rows.append(dict(gene=c[0].split()[0],acc=c[1],F=s[0],R=s[1],L=int(re.findall(r'\d+',c[4])[0]),eff=c[6]))
# corrections (reasons in the evidence file)
CORR={'Hes1':dict(acc='NM_024360.4',F='CACCGGACAAACCAAAGACAGC',R='GAGGTGCTTCACTGTCATTTCC',note='r_Hes1_Fw_NEW / r_Hes1_Rv_NEW (13 Feb 2020 primer sheet); printed pair was the canine K9_Hes1 design'),
      'Tcf7':dict(acc='NM_001427430.1',F='GAGCACACTTCGCAGAGACTTTTC',R='GTAGTGCATGAGCAGATTGAAGGC',note='r_Tcf7_Fw_NEW / r_Tcf7_Rv_NEW (13 Feb 2020 primer sheet)'),
      'Ryr3':dict(acc='XM_087613279.1',note='accession updated; primers unchanged (previous accession XM_039106806.1 no longer matches)')}
EXCL={'Hey1':'identity unresolved (printed pair is a canine design; melt identical to Hes1)','Ryr1':'non-specific melt (75.1 \u00b0C); 12 of 12 wells flagged by rule','Wnt2':'17% of wells at Cq \u2265 35 or without melt peak','Wnt2b':'46% of wells at Cq \u2265 35 or without melt peak'}
# calibration from well-matched genes
pts=[]
def ampseq(acc,F,R):
    h,seq=fetch(acc); mf,pf=find(seq,F); mr,pr=find(seq,rc(R))
    if mf is None or mr is None: return None,None,None,h
    return seq[pf:pr+len(R)],mf,mr,h
def pred(a):
    gc=100*sum(1 for ch in a if ch in 'GC')/len(a); return 81.5+0.41*gc-675/len(a)-21.6
base={}
for r in rows:
    g=r['gene']; c=CORR.get(g,{})
    acc=c.get('acc',r['acc']); F=c.get('F',r['F']); R=c.get('R',r['R'])
    a,mf,mr,h=ampseq(acc,F,R)
    if g=='Hey1': a,mf,mr,h=ampseq('NM_001191845.2',F,R)
    base[g]=dict(gene=g,acc=acc,F=F,R=R,amp=len(a) if a else None,Fmm=mf,Rmm=mr,hdr=h[1:90],pred=pred(a) if a else None,tmF=round(tm(F),2),tmR=round(tm(R),2),eff=r['eff'],old_acc=r['acc'],old_L=r['L'])
key={'Glut2':'Glut2','Ins1':'Ins1'}
cal=[(v['pred'],summ[g]['tm']) for g,v in base.items() if g in summ and v['pred'] and g not in('Hes1','Hey1','Ryr3','Ryr1','Tcf7') ]
X=np.array([c[0] for c in cal]); Y=np.array([c[1] for c in cal]); b,a0=np.polyfit(X,Y,1); sd=(Y-(a0+b*X)).std(ddof=2)
print('calibration genes',len(cal),'SD',round(sd,2))
out=[]
for g,v in base.items():
    s=summ.get(g)
    p=a0+b*v['pred'] if v['pred'] else None
    out.append(dict(gene=g,reference=v['acc'],species=('Rattus norvegicus' if 'Rattus' in v['hdr'] else 'Canis lupus' if 'Canis' in v['hdr'] else '?'),
      forward=v['F'],reverse=v['R'],Tm_F=v['tmF'],Tm_R=v['tmR'],F_mismatch=v['Fmm'],R_mismatch=v['Rmm'],amplicon_bp=v['amp'],table_bp_before=v['old_L'],
      melt_mean=(round(s['tm'],2) if s and s['tm'] else None),melt_sd=(round(s['sd'],2) if s and s['sd'] is not None else None),wells=(s['n'] if s else None),
      wells_late_or_no_melt=(max(s['late'],s['nomelt']) if s else None),predicted_melt=(round(p,1) if p else None),
      melt_diff=(round(s['tm']-p,1) if s and p else None),melt_diff_SD=(round((s['tm']-p)/sd,1) if s and p else None),
      correction=CORR.get(g,{}).get('note',''),analysis=('excluded: '+EXCL[g]) if g in EXCL else 'included'))
os.makedirs(r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Analysis_outputs\primer_check',exist_ok=True)
OUT=r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Analysis_outputs\primer_check'
w=csv.DictWriter(open(OUT+r'\primer_evidence.csv','w',newline='',encoding='utf-8-sig'),fieldnames=list(out[0]))
w.writeheader(); [w.writerow(o) for o in out]
json.dump(dict(rows=out,cal=dict(a=a0,b=b,sd=sd,n=len(cal))),open('evidence.json','w',encoding='utf8'),ensure_ascii=False)
for o in out: print(f"{o['gene']:8s} {o['reference']:15s} {o['species'][:6]:6s} mm={o['F_mismatch']}/{o['R_mismatch']} amp={o['amplicon_bp']} (was {o['table_bp_before']}) melt={o['melt_mean']} pred={o['predicted_melt']} diff={o['melt_diff']} ({o['melt_diff_SD']} SD) late/nomelt={o['wells_late_or_no_melt']}/{o['wells']} {o['analysis'][:30]}")
