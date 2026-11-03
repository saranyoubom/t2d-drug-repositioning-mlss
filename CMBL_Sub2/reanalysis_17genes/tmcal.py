import sys,re,json,math,statistics as st
sys.stdout.reconfigure(encoding='utf8')
exec(open('insilico.py',encoding='utf8').read().split("d=docx.Document")[0])
import numpy as np
summ=json.load(open('cfx_summary.json',encoding='utf8'))
tab=json.load(open('insilico_table.json',encoding='utf8'))
import docx
d=docx.Document(r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Supplementary\Additional file 1.docx')
prim={}
for r in d.tables[0].rows[2:]:
    c=[x.text.strip() for x in r.cells]
    if c[0]==c[1]: continue
    s=c[3].split()
    if len(s)>=2: prim[c[0].split()[0]]=(c[1],s[0],s[1])
def amp(acc,F,R):
    hdr,seq=fetch(acc); mf,pf=find(seq,F); mr,pr=find(seq,rc(R))
    if mf is None or mr is None: return None
    return seq[pf:pr+len(R)]
def pred(a):
    gc=100*sum(1 for ch in a if ch in 'GC')/len(a); return 81.5+0.41*gc-675/len(a)-21.6+0  # 0.05 M Na+
pts=[]
for g,(acc,F,R) in prim.items():
    key={'Glut2':'Glut2','Ins1':'Ins1'}.get(g,g)
    if g in('Hes1','Hey1','Ryr3'): continue
    a=amp(acc,F,R)
    if not a or key not in summ or summ[key]['tm'] is None: continue
    if key=='Ryr1': continue
    pts.append((g,len(a),round(pred(a),2),summ[key]['tm']))
X=np.array([p[2] for p in pts]); Y=np.array([p[3] for p in pts])
b,a0=np.polyfit(X,Y,1); res=Y-(a0+b*X); sd=res.std(ddof=2)
print(f'calibration over {len(pts)} genes: Tm_meas = {a0:.2f} + {b:.2f} * Tm_pred ; residual SD {sd:.2f} C ; r = {np.corrcoef(X,Y)[0,1]:.2f}')
for p,r in zip(pts,res): print(f'   {p[0]:8s} N={p[1]:3d} pred={p[2]:6.2f} meas={p[3]:6.2f} resid={r:+.2f}')
cand={
 'Hes1 printed(canine pair) on rat NM_024360':('NM_024360','GAGAAGGCGGACATTCTGGA','ACCTCGTTCATACACTCGCTG',summ['Hes1']['tm']),
 'Hes1 r_NEW on rat NM_024360':('NM_024360','CACCGGACAAACCAAAGACAGC','GAGGTGCTTCACTGTCATTTCC',summ['Hes1']['tm']),
 'Hes1 r_old (r_Hes1_Fw/Rv)':('NM_024360','CTGGAGAAGGCAGACATTCTGG','GTCACCTCGTTCATGCACTCG',summ['Hes1']['tm']),
 'Tcf7 printed':('NM_001427430','CGTCCTTGATGCTGGGATCTG','GATGATTGGCTTCTTCGCCTCC',summ['Tcf7']['tm']),
 'Tcf7 r_NEW':('NM_001427430','GAGCACACTTCGCAGAGACTTTTC','GTAGTGCATGAGCAGATTGAAGGC',summ['Tcf7']['tm']),
 'Ryr3 printed on XM_087613279':('XM_087613279','CAACCTACTACCCCGAAGGAAC','GAGTGGGCTCTGCTGTTAGG',summ['Ryr3']['tm']),
 'Hey1 printed(canine) on rat NM_001191845':('NM_001191845','ACCTGAAAATGCTGCACACG','GCTGGGAGGCGTAGTTGTTA',summ['Hey1']['tm']),
}
print()
for k,(acc,F,R,meas) in cand.items():
    hdr,seq=fetch(acc)
    # allow mismatches for the canine pair
    mf,pf=find(seq,F); mr,pr=find(seq,rc(R))
    a=seq[pf:pr+len(R)]; p=a0+b*pred(a)
    print(f'{k:44s} N={len(a):3d} predicted Tm={p:5.1f}  measured={meas:5.1f}  diff={meas-p:+.1f}  ({(meas-p)/sd:+.1f} SD)')
