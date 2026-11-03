import json,re,sys,math,collections,csv
import numpy as np, docx
sys.stdout.reconfigure(encoding='utf8')
rows=json.load(open('cfx_long.json',encoding='utf8'))
d=docx.Document(r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Supplementary\Additional file 1_CMBL_working.docx')
EFF={}
for r in d.tables[0].rows[2:]:
    c=[x.text.strip() for x in r.cells]
    if c[0]==c[1]: continue
    m=re.findall(r'\d+\.\d+|\d+',c[6])
    if m: EFF[c[0].split()[0]]=float(m[0])/100
EXCL={'Hey1','Wnt2','Wnt2b','Ryr1'}
def grp(s):
    s=str(s).lower()
    if 'dapt' in s: return 'DAPT'
    if 'dkk' in s: return 'DKK'
    if 'ctrl' in s or 'control' in s: return 'CTRL'
    return None
Cq=collections.defaultdict(list)
for r in rows:
    g=grp(r['sample']);
    if g and isinstance(r['cq'],(int,float)): Cq[(r['run'][:22],r['gene'],g)].append(r['cq'])
plates=sorted({k[0] for k in Cq}); genes=sorted({k[1] for k in Cq}-{'Gapdh'}-EXCL)
def fc(plate,gene,treat,E_t,E_g):
    need=[(plate,gene,treat),(plate,gene,'CTRL'),(plate,'Gapdh',treat),(plate,'Gapdh','CTRL')]
    if any(len(Cq.get(k,[]))==0 for k in need): return None
    try:
        dt=np.mean(Cq[(plate,gene,treat)])-np.mean(Cq[(plate,gene,'CTRL')])
        dg=np.mean(Cq[(plate,'Gapdh',treat)])-np.mean(Cq[(plate,'Gapdh','CTRL')])
    except Exception: return None
    return (1+E_t)**(-dt)/((1+E_g)**(-dg))
out=[]
print(f"Gapdh efficiency used: {EFF['Gapdh']:.3f}")
print(f"{'gene':8s} {'E_t':>5s} {'plates':>6s} | DAPT: FC(E=1) FC(E_amp) ratio | DKK-1: FC(E=1) FC(E_amp) ratio | flag")
for g in genes:
    rec={'gene':g,'E':EFF.get(g)}
    for t in ('DAPT','DKK'):
        a=[fc(p,g,t,1.0,1.0) for p in plates]; a=[x for x in a if x]
        b=[fc(p,g,t,EFF[g],EFF['Gapdh']) for p in plates]; b=[x for x in b if x]
        rec[t]=(float(np.exp(np.mean(np.log(a)))),float(np.exp(np.mean(np.log(b)))),len(a),[round(x,3) for x in a],[round(x,3) for x in b])
    flags=[]
    for t in ('DAPT','DKK'):
        A,B=rec[t][0],rec[t][1]
        if (A-1)*(B-1)<0: flags.append(f'{t}: direction flips')
        elif abs(math.log2(B/A))>0.585: flags.append(f'{t}: >1.5x shift')
        elif abs(math.log2(B/A))>0.32: flags.append(f'{t}: >1.25x shift')
    rec['flag']='; '.join(flags); out.append(rec)
    print(f"{g:8s} {EFF.get(g,0):5.2f} {rec['DAPT'][2]:6d} | {rec['DAPT'][0]:7.2f} {rec['DAPT'][1]:8.2f} {rec['DAPT'][1]/rec['DAPT'][0]:6.2f} | {rec['DKK'][0]:7.2f} {rec['DKK'][1]:8.2f} {rec['DKK'][1]/rec['DKK'][0]:6.2f} | {rec['flag']}")
json.dump(out,open('effsens.json','w',encoding='utf8'),ensure_ascii=False)
OUT=r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Analysis_outputs\primer_check'
w=csv.writer(open(OUT+r'\efficiency_sensitivity_group_level.csv','w',newline='',encoding='utf-8-sig'))
w.writerow(['gene','mean_efficiency','plates','DAPT_FC_E100','DAPT_FC_Eamplicon','DAPT_ratio','DKK1_FC_E100','DKK1_FC_Eamplicon','DKK1_ratio','flag'])
for r in out: w.writerow([r['gene'],r['E'],r['DAPT'][2],round(r['DAPT'][0],3),round(r['DAPT'][1],3),round(r['DAPT'][1]/r['DAPT'][0],3),round(r['DKK'][0],3),round(r['DKK'][1],3),round(r['DKK'][1]/r['DKK'][0],3),r['flag']])
