import sys,re,json,time,urllib.request,urllib.parse,docx
sys.stdout.reconfigure(encoding='utf8')
UA={'User-Agent':'primer-check/1.0 (mailto:YOUR_EMAIL_HERE)'}
def get(u):
    for _ in range(4):
        try: return urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=40).read().decode('utf8')
        except Exception as e: err=str(e); time.sleep(1.5)
    return ''
def rc(s): return s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1].upper()
def find(seq,p,maxmm=2):
    p=p.upper(); L=len(p); best=(99,-1)
    i=seq.find(p)
    if i>=0: return (0,i)
    for i in range(len(seq)-L+1):
        mm=0
        for a,b in zip(seq[i:i+L],p):
            if a!=b:
                mm+=1
                if mm>=best[0] or mm>maxmm: break
        else:
            if mm<best[0]: best=(mm,i)
    return best if best[0]<=maxmm else (None,-1)
cache={}
def fetch(acc):
    if acc in cache: return cache[acc]
    fa=get(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={acc}')
    seq=''.join(fa.split('\n')[1:]).upper(); hdr=fa.split('\n')[0]
    cache[acc]=(hdr,seq); time.sleep(0.4); return cache[acc]
def check(gene,acc,F,R,tablelen=None):
    hdr,seq=fetch(acc)
    mf,pf=find(seq,F); mr,pr=find(seq,rc(R))
    if mf is None or mr is None:
        # try the other orientation
        mf2,pf2=find(seq,rc(F)); mr2,pr2=find(seq,R)
        if mf2 is not None and mr2 is not None: return dict(gene=gene,acc=acc,hdr=hdr[:110],result='REVERSED ORIENTATION',F_mm=mf2,R_mm=mr2)
        return dict(gene=gene,acc=acc,hdr=hdr[:110],result='NO MATCH',F_mm=mf,R_mm=mr,seqlen=len(seq))
    amp=pr+len(R)-pf
    return dict(gene=gene,acc=acc,hdr=hdr[:110],result='MATCH' if mf==0 and mr==0 else f'MATCH with {mf}+{mr} mismatches',F_mm=mf,R_mm=mr,amplicon=amp,tablelen=tablelen,len_ok=(tablelen is None or amp==tablelen),F_pos=pf+1,R_pos=pr+len(R))
d=docx.Document(r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Supplementary\Additional file 1.docx')
tb=d.tables[0]; rows=[]
for r in tb.rows[2:]:
    c=[x.text.strip() for x in r.cells]
    if c[0]==c[1]: continue
    seqs=c[3].split();
    if len(seqs)<2: continue
    rows.append((c[0].split(' ')[0],c[1],seqs[0],seqs[1],int(re.findall(r'\d+',c[4])[0])))
out=[]
for g,acc,F,R,L in rows:
    out.append(check(g,acc,F,R,L));
# candidates
NEW={'Hes1_NEW':(F:='CACCGGACAAACCAAAGACAGC','GAGGTGCTTCACTGTCATTTCC'),'Tcf7_NEW':('GAGCACACTTCGCAGAGACTTTTC','GTAGTGCATGAGCAGATTGAAGGC'),
     'Hes1_r_old':('CTGGAGAAGGCAGACATTCTGG','GTCACCTCGTTCATGCACTCG'),'Tcf7_r_Fw+table_Rv':('CGTCCTTGATGCTGGGATCTG','GATGATTGGCTTCTTCGCCTCC')}
json.dump(out,open('insilico_table.json','w',encoding='utf8'),ensure_ascii=False)
for o in out:
    print(f"{o['gene']:8s} {o['acc']:16s} {o['result']:30s} amp={o.get('amplicon')} table={o.get('tablelen')} | {o['hdr'][:85]}")
