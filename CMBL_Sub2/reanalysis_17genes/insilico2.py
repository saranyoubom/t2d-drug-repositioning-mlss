import sys,re,json,time,urllib.request,urllib.parse
sys.stdout.reconfigure(encoding='utf8')
exec(open('insilico.py',encoding='utf8').read().split("d=docx.Document")[0])
def esearch(term,n=8):
    j=json.loads(get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmode=json&retmax=%d&term=%s'%(n,urllib.parse.quote(term))))
    ids=j['esearchresult']['idlist']; time.sleep(0.4)
    if not ids: return []
    s=json.loads(get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id='+','.join(ids)))['result']; time.sleep(0.4)
    return [(s[i]['caption'],s[i]['title'][:90],s[i].get('slen'),s[i].get('statusflag','')) for i in ids]
P={
 'Hes1_table(canine)':('GAGAAGGCGGACATTCTGGA','ACCTCGTTCATACACTCGCTG'),
 'Hes1_r_NEW':('CACCGGACAAACCAAAGACAGC','GAGGTGCTTCACTGTCATTTCC'),
 'Hes1_r_old':('CTGGAGAAGGCAGACATTCTGG','GTCACCTCGTTCATGCACTCG'),
 'Hey1_table(canine)':('ACCTGAAAATGCTGCACACG','GCTGGGAGGCGTAGTTGTTA'),
 'Tcf7_table':('CGTCCTTGATGCTGGGATCTG','GATGATTGGCTTCTTCGCCTCC'),
 'Tcf7_r_NEW':('GAGCACACTTCGCAGAGACTTTTC','GTAGTGCATGAGCAGATTGAAGGC'),
 'Ryr3_table':('CAACCTACTACCCCGAAGGAAC','GAGTGGGCTCTGCTGTTAGG'),
}
G={'Hes1':'Hes1','Hey1':'Hey1','Tcf7':'Tcf7','Ryr3':'Ryr3'}
cands={}
for g in G:
    r=esearch(f'{g}[Gene Name] AND Rattus norvegicus[Organism] AND (biomol_mrna[PROP]) AND srcdb_refseq[PROP]')
    cands[g]=r; print('##',g)
    for c in r: print('   ',c)
res=[]
for name,(F,R) in P.items():
    g=name.split('_')[0]
    for cap,title,slen,flag in cands[g]:
        if not cap.startswith(('NM_','XM_')): continue
        o=check(name,cap,F,R); 
        res.append({'primer':name,'acc':cap,'title':title,'result':o['result'],'amp':o.get('amplicon'),'F_mm':o.get('F_mm'),'R_mm':o.get('R_mm')})
        print(f"{name:20s} vs {cap:16s} {o['result']:28s} amp={o.get('amplicon')} F_mm={o.get('F_mm')} R_mm={o.get('R_mm')}")
json.dump({'cands':cands,'res':res},open('insilico2.json','w',encoding='utf8'),ensure_ascii=False)
