import math,re,sys,docx
sys.stdout.reconfigure(encoding='utf8')
NN={'AA':(-7.9,-22.2),'TT':(-7.9,-22.2),'AT':(-7.2,-20.4),'TA':(-7.2,-21.3),'CA':(-8.5,-22.7),'TG':(-8.5,-22.7),'GT':(-8.4,-22.4),'AC':(-8.4,-22.4),
    'CT':(-7.8,-21.0),'AG':(-7.8,-21.0),'GA':(-8.2,-22.2),'TC':(-8.2,-22.2),'CG':(-10.6,-27.2),'GC':(-9.8,-24.4),'GG':(-8.0,-19.9),'CC':(-8.0,-19.9)}
def tm(seq,mono=50e-3,div=1.5e-3,dntp=0.6e-3,conc=50e-9):
    s=seq.upper(); dH=0.0; dS=0.0
    for a,b in zip(s,s[1:]): h,e=NN[a+b]; dH+=h; dS+=e
    for end in (s[0],s[-1]):
        if end in 'GC': dH+=0.1; dS+=-2.8
        else: dH+=2.3; dS+=4.1
    na=(mono*1000+120*math.sqrt(max(0.0,div*1000-dntp*1000)))/1000
    dS+=0.368*(len(s)-1)*math.log(na)
    return dH*1000/(dS+1.987*math.log(conc/4))-273.15
if __name__=='__main__':
    d=docx.Document(r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Supplementary\Additional file 1.docx')
    errs=[]
    for r in d.tables[0].rows[2:]:
        c=[x.text.strip() for x in r.cells]
        if c[0]==c[1]: continue
        s=c[3].split(); tms=[float(x) for x in re.findall(r'\d+\.\d+',c[5])]
        if len(s)<2 or len(tms)<2: continue
        for seq,t in zip(s[:2],tms[:2]):
            errs.append((c[0].split()[0],round(tm(seq),2),t,round(tm(seq)-t,2)))
    import statistics as st
    e=[x[3] for x in errs]; print('n',len(e),'mean err',round(st.mean(e),2),'SD',round(st.pstdev(e),2),'max abs',max(abs(x) for x in e))
    for x in errs[:8]: print(x)
    for name,seq in (('r_Hes1_Fw_NEW','CACCGGACAAACCAAAGACAGC'),('r_Hes1_Rv_NEW','GAGGTGCTTCACTGTCATTTCC'),('r_Tcf7_Fw_NEW','GAGCACACTTCGCAGAGACTTTTC'),('r_Tcf7_Rv_NEW','GTAGTGCATGAGCAGATTGAAGGC')):
        print(name,round(tm(seq),2),len(seq))
