import zipfile,io,openpyxl,glob,re,os
EXP='C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/SBA/Rev2_250626/Export_qPCR'
CANON={'[content_types].xml':'[Content_Types].xml','xl/sharedstrings.xml':'xl/sharedStrings.xml','docprops/core.xml':'docProps/core.xml','docprops/app.xml':'docProps/app.xml'}
def fix_name(n):
    low=n.lower()
    if low in CANON: return CANON[low]
    if low.startswith('xl/worksheets/_rels/') and low.endswith('.xml.rels'): return 'xl/worksheets/_rels/'+n.split('/')[-1]
    if low.startswith('xl/worksheets/') and low.endswith('.xml'): return 'xl/worksheets/'+n.split('/')[-1]
    return n
def load(path):
    z=zipfile.ZipFile(path); buf=io.BytesIO(); zo=zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED)
    for n in z.namelist(): zo.writestr(fix_name(n),z.read(n))
    zo.close(); buf.seek(0)
    return openpyxl.load_workbook(buf,data_only=True)
def runs():
    return sorted({re.sub(r' -  .*','',os.path.basename(p)) for p in glob.glob(EXP+'/*.xlsx')})
def sheet_rows(run,kind):
    wb=load(EXP+'/'+run+' -  '+kind+'.xlsx'); return {ws.title:[r for r in ws.iter_rows(values_only=True)] for ws in wb}
