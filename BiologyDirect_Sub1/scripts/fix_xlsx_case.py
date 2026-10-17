import zipfile, os, sys

CANON = {
    '[content_types].xml': '[Content_Types].xml',
    '_rels/.rels': '_rels/.rels',
    'xl/workbook.xml': 'xl/workbook.xml',
    'xl/_rels/workbook.xml.rels': 'xl/_rels/workbook.xml.rels',
    'xl/styles.xml': 'xl/styles.xml',
    'xl/sharedstrings.xml': 'xl/sharedStrings.xml',
    'xl/theme/theme1.xml': 'xl/theme/theme1.xml',
    'docprops/core.xml': 'docProps/core.xml',
    'docprops/app.xml': 'docProps/app.xml',
}

def fix_name(name):
    low = name.lower()
    if low in CANON:
        return CANON[low]
    # worksheet parts and their rels: preserve numeric suffix, fix directory case
    if low.startswith('xl/worksheets/_rels/') and low.endswith('.xml.rels'):
        return 'xl/worksheets/_rels/' + name.split('/')[-1]
    if low.startswith('xl/worksheets/') and low.endswith('.xml'):
        return 'xl/worksheets/' + name.split('/')[-1]
    return name

def fix_file(src, dst):
    zin = zipfile.ZipFile(src, 'r')
    zout = zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED)
    for item in zin.infolist():
        data = zin.read(item.filename)
        zout.writestr(fix_name(item.filename), data)
    zout.close()
    zin.close()

if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    fix_file(src, dst)
    print('fixed ->', dst)
