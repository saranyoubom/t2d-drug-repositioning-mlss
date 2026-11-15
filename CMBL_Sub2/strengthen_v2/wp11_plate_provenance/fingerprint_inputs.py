"""fingerprint_inputs.py - WP11 (8 Nov 2026): sha256 of every raw input that the WP11 scripts read, so the results can be tied to the exact files.
Paths are relative to the project folder (CMBL_Sub2_291026). Usage: python -I fingerprint_inputs.py   (writes INPUT_FINGERPRINTS.txt)"""
import os, hashlib, time
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
RT = os.path.abspath(os.path.join(P, '..', '..', 'Undergrad_Glucose uptake signaling', 'RT-PCR'))
SBA = os.path.abspath(os.path.join(P, '..', 'SBA', 'Rev2_250626', 'Export_qPCR'))
PLATES = ['2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1', '2020-Jun-03_Gr2_13genes-Ctrl-DAPT-DKK1', '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2']
files = []
for ver in ('Before', 'After'):
    for pl in PLATES:
        for kind in ('Quantification Cq Results.xlsx', 'Quantification Amplification Results.xlsx'):
            files.append((os.path.join(P, 'Exported_data', ver, pl + ' -  ' + kind), 'user export of 8 Nov 2026, threshold %s' % ('as found' if ver == 'Before' else '150')))
for pl in PLATES:       # the saved .pcrd of the second 24 June plate lost the ' 2' of its name (its exports keep it)
    files.append((os.path.join(P, 'Exported_data', 'After', pl.replace(' 2', '') + '_Threshold150.pcrd'), 'user file saved with threshold 150 (8 Nov 2026)'))
files += [(os.path.join(P, 'Database', 'STRING12.5', 'string_interactions_21_nodes.tsv'), "user's STRING 12.5 web export, 8 Nov 2026"),
          (os.path.join(P, 'Database', 'STRING12.5', 'STRING12.5.txt'), 'network URL of that export')]
for f, d in [('2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx', '26 Aug instrument export (2025 re-export of the .pcrd, wells in plate order)'),
             ('2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx', "26 Aug workbook with the user's calculation columns"),
             ('2020-11-09_ITPR, RYR1-3 -  Quantification Cq Results.xlsx', '9 Nov instrument export')]:
    files.append((os.path.join(RT, f), d))
for f, d in [('2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx', '24 June plate 1, 2025 re-export'), ('2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 -  Quantification Amplification Results.xlsx', '24 June plate 1, 2025 re-export (curves)')]:
    files.append((os.path.join(SBA, f), d))
for f in ('ins1_sheet_cq.csv', 'gene_run_map.csv'):
    files.append((os.path.join(P, 'CMBL_Sub2', 'strengthen_v2', 'wp1_ins1_structure', f), 'WP1 matrix rows'))
out = ['Inputs read by the WP11 scripts; sha256, then description; fingerprinted %s' % time.strftime('%Y-%m-%d %H:%M:%S %z')]
for path, desc in files:
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    out.append('%s *%s' % (h, os.path.relpath(path, P).replace(os.sep, '/')))
    out.append(desc)
open(os.path.join(HERE, 'INPUT_FINGERPRINTS.txt'), 'w', encoding='utf8').write('\n'.join(out) + '\n')
print('%d files fingerprinted' % len(files))
