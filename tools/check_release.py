from pathlib import Path
import json
import re

ROOT=Path(__file__).resolve().parents[1]
required=[
    ROOT/'README.md',ROOT/'CITATION.cff',ROOT/'docs/index.html',ROOT/'docs/styles.css',ROOT/'docs/app.js',
    ROOT/'docs/assets/study-data.json',ROOT/'manuscript/JCF_MTI_Brief_Report.md',ROOT/'manuscript/Table1_Multipanel.md'
]
missing=[str(p.relative_to(ROOT)) for p in required if not p.exists()]
if missing: raise SystemExit('Missing required files: '+', '.join(missing))

# Core release files and headline values must remain synchronized.
D=json.loads((ROOT/'docs/assets/study-data.json').read_text())
checks={
    'tdi_n':(D['meta']['tdi_n'],270),
    'tdi_events':(D['meta']['tdi_events'],27),
    'mti_auroc':(D['headline']['mti_auroc'],0.695),
    'combo_auroc':(D['headline']['combo_auroc'],0.712),
    'hr_lvef':(D['headline']['hr_lvef'],2.20),
    'preserved_hr':(D['headline']['preserved_hr'],2.22),
    'one_year_n':(D['added_value']['n'],172),
    'one_year_events':(D['added_value']['events'],16),
}
for name,(got,want) in checks.items():
    if got!=want: raise SystemExit(f'{name}: expected {want}, found {got}')

print('MTI release check passed.')
