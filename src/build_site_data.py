from pathlib import Path
import json, sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT.parent / 'bc_cardiotox_raw' / 'BC_cardiotox_clinical_and_functional_variables.csv'
if not RAW.exists():
    RAW = Path('/mnt/data/bc_cardiotox_raw/BC_cardiotox_clinical_and_functional_variables.csv')
PRED = ROOT / 'results' / 'crossfitted_predictions.csv'
SITE = ROOT / 'site' / 'assets'
SITE.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT/'src'))
from run_analysis import extract_tdi_features, make_mti, COMPONENTS, COMPONENT_LABELS

df = pd.read_csv(RAW, sep=';', decimal=',')
tdi_cols = [c for c in df.columns if c.startswith('t ')]
tdi = df.loc[df[tdi_cols].notna().all(axis=1)].copy().reset_index(drop=True)

# Feature extraction and full-cohort descriptive MTI. This score is for visualization only;
# inferential/predictive results use training-fold standardization in run_analysis.py.
features = pd.DataFrame([extract_tdi_features(row) for row in tdi[tdi_cols].to_numpy(float)])
analysis = pd.concat([tdi.drop(columns=tdi_cols).reset_index(drop=True), features], axis=1)
analysis['MTI_visual'] = make_mti(analysis, reference=analysis).to_numpy()

# Downsample aggregated waveforms to keep the website light.
phase_full = np.linspace(0, 1, len(tdi_cols))
idx = np.unique(np.linspace(0, len(tdi_cols)-1, 201).round().astype(int))
phase = phase_full[idx]
X = tdi[tdi_cols].to_numpy(float)[:, idx]

def summary_wave(mask):
    xx=X[np.asarray(mask)]
    return {
        'median': np.median(xx,axis=0).round(4).tolist(),
        'q1': np.quantile(xx,.25,axis=0).round(4).tolist(),
        'q3': np.quantile(xx,.75,axis=0).round(4).tolist(),
    }

# MTI decile waveform atlas using descriptive full-cohort score.
ranks = pd.qcut(analysis['MTI_visual'].rank(method='first'), 10, labels=False)
deciles=[]
for d in range(10):
    mask=(ranks==d).to_numpy()
    deciles.append({
        'decile': d+1,
        'n': int(mask.sum()),
        'median_mti': round(float(analysis.loc[mask,'MTI_visual'].median()),3),
        'event_rate': round(float(analysis.loc[mask,'CTRCD'].mean()),4),
        'waveform': np.median(X[mask],axis=0).round(4).tolist(),
    })

pred = pd.read_csv(PRED)
y=pred['CTRCD'].to_numpy(int)
risk=pred['MTI+LVEF'].to_numpy(float)
order=np.argsort(-risk)
rows=[]
for frac in np.arange(.05,1.001,.05):
    k=max(1,int(np.ceil(frac*len(pred))))
    ids=order[:k]
    ev=int(y[ids].sum())
    rows.append({
        'fraction': round(float(frac),2),
        'selected_n': int(k),
        'events': ev,
        'event_capture': round(ev/max(int(y.sum()),1),4),
        'ppv': round(ev/k,4),
    })

# Baseline summaries for website cards.
def mediqr(s):
    s=pd.Series(s).dropna().astype(float)
    return {'median':round(float(s.median()),2),'q1':round(float(s.quantile(.25)),2),'q3':round(float(s.quantile(.75)),2)}

components=[]
for c in COMPONENTS:
    a=analysis.loc[analysis.CTRCD==1,c]
    b=analysis.loc[analysis.CTRCD==0,c]
    components.append({
        'label': COMPONENT_LABELS[c],
        'event': mediqr(a),
        'nonevent': mediqr(b),
    })

payload={
    'meta': {
        'clinical_n':531,'clinical_events':54,'tdi_n':270,'tdi_events':27,
        'median_followup_days':437,'median_event_time_days':316,
        'dataset':'BC_cardiotox','license':'CC BY 4.0',
        'dataset_doi':'10.6084/m9.figshare.22650748.v4'
    },
    'headline': {
        'mti_auc':0.695,'mti_auc_ci':[0.590,0.794],
        'combo_auc':0.712,'combo_auc_ci':[0.608,0.807],
        'mti_hr_lvef':2.20,'mti_hr_lvef_ci':[1.36,3.55],
        'preserved_lvef_hr':2.22,'preserved_lvef_ci':[1.41,3.51],
        'top20_capture':13/27,'top20_events':13,'top20_total_events':27,'top20_ppv':13/54,
    },
    'phase': phase.round(4).tolist(),
    'waveforms': {
        'CTRCD': summary_wave(tdi.CTRCD.to_numpy()==1),
        'No CTRCD': summary_wave(tdi.CTRCD.to_numpy()==0),
        'deciles':deciles,
    },
    'components':components,
    'risk_capture':rows,
    'model_performance':[
        {'model':'MTI','auc':0.695,'cindex':0.677,'auprc':0.192},
        {'model':'LVEF','auc':0.657,'cindex':0.678,'auprc':0.166},
        {'model':'Published vE′','auc':0.663,'cindex':0.660,'auprc':0.150},
        {'model':'MTI + LVEF','auc':0.712,'cindex':0.715,'auprc':0.211},
    ],
    'cox':[
        {'label':'MTI unadjusted','hr':2.37,'lo':1.49,'hi':3.78},
        {'label':'Adjusted for LVEF','hr':2.20,'lo':1.36,'hi':3.55},
        {'label':'Age + LVEF','hr':2.19,'lo':1.34,'hi':3.58},
        {'label':'Treatment-adjusted','hr':1.95,'lo':1.19,'hi':3.18},
        {'label':'Baseline LVEF ≥50%','hr':2.22,'lo':1.41,'hi':3.51},
    ],
    'legacy': {
        'name':'Systolic Completion Instability (SCI)',
        'modality':'2D apical 4-chamber cavity motion',
        'concept':'Temporal completion of LV cavity narrowing',
        'components':['Systolic completion debt','Delayed completion burden','Completion heterogeneity','Delayed segment fraction'],
        'status':'Legacy conceptual prototype. EchoRisk/CARDIOCARE results are not used in the current public-data study.'
    }
}
(SITE/'study-data.json').write_text(json.dumps(payload,separators=(',',':')))
print('wrote',SITE/'study-data.json', 'bytes', (SITE/'study-data.json').stat().st_size)
