from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'BC_cardiotox_clinical_and_functional_variables.csv'
PRED = ROOT / 'results' / 'crossfitted_predictions.csv'
OUT = ROOT / 'docs' / 'assets'
OUT.mkdir(parents=True, exist_ok=True)

if not RAW.exists():
    raise FileNotFoundError('Download BC_cardiotox_clinical_and_functional_variables.csv into data/.')
if not PRED.exists():
    raise FileNotFoundError('Run src/run_analysis.py before building website data.')

sys.path.insert(0, str(ROOT/'src'))
from run_analysis import extract_tdi_features, make_mti, COMPONENTS, COMPONENT_LABELS

full = pd.read_csv(RAW, sep=';', decimal=',')
tdi_cols = [c for c in full.columns if c.startswith('t ')]
tdi = full.loc[full[tdi_cols].notna().all(axis=1)].copy().reset_index(drop=True)
Xfull = tdi[tdi_cols].to_numpy(float)

features = pd.DataFrame([extract_tdi_features(w) for w in Xfull])
analysis = pd.concat([tdi.drop(columns=tdi_cols).reset_index(drop=True), features], axis=1)
analysis['MTI_visual'] = make_mti(analysis, reference=analysis).to_numpy()
analysis['MTI_visual'] = (analysis['MTI_visual']-analysis['MTI_visual'].mean())/analysis['MTI_visual'].std(ddof=0)

idx = np.unique(np.linspace(0, len(tdi_cols)-1, 241).round().astype(int))
phase = np.linspace(0,1,len(tdi_cols))[idx]
X = Xfull[:,idx]

def wave_summary(mask):
    xx = X[np.asarray(mask)]
    return {
        'mean': np.mean(xx,axis=0).round(4).tolist(),
        'sd': np.std(xx,axis=0,ddof=1).round(4).tolist(),
        'median': np.median(xx,axis=0).round(4).tolist(),
        'q1': np.quantile(xx,.25,axis=0).round(4).tolist(),
        'q3': np.quantile(xx,.75,axis=0).round(4).tolist(),
    }

ranks = pd.qcut(analysis['MTI_visual'].rank(method='first'), 10, labels=False)
deciles=[]
for d in range(10):
    mask=(ranks==d).to_numpy()
    deciles.append({
        'decile':d+1,
        'n':int(mask.sum()),
        'median_mti':round(float(analysis.loc[mask,'MTI_visual'].median()),3),
        'event_rate':round(float(analysis.loc[mask,'CTRCD'].mean()),4),
        'waveform':np.mean(X[mask],axis=0).round(4).tolist(),
    })

components=[]
interpretation={
    'transition_slope_magnitude':'Slower resolution of systolic velocity',
    'transition_curvature':'Blunted systolic-to-relaxation bend',
    'early_relaxation_latency':'Later E′ after systolic termination',
    'e_prime_magnitude':'Smaller early-diastolic velocity',
}
for c in COMPONENTS:
    ev=analysis.loc[analysis.CTRCD.eq(1),c]
    ne=analysis.loc[analysis.CTRCD.eq(0),c]
    components.append({
        'label': COMPONENT_LABELS[c].replace('E-prime',"E′"),
        'interpretation': interpretation[c],
        'event_median':round(float(ev.median()),3),
        'nonevent_median':round(float(ne.median()),3),
    })

pred=pd.read_csv(PRED)
y=pred['CTRCD'].to_numpy(int)
risk=pred['MTI+LVEF'].to_numpy(float)
order=np.argsort(-risk)
risk_capture=[]
for frac in np.arange(.05,1.001,.05):
    k=max(1,int(np.ceil(frac*len(pred))))
    ids=order[:k]
    events=int(y[ids].sum())
    risk_capture.append({
        'fraction':round(float(frac),2),
        'selected_n':int(k),
        'events':events,
        'capture':round(events/int(y.sum()),4),
        'ppv':round(events/k,4),
    })

models=[]
for model,label in [('MTI','MTI'),('LVEF','LVEF'),('vEprime',"Published E′ metric"),('MTI+LVEF','MTI + LVEF'),('Clinical','Clinical model'),('Clinical+MTI','Clinical + MTI')]:
    models.append({
        'model':label,
        'auroc':round(float(roc_auc_score(y,pred[model])),4),
        'auprc':round(float(average_precision_score(y,pred[model])),4),
    })

payload={
  'meta':{
    'clinical_n':531,'clinical_events':54,'tdi_n':270,'tdi_events':27,
    'median_followup_days':437,'median_event_time_days':316,
    'dataset':'BC_cardiotox','dataset_doi':'10.6084/m9.figshare.22650748.v4',
    'site_version':'2026-10-25'
  },
  'headline':{
    'mti_auroc':0.695,'mti_auroc_ci':[0.590,0.794],
    'mti_auprc':0.192,'mti_cindex':0.677,'mti_cindex_ci':[0.570,0.784],
    'combo_auroc':0.712,'combo_auroc_ci':[0.608,0.807],
    'combo_cindex':0.715,'combo_cindex_ci':[0.611,0.810],
    'delta_c':0.037,'delta_c_ci':[-0.059,0.137],
    'hr_lvef':2.20,'hr_lvef_ci':[1.36,3.55],
    'preserved_hr':2.22,'preserved_hr_ci':[1.41,3.51],
    'preserved_high_events':15,'preserved_high_n':83,
    'preserved_low_events':11,'preserved_low_n':178,
    'preserved_tertile_hr':3.38,'preserved_tertile_hr_ci':[1.55,7.38]
  },
  'phase':phase.round(4).tolist(),
  'waveforms':{
    'No CTRCD':wave_summary(tdi.CTRCD.to_numpy()==0),
    'CTRCD':wave_summary(tdi.CTRCD.to_numpy()==1),
    'deciles':deciles,
  },
  'components':components,
  'risk_capture':risk_capture,
  'models':models,
  'cox':[
    {'label':'MTI unadjusted','hr':2.37,'lo':1.49,'hi':3.78,'p':'<0.001'},
    {'label':'Adjusted for LVEF','hr':2.20,'lo':1.36,'hi':3.55,'p':'0.001'},
    {'label':'Age + LVEF','hr':2.19,'lo':1.34,'hi':3.58,'p':'0.002'},
    {'label':'Age + LVEF + treatment','hr':1.95,'lo':1.19,'hi':3.18,'p':'0.008'},
    {'label':'Baseline LVEF ≥50%','hr':2.22,'lo':1.41,'hi':3.51,'p':'<0.001'},
  ],
  'benchmark':[
    {'model':'MTI + LVEF','auroc':0.705,'auprc':0.213},
    {'model':'Clinical + MTI','auroc':0.704,'auprc':0.235},
    {'model':'MTI','auroc':0.684,'auprc':0.186},
    {'model':'Clinical model','auroc':0.659,'auprc':0.181},
    {'model':'Fractional shortening','auroc':0.658,'auprc':0.156},
    {'model':"Published E′ metric",'auroc':0.648,'auprc':0.141},
    {'model':'LVEF','auroc':0.645,'auprc':0.156},
    {'model':'Age + LVEF','auroc':0.633,'auprc':0.152},
    {'model':"E′ amplitude",'auroc':0.613,'auprc':0.138},
    {'model':'Expanded TDI landmarks','auroc':0.601,'auprc':0.130},
    {'model':"E′ timing",'auroc':0.592,'auprc':0.147},
    {'model':"S′ + E′ + A′ landmarks",'auroc':0.592,'auprc':0.135},
    {'model':'FPCA-5','auroc':0.619,'auprc':0.138},
    {'model':'Waveform ridge-25','auroc':0.533,'auprc':0.106},
  ],
  'added_value':{
    'n':172,'events':16,
    'lvef_auc':0.652,'lvef_mti_auc':0.698,
    'clinical_auc':0.618,'clinical_mti_auc':0.653,
    'clinical_auprc':0.144,'clinical_mti_auprc':0.175,
    'nri':0.378,'nri_ci':[-0.122,0.869],
    'idi':0.015,'idi_ci':[-0.003,0.034],
    'net_benefit_lvef_10':0.014,'net_benefit_combo_10':0.032,
  },
  'functional':{
    'global_p':0.047,'lvef_adjusted_p':0.056,'cluster_phase':[0.499,0.600],'cluster_fwer_p':0.094
  }
}

(OUT/'study-data.json').write_text(json.dumps(payload,separators=(',',':')))
print('wrote',OUT/'study-data.json')
