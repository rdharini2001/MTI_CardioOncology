from __future__ import annotations
import importlib.util, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'tables'; OUT.mkdir(parents=True,exist_ok=True)
SEED=20261023; REPEATS=100; HORIZON=365.0; BOOT=5000
spec=importlib.util.spec_from_file_location('core',ROOT/'src'/'run_analysis.py'); core=importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
raw=pd.read_csv(ROOT/'data'/'BC_cardiotox_clinical_and_functional_variables.csv',sep=';',decimal=',')
tcols=[c for c in raw if c.startswith('t ')]
tdi=raw.loc[raw[tcols].notna().all(axis=1)].copy().reset_index(drop=True)
feats=pd.DataFrame([core.extract_tdi_features(w) for w in tdi[tcols].to_numpy(float)])
a=pd.concat([tdi.drop(columns=tcols).reset_index(drop=True),feats],axis=1)
known=((a.CTRCD==1)&(a.time<=HORIZON)) | (a.time>HORIZON)
d=a.loc[known].copy().reset_index(drop=True)
y=((d.CTRCD==1)&(d.time<=HORIZON)).astype(int).to_numpy()
models=['LVEF','MTI+LVEF','Clinical','Clinical+MTI']
acc={m:np.zeros(len(d),float) for m in models}; cnt=np.zeros(len(d),int)

def fit_prob(Xtr,ytr,Xte,C=1.0):
    Xtr=np.asarray(Xtr,float); Xte=np.asarray(Xte,float)
    med=np.nanmedian(Xtr,axis=0); Xtr=np.where(np.isnan(Xtr),med,Xtr); Xte=np.where(np.isnan(Xte),med,Xte)
    mu=Xtr.mean(0); sd=Xtr.std(0); sd[sd==0]=1
    A=(Xtr-mu)/sd; B=(Xte-mu)/sd
    lr=LogisticRegression(penalty='l2',C=C,solver='liblinear',max_iter=3000,random_state=SEED)
    lr.fit(A,ytr)
    return lr.predict_proba(B)[:,1]

for r in range(REPEATS):
    cv=StratifiedKFold(5,shuffle=True,random_state=SEED+r)
    for tr,te in cv.split(np.zeros(len(d)),y):
        A=d.iloc[tr].copy(); B=d.iloc[te].copy(); yy=y[tr]
        mti_tr=core.make_mti(A,reference=A).to_numpy(); mti_te=core.make_mti(B,reference=A).to_numpy()
        acc['LVEF'][te]+=fit_prob(A[['LVEF']],yy,B[['LVEF']])
        acc['MTI+LVEF'][te]+=fit_prob(np.c_[mti_tr,A[['LVEF']].to_numpy()],yy,np.c_[mti_te,B[['LVEF']].to_numpy()])
        acc['Clinical'][te]+=fit_prob(A[['age','LVEF','AC','antiHER2']],yy,B[['age','LVEF','AC','antiHER2']])
        acc['Clinical+MTI'][te]+=fit_prob(np.c_[A[['age','LVEF','AC','antiHER2']].to_numpy(),mti_tr],yy,np.c_[B[['age','LVEF','AC','antiHER2']].to_numpy(),mti_te])
        cnt[te]+=1
pred=pd.DataFrame({'row_id':d.index,'CTRCD_1y':y,'time':d.time})
for m in models: pred[m]=np.clip(acc[m]/cnt,1e-6,1-1e-6)
pred.to_csv(OUT/'one_year_crossfitted_risk.csv',index=False)

def nri_idi(old,new,y):
    old=np.asarray(old,float); new=np.asarray(new,float); y=np.asarray(y,int); delta=new-old; eps=1e-12
    e=y==1; n=y==0
    eu=np.mean(delta[e]>eps); ed=np.mean(delta[e]<-eps); nd=np.mean(delta[n]<-eps); nu=np.mean(delta[n]>eps)
    nri_e=eu-ed; nri_ne=nd-nu; nri=nri_e+nri_ne
    ds_old=np.mean(old[e])-np.mean(old[n]); ds_new=np.mean(new[e])-np.mean(new[n]); idi=ds_new-ds_old
    return {'NRI_event':nri_e,'NRI_nonevent':nri_ne,'NRI_continuous':nri,'IDI':idi,
            'event_up':eu,'event_down':ed,'nonevent_down':nd,'nonevent_up':nu,'slope_old':ds_old,'slope_new':ds_new}

rng=np.random.default_rng(SEED+91); rows=[]
for oldn,newn in [('LVEF','MTI+LVEF'),('Clinical','Clinical+MTI')]:
    point=nri_idi(pred[oldn],pred[newn],y); boot={k:[] for k in ['NRI_event','NRI_nonevent','NRI_continuous','IDI']}
    pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
    for b in range(BOOT):
        idx=np.r_[rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)]
        z=nri_idi(pred[oldn].to_numpy()[idx],pred[newn].to_numpy()[idx],y[idx])
        for k in boot: boot[k].append(z[k])
    row={'comparison':f'{oldn} -> {newn}',**point}
    for k,v in boot.items():
        vv=np.asarray(v,float); row[k+'_low']=np.quantile(vv,.025); row[k+'_high']=np.quantile(vv,.975); row[k+'_bootstrap_p']=2*min(np.mean(vv<=0),np.mean(vv>=0))
    rows.append(row)
nri=pd.DataFrame(rows); nri.to_csv(OUT/'one_year_reclassification_nri_idi.csv',index=False)

# performance + Brier
perf=[]
for m in models:
    p=pred[m].to_numpy(); perf.append({'model':m,'n':len(y),'events':int(y.sum()),'AUROC_1y':roc_auc_score(y,p),'AUPRC_1y':average_precision_score(y,p),'Brier_1y':brier_score_loss(y,p)})
pd.DataFrame(perf).to_csv(OUT/'one_year_performance_calibration.csv',index=False)

# calibration quartiles
cal=[]
for m in models:
    p=pred[m].to_numpy(); q=np.quantile(p,[0,.25,.5,.75,1]); q[0]-=1e-12; q[-1]+=1e-12; g=np.digitize(p,q[1:-1],right=True)
    for k in range(4):
        ix=g==k; cal.append({'model':m,'quartile':k+1,'n':int(ix.sum()),'predicted_mean':np.mean(p[ix]),'observed_rate':np.mean(y[ix])})
pd.DataFrame(cal).to_csv(OUT/'one_year_calibration_quartiles.csv',index=False)

# decision curve
thr=np.linspace(.02,.30,57); prev=y.mean(); dr=[]
for m in models:
    p=pred[m].to_numpy()
    for pt in thr:
        hi=p>=pt; tp=np.sum(hi&(y==1))/len(y); fp=np.sum(hi&(y==0))/len(y); dr.append({'model':m,'threshold':pt,'net_benefit':tp-fp*pt/(1-pt)})
for pt in thr:
    dr.append({'model':'Treat all','threshold':pt,'net_benefit':prev-(1-prev)*pt/(1-pt)}); dr.append({'model':'Treat none','threshold':pt,'net_benefit':0.0})
pd.DataFrame(dr).to_csv(OUT/'one_year_decision_curve.csv',index=False)

# rank quartile reclassification for clinical vs clinical+MTI; exploratory, thresholds from baseline model distribution
old=pred['Clinical'].to_numpy(); new=pred['Clinical+MTI'].to_numpy(); cuts=np.quantile(old,[.25,.5,.75]); oc=np.digitize(old,cuts); nc=np.digitize(new,cuts)
mat=[]; direct=[]
for status,label in [(1,'1-year CTRCD'),(0,'1-year event-free')]:
    m=y==status; mm=np.zeros((4,4),int)
    for o,nw in zip(oc[m],nc[m]): mm[o,nw]+=1
    for i in range(4):
        for j in range(4): mat.append({'status':label,'old_quartile':i+1,'new_quartile':j+1,'n':int(mm[i,j])})
    delta=nc[m]-oc[m]; direct.append({'status':label,'n':int(m.sum()),'up_pct':100*np.mean(delta>0),'same_pct':100*np.mean(delta==0),'down_pct':100*np.mean(delta<0)})
pd.DataFrame(mat).to_csv(OUT/'one_year_reclassification_quartile_matrix.csv',index=False); pd.DataFrame(direct).to_csv(OUT/'one_year_reclassification_direction.csv',index=False)

# bootstrap AUROC difference / Brier difference for added value clinical+MTI vs clinical
add=[]
for oldn,newn in [('LVEF','MTI+LVEF'),('Clinical','Clinical+MTI')]:
    pa=pred[oldn].to_numpy(); pb=pred[newn].to_numpy();
    dau=roc_auc_score(y,pb)-roc_auc_score(y,pa); db=brier_score_loss(y,pb)-brier_score_loss(y,pa)
    va=[]; vb=[]
    for b in range(BOOT):
        idx=np.r_[rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)]
        va.append(roc_auc_score(y[idx],pb[idx])-roc_auc_score(y[idx],pa[idx])); vb.append(brier_score_loss(y[idx],pb[idx])-brier_score_loss(y[idx],pa[idx]))
    add.append({'comparison':f'{oldn} -> {newn}','delta_AUROC':dau,'delta_AUROC_low':np.quantile(va,.025),'delta_AUROC_high':np.quantile(va,.975),'delta_Brier':db,'delta_Brier_low':np.quantile(vb,.025),'delta_Brier_high':np.quantile(vb,.975)})
pd.DataFrame(add).to_csv(OUT/'one_year_added_value.csv',index=False)

print('1-year evaluable n/events:',len(y),int(y.sum()))
print('\nPerformance\n',pd.DataFrame(perf).to_string(index=False))
print('\nNRI/IDI\n',nri.to_string(index=False))
print('\nReclassification direction\n',pd.DataFrame(direct).to_string(index=False))
print('\nAdded value\n',pd.DataFrame(add).to_string(index=False))
