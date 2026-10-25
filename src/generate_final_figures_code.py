from __future__ import annotations
import importlib.util, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from scipy.signal import savgol_filter
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score
from statsmodels.duration.hazard_regression import PHReg

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'final_figures_code'; OUT.mkdir(parents=True,exist_ok=True)
TABLES=ROOT/'results'/'tables'
SEED=20261023

# -------- style --------
NAVY='#102A43'; BLUE='#1D5D83'; TEAL='#0B7A83'; CYAN='#54A9B7'; ORANGE='#E55A2A'; ROSE='#D95B73'
GRAY='#7A8A99'; MID='#B8C4CC'; LIGHT='#F4F7F9'; PALE_BLUE='#EAF3F7'; PALE_ORANGE='#FCEFEA'; GREEN='#3C8D7B'; BLACK='#17212B'; WHITE='#FFFFFF'
CMAP=LinearSegmentedColormap.from_list('mti_div',[TEAL, '#F7F7F7', ORANGE])

mpl.rcParams.update({
    'font.family':'Liberation Sans','font.size':9,'axes.titlesize':10.5,'axes.labelsize':9.2,
    'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.linewidth':0.8,
    'axes.edgecolor':'#7F8C96','axes.labelcolor':BLACK,'text.color':BLACK,
    'xtick.color':'#586875','ytick.color':'#586875','figure.facecolor':'white','axes.facecolor':'white',
    'savefig.facecolor':'white','pdf.fonttype':42,'ps.fonttype':42
})

def clean(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#8A98A5'); ax.spines['bottom'].set_color('#8A98A5')
    ax.grid(False)

def panel_label(ax, label, title=None):
    ax.text(-0.07,1.06,label,transform=ax.transAxes,fontsize=16,fontweight='bold',color=NAVY,va='top')
    if title: ax.text(0.015,1.045,title,transform=ax.transAxes,fontsize=10.5,fontweight='bold',color=NAVY,va='top')

def rounded_panel(fig, ax, pad=0.008):
    bb=ax.get_position(); x0=bb.x0-pad; y0=bb.y0-pad; w=bb.width+2*pad; h=bb.height+2*pad
    p=FancyBboxPatch((x0,y0),w,h,transform=fig.transFigure,boxstyle='round,pad=0.004,rounding_size=0.008',
                     facecolor='white',edgecolor='#D5E1E8',linewidth=0.9,zorder=-10)
    fig.patches.append(p)

def save(fig,name):
    for ext,dpi in [('png',400),('tiff',400),('svg',None)]:
        p=OUT/f'{name}.{ext}'
        kw={'bbox_inches':'tight','pad_inches':0.08}
        if dpi: kw['dpi']=dpi
        if ext=='tiff': kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(p,**kw)
    plt.close(fig)

# -------- data --------
spec=importlib.util.spec_from_file_location('core',ROOT/'src'/'run_analysis.py'); core=importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
raw=pd.read_csv(ROOT/'data'/'BC_cardiotox_clinical_and_functional_variables.csv',sep=';',decimal=',')
tcols=[c for c in raw if c.startswith('t ')]
tdi=raw.loc[raw[tcols].notna().all(axis=1)].copy().reset_index(drop=True)
waves=tdi[tcols].to_numpy(float); smooth=savgol_filter(waves,31,3,axis=1)
feat=pd.DataFrame([core.extract_tdi_features(w) for w in waves])
a=pd.concat([tdi.drop(columns=tcols).reset_index(drop=True),feat],axis=1)
a['MTI_raw']=core.make_mti(a); a['MTI']=(a.MTI_raw-a.MTI_raw.mean())/a.MTI_raw.std(ddof=0)
pred=pd.read_csv(ROOT/'results'/'crossfitted_predictions.csv')
bench=pd.read_csv(TABLES/'reviewer_benchmark_logistic_performance.csv')
cox=pd.read_csv(TABLES/'cox_models.csv')
comp=pd.read_csv(TABLES/'mti_components.csv')
reclass=pd.read_csv(TABLES/'one_year_reclassification_nri_idi.csv')
risk1=pd.read_csv(TABLES/'one_year_crossfitted_risk.csv')
cal=pd.read_csv(TABLES/'one_year_calibration_quartiles.csv')
dca=pd.read_csv(TABLES/'one_year_decision_curve.csv')
mat=pd.read_csv(TABLES/'one_year_reclassification_quartile_matrix.csv')
perf1=pd.read_csv(TABLES/'one_year_performance_calibration.csv')
add1=pd.read_csv(TABLES/'one_year_added_value.csv')

t=np.linspace(0,1,waves.shape[1]); y=a.CTRCD.to_numpy(int)

# ---------- Figure 1: phenotype ----------
fig=plt.figure(figsize=(14,10.3),constrained_layout=False)
gs=GridSpec(3,12,figure=fig,height_ratios=[1.0,1.1,0.95],hspace=0.52,wspace=0.55,left=0.055,right=0.985,top=0.93,bottom=0.07)
fig.suptitle('Pretreatment myocardial mechanical-transition phenotype and MTI construction',x=0.055,y=0.975,ha='left',fontsize=18,fontweight='bold',color=NAVY)
fig.text(0.055,0.948,'Code-generated from the public BC_cardiotox Tissue Doppler cohort',ha='left',fontsize=9.5,color=GRAY)

# A cohort flow
ax=fig.add_subplot(gs[0,0:4]); ax.axis('off'); panel_label(ax,'A','Cohort and analysis')
boxes=[(0.02,0.52,0.25,0.25,'531','women in public\nclinical cohort',PALE_BLUE),(0.37,0.52,0.25,0.25,'270','with complete\npretreatment TDI',PALE_BLUE),(0.72,0.52,0.25,0.25,'27','subsequent\nCTRCD events',PALE_ORANGE)]
for x,y0,w,h,n,txt,fc in boxes:
    ax.add_patch(FancyBboxPatch((x,y0),w,h,boxstyle='round,pad=0.02,rounding_size=0.03',fc=fc,ec='#C7D9E3',lw=1.2,transform=ax.transAxes))
    ax.text(x+w/2,y0+0.16,n,ha='center',va='center',fontsize=22,fontweight='bold',color=NAVY,transform=ax.transAxes)
    ax.text(x+w/2,y0+0.055,txt,ha='center',va='center',fontsize=8.7,color=BLACK,transform=ax.transAxes)
for x1,x2 in [(0.27,0.37),(0.62,0.72)]:
    ax.add_patch(FancyArrowPatch((x1,0.645),(x2,0.645),arrowstyle='-|>',mutation_scale=12,color=GRAY,lw=1.3,transform=ax.transAxes))
ax.text(0.03,0.29,'Primary endpoint',fontsize=8.5,fontweight='bold',color=NAVY,transform=ax.transAxes)
ax.text(0.03,0.17,'Cancer therapy-related cardiac dysfunction (CTRCD):\nLVEF <50% with an absolute decline ≥10 percentage points.',fontsize=8.6,color=BLACK,transform=ax.transAxes,linespacing=1.35)
rounded_panel(fig,ax)

# B waveform overlay
ax=fig.add_subplot(gs[0,4:9]); panel_label(ax,'B','Pretreatment TDI waveform')
evt=smooth[y==1]; ctl=smooth[y==0]
for arr,color,label in [(ctl,TEAL,'No CTRCD'),(evt,ORANGE,'CTRCD')]:
    med=np.median(arr,axis=0); lo=np.quantile(arr,.25,axis=0); hi=np.quantile(arr,.75,axis=0)
    ax.fill_between(t,lo,hi,color=color,alpha=.11,lw=0); ax.plot(t,med,color=color,lw=2.2,label=label)
ax.axhline(0,color=MID,lw=.9,ls='--'); ax.axvspan(.42,.62,color=ORANGE,alpha=.055)
ax.text(.52,ax.get_ylim()[1]*.83,'systolic → early relaxation',ha='center',fontsize=8,color=ROSE)
ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Myocardial velocity (cm/s)'); ax.legend(frameon=False,ncol=2,loc='lower left')
clean(ax); rounded_panel(fig,ax)

# C four components as clean schematics
axc=fig.add_subplot(gs[0,9:12]); axc.axis('off'); panel_label(axc,'C','Four interpretable components')
items=[('1','Transition slope deficit','slower systolic deceleration',BLUE),('2','Transition curvature deficit','blunted transition bend',ROSE),('3','Early-relaxation latency','later E′ after systole',TEAL),('4','E′ attenuation','smaller early-relaxation velocity',ORANGE)]
for i,(num,title,sub,col) in enumerate(items):
    yy=.80-i*.21
    axc.add_patch(FancyBboxPatch((.02,yy-.12),.96,.16,boxstyle='round,pad=.012,rounding_size=.02',fc=LIGHT,ec='#D6E1E7',lw=.8,transform=axc.transAxes))
    axc.add_patch(FancyBboxPatch((.04,yy-.085),.08,.09,boxstyle='round,pad=.01,rounding_size=.02',fc=col,ec='none',transform=axc.transAxes))
    axc.text(.08,yy-.04,num,ha='center',va='center',color='white',fontweight='bold',fontsize=10,transform=axc.transAxes)
    axc.text(.15,yy-.015,title,fontweight='bold',fontsize=8.7,color=NAVY,transform=axc.transAxes)
    axc.text(.15,yy-.075,sub,fontsize=8.1,color=GRAY,transform=axc.transAxes)
axc.text(.02,.02,'MTI = mean of 4 training-standardized components',fontsize=8.2,fontweight='bold',color=NAVY,transform=axc.transAxes)
rounded_panel(fig,axc)

# D atlas
ax=fig.add_subplot(gs[1,0:7]); panel_label(ax,'D','Waveform atlas ordered by increasing MTI')
z=(smooth-smooth.mean(1,keepdims=True))/smooth.std(1,keepdims=True)
ordr=np.argsort(a.MTI.to_numpy())
im=ax.imshow(z[ordr],aspect='auto',extent=[0,1,len(a),0],cmap=CMAP,vmin=-2.4,vmax=2.4,interpolation='nearest')
ax.axvline(.5,color='white',lw=1.0,ls='--',alpha=.85); ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Women ordered by MTI')
cbar=fig.colorbar(im,ax=ax,pad=.010,fraction=.020); cbar.set_label('z velocity',fontsize=7.5); cbar.ax.tick_params(labelsize=7)
clean(ax); rounded_panel(fig,ax)

# E phase difference
ax=fig.add_subplot(gs[1,7:12]); panel_label(ax,'E','Phase-resolved event minus non-event difference')
rng=np.random.default_rng(SEED); B=2500; dif=np.empty((B,len(t)),np.float32)
for b in range(B):
    ee=evt[rng.integers(0,len(evt),len(evt))].mean(0); cc=ctl[rng.integers(0,len(ctl),len(ctl))].mean(0); dif[b]=ee-cc
delta=evt.mean(0)-ctl.mean(0); lo=np.quantile(dif,.025,axis=0); hi=np.quantile(dif,.975,axis=0)
ax.fill_between(t,lo,hi,color=ORANGE,alpha=.16,lw=0); ax.plot(t,delta,color=ORANGE,lw=2.2); ax.axhline(0,color=MID,ls='--',lw=.9); ax.axvspan(.499,.600,color=ORANGE,alpha=.08)
ax.text(.03,.06,'Global functional permutation P = 0.047\nLVEF-adjusted P = 0.056\nBest cluster FWER P = 0.094',transform=ax.transAxes,fontsize=8.3,color=NAVY,va='bottom')
ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Δ velocity (cm/s)'); clean(ax); rounded_panel(fig,ax)

# F component distributions
ax=fig.add_subplot(gs[2,0:12]); panel_label(ax,'F','Individual MTI components are modest but directionally concordant')
xs=np.arange(4); vals=comp.event_higher_oriented_AUROC.to_numpy(); lows=comp.AUROC_CI_low.to_numpy(); highs=comp.AUROC_CI_high.to_numpy()
ax.errorbar(xs,vals,yerr=[vals-lows,highs-vals],fmt='o',ms=8,capsize=4,color=NAVY,ecolor='#8CA4B3',lw=1.4)
ax.axhline(.5,color=MID,ls='--',lw=1); ax.set_ylim(.45,.82); ax.set_ylabel('Event-oriented AUROC')
labels=['Slope\ndeficit','Curvature\ndeficit','Early-relaxation\nlatency','E′\nattenuation']; ax.set_xticks(xs,labels)
for i,v in enumerate(vals): ax.text(i,v+.025,f'{v:.3f}\nP={comp.Mann_Whitney_p.iloc[i]:.3f}',ha='center',fontsize=8,color=NAVY)
clean(ax); rounded_panel(fig,ax)
save(fig,'Figure1_MTI_mechanism_CODE')

# ---------- Figure 2: clinical validation ----------
fig=plt.figure(figsize=(14,10.4)); gs=GridSpec(3,12,figure=fig,hspace=.56,wspace=.78,left=.055,right=.985,top=.92,bottom=.065)
fig.suptitle('Clinical validation and prognostic performance of Mechanical Transition Instability',x=.055,y=.97,ha='left',fontsize=18,fontweight='bold',color=NAVY)
models=[('MTI',ORANGE,'MTI'),('LVEF',BLUE,'LVEF'),('vEprime',CYAN,'Published vE′'),('MTI+LVEF',TEAL,'MTI + LVEF')]
# A ROC
ax=fig.add_subplot(gs[0,0:6]); panel_label(ax,'A','Cross-fitted discrimination')
for col,c,lab in models:
    fpr,tpr,_=roc_curve(pred.CTRCD,pred[col]); auc=roc_auc_score(pred.CTRCD,pred[col]); ax.plot(fpr,tpr,lw=2,color=c,label=f'{lab}  {auc:.3f}')
ax.plot([0,1],[0,1],ls='--',color=MID,lw=1); ax.set_xlabel('1 − specificity'); ax.set_ylabel('Sensitivity'); ax.legend(title='AUROC',frameon=False,loc='lower right'); clean(ax); rounded_panel(fig,ax)
# B PR
ax=fig.add_subplot(gs[0,6:12]); panel_label(ax,'B','Precision-recall performance')
for col,c,lab in models:
    pr,rc,_=precision_recall_curve(pred.CTRCD,pred[col]); ap=average_precision_score(pred.CTRCD,pred[col]); ax.plot(rc,pr,lw=2,color=c,label=f'{lab}  {ap:.3f}')
ax.axhline(pred.CTRCD.mean(),ls='--',color=MID,lw=1,label=f'Prevalence {pred.CTRCD.mean():.2f}'); ax.set_xlabel('Recall'); ax.set_ylabel('Precision'); ax.legend(title='AUPRC',frameon=False,loc='upper right'); clean(ax); rounded_panel(fig,ax)
# C KM
ax=fig.add_subplot(gs[1,0:5]); panel_label(ax,'C','CTRCD-free survival across MTI tertiles')
ter=pd.qcut(a.MTI,3,labels=['Low','Middle','High']); colors={'Low':TEAL,'Middle':GRAY,'High':ORANGE}
def km_xy(time,event):
    time=np.asarray(time,float); event=np.asarray(event,int); evts=np.sort(np.unique(time[event==1])); xs=[0]; ys=[1.]; s=1.
    for tt in evts:
        n=np.sum(time>=tt); d=np.sum((time==tt)&(event==1)); xs.extend([tt,tt]); ys.extend([s,s*(1-d/n)]); s=ys[-1]
    xs.append(time.max()); ys.append(s); return np.array(xs),np.array(ys)
for g in ['Low','Middle','High']:
    m=(ter.astype(str)==g).to_numpy(); xx,yy=km_xy(a.time[m],a.CTRCD[m]); ax.step(xx,yy,where='post',lw=2.2,color=colors[g],label=f'{g}: {int(a.CTRCD[m].sum())}/{int(m.sum())} events')
ax.set_xlabel('Days from treatment initiation'); ax.set_ylabel('Freedom from CTRCD'); ax.set_ylim(.70,1.01); ax.legend(frameon=False,loc='lower left'); clean(ax); rounded_panel(fig,ax)
# D forest
ax=fig.add_subplot(gs[1,5:12]); panel_label(ax,'D','Association with time to CTRCD')
rows=[('Unadjusted','MTI unadjusted'),('Adjusted for LVEF','MTI + LVEF'),('Age + LVEF','MTI + age + LVEF'),('Treatment-adjusted','MTI + age + LVEF + treatment')]
vals=[]
for lab,mod in rows:
    r=cox[(cox.model==mod)&(cox.term=='MTI')].iloc[0]; vals.append((lab,r.HR_per_1SD,r.CI_low,r.CI_high,r.p_value))
pres=pd.read_csv(TABLES/'cox_preserved_lvef_sensitivity.csv'); r=pres[pres.term=='MTI'].iloc[0]; vals.append(('Baseline LVEF ≥50%',r.HR_per_1SD,r.CI_low,r.CI_high,r.p_value))
ypos=np.arange(len(vals))[::-1]
for yp,(lab,hr,lo,hi,pv) in zip(ypos,vals):
    ax.errorbar(hr,yp,xerr=[[hr-lo],[hi-hr]],fmt='o',color=TEAL,ecolor=NAVY,capsize=3,ms=6); ax.text(4.2,yp,f'{hr:.2f} ({lo:.2f}–{hi:.2f})   P={pv:.3g}',va='center',fontsize=8.4)
ax.axvline(1,color=MID,ls='--',lw=1); ax.set_xscale('log'); ax.set_xlim(.75,7.2); ax.set_yticks(ypos,[v[0] for v in vals]); ax.set_xlabel('Hazard ratio per 1-SD higher MTI (log scale)'); clean(ax); rounded_panel(fig,ax)
# E risk concentration
ax=fig.add_subplot(gs[2,0:5]); panel_label(ax,'E','Risk concentration for surveillance triage')
n=len(pred); total=int(pred.CTRCD.sum()); fr=np.linspace(.05,1,40)
for col,c,lab in [('MTI',ORANGE,'MTI'),('LVEF',BLUE,'LVEF'),('MTI+LVEF',TEAL,'MTI + LVEF')]:
    ordx=np.argsort(-pred[col].to_numpy()); cap=[]
    for f in fr:
        k=max(1,int(np.ceil(n*f))); cap.append(pred.CTRCD.to_numpy()[ordx[:k]].sum()/total)
    ax.plot(fr*100,np.array(cap)*100,lw=2.1,color=c,label=lab)
ax.plot([0,100],[0,100],ls='--',color=MID,lw=1,label='Random selection'); ax.scatter([20],[48.1],s=45,color=TEAL,zorder=4); ax.annotate('Top 20% contains 13/27 events\nPPV 24.1% vs 10.0% prevalence',xy=(20,48.1),xytext=(34,29),arrowprops=dict(arrowstyle='->',color=GRAY),fontsize=8.2,color=NAVY)
ax.set_xlabel('Cohort selected for intensified surveillance (%)'); ax.set_ylabel('Future CTRCD events captured (%)'); ax.set_xlim(0,100); ax.set_ylim(0,103); ax.legend(frameon=False,loc='upper left',fontsize=7.3); clean(ax); rounded_panel(fig,ax)
# F benchmark
ax=fig.add_subplot(gs[2,6:12]); panel_label(ax,'F','Reviewer-facing benchmark comparison')
show=['MTI+LVEF','Clinical+MTI','MTI','Clinical','Fractional shortening','Published vEprime','LVEF','Expanded TDI landmarks','S+E+A landmarks','Sprime amplitude']
b=bench.set_index('model').loc[show].iloc[::-1]; yy=np.arange(len(b)); cols=[TEAL if 'MTI' in s else BLUE for s in b.index]
for i,(name,row) in enumerate(b.iterrows()): ax.errorbar(row.AUROC,i,xerr=[[row.AUROC-row.AUROC_low],[row.AUROC_high-row.AUROC]],fmt='o',ms=6,capsize=3,color=cols[i],ecolor='#8FA5B2')
ax.axvline(.5,color=MID,ls='--',lw=1); ax.set_yticks(yy,[x.replace('vEprime','vE′') for x in b.index]); ax.set_xlabel('Repeated cross-fitted AUROC'); ax.set_xlim(.42,.84); clean(ax); rounded_panel(fig,ax)
save(fig,'Figure2_MTI_clinical_CODE')

# ---------- Figure 3: reclassification and added value ----------
fig=plt.figure(figsize=(14,10.5)); gs=GridSpec(3,12,figure=fig,hspace=.62,wspace=.70,left=.055,right=.985,top=.92,bottom=.065)
fig.suptitle('Incremental clinical utility, reclassification, and robustness',x=.055,y=.97,ha='left',fontsize=18,fontweight='bold',color=NAVY)
fig.text(.055,.944,'Exploratory 1-year reclassification analyses use 172 women with known 1-year status (16 CTRCD events).',fontsize=9.2,color=GRAY)
# A NRI/IDI
ax=fig.add_subplot(gs[0,0:6]); panel_label(ax,'A','Continuous 1-year reclassification after adding MTI')
metrics=['NRI_continuous','NRI_event','NRI_nonevent','IDI']; labels=['Continuous NRI','Event NRI','Non-event NRI','IDI'];
row=reclass[reclass.comparison=='Clinical -> Clinical+MTI'].iloc[0]; xx=np.arange(4)
for i,(m,lab) in enumerate(zip(metrics,labels)):
    v=row[m]; lo=row[m+'_low']; hi=row[m+'_high']; ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',ms=7,capsize=4,color=TEAL if m!='IDI' else ORANGE,ecolor=NAVY)
    ax.text(.90,i,f'{v:+.3f}  [{lo:+.3f}, {hi:+.3f}]',va='center',fontsize=7.8,color=NAVY)
ax.axvline(0,color=MID,ls='--',lw=1); ax.set_yticks(xx,labels); ax.set_xlim(-.35,1.28); ax.set_xlabel('Improvement after adding MTI'); clean(ax); rounded_panel(fig,ax)
# B DCA
ax=fig.add_subplot(gs[0,6:12]); panel_label(ax,'B','Decision-curve analysis at 1 year')
for m,c,lab in [('Clinical',BLUE,'Clinical model'),('Clinical+MTI',TEAL,'Clinical + MTI'),('LVEF',GRAY,'LVEF'),('MTI+LVEF',ORANGE,'MTI + LVEF')]:
    z=dca[dca.model==m]; ax.plot(z.threshold*100,z.net_benefit,lw=2,color=c,label=lab)
for m,c,ls in [('Treat all',MID,'--'),('Treat none','#D5DDE2',':')]:
    z=dca[dca.model==m]; ax.plot(z.threshold*100,z.net_benefit,lw=1.2,color=c,ls=ls,label=m)
ax.axvspan(8,15,color=TEAL,alpha=.045); ax.annotate('At 10% threshold:\nMTI+LVEF 0.032 vs LVEF 0.014',xy=(10,.0323),xytext=(16,.062),arrowprops=dict(arrowstyle='->',color=GRAY,lw=.9),fontsize=7.6,color=NAVY); ax.set_xlabel('1-year risk threshold (%)'); ax.set_ylabel('Net benefit'); ax.set_xlim(2,30); ax.set_ylim(-.025,.09); ax.legend(frameon=False,ncol=2,loc='upper right'); clean(ax); rounded_panel(fig,ax)
# C matrices events/non-events
sub=GridSpecFromSubplotSpec(1,2,subplot_spec=gs[1,0:7],wspace=.28)
axs=[]
for j,status in enumerate(['1-year CTRCD','1-year event-free']):
    ax=fig.add_subplot(sub[0,j]); axs.append(ax)
    M=mat[mat.status==status].pivot(index='old_quartile',columns='new_quartile',values='n').reindex(index=[1,2,3,4],columns=[1,2,3,4]).fillna(0).to_numpy(int)
    im=ax.imshow(M,cmap='Blues',vmin=0,vmax=max(1,M.max()))
    for r0 in range(4):
        for c0 in range(4):
            val=M[r0,c0]; ax.text(c0,r0,str(val),ha='center',va='center',fontsize=9,fontweight='bold' if val else 'normal',color='white' if val>M.max()/2 else NAVY)
    ax.set_xticks(range(4),['Q1','Q2','Q3','Q4']); ax.set_yticks(range(4),['Q1','Q2','Q3','Q4']); ax.set_xlabel('Clinical + MTI risk quartile');
    if j==0: ax.set_ylabel('Clinical-model risk quartile')
    ax.set_title(status,fontsize=9.5,fontweight='bold',color=NAVY)
    for sp in ax.spines.values(): sp.set_visible(False)
panel_label(axs[0],'C','Exploratory risk-quartile reclassification')
# overlay note
for ax in axs: rounded_panel(fig,ax,pad=.005)
# D calibration
ax=fig.add_subplot(gs[1,7:12]); panel_label(ax,'D','Cross-fitted 1-year calibration')
for m,c,lab in [('LVEF',BLUE,'LVEF'),('MTI+LVEF',ORANGE,'MTI + LVEF'),('Clinical',GRAY,'Clinical'),('Clinical+MTI',TEAL,'Clinical + MTI')]:
    z=cal[cal.model==m]; ax.plot(z.predicted_mean,z.observed_rate,'o-',lw=1.8,ms=5,color=c,label=lab)
ax.plot([0,.3],[0,.3],ls='--',color=MID,lw=1); ax.set_xlim(0,.3); ax.set_ylim(0,.35); ax.set_xlabel('Mean predicted 1-year risk'); ax.set_ylabel('Observed 1-year event rate'); ax.legend(frameon=False,fontsize=7.5); clean(ax); rounded_panel(fig,ax)
# E benchmark expanded
ax=fig.add_subplot(gs[2,0:7]); panel_label(ax,'E','MTI remains competitive against increasingly flexible baselines')
show2=['MTI+LVEF','Clinical+MTI','MTI','Clinical','Fractional shortening','Published vEprime','Eprime amplitude','Eprime timing','Expanded TDI landmarks','S+E+A landmarks']
bb=bench.set_index('model').loc[show2].sort_values('AUROC'); yy=np.arange(len(bb))
for i,(name,row) in enumerate(bb.iterrows()):
    col=TEAL if 'MTI' in name else BLUE; ax.errorbar(row.AUROC,i,xerr=[[row.AUROC-row.AUROC_low],[row.AUROC_high-row.AUROC]],fmt='o',ms=6,capsize=3,color=col,ecolor='#91A6B2')
ax.axvline(.5,color=MID,ls='--',lw=1); ax.set_yticks(yy,[s.replace('vEprime','vE′') for s in bb.index]); ax.set_xlabel('Repeated cross-fitted AUROC'); ax.set_xlim(.45,.83); clean(ax); rounded_panel(fig,ax)
# F preserved LVEF risk enrichment
ax=fig.add_subplot(gs[2,7:12]); panel_label(ax,'F','High MTI identifies risk despite preserved baseline LVEF')
q=a.MTI.quantile(2/3); subp=a[a.LVEF>=50].copy(); high=subp.MTI>=q; rates=[subp.loc[~high,'CTRCD'].mean()*100,subp.loc[high,'CTRCD'].mean()*100]; ns=[(~high).sum(),high.sum()]; evs=[subp.loc[~high,'CTRCD'].sum(),subp.loc[high,'CTRCD'].sum()]
ax.bar([0,1],rates,width=.55,color=[PALE_BLUE,PALE_ORANGE],edgecolor=[BLUE,ORANGE],linewidth=1.5)
for i,r0 in enumerate(rates): ax.text(i,r0+1.1,f'{int(evs[i])}/{int(ns[i])}\n{r0:.1f}%',ha='center',fontweight='bold',color=NAVY)
# adjusted HR high tertile vs rest
subp['high']=(high).astype(int); fit=PHReg(subp.time.astype(float),subp[['high','LVEF']].astype(float),status=subp.CTRCD.astype(int),ties='breslow').fit(disp=0); hr=np.exp(fit.params[0]); ci=np.exp(fit.conf_int()[0]); pv=fit.pvalues[0]
ax.text(.5,.88,f'LVEF-adjusted HR {hr:.2f}\n95% CI {ci[0]:.2f}–{ci[1]:.2f}; P={pv:.003f}',transform=ax.transAxes,ha='center',va='top',fontsize=9,color=NAVY,bbox=dict(boxstyle='round,pad=.35',fc=LIGHT,ec='#D9E3E8'))
ax.set_xticks([0,1],['Lower 2 MTI tertiles','Highest MTI tertile']); ax.set_ylabel('CTRCD event rate (%)'); ax.set_ylim(0,25); clean(ax); rounded_panel(fig,ax)
save(fig,'Figure3_MTI_reclassification_CODE')

# concise results summary for manuscript/reviewer response
rcli=reclass[reclass.comparison=='Clinical -> Clinical+MTI'].iloc[0]; rlvef=reclass[reclass.comparison=='LVEF -> MTI+LVEF'].iloc[0]
summary=f'''# Added-value and reclassification analyses\n\n## Fixed 1-year analysis cohort\n- 172 women had known 1-year status; 16 developed CTRCD by 1 year.\n- Predictions were generated with 100 repeated stratified 5-fold cross-validation. MTI standardization, imputation, scaling, and model fitting occurred in training folds only.\n\n## 1-year discrimination and calibration\n{perf1.to_markdown(index=False)}\n\n## Continuous reclassification\n- LVEF -> MTI+LVEF: continuous NRI {rlvef.NRI_continuous:+.3f} (95% bootstrap CI {rlvef.NRI_continuous_low:+.3f} to {rlvef.NRI_continuous_high:+.3f}); IDI {rlvef.IDI:+.3f} ({rlvef.IDI_low:+.3f} to {rlvef.IDI_high:+.3f}).\n- Clinical -> Clinical+MTI: continuous NRI {rcli.NRI_continuous:+.3f} (95% bootstrap CI {rcli.NRI_continuous_low:+.3f} to {rcli.NRI_continuous_high:+.3f}); IDI {rcli.IDI:+.3f} ({rcli.IDI_low:+.3f} to {rcli.IDI_high:+.3f}).\n- These estimates are directionally favorable but imprecise and are exploratory.\n\n## Decision-curve interpretation\n- At a 10% 1-year risk threshold, MTI+LVEF net benefit was 0.0323 versus 0.0142 for LVEF alone. Clinical+MTI net benefit was 0.0239 versus 0.0084 for the clinical model.\n- These threshold analyses are exploratory and should not be treated as validated surveillance cutoffs.\n\n## Preserved-LVEF phenotype\n- Among 261 women with baseline LVEF >=50%, the highest MTI tertile had 15/83 events (18.1%) versus 11/178 (6.2%) in the lower two tertiles.\n- After adjustment for baseline LVEF, the highest MTI tertile had HR {hr:.2f} (95% CI {ci[0]:.2f}-{ci[1]:.2f}; P={pv:.4f}) versus the lower two tertiles. This tertile analysis is exploratory and the threshold is cohort-relative.\n'''
(ROOT/'results'/'ADDED_VALUE_RECLASSIFICATION_RESULTS.md').write_text(summary)
print(summary)
print('Saved figures to',OUT)
