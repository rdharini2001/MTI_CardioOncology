from __future__ import annotations
from pathlib import Path
import importlib.util
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from scipy.signal import savgol_filter
from sklearn.metrics import roc_curve, roc_auc_score
from statsmodels.duration.survfunc import SurvfuncRight

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'individual_cardiomorph_code'
OUT.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / 'results' / 'tables'
SEED = 20261025

# Shared CardioMorph / Nature-like palette
NAVY = '#183B56'
DEEP = '#0B5366'
TEAL = '#0F7C82'
CYAN = '#65B8C2'
ORANGE = '#D95D2B'
ORANGE2 = '#E97A46'
SLATE = '#6F8191'
MID = '#A9B5BE'
LIGHT = '#F5F6F6'
HEADER = '#ECEEEE'
PALE_TEAL = '#E7F2F3'
PALE_ORANGE = '#FAEEE8'
WHITE = '#FFFFFF'
BLACK = '#202A33'
CMAP = LinearSegmentedColormap.from_list('mti', [TEAL, '#F7F7F4', ORANGE2])

mpl.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9.5,
    'axes.linewidth': 1.0,
    'axes.edgecolor': '#536675',
    'axes.labelcolor': BLACK,
    'xtick.color': '#455A64',
    'ytick.color': '#455A64',
    'figure.facecolor': WHITE,
    'axes.facecolor': WHITE,
    'savefig.facecolor': WHITE,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})


def setup_panel(title: str, figsize=(9.2, 6.9), header_h=0.13):
    fig = plt.figure(figsize=figsize)
    # CardioMorph-like dashed outer boundary and header strip
    fig.patches.extend([
        Rectangle((0.02, 0.025), 0.96, 0.95, transform=fig.transFigure,
                  fill=False, edgecolor='#C4C8CA', linewidth=1.4,
                  linestyle=(0, (5, 4)), zorder=-20),
        Rectangle((0.02, 0.975-header_h), 0.96, header_h, transform=fig.transFigure,
                  facecolor=HEADER, edgecolor='none', zorder=-25)
    ])
    fig.text(0.045, 0.935, title, fontsize=18, fontweight='bold', color=BLACK,
             va='center', ha='left')
    ax = fig.add_axes([0.10, 0.12, 0.84, 0.70])
    return fig, ax


def clean(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#536675')
    ax.spines['bottom'].set_color('#536675')
    ax.tick_params(length=4, width=0.9)
    ax.grid(False)


def save(fig, stem):
    for ext, dpi in [('png', 300), ('svg', None), ('tiff', 300)]:
        p = OUT / f'{stem}.{ext}'
        kwargs = dict(bbox_inches='tight', pad_inches=0.05)
        if dpi:
            kwargs['dpi'] = dpi
        if ext == 'tiff':
            kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}
        fig.savefig(p, **kwargs)
    plt.close(fig)


def draw_person(ax, x, y, scale=1.0, color=TEAL, female=True, alpha=1.0):
    """Minimal journal-style person glyph drawn only with Matplotlib primitives."""
    ax.add_patch(Circle((x, y+0.055*scale), 0.024*scale, fc=color, ec='none', alpha=alpha))
    if female:
        ax.add_patch(plt.Polygon([[x, y+0.025*scale], [x-0.035*scale, y-0.045*scale], [x+0.035*scale, y-0.045*scale]],
                                 closed=True, fc=color, ec='none', alpha=alpha))
    else:
        ax.add_patch(Rectangle((x-0.022*scale, y-0.045*scale), 0.044*scale, 0.070*scale, fc=color, ec='none', alpha=alpha))
    ax.plot([x-0.012*scale, x-0.012*scale], [y-0.045*scale, y-0.095*scale], color=color, lw=2.4*scale, solid_capstyle='round', alpha=alpha)
    ax.plot([x+0.012*scale, x+0.012*scale], [y-0.045*scale, y-0.095*scale], color=color, lw=2.4*scale, solid_capstyle='round', alpha=alpha)


# ---------- data ----------
spec = importlib.util.spec_from_file_location('core', ROOT/'src'/'run_analysis.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
raw = pd.read_csv(ROOT/'data'/'BC_cardiotox_clinical_and_functional_variables.csv', sep=';', decimal=',')
tcols = [c for c in raw.columns if c.startswith('t ')]
tdi = raw.loc[raw[tcols].notna().all(axis=1)].copy().reset_index(drop=True)
waves = tdi[tcols].to_numpy(float)
smooth = savgol_filter(waves, 31, 3, axis=1)
features = pd.DataFrame([core.extract_tdi_features(w) for w in waves])
a = pd.concat([tdi.drop(columns=tcols).reset_index(drop=True), features], axis=1)
a['MTI_raw'] = core.make_mti(a)
a['MTI'] = (a['MTI_raw'] - a['MTI_raw'].mean()) / a['MTI_raw'].std(ddof=0)
t = np.linspace(0, 1, waves.shape[1])
y = a['CTRCD'].to_numpy(int)
pred = pd.read_csv(ROOT/'results'/'crossfitted_predictions.csv')
cox = pd.read_csv(TABLES/'cox_models.csv')
pres = pd.read_csv(TABLES/'cox_preserved_lvef_sensitivity.csv')
reclass = pd.read_csv(TABLES/'one_year_reclassification_nri_idi.csv')
dca = pd.read_csv(TABLES/'one_year_decision_curve.csv')
mat = pd.read_csv(TABLES/'one_year_reclassification_quartile_matrix.csv')

# ---------- 01 Cohort flow ----------
fig, ax = setup_panel('BC_cardiotox cohort and analytic sample', figsize=(10.0, 5.6))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
boxes = [
    (0.02, 0.36, 0.27, 0.45, PALE_TEAL, TEAL, '531', 'women with HER2-positive\nbreast cancer', '54 CTRCD events in source cohort'),
    (0.365, 0.36, 0.27, 0.45, PALE_TEAL, TEAL, '270', 'with complete pretreatment\nTDI waveforms', 'analysis cohort'),
    (0.71, 0.36, 0.27, 0.45, PALE_ORANGE, ORANGE, '27', 'subsequent CTRCD events', 'within TDI analysis cohort'),
]
for j,(x0,y0,w,h,fc,ec,n,main,sub) in enumerate(boxes):
    ax.add_patch(FancyBboxPatch((x0,y0),w,h, boxstyle='round,pad=0.018,rounding_size=0.03',
                                fc=fc, ec=ec, lw=1.3))
    # people icons
    if j == 0:
        for k,xx in enumerate(np.linspace(x0+0.055,x0+w-0.055,5)):
            draw_person(ax, xx, y0+0.31, scale=0.78, color=TEAL if k==2 else '#87A7B3', alpha=1.0)
    elif j == 1:
        for k,xx in enumerate(np.linspace(x0+0.08,x0+w-0.08,3)):
            draw_person(ax, xx, y0+0.31, scale=0.82, color=TEAL if k==1 else '#87A7B3')
    else:
        draw_person(ax, x0+w/2, y0+0.31, scale=0.95, color=ORANGE)
    ax.text(x0+w/2, y0+0.205, n, ha='center', va='center', fontsize=29, fontweight='bold', color=ec)
    ax.text(x0+w/2, y0+0.11, main, ha='center', va='center', fontsize=11.5, color=BLACK, linespacing=1.25)
    ax.text(x0+w/2, y0+0.035, sub, ha='center', va='bottom', fontsize=8.6, color=SLATE)
for x1,x2 in [(0.29,0.365),(0.635,0.71)]:
    ax.add_patch(FancyArrowPatch((x1+0.008,0.585),(x2-0.008,0.585), arrowstyle='-|>',
                                 mutation_scale=18, color=SLATE, lw=1.8))
ax.add_patch(FancyBboxPatch((0.06,0.07),0.88,0.18,boxstyle='round,pad=0.012,rounding_size=0.025',
                            fc='#F2F4F5',ec='#C7CCCF',lw=1.0))
ax.text(0.50,0.18,'Primary endpoint: cancer therapy-related cardiac dysfunction (CTRCD)',
        ha='center',va='center',fontsize=11.2,fontweight='bold',color=BLACK)
ax.text(0.50,0.11,'Absolute LVEF decrease ≥10 percentage points to <50%',
        ha='center',va='center',fontsize=10.5,color=NAVY)
save(fig,'01_cohort_flow')

# ---------- 02 Pretreatment waveform mean ± SD ----------
fig, ax = setup_panel('Pretreatment TDI waveform phenotype')
ctl = smooth[y==0]; evt = smooth[y==1]
for arr,c,label in [(ctl,TEAL,f'No CTRCD (n={len(ctl)})'),(evt,ORANGE,f'CTRCD (n={len(evt)})')]:
    m = arr.mean(axis=0); sd = arr.std(axis=0,ddof=1)
    ax.fill_between(t,m-sd,m+sd,color=c,alpha=0.13,lw=0)
    ax.plot(t,m,color=c,lw=2.5,label=label)
ax.axhline(0,color=MID,lw=1,ls='--')
# actual median zero-crossing and E' region from all patients
zero_t = float(np.median(a['systolic_zero_time']))
e_t = float(np.median(a['e_prime_time']))
ax.axvline(zero_t,color=SLATE,lw=1.1,ls='--')
ax.axvspan(max(zero_t,0.46), min(e_t+0.07,0.72), color=ORANGE, alpha=0.07)
mean_ctl=ctl.mean(axis=0); iS=np.argmax(mean_ctl[(t>=.02)&(t<=.40)])+np.where((t>=.02)&(t<=.40))[0][0]
iE=np.argmin(mean_ctl[(t>=.45)&(t<=.78)])+np.where((t>=.45)&(t<=.78))[0][0]
ax.scatter([t[iS]],[mean_ctl[iS]],s=48,color=TEAL,edgecolor='white',zorder=5)
ax.text(t[iS],mean_ctl[iS]+1.0,'S',ha='center',fontweight='bold',color=NAVY,fontsize=12)
ax.scatter([t[iE]],[mean_ctl[iE]],s=48,color=TEAL,edgecolor='white',zorder=5)
ax.text(t[iE],mean_ctl[iE]-1.25,"E′",ha='center',fontweight='bold',color=NAVY,fontsize=12)
ax.annotate('Systole', xy=(0.19,0.98), xytext=(0.19,1.08), xycoords='axes fraction', textcoords='axes fraction',
            ha='center',fontsize=10,color=SLATE,annotation_clip=False)
ax.annotate('Early relaxation', xy=(0.67,0.98), xytext=(0.67,1.08), xycoords='axes fraction', textcoords='axes fraction',
            ha='center',fontsize=10,color=SLATE,annotation_clip=False)
ax.text(0.02,0.03,'Mean waveform ± 1 SD',transform=ax.transAxes,fontsize=9,color=SLATE)
ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Myocardial velocity (cm/s)')
ax.legend(frameon=False,loc='upper right'); clean(ax)
save(fig,'02_pretreatment_tdi_waveforms')

# ---------- 03 MTI construction actual representative waveform ----------
fig, ax = setup_panel('MTI construction from interpretable transition features')
mean_wave = smooth.mean(axis=0)
dist = np.sqrt(((smooth-mean_wave)**2).mean(axis=1))
cand = (a['MTI'].abs()<0.5) & (a['s_peak_velocity']>5) & (a['e_prime_velocity']<-3) & a['systolic_zero_time'].between(.35,.55) & a['e_prime_time'].between(.50,.72)
ids = np.where(cand.to_numpy())[0]
score = dist[ids] + 0.8*np.abs(a.loc[ids,'MTI'].to_numpy())
rep_idx = int(ids[np.argmin(score)])
rep_wave = smooth[rep_idx]
rep = features.iloc[rep_idx]
ax.plot(t,rep_wave,color=DEEP,lw=2.5)
ax.axhline(0,color=MID,lw=1,ls='--')
z0=float(rep.systolic_zero_time); e0=float(rep.e_prime_time); s0=float(rep.s_peak_time)
ax.axvline(z0,color=SLATE,lw=1.0,ls='--')
ax.axvspan(z0,e0,color=ORANGE,alpha=.055)
# slope annotation
mid_s=(s0+z0)/2; ymid=np.interp(mid_s,t,rep_wave)
ax.annotate('1  Transition slope deficit\nslower systolic deceleration',xy=(mid_s,ymid),xytext=(0.06,0.87),textcoords='axes fraction',
            arrowprops=dict(arrowstyle='-|>',lw=1.4,color=TEAL),bbox=dict(boxstyle='round,pad=.4',fc=PALE_TEAL,ec=TEAL),
            fontsize=9.5,color=BLACK)
# curvature around zero crossing
ax.annotate('2  Transition curvature deficit\nblunted transition bend',xy=(z0,np.interp(z0,t,rep_wave)),xytext=(0.51,0.79),textcoords='axes fraction',
            arrowprops=dict(arrowstyle='-|>',lw=1.4,color=ORANGE),bbox=dict(boxstyle='round,pad=.4',fc=PALE_ORANGE,ec=ORANGE),
            fontsize=9.5,color=BLACK)
# latency bracket
br_y=0.8
ax.annotate('',xy=(e0,br_y),xytext=(z0,br_y),arrowprops=dict(arrowstyle='<->',lw=1.6,color=TEAL))
ax.annotate('3  Early-relaxation latency\nzero crossing → E′',xy=((z0+e0)/2,br_y),xytext=(0.60,0.56),textcoords='axes fraction',
            arrowprops=dict(arrowstyle='-|>',lw=1.2,color=TEAL),bbox=dict(boxstyle='round,pad=.4',fc=PALE_TEAL,ec=TEAL),fontsize=9.5)
# E attenuation
ye=np.interp(e0,t,rep_wave); ax.scatter([e0],[ye],s=58,color=ORANGE,edgecolor='white',zorder=4)
ax.annotate("4  E′ attenuation\nsmaller early-relaxation velocity",xy=(e0,ye),xytext=(0.68,0.27),textcoords='axes fraction',
            arrowprops=dict(arrowstyle='-|>',lw=1.4,color=ORANGE),bbox=dict(boxstyle='round,pad=.4',fc=PALE_ORANGE,ec=ORANGE),fontsize=9.5)
ax.text(e0,ye-0.9,"E′",ha='center',fontsize=11,fontweight='bold',color=NAVY)
ax.text(0.02,0.03,f'Representative waveform: woman closest to cohort median MTI (MTI={a.loc[rep_idx,"MTI"]:.2f} SD)',transform=ax.transAxes,fontsize=8.8,color=SLATE)
ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Myocardial velocity (cm/s)'); clean(ax)
save(fig,'03_mti_construction')

# ---------- 04 Real atlas ----------
fig, ax = setup_panel('Mechanical phenotype across patients')
z=(smooth-smooth.mean(axis=1,keepdims=True))/smooth.std(axis=1,keepdims=True)
ordr=np.argsort(a['MTI'].to_numpy())
im=ax.imshow(z[ordr],aspect='auto',extent=[0,1,len(a),0],cmap=CMAP,vmin=-2.3,vmax=2.3,interpolation='nearest')
ax.axvline(np.median(a['systolic_zero_time']),color=NAVY,lw=1.0,ls='--')
ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Women ordered by increasing MTI')
cb=fig.colorbar(im,ax=ax,pad=0.015,fraction=0.035); cb.set_label('Within-patient standardized velocity',fontsize=10); cb.ax.tick_params(labelsize=9)
clean(ax)
save(fig,'04_mti_ordered_waveform_atlas')

# ---------- 05 Phase-resolved difference + bootstrap ----------
fig, ax = setup_panel('Phase-resolved difference in pretreatment myocardial motion')
rng=np.random.default_rng(SEED); B=1000
delta=evt.mean(axis=0)-ctl.mean(axis=0)
boot=np.empty((B,len(t)),dtype=np.float32)
for b in range(B):
    be=evt[rng.integers(0,len(evt),len(evt))].mean(axis=0)
    bc=ctl[rng.integers(0,len(ctl),len(ctl))].mean(axis=0)
    boot[b]=be-bc
lo=np.quantile(boot,.025,axis=0); hi=np.quantile(boot,.975,axis=0)
ax.fill_between(t,lo,hi,color=ORANGE,alpha=.16,lw=0,label='Bootstrap 95% CI')
ax.plot(t,delta,color=ORANGE,lw=2.5)
ax.axhline(0,color=MID,lw=1,ls='--'); ax.axvspan(.499,.600,color=ORANGE,alpha=.08)
ax.text(0.025,0.06,'Global functional permutation: P=0.047\nAfter residualizing LVEF: P=0.056\nBest contiguous cluster: phase 0.499–0.600, cluster-FWER P=0.094',
        transform=ax.transAxes,fontsize=10.0,color=NAVY,va='bottom')
ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Δ myocardial velocity\n(CTRCD − no CTRCD), cm/s'); clean(ax)
save(fig,'05_phase_resolved_difference')

# ---------- 06 Cross-fitted ROC ----------
fig, ax = setup_panel('Cross-fitted discrimination')
models=[('MTI',ORANGE,'MTI'),('LVEF','#245B79','LVEF'),('vEprime',CYAN,"Published vE′"),('MTI+LVEF',TEAL,'MTI + LVEF')]
for col,c,lab in models:
    fpr,tpr,_=roc_curve(pred['CTRCD'],pred[col]); auc=roc_auc_score(pred['CTRCD'],pred[col])
    ax.plot(fpr,tpr,color=c,lw=2.2,label=f'{lab}   {auc:.3f}')
ax.plot([0,1],[0,1],color=MID,lw=1.2,ls='--')
ax.set_xlim(0,1); ax.set_ylim(0,1.02); ax.set_xlabel('1 − specificity'); ax.set_ylabel('Sensitivity')
leg=ax.legend(title='AUROC',frameon=False,loc='lower right'); leg._legend_box.align='left'; clean(ax)
save(fig,'06_crossfitted_roc')

# ---------- 07 KM with actual CI ----------
fig, ax = setup_panel('CTRCD-free survival across MTI tertiles')
ter=pd.qcut(a['MTI'],3,labels=['Low','Middle','High'])
colors={'Low':TEAL,'Middle':SLATE,'High':ORANGE}
for g in ['Low','Middle','High']:
    m=(ter.astype(str)==g).to_numpy(); times=a.loc[m,'time'].to_numpy(float); ev=a.loc[m,'CTRCD'].to_numpy(int)
    sf=SurvfuncRight(times,ev)
    st=np.r_[0,sf.surv_times,times.max()]; sp=np.r_[1,sf.surv_prob,sf.surv_prob[-1]]; se=np.r_[0,sf.surv_prob_se,sf.surv_prob_se[-1]]
    lo=np.clip(sp-1.96*se,0,1); hi=np.clip(sp+1.96*se,0,1)
    ax.step(st,sp,where='post',lw=2.3,color=colors[g],label=f'{g}: {int(ev.sum())}/{int(m.sum())} events')
    ax.fill_between(st,lo,hi,step='post',color=colors[g],alpha=.11,lw=0)
    # censor marks at a sparse subset of censor times using last survival value
    cen=np.sort(times[ev==0]);
    if len(cen):
        sel=cen[np.linspace(0,len(cen)-1,min(9,len(cen))).astype(int)]
        vals=np.array([sp[np.searchsorted(st,x,side='right')-1] for x in sel])
        ax.plot(sel,vals,linestyle='none',marker='|',markersize=9,markeredgewidth=1.4,color=colors[g])
ax.set_xlabel('Days from treatment initiation'); ax.set_ylabel('Freedom from CTRCD'); ax.set_ylim(.70,1.01)
ax.legend(frameon=False,loc='lower left'); clean(ax)
save(fig,'07_kaplan_meier_mti_tertiles')

# ---------- 08 Cox forest ----------
fig, ax = setup_panel('Association with time to CTRCD', figsize=(9.6,6.3))
rows=[('Unadjusted','MTI unadjusted'),('Adjusted for LVEF','MTI + LVEF'),('Adjusted for age + LVEF','MTI + age + LVEF'),('Treatment-adjusted','MTI + age + LVEF + treatment')]
vals=[]
for lab,mod in rows:
    r=cox[(cox['model']==mod)&(cox['term']=='MTI')].iloc[0]
    vals.append((lab,float(r.HR_per_1SD),float(r.CI_low),float(r.CI_high),float(r.p_value)))
r=pres[pres['term']=='MTI'].iloc[0]
vals.append(('Baseline LVEF ≥50%',float(r.HR_per_1SD),float(r.CI_low),float(r.CI_high),float(r.p_value)))
ypos=np.arange(len(vals))[::-1]
for yp,(lab,hr,lo_,hi_,pv) in zip(ypos,vals):
    ax.errorbar(hr,yp,xerr=[[hr-lo_],[hi_-hr]],fmt='o',ms=7,capsize=4,color=TEAL,ecolor=NAVY,lw=1.4)
    ax.text(4.2,yp,f'{hr:.2f} ({lo_:.2f}–{hi_:.2f})',va='center',fontsize=9.2,color=BLACK)
ax.axvline(1,color=MID,lw=1.1,ls='--'); ax.set_xscale('log'); ax.set_xlim(.75,6.7)
ax.set_yticks(ypos,[v[0] for v in vals]); ax.set_xlabel('Hazard ratio per 1-SD higher MTI (log scale)'); ax.set_ylim(-.6,len(vals)-.35)
ax.text(.98,.98,'HR (95% CI)',transform=ax.transAxes,ha='right',va='top',fontsize=9,fontweight='bold',color=NAVY)
clean(ax)
save(fig,'08_cox_forest')

# ---------- 09 Surveillance / event capture ----------
fig, ax = setup_panel('Risk concentration for surveillance triage')
n=len(pred); total=int(pred['CTRCD'].sum()); fr=np.linspace(.05,1,40)
for col,c,lab in [('MTI',ORANGE,'MTI'),('LVEF','#245B79','LVEF'),('MTI+LVEF',TEAL,'MTI + LVEF')]:
    order=np.argsort(-pred[col].to_numpy()); cap=[]
    for f in fr:
        k=max(1,int(np.ceil(n*f))); cap.append(pred['CTRCD'].to_numpy()[order[:k]].sum()/total)
    ax.plot(fr*100,np.array(cap)*100,color=c,lw=2.2,label=lab)
ax.plot([0,100],[0,100],color=MID,lw=1.1,ls='--',label='Random selection')
# exact top 20 and bootstrap CI for MTI+LVEF
k=int(np.ceil(n*.20)); order=np.argsort(-pred['MTI+LVEF'].to_numpy()); captured=int(pred['CTRCD'].to_numpy()[order[:k]].sum()); capture=captured/total
ppv=captured/k
rng=np.random.default_rng(SEED+1); vals=[]
for b in range(3000):
    idx=rng.integers(0,n,n); d=pred.iloc[idx].reset_index(drop=True); tot=d.CTRCD.sum()
    if tot==0: continue
    kk=max(1,int(np.ceil(len(d)*.20))); oo=np.argsort(-d['MTI+LVEF'].to_numpy()); vals.append(d.CTRCD.to_numpy()[oo[:kk]].sum()/tot)
ci=np.quantile(vals,[.025,.975])
ax.scatter([20],[capture*100],s=65,color=TEAL,edgecolor='white',zorder=4)
ax.annotate(f'Top 20% captures {captured}/{total} events ({capture*100:.1f}%)\n95% CI {ci[0]*100:.0f}–{ci[1]*100:.0f}% · PPV {ppv*100:.1f}%',
            xy=(20,capture*100),xytext=(36,28),arrowprops=dict(arrowstyle='-|>',color=SLATE,lw=1.2),fontsize=10,color=NAVY)
ax.set_xlabel('Cohort selected for intensified surveillance (%)'); ax.set_ylabel('Future CTRCD events captured (%)'); ax.set_xlim(0,100); ax.set_ylim(0,103)
ax.legend(frameon=False,loc='upper left'); clean(ax)
save(fig,'09_surveillance_event_capture')

# ---------- 10 Reclassification + decision curve ----------
fig = plt.figure(figsize=(11.0,6.8))
fig.patches.extend([
    Rectangle((0.02,0.025),0.96,0.95,transform=fig.transFigure,fill=False,edgecolor='#C4C8CA',linewidth=1.4,linestyle=(0,(5,4)),zorder=-20),
    Rectangle((0.02,0.845),0.96,0.13,transform=fig.transFigure,facecolor=HEADER,edgecolor='none',zorder=-25)
])
fig.text(0.045,0.91,'Incremental 1-year risk reclassification and clinical utility',fontsize=18,fontweight='bold',color=BLACK,va='center')
ax1=fig.add_axes([0.08,0.14,0.38,0.62]); ax2=fig.add_axes([0.55,0.14,0.39,0.62])
# NRI/IDI forest
row=reclass[reclass['comparison']=='Clinical -> Clinical+MTI'].iloc[0]
metrics=[('NRI_continuous','Continuous NRI'),('NRI_event','Event NRI'),('NRI_nonevent','Non-event NRI'),('IDI','IDI')]
yy=np.arange(4)[::-1]
for yp,(m,label) in zip(yy,metrics):
    v=float(row[m]); lo_=float(row[m+'_low']); hi_=float(row[m+'_high'])
    c=TEAL if m!='IDI' else ORANGE
    ax1.errorbar(v,yp,xerr=[[v-lo_],[hi_-v]],fmt='o',ms=7,capsize=4,color=c,ecolor=NAVY,lw=1.3)
    ax1.text(.94,yp,f'{v:+.3f}',va='center',fontsize=9,color=BLACK)
ax1.axvline(0,color=MID,lw=1,ls='--'); ax1.set_yticks(yy,[m[1] for m in metrics]); ax1.set_xlim(-.4,1.15); ax1.set_xlabel('Change after adding MTI to clinical model')
ax1.set_title('Continuous reclassification',loc='left',fontweight='bold',fontsize=12,color=NAVY); clean(ax1)
# decision curve
for m,c,lab in [('Clinical','#245B79','Clinical'),('Clinical+MTI',TEAL,'Clinical + MTI'),('LVEF',SLATE,'LVEF'),('MTI+LVEF',ORANGE,'MTI + LVEF')]:
    z=dca[dca['model']==m]
    ax2.plot(z['threshold']*100,z['net_benefit'],color=c,lw=2,label=lab)
for m,c,ls in [('Treat all',MID,'--'),('Treat none','#D4D9DC',':')]:
    z=dca[dca['model']==m]
    ax2.plot(z['threshold']*100,z['net_benefit'],color=c,lw=1.1,ls=ls,label=m)
ax2.axvspan(8,15,color=TEAL,alpha=.05)
ax2.set_xlim(2,30); ax2.set_ylim(-.025,.09); ax2.set_xlabel('1-year risk threshold (%)'); ax2.set_ylabel('Net benefit')
ax2.set_title('Decision-curve analysis',loc='left',fontweight='bold',fontsize=12,color=NAVY); ax2.legend(frameon=False,ncol=2,fontsize=8,loc='upper right'); clean(ax2)
fig.text(0.08,0.07,'Exploratory fixed-horizon analysis among women with known 1-year status.',fontsize=8.8,color=SLATE)
save(fig,'10_reclassification_and_decision_curve')

print(f'Wrote files to {OUT}')
for p in sorted(OUT.glob('*.png')):
    print(p.name)
