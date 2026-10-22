from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.ndimage import label as nd_label
from scipy.signal import savgol_filter
from scipy.stats import chi2, ttest_ind
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
from statsmodels.duration.hazard_regression import PHReg
import patsy
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
TABLES = RESULTS / 'tables'
FIGDIR = RESULTS / 'nature_figures'
PANELDIR = FIGDIR / 'panels'
for p in (TABLES, FIGDIR, PANELDIR):
    p.mkdir(parents=True, exist_ok=True)

# Local modules
spec = importlib.util.spec_from_file_location('run_analysis', ROOT / 'src' / 'run_analysis.py')
ra = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ra)
from nature_theme import (
    apply_nature_style, NAVY, TEAL, TEAL_LIGHT, ORANGE, ORANGE_LIGHT,
    SLATE, GRAY, LIGHT_GRAY, PALE_BLUE, PALE_ORANGE, WHITE, BLACK,
    MODEL_COLORS,
)

SEED = 20261022
apply_nature_style()


def save(fig, name: str, dpi=420):
    png = PANELDIR / f'{name}.png'
    svg = PANELDIR / f'{name}.svg'
    fig.savefig(png, dpi=dpi, bbox_inches='tight', pad_inches=0.04)
    fig.savefig(svg, bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
    return png


def load_data():
    raw_path = ROOT / 'data' / 'BC_cardiotox_clinical_and_functional_variables.csv'
    raw = pd.read_csv(raw_path, sep=';', decimal=',')
    tdi_cols = [c for c in raw.columns if c.startswith('t ')]
    tdi = raw.loc[raw[tdi_cols].notna().all(axis=1)].copy().reset_index(drop=True)
    waves = tdi[tdi_cols].to_numpy(float)
    smooth = savgol_filter(waves, 31, 3, axis=1)
    feats = pd.DataFrame([ra.extract_tdi_features(w) for w in waves])
    analysis = pd.concat([tdi.drop(columns=tdi_cols).reset_index(drop=True), feats], axis=1)
    analysis['fractional_shortening'] = (analysis['LVDd'] - analysis['LVSd']) / analysis['LVDd'] * 100
    analysis['MTI_raw'] = ra.make_mti(analysis)
    analysis['MTI'] = (analysis['MTI_raw'] - analysis['MTI_raw'].mean()) / analysis['MTI_raw'].std(ddof=0)
    pred = pd.read_csv(RESULTS / 'crossfitted_predictions.csv')
    return raw, tdi, analysis, waves, smooth, pred


def bh_adjust(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def harrell_c(time, event, risk):
    time = np.asarray(time, float); event = np.asarray(event, int); risk = np.asarray(risk, float)
    total = concord = ties = 0
    for i in range(len(time)):
        if event[i] != 1:
            continue
        js = np.where(time > time[i])[0]
        if len(js) == 0:
            continue
        total += len(js)
        concord += int(np.sum(risk[i] > risk[js]))
        ties += int(np.sum(risk[i] == risk[js]))
    return (concord + 0.5 * ties) / total if total else np.nan


def strat_boot_idx(y, rng):
    y = np.asarray(y, int)
    pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    return np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]


def compute_performance(pred, B=5000):
    y = pred['CTRCD'].to_numpy(int); time = pred['time'].to_numpy(float)
    models = ['MTI', 'LVEF', 'vEprime', 'MTI+LVEF']
    point = []
    for m in models:
        risk = pred[m].to_numpy(float)
        point.append({
            'model': m,
            'AUROC': roc_auc_score(y, risk),
            'AUPRC': average_precision_score(y, risk),
            'C_index': harrell_c(time, y, risk),
        })
    point = pd.DataFrame(point).set_index('model')
    rng = np.random.default_rng(SEED + 10)
    store = {m: {'AUROC': [], 'AUPRC': [], 'C_index': []} for m in models}
    diff = {'AUROC': [], 'C_index': []}
    for _ in range(B):
        idx = strat_boot_idx(y, rng)
        yy = y[idx]; tt = time[idx]
        for m in models:
            r = pred[m].to_numpy(float)[idx]
            store[m]['AUROC'].append(roc_auc_score(yy, r))
            store[m]['AUPRC'].append(average_precision_score(yy, r))
            store[m]['C_index'].append(harrell_c(tt, yy, r))
        diff['AUROC'].append(store['MTI+LVEF']['AUROC'][-1] - store['LVEF']['AUROC'][-1])
        diff['C_index'].append(store['MTI+LVEF']['C_index'][-1] - store['LVEF']['C_index'][-1])
    rows = []
    for m in models:
        row = {'model': m}
        for metric in ['AUROC','AUPRC','C_index']:
            vals = np.asarray(store[m][metric])
            row[metric] = float(point.loc[m, metric])
            row[f'{metric}_low'] = float(np.quantile(vals, .025))
            row[f'{metric}_high'] = float(np.quantile(vals, .975))
        rows.append(row)
    perf = pd.DataFrame(rows)
    perf.to_csv(TABLES / 'nature_performance_bootstrap.csv', index=False)
    diff_df = pd.DataFrame([
        {'comparison':'MTI+LVEF vs LVEF','metric':'AUROC', 'difference':float(point.loc['MTI+LVEF','AUROC']-point.loc['LVEF','AUROC']),
         'CI_low':float(np.quantile(diff['AUROC'],.025)), 'CI_high':float(np.quantile(diff['AUROC'],.975))},
        {'comparison':'MTI+LVEF vs LVEF','metric':'C_index','difference':float(point.loc['MTI+LVEF','C_index']-point.loc['LVEF','C_index']),
         'CI_low':float(np.quantile(diff['C_index'],.025)), 'CI_high':float(np.quantile(diff['C_index'],.975))},
    ])
    diff_df.to_csv(TABLES / 'paired_performance_differences.csv', index=False)
    return perf, diff_df


def functional_waveform_test(smooth, analysis, permutations=5000, boot=2500):
    y = analysis['CTRCD'].to_numpy(int)
    t = np.linspace(0, 1, smooth.shape[1])

    # Residualize the waveform on baseline LVEF without using outcome labels.
    lvef = analysis['LVEF'].to_numpy(float)
    z = (lvef - np.nanmean(lvef)) / np.nanstd(lvef)
    C = np.column_stack([np.ones(len(z)), z])
    beta = np.linalg.lstsq(C, smooth, rcond=None)[0]
    resid = smooth - C @ beta

    n = len(y); n1 = int(y.sum()); n0 = n - n1

    def welch_from_mask(A, mask):
        # mask is B x n with exactly n1 event positions per permutation.
        total = A.sum(axis=0); total2 = (A*A).sum(axis=0)
        s1 = mask @ A; ss1 = mask @ (A*A)
        s0 = total[None, :] - s1; ss0 = total2[None, :] - ss1
        m1 = s1 / n1; m0 = s0 / n0
        v1 = (ss1 - (s1*s1)/n1) / max(n1-1, 1)
        v0 = (ss0 - (s0*s0)/n0) / max(n0-1, 1)
        den = np.sqrt(np.maximum(v1/n1 + v0/n0, 1e-12))
        return (m1-m0)/den

    def welch_obs(A, labels):
        return welch_from_mask(A, labels[None, :].astype(float))[0]

    def clusters(tv, threshold=2.0):
        mask = np.abs(tv) >= threshold
        labs, nc = nd_label(mask.astype(int))
        out = []
        for k in range(1, nc+1):
            idx = np.where(labs == k)[0]
            out.append((idx, float(np.sum(np.abs(tv[idx])))))
        return out

    obs_t_raw = welch_obs(smooth, y)
    obs_t_adj = welch_obs(resid, y)
    raw_stat = float(np.mean(obs_t_raw**2))
    adj_stat = float(np.mean(obs_t_adj**2))
    obs_clusters = clusters(obs_t_adj)

    rng = np.random.default_rng(SEED + 20)
    # Generate fixed-size event masks. This exactly preserves the observed class count.
    mask = np.zeros((permutations, n), dtype=np.float32)
    for b in range(permutations):
        mask[b, rng.choice(n, n1, replace=False)] = 1.0

    raw_perm = np.empty(permutations, float)
    adj_perm = np.empty(permutations, float)
    max_cluster = np.zeros(permutations, float)
    batch = 250
    for lo_b in range(0, permutations, batch):
        hi_b = min(permutations, lo_b+batch)
        M = mask[lo_b:hi_b]
        tr = welch_from_mask(smooth, M)
        ta = welch_from_mask(resid, M)
        raw_perm[lo_b:hi_b] = np.mean(tr*tr, axis=1)
        adj_perm[lo_b:hi_b] = np.mean(ta*ta, axis=1)
        for j in range(hi_b-lo_b):
            cc = clusters(ta[j])
            max_cluster[lo_b+j] = max((mass for _, mass in cc), default=0.0)

    global_raw_p = (1 + np.sum(raw_perm >= raw_stat)) / (permutations + 1)
    global_adj_p = (1 + np.sum(adj_perm >= adj_stat)) / (permutations + 1)

    # Stratified bootstrap of the event minus non-event mean waveform difference.
    evt = smooth[y == 1]; ctl = smooth[y == 0]
    rng2 = np.random.default_rng(SEED + 21)
    diffs = np.empty((boot, smooth.shape[1]), dtype=np.float32)
    for b in range(boot):
        e = evt[rng2.integers(0, len(evt), len(evt))].mean(0)
        c = ctl[rng2.integers(0, len(ctl), len(ctl))].mean(0)
        diffs[b] = e - c
    delta = evt.mean(0) - ctl.mean(0)
    dlow = np.quantile(diffs, .025, axis=0); dhigh = np.quantile(diffs, .975, axis=0)

    cluster_id = np.zeros(len(t), int); cluster_p = np.ones(len(t), float)
    crows = []
    for j, (idx, mass) in enumerate(obs_clusters, start=1):
        p_cluster = (1 + np.sum(max_cluster >= mass)) / (permutations + 1)
        cluster_id[idx] = j; cluster_p[idx] = p_cluster
        crows.append({'cluster':j,'phase_start':float(t[idx[0]]),'phase_end':float(t[idx[-1]]),'mass':mass,'FWER_p':float(p_cluster),'mean_t':float(obs_t_adj[idx].mean())})
    pd.DataFrame(crows).to_csv(TABLES / 'phase_cluster_permutation.csv', index=False)
    pd.DataFrame({
        'phase':t, 'delta_velocity_event_minus_nonevent':delta,
        'delta_low':dlow,'delta_high':dhigh,
        't_unadjusted':obs_t_raw,'t_LVEF_adjusted':obs_t_adj,
        'cluster_id':cluster_id,'cluster_FWER_p':cluster_p,
    }).to_csv(TABLES / 'phase_functional_analysis.csv', index=False)
    global_df = pd.DataFrame([
        {'test':'Global waveform difference, unadjusted','statistic_mean_t2':raw_stat,'permutation_p':global_raw_p},
        {'test':'Global waveform difference, residualized for baseline LVEF','statistic_mean_t2':adj_stat,'permutation_p':global_adj_p},
    ])
    global_df.to_csv(TABLES / 'functional_global_permutation.csv', index=False)
    return t, delta, dlow, dhigh, obs_t_adj, crows, global_df

def risk_concentration(pred, B=5000):
    y = pred['CTRCD'].to_numpy(int); n=len(y); total=y.sum()
    models = ['MTI','LVEF','MTI+LVEF']
    fractions = np.linspace(.05, 1.0, 40)
    curve_rows=[]
    for m in models:
        order=np.argsort(-pred[m].to_numpy(float))
        for f in fractions:
            k=max(1,int(np.ceil(n*f)))
            ev=int(y[order[:k]].sum())
            curve_rows.append({'model':m,'fraction_selected':f,'event_capture':ev/total,'PPV':ev/k,'enrichment':(ev/k)/(total/n)})
    curve=pd.DataFrame(curve_rows)
    curve.to_csv(TABLES/'risk_concentration_curve.csv',index=False)

    rng=np.random.default_rng(SEED+30)
    summary=[]
    for f in [.20, 1/3, .50]:
        for m in models:
            vals=[]
            score=pred[m].to_numpy(float)
            for _ in range(B):
                idx=strat_boot_idx(y,rng); yy=y[idx]; ss=score[idx]; nn=len(idx); ee=yy.sum(); k=max(1,int(np.ceil(nn*f)))
                o=np.argsort(-ss); ev=yy[o[:k]].sum(); ppv=ev/k
                vals.append((ev/ee,ppv,ppv/(ee/nn)))
            vals=np.asarray(vals)
            point=curve[(curve.model==m)].iloc[(curve[curve.model==m].fraction_selected-f).abs().argmin()]
            # exact point at requested fraction
            k=max(1,int(np.ceil(n*f))); o=np.argsort(-score); ev=y[o[:k]].sum(); ppv=ev/k
            summary.append({'model':m,'fraction_selected':f,'n_selected':k,'events_captured':int(ev),
                'event_capture':ev/total,'event_capture_low':np.quantile(vals[:,0],.025),'event_capture_high':np.quantile(vals[:,0],.975),
                'PPV':ppv,'PPV_low':np.quantile(vals[:,1],.025),'PPV_high':np.quantile(vals[:,1],.975),
                'enrichment_vs_prevalence':ppv/(total/n),'enrichment_low':np.quantile(vals[:,2],.025),'enrichment_high':np.quantile(vals[:,2],.975)})
    summ=pd.DataFrame(summary); summ.to_csv(TABLES/'risk_concentration_summary.csv',index=False)
    return curve,summ


def spline_dose_response(analysis):
    d=analysis[['time','CTRCD','MTI','LVEF']].dropna().copy()
    d['LVEF_z']=(d.LVEF-d.LVEF.mean())/d.LVEF.std(ddof=0)
    lin=PHReg(d.time,np.column_stack([d.MTI,d.LVEF_z]),status=d.CTRCD,ties='breslow').fit(disp=0)
    dm=patsy.dmatrix('cr(MTI, df=3)-1',d,return_type='dataframe')
    info=dm.design_info
    X=np.column_stack([np.asarray(dm),d.LVEF_z])
    spl=PHReg(d.time,X,status=d.CTRCD,ties='breslow').fit(disp=0)
    LR=2*(spl.llf-lin.llf); p_nonlin=float(chi2.sf(LR,2))
    grid=np.linspace(d.MTI.quantile(.02),d.MTI.quantile(.98),180)
    Bg=np.asarray(patsy.build_design_matrices([info],{'MTI':grid})[0])
    B0=np.asarray(patsy.build_design_matrices([info],{'MTI':[0.0]})[0])[0]
    Xg=np.column_stack([Bg,np.zeros(len(grid))]); x0=np.r_[B0,0.0]
    delta=Xg-x0; beta=np.asarray(spl.params); cov=np.asarray(spl.cov_params())
    loghr=delta@beta; se=np.sqrt(np.einsum('ij,jk,ik->i',delta,cov,delta))
    out=pd.DataFrame({'MTI':grid,'HR_vs_MTI0':np.exp(loghr),'CI_low':np.exp(loghr-1.96*se),'CI_high':np.exp(loghr+1.96*se)})
    out['nonlinearity_p']=p_nonlin
    out.to_csv(TABLES/'mti_spline_dose_response.csv',index=False)

    # Continuous interaction with LVEF.
    mz=(d.MTI-d.MTI.mean())/d.MTI.std(ddof=0); lz=d.LVEF_z
    inter=PHReg(d.time,np.column_stack([mz,lz,mz*lz]),status=d.CTRCD,ties='breslow').fit(disp=0)
    intertab=pd.DataFrame({'term':['MTI','LVEF','MTI_x_LVEF'],'logHR':inter.params,'SE':inter.bse,'p_value':inter.pvalues,'HR':np.exp(inter.params)})
    intertab.to_csv(TABLES/'mti_lvef_interaction.csv',index=False)
    return out,p_nonlin,intertab


def fpca_analysis(smooth,analysis,ncomp=5):
    X=smooth-smooth.mean(0)
    pca=PCA(n_components=ncomp,svd_solver='randomized',random_state=SEED).fit(X)
    scores=pca.transform(X)
    rows=[]; lz=(analysis.LVEF-analysis.LVEF.mean())/analysis.LVEF.std(ddof=0)
    for k in range(ncomp):
        z=(scores[:,k]-scores[:,k].mean())/scores[:,k].std(ddof=0)
        fit=PHReg(analysis.time,np.column_stack([z,lz]),status=analysis.CTRCD,ties='breslow').fit(disp=0)
        b,se,p=fit.params[0],fit.bse[0],fit.pvalues[0]
        rows.append({'PC':k+1,'variance_explained':pca.explained_variance_ratio_[k],'HR_per_1SD':np.exp(b),'CI_low':np.exp(b-1.96*se),'CI_high':np.exp(b+1.96*se),'p_value':p,'corr_with_MTI':np.corrcoef(z,analysis.MTI)[0,1]})
    tab=pd.DataFrame(rows); tab['BH_q']=bh_adjust(tab.p_value)
    tab.to_csv(TABLES/'functional_PCA_survival.csv',index=False)
    np.savez(RESULTS/'functional_pca.npz',mean=pca.mean_,components=pca.components_,explained=pca.explained_variance_ratio_,scores=scores)
    return pca,scores,tab


def tdi_selection_audit(raw):
    tdi_cols=[c for c in raw.columns if c.startswith('t ')]
    avail=raw[tdi_cols].notna().all(axis=1)
    rows=[]
    cont=['age','LVEF','heart_rate','time']
    binary=['CTRCD','AC','antiHER2','HTA','DL','DM','RTprev']
    for c in cont:
        a=raw.loc[avail,c].dropna(); b=raw.loc[~avail,c].dropna(); sp=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        rows.append({'variable':c,'SMD':(a.mean()-b.mean())/sp if sp else np.nan,'TDI':a.mean(),'No_TDI':b.mean()})
    for c in binary:
        a=raw.loc[avail,c].dropna(); b=raw.loc[~avail,c].dropna(); p1=a.mean(); p0=b.mean(); pb=(p1+p0)/2; den=np.sqrt(pb*(1-pb)) if 0<pb<1 else np.nan
        rows.append({'variable':c,'SMD':(p1-p0)/den if den else np.nan,'TDI':p1,'No_TDI':p0})
    tab=pd.DataFrame(rows); tab.to_csv(TABLES/'tdi_selection_audit.csv',index=False); return tab


def kaplan_meier(time,event):
    order=np.argsort(time); time=np.asarray(time)[order]; event=np.asarray(event)[order]
    unique=np.unique(time[event==1]); xs=[0.0]; ys=[1.0]; surv=1.0
    for ti in unique:
        nr=np.sum(time>=ti); de=np.sum((time==ti)&(event==1)); xs.extend([ti,ti]); ys.extend([ys[-1],surv*(1-de/nr)]); surv=ys[-1]
    xs.append(float(time.max())); ys.append(surv)
    return np.asarray(xs),np.asarray(ys)


def make_panels(raw,analysis,smooth,pred,perf,diffdf,phase,risk,spline,pca_data,sel):
    t=np.linspace(0,1,smooth.shape[1]); y=analysis.CTRCD.to_numpy(int)
    panel_paths=[]

    # Panel: cohort flow diagram, generated entirely in matplotlib.
    fig=plt.figure(figsize=(7.2,2.5)); ax=fig.add_axes([0,0,1,1]); ax.set_axis_off()
    boxes=[(0.03,0.23,0.24,0.56,'Public BC_cardiotox\nclinical cohort','n = 531\n54 CTRCD events'),
           (0.38,0.23,0.24,0.56,'Pretreatment TDI\navailable','n = 270\n27 CTRCD events'),
           (0.73,0.23,0.24,0.56,'MTI analysis','4 components\nTime-to-CTRCD')]
    for i,(x,y0,w,h,title,sub) in enumerate(boxes):
        fc=PALE_BLUE if i!=1 else '#F3F7F8'
        ax.add_patch(FancyBboxPatch((x,y0),w,h,boxstyle='round,pad=0.014,rounding_size=0.02',fc=fc,ec=LIGHT_GRAY,lw=1.1))
        ax.text(x+w/2,y0+h*.68,title,ha='center',va='center',fontsize=11,weight='bold',color=NAVY)
        ax.text(x+w/2,y0+h*.33,sub,ha='center',va='center',fontsize=9.4,color=SLATE,linespacing=1.4)
    for x1,x2 in [(0.27,0.38),(0.62,0.73)]:
        ax.add_patch(FancyArrowPatch((x1+0.01,.51),(x2-0.01,.51),arrowstyle='-|>',mutation_scale=15,lw=1.4,color=TEAL))
    ax.text(.5,.94,'Study cohort and analysis flow',ha='center',va='center',fontsize=12,weight='bold',color=BLACK)
    panel_paths.append(save(fig,'A_cohort_flow'))

    # Group waveforms.
    fig=plt.figure(figsize=(7.2,4.1)); ax=fig.add_axes([.12,.15,.84,.77])
    for val,label,color in [(0,'No CTRCD',TEAL),(1,'CTRCD',ORANGE)]:
        W=smooth[y==val]; med=np.median(W,axis=0); q1=np.quantile(W,.25,axis=0); q3=np.quantile(W,.75,axis=0)
        ax.plot(t,med,color=color,lw=2.2,label=label); ax.fill_between(t,q1,q3,color=color,alpha=.13,lw=0)
    ez=float(analysis.e_prime_time.median()); sz=float(analysis.systolic_zero_time.median())
    ax.axvline(sz,color=GRAY,lw=1,ls='--'); ax.axvspan(max(.45,ez-.055),min(.78,ez+.055),color=PALE_ORANGE,alpha=.8,zorder=0)
    ax.axhline(0,color=GRAY,lw=.8); ax.set(xlabel='Normalized cardiac cycle',ylabel='Myocardial velocity (cm/s)',title='Outcome-associated pretreatment TDI waveforms')
    ax.legend(loc='upper right'); ax.spines[['top','right']].set_visible(False)
    ax.text(.18,ax.get_ylim()[1]*.88,'S',color=NAVY,weight='bold'); ax.text(ez,ax.get_ylim()[0]*.92,"E′",color=NAVY,weight='bold',ha='center')
    panel_paths.append(save(fig,'B_group_waveforms'))

    # Actual-data waveform atlas sorted by MTI.
    fig=plt.figure(figsize=(7.2,5.0)); ax=fig.add_axes([.14,.12,.78,.80])
    order=np.argsort(analysis.MTI.to_numpy())
    Z=smooth[order]
    # row-normalize only for visual shape atlas while preserving phase morphology
    Z=(Z-Z.mean(1,keepdims=True))/(Z.std(1,keepdims=True)+1e-9)
    cmap=LinearSegmentedColormap.from_list('nature_div',[TEAL_LIGHT,WHITE,ORANGE_LIGHT],N=256)
    im=ax.imshow(Z,aspect='auto',interpolation='nearest',extent=[0,1,len(Z),0],cmap=cmap,vmin=-2.5,vmax=2.5)
    ax.set(xlabel='Normalized cardiac cycle',ylabel='Patients ordered by increasing MTI',title='Waveform atlas ordered by Mechanical Transition Instability')
    ax.axvline(float(analysis.systolic_zero_time.median()),color=NAVY,lw=1,ls='--')
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.02); cb.set_label('Within-patient standardized velocity')
    ax.spines[['top','right']].set_visible(False)
    panel_paths.append(save(fig,'C_waveform_atlas'))

    # Phase-resolved waveform difference with rigorous functional test annotation.
    tt,delta,lo,hi,tadj,crows,gdf=phase
    fig=plt.figure(figsize=(7.2,4.2)); ax=fig.add_axes([.12,.16,.84,.74])
    ax.fill_between(tt,lo,hi,color=ORANGE,alpha=.15,lw=0); ax.plot(tt,delta,color=ORANGE,lw=2.2)
    ax.axhline(0,color=GRAY,lw=.9,ls='--'); ax.axvspan(.50,.60,color=PALE_ORANGE,alpha=.65,zorder=0)
    ax.set(xlabel='Normalized cardiac cycle',ylabel='Δ myocardial velocity (CTRCD − no CTRCD), cm/s',title='Phase-resolved difference in pretreatment myocardial motion')
    ax.spines[['top','right']].set_visible(False)
    pr=float(gdf.iloc[0].permutation_p); pa=float(gdf.iloc[1].permutation_p)
    ax.text(.02,.04,f'Global functional permutation: P={pr:.3f}\nAfter residualizing LVEF: P={pa:.3f}',transform=ax.transAxes,fontsize=8.4,color=SLATE,va='bottom')
    if crows:
        best=min(crows,key=lambda x:x['FWER_p'])
        ax.text(.98,.04,f"Best contiguous cluster: phase {best['phase_start']:.3f}–{best['phase_end']:.3f}\ncluster-FWER P={best['FWER_p']:.3f}",transform=ax.transAxes,ha='right',va='bottom',fontsize=8.4,color=SLATE)
    panel_paths.append(save(fig,'D_phase_difference'))

    # MTI four-component construction from actual median waveform.
    fig=plt.figure(figsize=(7.2,4.2)); ax=fig.add_axes([.10,.15,.86,.70])
    fig.suptitle('MTI construction from four interpretable transition features',y=.97,fontsize=11.2,fontweight='bold',color=BLACK)
    med=np.median(smooth,axis=0); i0=int(np.argmin(abs(t-analysis.systolic_zero_time.median()))); ie=int(np.argmin(abs(t-analysis.e_prime_time.median())))
    ax.plot(t,med,color=TEAL,lw=2.3); ax.axhline(0,color=GRAY,lw=.8); ax.axvline(t[i0],color=NAVY,ls='--',lw=1)
    ax.scatter([t[ie]],[med[ie]],s=35,color=ORANGE,zorder=4)
    ann_kw=dict(arrowprops=dict(arrowstyle='->',color=SLATE,lw=.9),fontsize=8.2,fontweight='bold')
    ax.annotate('1  Transition slope deficit',xy=(.36,np.interp(.36,t,med)),xycoords='data',xytext=(.04,.96),textcoords='axes fraction',color=NAVY,ha='left',va='top',**ann_kw)
    ax.annotate('2  Transition curvature deficit',xy=(t[i0]-.025,1),xycoords='data',xytext=(.38,.96),textcoords='axes fraction',color=ORANGE,ha='left',va='top',**ann_kw)
    ax.annotate('3  Early-relaxation latency',xy=((t[i0]+t[ie])/2,-1),xycoords='data',xytext=(.62,.78),textcoords='axes fraction',color=NAVY,ha='left',va='top',**ann_kw)
    ax.annotate("4  E′ attenuation",xy=(t[ie],med[ie]),xycoords='data',xytext=(.72,.22),textcoords='axes fraction',color=ORANGE,ha='left',va='top',**ann_kw)
    ax.text(t[i0],.5,'Systolic zero',rotation=90,ha='right',va='bottom',fontsize=7.6,color=SLATE)
    ax.text(t[ie],med[ie]-1.1,"E′ peak",ha='center',va='top',fontsize=7.6,color=SLATE)
    ax.set(xlabel='Normalized cardiac cycle',ylabel='Myocardial velocity (cm/s)')
    ax.spines[['top','right']].set_visible(False)
    panel_paths.append(save(fig,'E_MTI_components'))

    # ROC.
    fig=plt.figure(figsize=(5.4,4.7)); ax=fig.add_axes([.14,.14,.82,.78])
    labels={'MTI':'MTI','LVEF':'LVEF','vEprime':"Published vE′",'MTI+LVEF':'MTI + LVEF'}
    for m in ['MTI','LVEF','vEprime','MTI+LVEF']:
        fpr,tpr,_=roc_curve(pred.CTRCD,pred[m]); auc=roc_auc_score(pred.CTRCD,pred[m]); ax.plot(fpr,tpr,color=MODEL_COLORS[m],label=f'{labels[m]}  {auc:.3f}',lw=2)
    ax.plot([0,1],[0,1],color=GRAY,ls='--',lw=1); ax.set(xlabel='False-positive rate',ylabel='True-positive rate',title='Cross-fitted discrimination'); ax.legend(title='AUROC',loc='lower right'); ax.spines[['top','right']].set_visible(False)
    panel_paths.append(save(fig,'F_ROC'))

    # Performance forest, AUROC only, with exact bootstrap CIs.
    fig=plt.figure(figsize=(5.6,3.8)); ax=fig.add_axes([.28,.17,.68,.73])
    pp=perf.set_index('model').loc[['MTI','LVEF','vEprime','MTI+LVEF']].reset_index(); yy=np.arange(len(pp))[::-1]
    for ypos,row in zip(yy,pp.itertuples()):
        ax.errorbar(row.AUROC,ypos,xerr=[[row.AUROC-row.AUROC_low],[row.AUROC_high-row.AUROC]],fmt='o',ms=6,capsize=3,color=MODEL_COLORS[row.model],ecolor=MODEL_COLORS[row.model])
        ax.text(row.AUROC_high+.015,ypos,f'{row.AUROC:.3f}',va='center',fontsize=8.5,color=BLACK)
    ax.set_yticks(yy,['MTI','LVEF',"Published vE′",'MTI + LVEF']); ax.set_xlim(.45,.88); ax.set_xlabel('AUROC (95% bootstrap CI)'); ax.set_title('Cross-fitted model comparison'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
    panel_paths.append(save(fig,'G_performance_forest'))

    # Kaplan-Meier MTI tertiles.
    q1,q2=analysis.MTI.quantile([1/3,2/3]); groups=pd.cut(analysis.MTI,[-np.inf,q1,q2,np.inf],labels=['Low','Middle','High'])
    fig=plt.figure(figsize=(5.5,4.2)); ax=fig.add_axes([.14,.15,.82,.77])
    for g,c in [('Low',TEAL),('Middle',GRAY),('High',ORANGE)]:
        m=groups.astype(str).eq(g); xk,yk=kaplan_meier(analysis.loc[m,'time'],analysis.loc[m,'CTRCD']); e=int(analysis.loc[m,'CTRCD'].sum()); n=int(m.sum()); ax.step(xk,yk,where='post',color=c,lw=2,label=f'{g} MTI  ({e}/{n} events)')
    ax.set_ylim(.70,1.01); ax.set(xlabel='Days from treatment initiation',ylabel='Freedom from CTRCD',title='Risk gradient across pretreatment MTI tertiles'); ax.legend(loc='lower left'); ax.spines[['top','right']].set_visible(False)
    panel_paths.append(save(fig,'H_KM_tertiles'))

    # Cox forest.
    cox=pd.read_csv(TABLES/'cox_models.csv'); pres=pd.read_csv(TABLES/'cox_preserved_lvef_sensitivity.csv')
    rows=[]
    for model,label in [('MTI unadjusted','Unadjusted'),('MTI + LVEF','Adjusted for LVEF'),('MTI + age + LVEF','Adjusted for age + LVEF'),('MTI + age + LVEF + treatment','Treatment-adjusted')]:
        r=cox[(cox.model==model)&(cox.term=='MTI')].iloc[0]; rows.append((label,r.HR_per_1SD,r.CI_low,r.CI_high,r.p_value))
    r=pres[pres.term=='MTI'].iloc[0]; rows.append(('Baseline LVEF ≥50%',r.HR_per_1SD,r.CI_low,r.CI_high,r.p_value))
    fig=plt.figure(figsize=(5.8,4.1)); ax=fig.add_axes([.36,.16,.60,.76]); yy=np.arange(len(rows))[::-1]
    for ypos,(lab,hr,lo_,hi_,p_) in zip(yy,rows):
        ax.errorbar(hr,ypos,xerr=[[hr-lo_],[hi_-hr]],fmt='o',ms=6,capsize=3,color=TEAL,ecolor=SLATE)
        ax.text(hi_*1.04,ypos,f'{hr:.2f} ({lo_:.2f}–{hi_:.2f})',va='center',fontsize=8)
    ax.axvline(1,color=GRAY,ls='--',lw=1); ax.set_xscale('log'); ax.set_yticks(yy,[r[0] for r in rows]); ax.set_xlabel('Hazard ratio per 1-SD higher MTI'); ax.set_title('Association with time to CTRCD'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
    panel_paths.append(save(fig,'I_Cox_forest'))

    # Risk concentration / surveillance triage curve.
    curve,summ=risk
    fig=plt.figure(figsize=(5.6,4.4)); ax=fig.add_axes([.14,.15,.82,.77])
    for m,label_ in [('MTI','MTI'),('LVEF','LVEF'),('MTI+LVEF','MTI + LVEF')]:
        d=curve[curve.model==m]; ax.plot(d.fraction_selected*100,d.event_capture*100,color=MODEL_COLORS[m],label=label_,lw=2.2)
    ax.plot([0,100],[0,100],color=GRAY,ls='--',lw=1,label='Random selection')
    ax.set(xlabel='Cohort selected for intensified surveillance (%)',ylabel='Future CTRCD events captured (%)',title='Risk concentration for surveillance triage'); ax.legend(loc='lower right'); ax.spines[['top','right']].set_visible(False)
    s20=summ[(summ.model=='MTI+LVEF')&(np.isclose(summ.fraction_selected,.20))].iloc[0]
    yerr=np.array([[100*(s20.event_capture-s20.event_capture_low)],[100*(s20.event_capture_high-s20.event_capture)]])
    ax.errorbar([20],[100*s20.event_capture],yerr=yerr,fmt='o',ms=5,capsize=3,color=TEAL,ecolor=TEAL,zorder=5)
    ax.annotate(f"Top 20% captures {100*s20.event_capture:.0f}% of events\n95% CI {100*s20.event_capture_low:.0f}–{100*s20.event_capture_high:.0f}%\nPPV {100*s20.PPV:.0f}% vs 10% prevalence",xy=(20,100*s20.event_capture),xytext=(35,27),arrowprops=dict(arrowstyle='->',color=SLATE),fontsize=8.0,color=SLATE)
    panel_paths.append(save(fig,'J_risk_concentration'))

    # Dose response spline.
    spline_df,p_nonlin,intertab=spline
    fig=plt.figure(figsize=(5.5,4.2)); ax=fig.add_axes([.14,.15,.82,.77])
    ax.fill_between(spline_df.MTI,spline_df.CI_low,spline_df.CI_high,color=TEAL,alpha=.14,lw=0); ax.plot(spline_df.MTI,spline_df.HR_vs_MTI0,color=TEAL,lw=2.3); ax.axhline(1,color=GRAY,ls='--',lw=1); ax.axvline(0,color=LIGHT_GRAY,lw=1)
    ax.set_yscale('log'); ax.set(xlabel='MTI (SD units)',ylabel='Relative hazard of CTRCD',title='Continuous MTI risk gradient after LVEF adjustment'); ax.spines[['top','right']].set_visible(False); ax.text(.03,.95,f'Nonlinearity P={p_nonlin:.3f}',transform=ax.transAxes,ha='left',va='top',fontsize=8.4,color=SLATE)
    panel_paths.append(save(fig,'K_spline_dose_response'))

    # Joint risk surface from continuous Cox model with interaction (descriptive association).
    mz=(analysis.MTI-analysis.MTI.mean())/analysis.MTI.std(ddof=0); lz=(analysis.LVEF-analysis.LVEF.mean())/analysis.LVEF.std(ddof=0)
    fit=PHReg(analysis.time,np.column_stack([mz,lz,mz*lz]),status=analysis.CTRCD,ties='breslow').fit(disp=0); beta=np.asarray(fit.params)
    gx=np.linspace(np.quantile(analysis.MTI,.03),np.quantile(analysis.MTI,.97),90); gy=np.linspace(np.quantile(analysis.LVEF,.03),np.quantile(analysis.LVEF,.97),90)
    XX,YY=np.meshgrid(gx,gy); MX=(XX-analysis.MTI.mean())/analysis.MTI.std(ddof=0); LY=(YY-analysis.LVEF.mean())/analysis.LVEF.std(ddof=0); lp=beta[0]*MX+beta[1]*LY+beta[2]*MX*LY; HR=np.exp(lp-np.nanmedian(lp))
    fig=plt.figure(figsize=(5.5,4.5)); ax=fig.add_axes([.15,.14,.78,.78]); levels=np.geomspace(max(.25,np.nanmin(HR)),min(4,np.nanmax(HR)),16); cf=ax.contourf(XX,YY,HR,levels=levels,cmap=LinearSegmentedColormap.from_list('risk',[PALE_BLUE,WHITE,PALE_ORANGE,ORANGE_LIGHT]),extend='both')
    ax.scatter(analysis.loc[y==0,'MTI'],analysis.loc[y==0,'LVEF'],s=8,c=SLATE,alpha=.30,lw=0); ax.scatter(analysis.loc[y==1,'MTI'],analysis.loc[y==1,'LVEF'],s=24,facecolor=ORANGE,edgecolor=WHITE,lw=.6,label='CTRCD')
    cb=fig.colorbar(cf,ax=ax,fraction=.045,pad=.02); cb.set_label('Relative hazard (descriptive Cox surface)'); ax.set(xlabel='MTI (SD units)',ylabel='Baseline LVEF (%)',title='Complementary continuous information from MTI and LVEF'); ax.legend(loc='upper right'); ax.spines[['top','right']].set_visible(False); pint=float(intertab.loc[intertab.term=='MTI_x_LVEF','p_value'].iloc[0]); ax.text(.03,.04,f'MTI×LVEF interaction P={pint:.3f}',transform=ax.transAxes,color=SLATE,fontsize=8.4)
    panel_paths.append(save(fig,'L_joint_risk_surface'))

    # FPCA orthogonal check.
    pca,scores,fpca_tab=pca_data
    fig=plt.figure(figsize=(5.6,4.2)); ax=fig.add_axes([.13,.15,.83,.76]); mean=pca.mean_; sdmode=np.sqrt(pca.explained_variance_[0])*pca.components_[0]
    ax.plot(t,mean,color=SLATE,lw=2,label='Mean waveform'); ax.plot(t,mean+sdmode,color=ORANGE,lw=2,label='Mean + 1 SD PC1'); ax.plot(t,mean-sdmode,color=TEAL,lw=2,label='Mean − 1 SD PC1'); ax.axhline(0,color=GRAY,lw=.8)
    r=fpca_tab.iloc[0]; ax.set(xlabel='Normalized cardiac cycle',ylabel='Myocardial velocity (cm/s)',title='Outcome-agnostic functional PCA identifies a related morphology axis'); ax.legend(loc='lower left'); ax.spines[['top','right']].set_visible(False); ax.text(.98,.95,f"PC1 explains {100*r.variance_explained:.1f}% variance\nLVEF-adjusted HR {r.HR_per_1SD:.2f}\nnominal P={r.p_value:.3f}; BH q={r.BH_q:.3f}",transform=ax.transAxes,ha='right',va='top',fontsize=8.2,color=SLATE)
    panel_paths.append(save(fig,'M_FPCA_validation'))

    # Selection audit love plot.
    labels_map={'age':'Age','LVEF':'LVEF','heart_rate':'Heart rate','time':'Follow-up time','CTRCD':'CTRCD','AC':'Anthracycline','antiHER2':'Anti-HER2','HTA':'Hypertension','DL':'Dyslipidemia','DM':'Diabetes','RTprev':'Prior radiotherapy'}
    ss=sel.copy(); ss['abs']=ss.SMD.abs(); ss=ss.sort_values('abs'); yy=np.arange(len(ss))
    fig=plt.figure(figsize=(5.5,4.6)); ax=fig.add_axes([.32,.13,.64,.80]); ax.scatter(ss.SMD,yy,s=38,color=[ORANGE if abs(x)>=.2 else TEAL for x in ss.SMD],zorder=3); ax.axvline(0,color=GRAY,lw=.8); ax.axvline(.2,color=LIGHT_GRAY,ls='--',lw=1); ax.axvline(-.2,color=LIGHT_GRAY,ls='--',lw=1); ax.set_yticks(yy,[labels_map.get(x,x) for x in ss.variable]); ax.set_xlabel('Standardized mean difference\n(TDI available − TDI unavailable)'); ax.set_title('TDI availability selection audit'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
    panel_paths.append(save(fig,'N_selection_audit'))

    return panel_paths


def stitch(paths, outname, ncols=2, bg='white'):
    imgs=[Image.open(p).convert('RGB') for p in paths]
    target_w=max(im.width for im in imgs)
    # normalize each panel to same width while preserving aspect ratio
    norm=[]
    for im in imgs:
        if im.width != target_w:
            h=round(im.height*target_w/im.width); im=im.resize((target_w,h),Image.Resampling.LANCZOS)
        norm.append(im)
    nrows=math.ceil(len(norm)/ncols); gap=36; row_heights=[]
    for r in range(nrows): row_heights.append(max(im.height for im in norm[r*ncols:(r+1)*ncols]))
    canvas=Image.new('RGB',(ncols*target_w+(ncols-1)*gap,sum(row_heights)+(nrows-1)*gap),bg)
    y0=0
    for r in range(nrows):
        x0=0
        for c in range(ncols):
            i=r*ncols+c
            if i>=len(norm): break
            canvas.paste(norm[i],(x0,y0)); x0+=target_w+gap
        y0+=row_heights[r]+gap
    out=FIGDIR/outname; canvas.save(out,dpi=(300,300),quality=95)
    tiff=out.with_suffix('.tiff'); canvas.save(tiff,dpi=(300,300),compression='tiff_lzw')
    return out,tiff


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--bootstrap',type=int,default=5000); ap.add_argument('--permutations',type=int,default=5000); args=ap.parse_args()
    raw,tdi,analysis,waves,smooth,pred=load_data()
    analysis.to_csv(RESULTS/'derived_analysis_dataset.csv',index=False)
    perf,diffdf=compute_performance(pred,B=args.bootstrap)
    phase=functional_waveform_test(smooth,analysis,permutations=args.permutations,boot=max(1000,args.bootstrap//2))
    risk=risk_concentration(pred,B=args.bootstrap)
    spline=spline_dose_response(analysis)
    pca_data=fpca_analysis(smooth,analysis)
    sel=tdi_selection_audit(raw)
    paths=make_panels(raw,analysis,smooth,pred,perf,diffdf,phase,risk,spline,pca_data,sel)

    # Main manuscript composites. Each constituent chart is generated independently above.
    main1_names=['A_cohort_flow','B_group_waveforms','E_MTI_components','C_waveform_atlas','D_phase_difference']
    main2_names=['F_ROC','G_performance_forest','H_KM_tertiles','I_Cox_forest','J_risk_concentration','K_spline_dose_response']
    extended_names=['L_joint_risk_surface','M_FPCA_validation','N_selection_audit']
    pathmap={p.stem:p for p in paths}
    stitch([pathmap[x] for x in main1_names], 'Figure1_Nature_MTI_mechanism.png', ncols=2)
    stitch([pathmap[x] for x in main2_names], 'Figure2_Nature_MTI_clinical.png', ncols=2)
    stitch([pathmap[x] for x in extended_names], 'Extended_Figure_Nature_robustness.png', ncols=2)

    manifest={'primary_figures':['Figure1_Nature_MTI_mechanism.png','Figure2_Nature_MTI_clinical.png'],
              'extended_figure':'Extended_Figure_Nature_robustness.png','panel_count':len(paths),
              'bootstrap_replicates':args.bootstrap,'permutations':args.permutations}
    (FIGDIR/'figure_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))

if __name__=='__main__':
    main()
