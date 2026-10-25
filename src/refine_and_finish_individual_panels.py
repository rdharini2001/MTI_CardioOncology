from __future__ import annotations
from pathlib import Path
import importlib.util
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, FancyArrowPatch
from scipy.signal import savgol_filter
from statsmodels.duration.survfunc import SurvfuncRight

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'/'individual_cardiomorph_code'; OUT.mkdir(parents=True,exist_ok=True); T=ROOT/'results'/'tables'
NAVY='#183B56'; DEEP='#0B5366'; TEAL='#0F7C82'; CYAN='#65B8C2'; ORANGE='#D95D2B'; SLATE='#6F8191'; MID='#A9B5BE'; HEADER='#ECEEEE'; PALE_TEAL='#E7F2F3'; PALE_ORANGE='#FAEEE8'; BLACK='#202A33'
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'axes.linewidth':1,'figure.facecolor':'white','savefig.facecolor':'white'})

def setup(title,figsize=(9.2,6.9)):
 fig=plt.figure(figsize=figsize); fig.patches.extend([Rectangle((.02,.025),.96,.95,transform=fig.transFigure,fill=False,edgecolor='#C4C8CA',lw=1.4,ls=(0,(5,4)),zorder=-20),Rectangle((.02,.845),.96,.13,transform=fig.transFigure,fc=HEADER,ec='none',zorder=-25)]); fig.text(.045,.91,title,fontsize=18,fontweight='bold',color=BLACK,va='center'); ax=fig.add_axes([.10,.12,.84,.70]); return fig,ax

def clean(ax):
 ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['left'].set_color('#536675'); ax.spines['bottom'].set_color('#536675'); ax.tick_params(length=4,width=.9); ax.grid(False)

def save(fig,stem):
 for ext,dpi in [('png',300),('svg',None),('tiff',300)]:
  kw={'bbox_inches':'tight','pad_inches':.05};
  if dpi: kw['dpi']=dpi
  if ext=='tiff': kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(OUT/f'{stem}.{ext}',**kw)
 plt.close(fig)

def person(ax,x,y,scale,color):
 ax.add_patch(Circle((x,y+.055*scale),.024*scale,fc=color,ec='none'))
 ax.add_patch(plt.Polygon([[x,y+.025*scale],[x-.035*scale,y-.045*scale],[x+.035*scale,y-.045*scale]],closed=True,fc=color,ec='none'))
 ax.plot([x-.012*scale,x-.012*scale],[y-.045*scale,y-.095*scale],color=color,lw=2.4*scale,solid_capstyle='round'); ax.plot([x+.012*scale,x+.012*scale],[y-.045*scale,y-.095*scale],color=color,lw=2.4*scale,solid_capstyle='round')

spec=importlib.util.spec_from_file_location('core',ROOT/'src'/'run_analysis.py'); core=importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
raw=pd.read_csv(ROOT/'data'/'BC_cardiotox_clinical_and_functional_variables.csv',sep=';',decimal=','); tc=[c for c in raw if c.startswith('t ')]; tdi=raw.loc[raw[tc].notna().all(1)].copy().reset_index(drop=True); waves=tdi[tc].to_numpy(float); smooth=savgol_filter(waves,31,3,axis=1); f=pd.DataFrame([core.extract_tdi_features(w) for w in waves]); a=pd.concat([tdi.drop(columns=tc),f],axis=1); a['MTI_raw']=core.make_mti(a); a['MTI']=(a.MTI_raw-a.MTI_raw.mean())/a.MTI_raw.std(ddof=0); t=np.linspace(0,1,waves.shape[1]); y=a.CTRCD.to_numpy(int); evt=smooth[y==1]; ctl=smooth[y==0]
pred=pd.read_csv(ROOT/'results'/'crossfitted_predictions.csv'); cox=pd.read_csv(T/'cox_models.csv'); pres=pd.read_csv(T/'cox_preserved_lvef_sensitivity.csv'); reclass=pd.read_csv(T/'one_year_reclassification_nri_idi.csv'); dca=pd.read_csv(T/'one_year_decision_curve.csv')

# 01
fig,ax=setup('BC_cardiotox cohort and analytic sample',(10,5.6)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
boxes=[(.02,.36,.27,.45,PALE_TEAL,TEAL,'531','women with HER2-positive\nbreast cancer','54 CTRCD events in source cohort'),(.365,.36,.27,.45,PALE_TEAL,TEAL,'270','with complete pretreatment\nTDI waveforms','analysis cohort'),(.71,.36,.27,.45,PALE_ORANGE,ORANGE,'27','subsequent CTRCD events','within TDI analysis cohort')]
for j,(x0,y0,w,h,fc,ec,n,main,sub) in enumerate(boxes):
 ax.add_patch(FancyBboxPatch((x0,y0),w,h,boxstyle='round,pad=.018,rounding_size=.03',fc=fc,ec=ec,lw=1.3))
 if j==0:
  for k,xx in enumerate(np.linspace(x0+.06,x0+w-.06,4)): person(ax,xx,y0+.31,.72,TEAL if k==1 else '#87A7B3')
 elif j==1:
  for k,xx in enumerate(np.linspace(x0+.085,x0+w-.085,3)): person(ax,xx,y0+.31,.76,TEAL if k==1 else '#87A7B3')
 else: person(ax,x0+w/2,y0+.31,.92,ORANGE)
 ax.text(x0+w/2,y0+.205,n,ha='center',va='center',fontsize=29,fontweight='bold',color=ec); ax.text(x0+w/2,y0+.11,main,ha='center',va='center',fontsize=11.5,color=BLACK,linespacing=1.25); ax.text(x0+w/2,y0+.035,sub,ha='center',va='bottom',fontsize=8.6,color=SLATE)
for x1,x2 in [(.29,.365),(.635,.71)]: ax.add_patch(FancyArrowPatch((x1+.008,.585),(x2-.008,.585),arrowstyle='-|>',mutation_scale=18,color=SLATE,lw=1.8))
ax.add_patch(FancyBboxPatch((.06,.07),.88,.18,boxstyle='round,pad=.012,rounding_size=.025',fc='#F2F4F5',ec='#C7CCCF',lw=1)); ax.text(.5,.18,'Primary endpoint: cancer therapy-related cardiac dysfunction (CTRCD)',ha='center',va='center',fontsize=11.2,fontweight='bold',color=BLACK); ax.text(.5,.11,'Absolute LVEF decrease ≥10 percentage points to <50%',ha='center',va='center',fontsize=10.5,color=NAVY); save(fig,'01_cohort_flow')

# 03 actual exemplar
fig,ax=setup('MTI construction from interpretable transition features'); mean=smooth.mean(0); dist=np.sqrt(((smooth-mean)**2).mean(1)); cand=(a.MTI.abs()<.5)&(a.s_peak_velocity>5)&(a.e_prime_velocity<-3)&a.systolic_zero_time.between(.35,.55)&a.e_prime_time.between(.50,.72); ids=np.where(cand)[0]; score=dist[ids]+.8*np.abs(a.MTI.iloc[ids].to_numpy()); idx=int(ids[np.argmin(score)]); w=smooth[idx]; r=f.iloc[idx]; ax.plot(t,w,color=DEEP,lw=2.6); ax.axhline(0,color=MID,lw=1,ls='--'); z=float(r.systolic_zero_time); e=float(r.e_prime_time); s=float(r.s_peak_time); ax.axvline(z,color=SLATE,lw=1,ls='--'); ax.axvspan(z,e,color=ORANGE,alpha=.055)
mid=(s+z)/2; ax.annotate('1  Transition slope deficit\nslower systolic deceleration',xy=(mid,np.interp(mid,t,w)),xytext=(.05,.87),textcoords='axes fraction',arrowprops=dict(arrowstyle='-|>',color=TEAL,lw=1.4),bbox=dict(boxstyle='round,pad=.4',fc=PALE_TEAL,ec=TEAL),fontsize=9.5); ax.annotate('2  Transition curvature deficit\nblunted transition bend',xy=(z,np.interp(z,t,w)),xytext=(.53,.80),textcoords='axes fraction',arrowprops=dict(arrowstyle='-|>',color=ORANGE,lw=1.4),bbox=dict(boxstyle='round,pad=.4',fc=PALE_ORANGE,ec=ORANGE),fontsize=9.5); by=.65; ax.annotate('',xy=(e,by),xytext=(z,by),arrowprops=dict(arrowstyle='<->',color=TEAL,lw=1.5)); ax.annotate('3  Early-relaxation latency\nzero crossing → E′',xy=((z+e)/2,by),xytext=(.60,.55),textcoords='axes fraction',arrowprops=dict(arrowstyle='-|>',color=TEAL,lw=1.2),bbox=dict(boxstyle='round,pad=.4',fc=PALE_TEAL,ec=TEAL),fontsize=9.5); ye=np.interp(e,t,w); ax.scatter([e],[ye],s=58,color=ORANGE,edgecolor='white',zorder=4); ax.annotate("4  E′ attenuation\nsmaller early-relaxation velocity",xy=(e,ye),xytext=(.68,.26),textcoords='axes fraction',arrowprops=dict(arrowstyle='-|>',color=ORANGE,lw=1.4),bbox=dict(boxstyle='round,pad=.4',fc=PALE_ORANGE,ec=ORANGE),fontsize=9.5); ax.text(e,ye-.9,"E′",ha='center',fontsize=11,fontweight='bold',color=NAVY); ax.text(.02,.03,f'Actual representative waveform: woman closest to cohort-average morphology among central MTI values (MTI={a.loc[idx,"MTI"]:.2f} SD)',transform=ax.transAxes,fontsize=8.3,color=SLATE); ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Myocardial velocity (cm/s)'); clean(ax); save(fig,'03_mti_construction')

# 05 actual bootstrap CI
fig,ax=setup('Phase-resolved difference in pretreatment myocardial motion'); delta=evt.mean(0)-ctl.mean(0); rng=np.random.default_rng(20261025); B=1000; boot=np.empty((B,len(t)),np.float32)
for b in range(B): boot[b]=evt[rng.integers(0,len(evt),len(evt))].mean(0)-ctl[rng.integers(0,len(ctl),len(ctl))].mean(0)
lo,hi=np.quantile(boot,[.025,.975],axis=0); ax.fill_between(t,lo,hi,color=ORANGE,alpha=.16,lw=0); ax.plot(t,delta,color=ORANGE,lw=2.5); ax.axhline(0,color=MID,lw=1,ls='--'); ax.axvspan(.499,.600,color=ORANGE,alpha=.08); ax.text(.025,.065,'Global functional permutation: P=0.047\nAfter residualizing LVEF: P=0.056\nBest contiguous cluster: phase 0.499–0.600\nCluster-FWER P=0.094',transform=ax.transAxes,fontsize=10,color=NAVY,va='bottom',linespacing=1.35); ax.set_xlabel('Normalized cardiac cycle'); ax.set_ylabel('Δ myocardial velocity\n(CTRCD − no CTRCD), cm/s'); clean(ax); save(fig,'05_phase_resolved_difference')

# 07 KM actual + extended CI
fig,ax=setup('CTRCD-free survival across MTI tertiles'); ter=pd.qcut(a.MTI,3,labels=['Low','Middle','High']); colors={'Low':TEAL,'Middle':SLATE,'High':ORANGE}
for g in ['Low','Middle','High']:
 m=(ter.astype(str)==g).to_numpy(); times=a.loc[m,'time'].to_numpy(float); ev=a.loc[m,'CTRCD'].to_numpy(int); sf=SurvfuncRight(times,ev); st=np.r_[0,sf.surv_times,times.max()]; sp=np.r_[1,sf.surv_prob,sf.surv_prob[-1]]; se=np.r_[0,sf.surv_prob_se,sf.surv_prob_se[-1]]; lo=np.clip(sp-1.96*se,0,1); hi=np.clip(sp+1.96*se,0,1); ax.step(st,sp,where='post',lw=2.3,color=colors[g],label=f'{g}: {int(ev.sum())}/{int(m.sum())} events'); ax.fill_between(st,lo,hi,step='post',color=colors[g],alpha=.11,lw=0); cen=np.sort(times[ev==0]); sel=cen[np.linspace(0,len(cen)-1,min(8,len(cen))).astype(int)]; vals=np.array([sp[np.searchsorted(st,x,side='right')-1] for x in sel]); ax.plot(sel,vals,ls='none',marker='|',markersize=9,markeredgewidth=1.4,color=colors[g])
ax.set_xlabel('Days from treatment initiation'); ax.set_ylabel('Freedom from CTRCD'); ax.set_ylim(.70,1.01); ax.legend(frameon=False,loc='lower left'); clean(ax); save(fig,'07_kaplan_meier_mti_tertiles')

# 09 surveillance with larger annotation
fig,ax=setup('Risk concentration for surveillance triage'); n=len(pred); total=int(pred.CTRCD.sum()); fr=np.linspace(.05,1,40)
for col,c,lab in [('MTI',ORANGE,'MTI'),('LVEF','#245B79','LVEF'),('MTI+LVEF',TEAL,'MTI + LVEF')]:
 order=np.argsort(-pred[col].to_numpy()); cap=[pred.CTRCD.to_numpy()[order[:max(1,int(np.ceil(n*f)))]].sum()/total for f in fr]; ax.plot(fr*100,np.array(cap)*100,color=c,lw=2.2,label=lab)
ax.plot([0,100],[0,100],color=MID,lw=1.1,ls='--',label='Random selection'); k=int(np.ceil(n*.20)); order=np.argsort(-pred['MTI+LVEF'].to_numpy()); captured=int(pred.CTRCD.to_numpy()[order[:k]].sum()); capture=captured/total; ppv=captured/k; rng=np.random.default_rng(20261026); bv=[]
for b in range(2500):
 idx=rng.integers(0,n,n); d=pred.iloc[idx].reset_index(drop=True); tot=d.CTRCD.sum();
 if tot==0: continue
 kk=max(1,int(np.ceil(len(d)*.20))); oo=np.argsort(-d['MTI+LVEF'].to_numpy()); bv.append(d.CTRCD.to_numpy()[oo[:kk]].sum()/tot)
ci=np.quantile(bv,[.025,.975]); ax.scatter([20],[capture*100],s=68,color=TEAL,edgecolor='white',zorder=4); ax.annotate(f'Top 20% captures {captured}/{total} events ({capture*100:.1f}%)\n95% CI {ci[0]*100:.0f}–{ci[1]*100:.0f}%\nPPV {ppv*100:.1f}% vs 10.0% prevalence',xy=(20,capture*100),xytext=(39,23),arrowprops=dict(arrowstyle='-|>',color=SLATE,lw=1.2),fontsize=10,color=NAVY); ax.set_xlabel('Cohort selected for intensified surveillance (%)'); ax.set_ylabel('Future CTRCD events captured (%)'); ax.set_xlim(0,100); ax.set_ylim(0,103); ax.legend(frameon=False,loc='upper left'); clean(ax); save(fig,'09_surveillance_event_capture')

# 10 reclass + DCA actual saved outputs
fig=plt.figure(figsize=(11,6.8)); fig.patches.extend([Rectangle((.02,.025),.96,.95,transform=fig.transFigure,fill=False,edgecolor='#C4C8CA',lw=1.4,ls=(0,(5,4)),zorder=-20),Rectangle((.02,.845),.96,.13,transform=fig.transFigure,fc=HEADER,ec='none',zorder=-25)]); fig.text(.045,.91,'Incremental 1-year reclassification and decision utility',fontsize=18,fontweight='bold',color=BLACK,va='center'); ax1=fig.add_axes([.08,.15,.38,.60]); ax2=fig.add_axes([.55,.15,.39,.60]); row=reclass[reclass.comparison=='Clinical -> Clinical+MTI'].iloc[0]; metrics=[('NRI_continuous','Continuous NRI'),('NRI_event','Event NRI'),('NRI_nonevent','Non-event NRI'),('IDI','IDI')]; yy=np.arange(4)[::-1]
for yp,(m,label) in zip(yy,metrics):
 v=float(row[m]); lo=float(row[m+'_low']); hi=float(row[m+'_high']); cc=TEAL if m!='IDI' else ORANGE; ax1.errorbar(v,yp,xerr=[[v-lo],[hi-v]],fmt='o',ms=7,capsize=4,color=cc,ecolor=NAVY,lw=1.3); ax1.text(.96,yp,f'{v:+.3f}',va='center',fontsize=9,color=BLACK)
ax1.axvline(0,color=MID,lw=1,ls='--'); ax1.set_yticks(yy,[x[1] for x in metrics]); ax1.set_xlim(-.4,1.18); ax1.set_xlabel('Change after adding MTI to clinical model'); ax1.set_title('Continuous reclassification',loc='left',fontweight='bold',fontsize=12,color=NAVY); clean(ax1)
for m,c,lab in [('Clinical','#245B79','Clinical'),('Clinical+MTI',TEAL,'Clinical + MTI'),('LVEF',SLATE,'LVEF'),('MTI+LVEF',ORANGE,'MTI + LVEF')]:
 z=dca[dca.model==m]; ax2.plot(z.threshold*100,z.net_benefit,color=c,lw=2,label=lab)
for m,c,ls in [('Treat all',MID,'--'),('Treat none','#D4D9DC',':')]:
 z=dca[dca.model==m]; ax2.plot(z.threshold*100,z.net_benefit,color=c,lw=1.1,ls=ls,label=m)
ax2.axvspan(8,15,color=TEAL,alpha=.05); ax2.set_xlim(2,30); ax2.set_ylim(-.025,.09); ax2.set_xlabel('1-year risk threshold (%)'); ax2.set_ylabel('Net benefit'); ax2.set_title('Decision-curve analysis',loc='left',fontweight='bold',fontsize=12,color=NAVY); ax2.legend(frameon=False,ncol=2,fontsize=8,loc='upper right'); clean(ax2); fig.text(.08,.07,'Exploratory analysis among 172 women with known 1-year status (16 CTRCD events).',fontsize=9,color=SLATE); save(fig,'10_reclassification_and_decision_curve')
print('done')
