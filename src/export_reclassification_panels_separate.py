from pathlib import Path
import pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'/'individual_cardiomorph_code'; T=ROOT/'results'/'tables'
NAVY='#183B56'; TEAL='#0F7C82'; ORANGE='#D95D2B'; SLATE='#6F8191'; MID='#A9B5BE'; HEADER='#ECEEEE'; BLACK='#202A33'
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'figure.facecolor':'white','savefig.facecolor':'white'})
def setup(title):
 fig=plt.figure(figsize=(9.2,6.9)); fig.patches.extend([Rectangle((.02,.025),.96,.95,transform=fig.transFigure,fill=False,edgecolor='#C4C8CA',lw=1.4,ls=(0,(5,4)),zorder=-20),Rectangle((.02,.845),.96,.13,transform=fig.transFigure,fc=HEADER,ec='none',zorder=-25)]); fig.text(.045,.91,title,fontsize=18,fontweight='bold',color=BLACK,va='center'); ax=fig.add_axes([.13,.15,.78,.60]); return fig,ax
def clean(ax):
 ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['left'].set_color('#536675'); ax.spines['bottom'].set_color('#536675'); ax.grid(False)
def save(fig,stem):
 for ext,dpi in [('png',300),('svg',None),('tiff',300)]:
  kw={'bbox_inches':'tight','pad_inches':.05};
  if dpi: kw['dpi']=dpi
  if ext=='tiff': kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(OUT/f'{stem}.{ext}',**kw)
 plt.close(fig)
reclass=pd.read_csv(T/'one_year_reclassification_nri_idi.csv'); row=reclass[reclass.comparison=='Clinical -> Clinical+MTI'].iloc[0]
fig,ax=setup('Continuous 1-year reclassification after adding MTI'); metrics=[('NRI_continuous','Continuous NRI'),('NRI_event','Event NRI'),('NRI_nonevent','Non-event NRI'),('IDI','IDI')]; yy=np.arange(4)[::-1]
for yp,(m,label) in zip(yy,metrics):
 v=float(row[m]); lo=float(row[m+'_low']); hi=float(row[m+'_high']); c=TEAL if m!='IDI' else ORANGE; ax.errorbar(v,yp,xerr=[[v-lo],[hi-v]],fmt='o',ms=9,capsize=5,color=c,ecolor=NAVY,lw=1.5); ax.text(1.02,yp,f'{v:+.3f}  [{lo:+.3f}, {hi:+.3f}]',va='center',fontsize=10,color=BLACK)
ax.axvline(0,color=MID,lw=1.1,ls='--'); ax.set_yticks(yy,[x[1] for x in metrics]); ax.set_xlim(-.4,1.55); ax.set_xlabel('Improvement after adding MTI to clinical model'); ax.text(.02,.03,'172 women with known 1-year status; 16 CTRCD events.\nConfidence intervals from patient bootstrap resampling.',transform=ax.transAxes,fontsize=9.5,color=SLATE); clean(ax); save(fig,'10a_reclassification_nri_idi')

dca=pd.read_csv(T/'one_year_decision_curve.csv'); fig,ax=setup('Decision-curve analysis at 1 year')
for m,c,lab in [('Clinical','#245B79','Clinical model'),('Clinical+MTI',TEAL,'Clinical + MTI'),('LVEF',SLATE,'LVEF'),('MTI+LVEF',ORANGE,'MTI + LVEF')]:
 z=dca[dca.model==m]; ax.plot(z.threshold*100,z.net_benefit,color=c,lw=2.3,label=lab)
for m,c,ls in [('Treat all',MID,'--'),('Treat none','#D4D9DC',':')]:
 z=dca[dca.model==m]; ax.plot(z.threshold*100,z.net_benefit,color=c,lw=1.3,ls=ls,label=m)
ax.axvspan(8,15,color=TEAL,alpha=.05); ax.set_xlim(2,30); ax.set_ylim(-.025,.09); ax.set_xlabel('1-year predicted-risk threshold (%)'); ax.set_ylabel('Net benefit'); ax.legend(frameon=False,ncol=2,loc='upper right'); ax.text(.02,.04,'Shaded region highlights 8%–15% thresholds.\nExploratory fixed-horizon analysis.',transform=ax.transAxes,fontsize=9.5,color=SLATE); clean(ax); save(fig,'10b_decision_curve')
print('done')
