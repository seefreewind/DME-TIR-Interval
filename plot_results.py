"""Two summary figures from final local TSVs; no individual data are exported."""
from pathlib import Path
import argparse,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();R=a.run;O=R/'figures';O.mkdir(exist_ok=True)
def rd(n):return list(csv.DictReader((R/'phase1/results'/n).open(),delimiter='\t'))
rows=rd('per_eye_strict.tsv');fig,axes=plt.subplots(1,2,figsize=(9,3.5),sharey=True)
for ax,t,n in zip(axes,['T','AC'],[401,157]):
 z=[r for r in rows if r['trial']==t];assert len(z)==n
 for k,col in [(8,'#0072B2'),(12,'#D55E00')]:
  x=np.sort([abs(float(r[f'D{k}_days'])) for r in z]);q=np.quantile(x,.9,method='linear')
  ax.step(np.r_[0,x],np.r_[0,np.arange(1,n+1)/n],where='post',c=col,label=f'{k}-week; P90 {q:.1f} d')
  ax.plot(q,.9,'o',color=col,ms=4)
 ax.set(title=f'Protocol {t} ({n} eyes)',xlabel='Absolute within-eye difference (days)',ylim=(0,1.02),xlim=(0,None));ax.legend(fontsize=8)
axes[0].set_ylabel('Empirical cumulative proportion');fig.tight_layout()
for ext in ['png','svg']:fig.savefig(O/f'absolute_difference_ecdf.{ext}',dpi=300)
plt.close(fig)
rs=rd('treatment_contrasts.tsv');fig,ax=plt.subplots(figsize=(7,3.4));i=0
for t in ['T','AC']:
 for k,col in [(8,'#0072B2'),(12,'#D55E00')]:
  z=[r for r in rs if r['trial']==t and r['analysis_set']=='strict' and r['unit']=='weeks' and r['metric']==f'contrast_change_{k}'];assert len(z)==1;r=z[0]
  v,lo,hi=[float(r[x]) for x in ['estimate','ci_low','ci_high']]
  ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color=col,capsize=3);ax.text(hi+.03,i,f'{v:.2f} [{lo:.2f}, {hi:.2f}]',va='center',fontsize=8);i+=1
ax.axvline(0,c='.4',ls=':',lw=1);ax.set(yticks=range(4),yticklabels=['T: 8 vs 4 weeks','T: 12 vs 4 weeks','AC: 8 vs 4 weeks','AC: 12 vs 4 weeks'],xlabel='Change in adjusted treatment contrast (weeks)');ax.invert_yaxis();ax.margins(x=.35);fig.tight_layout()
for ext in ['png','svg']:fig.savefig(O/f'treatment_contrast_changes.{ext}',dpi=300)
plt.close(fig)
