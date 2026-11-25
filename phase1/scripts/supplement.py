"""Frozen exploratory three-arm contrasts and paired phase strategy changes."""
import analyze as A
import numpy as np
import pathlib,csv,collections
P=A.P;O=P/'results'
def write(name,rows):
 with (O/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
buck=collections.defaultdict(list)
for r in A.rd(P/'normalized/visits.tsv'):
 for c in ['week','day','va','raw_actual_day']:r[c]=float(r[c]) if r[c] else None
 buck[r['trial'],r['participant_id'],r['eye']].append(r)
el={(r['trial'],r['participant_id'],r['eye']):r for r in A.rd(P/'eligibility/eye_eligibility_real.tsv')};rows=[]
for r in A.rd(P/'normalized/eyes.tsv'):
 key=r['trial'],r['participant_id'],r['eye']
 if el[key]['strict_8']!='True':continue
 r['assigned'],_=A.assign_planned(buck[key]);rows.append(r)
rows.sort(key=lambda r:(r['trial'],int(r['participant_id']),r['eye']))
def fit(rr,yy,t,historical=False):
 x=np.column_stack([np.ones(len(rr)),[r['group'] in ['Aflibercept','AC-A'] for r in rr],[float(r['baseline_va'])-69 for r in rr]])
 if t=='AC':x=np.column_stack([x,[r['bilateral']=='True' for r in rr]])
 if historical:x=np.column_stack([x,[r['group']=='Ranibizumab' for r in rr]])
 w=A.boot(rr).w;mat=np.einsum('bi,ij,ik->bjk',w,x,x,optimize=True);rhs=np.einsum('bi,ij,ik->bjk',w,x,yy,optimize=True)
 assert np.all(np.linalg.matrix_rank(mat)==x.shape[1]);return np.linalg.lstsq(x,yy,rcond=None)[0],np.linalg.solve(mat,rhs)
phaseres=[];explore=[]
for t in ['T','AC']:
 rr=[r for r in rows if r['trial']==t and (t=='AC' or r['group']!='Ranibizumab')]
 for k,phases in [(8,[4]),(12,[4,8])]:
  base=np.array([A.outcome(r,k) for r in rr])
  for ph in phases:
   cal=sorted({0,48,52}|{w for w in range(ph,48,k) if 0<w<48});yy=np.column_stack([base,[A.outcome(r,cal) for r in rr]]);pt,bs=fit(rr,yy,t)
   phaseres.append(A.statrow(dict(trial=t,calendar=k,phase=ph,reference_phase=0,eyes=len(rr),people=len({r['participant_id'] for r in rr})),'paired_strategy_contrast_shift',(pt[1,1]-pt[1,0])/7,(bs[:,1,1]-bs[:,1,0])/7,'weeks'))
rr=[r for r in rows if r['trial']=='T'];yy=np.array([[A.outcome(r,k) for k in [4,8,12]] for r in rr]);pt,bs=fit(rr,yy,'T',True)
for label,c in [('Aflibercept-Ranibizumab',np.array([0,1,0,-1])),('Bevacizumab-Ranibizumab',np.array([0,0,0,-1]))]:
 val=c@pt;dist=np.einsum('j,bjk->bk',c,bs)
 for j,k in enumerate([4,8,12]):
  meta=dict(trial='T',analysis_set='strict_historical_exploratory',eyes=401,contrast_direction=label)
  explore.append(A.statrow(meta,f'contrast_{k}',val[j]/7,dist[:,j]/7,'weeks'))
  if j:explore.append(A.statrow(meta,f'contrast_change_{k}',(val[j]-val[0])/7,(dist[:,j]-dist[:,0])/7,'weeks'))
write('phase_strategy_shift.tsv',phaseres);write('exploratory_T_ranibizumab.tsv',explore)
print('Frozen exploratory and paired-phase tables complete; all5000 models full rank')
