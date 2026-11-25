"""Authorized Phase 1. Exact threshold integration; participant cluster bootstrap.
Python3 + NumPy only. No input-file mutation and no manuscript generation.
"""
import pathlib,csv,json,collections,hashlib,sys,math,datetime
import numpy as np
P=pathlib.Path(__file__).resolve().parents[1]
ROOT=next(q for q in [P,*P.parents] if (q/'data/raw').is_dir())
sys.path.insert(0,str(ROOT/'phase0/scripts'))
from phase0_core import assign_planned,CALENDARS
B=5000;SEED=20261124;TOL=1e-8
OUT=None;CHECKS=[];BOOTLOG=[];CACHE={}
def rd(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rr):
 if not rr:return
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def check(name,expected,observed):
 passed=expected==observed if not isinstance(expected,float) else abs(expected-observed)<TOL
 CHECKS.append(dict(check=name,expected=expected,observed=observed,passed=passed));assert passed,(name,expected,observed)
def integral(pts,end=364,method='linear',inclusive=True):
 pts=sorted(pts);assert pts[0][0]==0 and all(b[0]>a[0] for a,b in zip(pts,pts[1:]))
 if pts[-1][0]<end:
  assert 0<=end-pts[-1][0]<=7
  pts=pts+[(end,pts[-1][1])]
 ok=lambda y:y>=69 if inclusive else y>69
 total=0.
 for (a,x),(b,y) in zip(pts,pts[1:]):
  if a>=end:break
  stop=min(b,end);dt=stop-a
  if method=='locf':total+=dt*ok(x);continue
  z=x+(y-x)*dt/(b-a)
  if x==z:total+=dt*ok(x)
  elif min(x,z)>=69:total+=dt
  elif max(x,z)<=69:pass
  else:
   q=(69-x)/(z-x);total+=dt*((1-q) if z>x else q)
 assert -TOL<=total<=end+TOL
 return total

def synthetic():
 check('constant69 inclusive',364.,integral([(0,69),(364,69)]))
 check('constant69 strict',0.,integral([(0,69),(364,69)],inclusive=False))
 check('dip complete days',28.,integral([(0,72),(28,66),(56,72)],end=56))
 check('dip hidden days',56.,integral([(0,72),(56,72)],end=56))
 check('irregular dates',25.,integral([(0,72),(20,66),(50,72)],end=50))
 check('truncate364',182.,integral([(0,60),(400,60+18*400/364)],end=364))
 check('tail7',364.,integral([(0,72),(357,72)]))
 check('LOCF',28.,integral([(0,72),(28,66),(56,72)],end=56,method='locf'))
 r=[dict(week=w,day=w*7,va=72,measurement_id=str(w)) for w in [0,8,12,16,24,32,40,48,52]]
 a,_=assign_planned(r);check('missing4 does not shift calendar8',True,4 not in a and set(CALENDARS[8])<=set(a))

class Boot:
 def __init__(self,rows):
  self.rows=rows;self.n=len(rows);self.ids=sorted({r['participant_id'] for r in rows},key=lambda x:int(x));self.m=len(self.ids)
  sig='|'.join(':'.join([r['trial'],r['participant_id'],r['eye'],r['group']]) for r in rows)
  self.key=hashlib.sha256(sig.encode()).hexdigest()
  seed=np.random.SeedSequence([SEED,int(self.key[:8],16)])
  rng=np.random.default_rng(seed);index={p:i for i,p in enumerate(self.ids)}
  counts=np.zeros((B,self.m),dtype=np.int16)
  groups=sorted({r['group'] for r in rows}) if rows[0]['trial']=='T' else ['ALL']
  for g in groups:
   ids=sorted({r['participant_id'] for r in rows if g=='ALL' or r['group']==g},key=lambda x:int(x))
   inds=np.array([index[p] for p in ids]);draw=rng.integers(0,len(ids),size=(B,len(ids)))
   for j in range(B):counts[j,inds]=np.bincount(draw[j],minlength=len(ids))
  self.w=counts[:,[index[r['participant_id']] for r in rows]].astype(float)
  assert np.all(counts.sum(axis=1)==self.m)
  # All eyes of a participant must have identical multiplicity in every replicate.
  for p in self.ids:
   ix=[j for j,r in enumerate(rows) if r['participant_id']==p]
   assert all(np.array_equal(self.w[:,ix[0]],self.w[:,k]) for k in ix)
  self.den=self.w.sum(axis=1)
  BOOTLOG.append(dict(trial=rows[0]['trial'],population_sha256=self.key,eyes=self.n,people=self.m,replicates=B,seed=SEED,derived_seed_int=int(self.key[:8],16),cluster_check=True))
 def ratio(self,num,den):
  n=self.w@np.asarray(num,float);d=self.w@np.asarray(den,float)
  return np.divide(n,d,out=np.full(B,np.nan),where=d>0)
 def mean(self,y):return self.ratio(y,np.ones(self.n))
 def quant(self,y,p):
  order=np.argsort(y,kind='stable');ys=np.asarray(y)[order];cs=np.cumsum(self.w[:,order],axis=1)
  h=(self.den-1)*p;lo=np.floor(h);hi=np.ceil(h)
  l=np.argmax(cs>lo[:,None],axis=1);u=np.argmax(cs>hi[:,None],axis=1)
  return ys[l]+(h-lo)*(ys[u]-ys[l])

def boot(rows):
 key=tuple((r['trial'],r['participant_id'],r['eye']) for r in rows)
 if key not in CACHE:CACHE[key]=Boot(rows)
 return CACHE[key]
def ci(xs):
 x=np.asarray(xs);good=x[np.isfinite(x)];assert len(good)>=.99*B,('bootstrap usable below99%',len(good))
 return *np.quantile(good,[.025,.975],method='linear'),len(good)
def statrow(meta,metric,value,reps,unit='days'):
 lo,hi,n=ci(reps)
 return dict(**meta,metric=metric,unit=unit,estimate=float(value),ci_low=float(lo),ci_high=float(hi),bootstrap_valid=n)

def outcome(e,calendar,method='linear',inclusive=True,end=364,nominal=False):
 a=e['assigned'];weeks=sorted(a) if calendar==4 else calendar if isinstance(calendar,list) else CALENDARS[calendar]
 points=[((w*7 if nominal else a[w]['day']),a[w]['va']) for w in weeks]
 return integral(points,end,method,inclusive)

def delta_stats(rows,y4,yk,meta):
 bb=boot(rows);d=np.asarray(yk)-np.asarray(y4);ad=np.abs(d)
 vals=[('mean_abs',np.mean(ad),bb.mean(ad)),('mean_signed',np.mean(d),bb.mean(d)),('p90_abs',np.quantile(ad,.9,method='linear'),bb.quant(ad,.9)),('positive_contribution',np.maximum(d,0).mean(),bb.mean(np.maximum(d,0))),('negative_contribution',np.maximum(-d,0).mean(),bb.mean(np.maximum(-d,0)))]
 ret=[]
 for name,v,bs in vals:
  ret.append(statrow(meta,name,v,bs))
  ret.append(statrow(meta,name,v/7,bs/7,'weeks'))
 cancellation=np.mean(ad)-abs(np.mean(d));bcancel=bb.mean(ad)-abs(bb.mean(d))
 ret.append(statrow(meta,'cancellation_days',cancellation,bcancel))
 frac=cancellation/np.mean(ad) if np.mean(ad)>0 else 0
 bs=np.divide(bcancel,bb.mean(ad),out=np.zeros(B),where=bb.mean(ad)>0)
 ret.append(statrow(meta,'cancellation_fraction',frac,bs,'proportion'))
 counts=dict(**meta,n_eyes=len(rows),n_people=len({r['participant_id'] for r in rows}),positive=int(np.sum(d>TOL)),zero=int(np.sum(abs(d)<=TOL)),negative=int(np.sum(d< -TOL)),mean_T4_days=float(np.mean(y4)),mean_sparse_days=float(np.mean(yk)))
 return ret,counts

def model(rows,ys,meta,historical=False):
 # Current explicit instruction: two-arm T primary. Historical three-arm is supplementary.
 if rows[0]['trial']=='T' and not historical:
  keep=[i for i,r in enumerate(rows) if r['group'] in ['Aflibercept','Bevacizumab']]
  rr=[rows[i] for i in keep];yy=np.asarray(ys)[keep]
 else:rr=rows;yy=np.asarray(ys)
 t=rr[0]['trial'];aa='Aflibercept' if t=='T' else 'AC-A'
 X=np.column_stack([np.ones(len(rr)),[r['group']==aa for r in rr],[float(r['baseline_va'])-69 for r in rr]])
 if t=='AC':X=np.column_stack([X,[r['bilateral']=='True' for r in rr]])
 if historical:X=np.column_stack([X,[r['group']=='Ranibizumab' for r in rr]])
 bb=boot(rr);w=bb.w
 mat=np.einsum('bi,ij,ik->bjk',w,X,X,optimize=True);rhs=np.einsum('bi,ij,ik->bjk',w,X,yy,optimize=True)
 good=np.linalg.matrix_rank(mat)==X.shape[1];check('model bootstrap rank '+str(meta),True,bool(np.mean(good)>=.99))
 co=np.full((B,yy.shape[1]),np.nan);co[good]=np.linalg.solve(mat[good],rhs[good])[:,1,:]
 point=np.linalg.lstsq(X,yy,rcond=None)[0][1,:]
 # Independently fit difference as outcome; OLS linearity must hold.
 direct=np.linalg.lstsq(X,yy[:,1:]-yy[:,[0]],rcond=None)[0][1,:]
 check('OLS delta identity '+str(meta),0.,float(np.max(abs(direct-(point[1:]-point[0])))))
 out=[]
 for j,k in enumerate([4,8,12][:yy.shape[1]]):
  mm=dict(**meta,model='T_three_arm_historical' if historical else 'T_two_arm274' if t=='T' and len(rr)==274 else 'T_two_arm' if t=='T' else 'AC_baseline_bilateral',eyes=len(rr),people=len({r['participant_id'] for r in rr}),contrast_direction='Aflibercept-Bevacizumab' if t=='T' else 'AC-A_minus_AC-B')
  out.append(statrow(mm,f'contrast_{k}',point[j]/7,co[:,j]/7,'weeks'))
  if j:out.append(statrow(mm,f'contrast_change_{k}',(point[j]-point[0])/7,(co[:,j]-co[:,0])/7,'weeks'))
 return out

def heldout(rows,k):
 points=[];intervals=[];eye_counts=[]
 for e in rows:
  a=e['assigned'];ret=CALENDARS[k];counts=np.zeros(4,dtype=int);prefix={c:e[c] for c in ['trial','participant_id','eye','group']}
  for w in sorted(a):
   if w in ret or not 0<=a[w]['day']<=364:continue
   d=a[w]['day'];br=[(l,r) for l,r in zip(ret,ret[1:]) if a[l]['day']<d<a[r]['day']];assert len(br)==1
   l,r=br[0];x=a[l];y=a[r];pred=x['va']+(y['va']-x['va'])*(d-x['day'])/(y['day']-x['day']);actual=a[w]['va']>=69;predicted=pred>=69
   idx=int(actual)*2+int(predicted);counts[idx]+=1
   points.append(dict(**prefix,calendar=k,week=w,actual_day=d,actual_va=a[w]['va'],predicted_va=pred,actual_attaining=actual,predicted_attaining=predicted,left_week=l,right_week=r))
  for l,r in zip(ret,ret[1:]):
   hidden=[w for w in a if w not in ret and a[l]['day']<a[w]['day']<a[r]['day'] and 0<=a[w]['day']<=364]
   if not hidden:continue
   both=a[l]['va']>=69 and a[r]['va']>=69;nbelow=sum(a[w]['va']<69 for w in hidden)
   intervals.append(dict(**prefix,calendar=k,left_week=l,right_week=r,hidden_points=len(hidden),both_endpoints_attain=both,hidden_below=nbelow,missed_dip_interval=both and nbelow>0))
  eye_counts.append(counts)
 return points,intervals,np.asarray(eye_counts)

def holdout_summary(rows,k,counts,intervals,aset):
 bb=boot(rows);ss=counts.sum(axis=0);n=counts.sum(axis=1);neg=counts[:,0]+counts[:,1];pos=counts[:,2]+counts[:,3];err=counts[:,1]+counts[:,2]
 meta=dict(trial=rows[0]['trial'],analysis_set=aset,calendar=k,eyes=len(rows),people=len({r['participant_id'] for r in rows}),points=int(sum(ss)))
 tab=[dict(**meta,actual_attaining=a,predicted_attaining=p,count=int(ss[a*2+p])) for a in [0,1] for p in [0,1]]
 metrics=[]
 for label,num,den in [('missed_dip',counts[:,1],neg),('missed_attainment',counts[:,2],pos),('overall_misclassification',err,n)]:
  pointval=float(sum(num)/sum(den)) if sum(den)>0 else math.nan
  if sum(den)>0:metrics.append(dict(**statrow(meta,label+'_point_pooled',pointval,bb.ratio(num,den),'proportion'),numerator=int(sum(num)),denominator=int(sum(den)),denominator_eyes=int(sum(den>0))))
  valid=den>0;rate=np.divide(num,den,out=np.zeros(len(rows)),where=valid)
  if sum(valid):metrics.append(dict(**statrow(meta,label+'_eye_equal',float(rate[valid].mean()),bb.ratio(rate,valid),'proportion'),numerator='',denominator=int(sum(valid)),denominator_eyes=int(sum(valid))))
 ic=collections.defaultdict(lambda:[0,0])
 for r in intervals:
  key=r['participant_id'],r['eye'];ic[key][0]+=r['missed_dip_interval'];ic[key][1]+=r['both_endpoints_attain']
 nums=[ic[r['participant_id'],r['eye']][0] for r in rows];dens=[ic[r['participant_id'],r['eye']][1] for r in rows]
 if sum(dens):metrics.append(dict(**statrow(meta,'missed_dip_interval_point_pooled',sum(nums)/sum(dens),bb.ratio(nums,dens),'proportion'),numerator=sum(nums),denominator=sum(dens),denominator_eyes=sum(x>0 for x in dens)))
 return tab,metrics

def main():
 global OUT
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=pathlib.Path,default=P/'results');args=ap.parse_args();OUT=args.out;OUT.mkdir(parents=True,exist_ok=False)
 config=json.loads((P/'RUN_CONFIG.json').read_text());assert config['PHASE1_AUTHORIZATION'] is True and config['MANUSCRIPT_AUTHORIZATION'] is False
 assert (P/'frozen_version_comparison.tsv').exists()
 synthetic()
 raw=rd(P/'normalized/visits.tsv');eyes=rd(P/'normalized/eyes.tsv');elig=rd(P/'eligibility/eye_eligibility_real.tsv');byelig={(r['trial'],r['participant_id'],r['eye']):r for r in elig};buck=collections.defaultdict(list)
 for r in raw:
  for c in ['week','day','va','raw_actual_day']:r[c]=float(r[c]) if r[c] else None
  buck[r['trial'],r['participant_id'],r['eye']].append(r)
 for e in eyes:
  key=e['trial'],e['participant_id'],e['eye'];e['assigned'],_=assign_planned(buck[key]);e['elig']=byelig[key]
 eyes.sort(key=lambda r:(r['trial'],int(r['participant_id']),r['eye']))
 summary=[];signs=[];contrasts=[];per=[];qa=[];sens=[];senscounts=[];senscontrast=[];phasecounts=[];phasediff=[];holdpoints=[];holdintervals=[];holdtables=[];holdmetrics=[];selection=[];overlap=[];excluded=[]
 for t in ['T','AC']:
  allrows=[e for e in eyes if e['trial']==t];strict=[e for e in allrows if e['elig']['strict_8']=='True'];check(t+' strict n',401 if t=='T' else 157,len(strict));check(t+' strict people',401 if t=='T' else 139,len({r['participant_id'] for r in strict}))
  ys=np.array([[outcome(e,k) for k in [4,8,12]] for e in strict])
  for i,e in enumerate(strict):
   dd=dict(**{c:e[c] for c in ['trial','participant_id','eye','group']},analysis_set='strict',baseline_va=e['baseline_va'],baseline_actual_day=e['baseline_actual_day'],actual48=e['assigned'][48]['day'],annual_day=e['assigned'][52]['day'])
   for j,k in enumerate([4,8,12]):dd[f'T{k}_days']=ys[i,j];dd[f'T{k}_weeks']=ys[i,j]/7
   for j,k in enumerate([8,12],1):dd[f'D{k}_days']=ys[i,j]-ys[i,0];dd[f'D{k}_weeks']=(ys[i,j]-ys[i,0])/7
   per.append(dd)
   a48=e['assigned'][48]['day'];until=np.array([outcome(e,k,end=a48) for k in [4,8,12]])
   for j,k in enumerate([8,12],1):
    err=(ys[i,j]-ys[i,0])-(until[j]-until[0]);assert abs(err)<TOL
    qa.append(dict(**{c:e[c] for c in ['trial','participant_id','eye']},calendar=k,actual48=a48,D_full_days=ys[i,j]-ys[i,0],D_until_actual48_days=until[j]-until[0],absolute_discrepancy=abs(err),passed=abs(err)<TOL))
  for j,k in enumerate([8,12],1):
   rr,cc=delta_stats(strict,ys[:,0],ys[:,j],dict(trial=t,analysis_set='strict',calendar=k));summary+=rr;signs.append(cc)
   hp,hi,hc=heldout(strict,k);check(t+f' holdout{k} points',(2406 if k==8 else 3208) if t=='T' else (942 if k==8 else 1256),len(hp));ht,hm=holdout_summary(strict,k,hc,hi,'strict');holdpoints+=hp;holdintervals+=hi;holdtables+=ht;holdmetrics+=hm
  contrasts+=model(strict,ys,dict(trial=t,analysis_set='strict'))
  if t=='T':contrasts+=model(strict,ys,dict(trial=t,analysis_set='strict_historical'),historical=True)
  # Inclusion and baseline distribution use all randomized eyes, original arms.
  for g in sorted({r['group'] for r in allrows}):
   ar=[r for r in allrows if r['group']==g]
   for included in [True,False]:
    rr=[r for r in ar if (r['elig']['strict_8']=='True')==included];base=[float(r['baseline_va']) for r in rr if r['baseline_va']!='']
    selection.append(dict(trial=t,group=g,included=included,randomized_eyes=len(ar),eyes=len(rr),people=len({r['participant_id'] for r in rr}),fraction=len(rr)/len(ar),baseline_n=len(base),baseline_mean=np.mean(base),baseline_sd=np.std(base,ddof=1),baseline_q25=np.quantile(base,.25),baseline_median=np.median(base),baseline_q75=np.quantile(base,.75),baseline_min=min(base),baseline_max=max(base)))
  # Isolated algorithm/calendar/population changes.
  scenarios=[('LOCF','algorithm',strict,dict(method='locf')),('strict_gt69','algorithm',strict,dict(inclusive=False)),('actual0_336','calendar_time_window',strict,dict(end=336)),('nominal0_48','calendar_time_axis',strict,dict(end=336,nominal=True)),('actual_day0_baseline','population',[e for e in strict if float(e['baseline_actual_day'])==0],{}),('no_annual_carry','population',[e for e in strict if e['elig']['boundary']!='tail_le7'],{}),('max_gap_le63','population',[e for e in strict if float(e['elig']['max_gap'])<=63],{})]
  if t=='AC':
   scenarios += [('AC_annual308_420','population',[e for e in strict if 308<=float(e['elig']['annual_day'])<=420],{}),('AC_gap_and_annual_union','population',[e for e in strict if float(e['elig']['max_gap'])<=63 and 308<=float(e['elig']['annual_day'])<=420],{})]
  for name,kind,rr,kw in scenarios:
   # fixed336 supported by common retained actual48/annual; nominal48 omits52.
   out=np.array([[outcome(e,(list(range(0,49,4)) if k==4 else list(range(0,49,k))) if name=='nominal0_48' else k,**kw) for k in [4,8,12]] for e in rr])
   for j,k in enumerate([8,12],1):
    ss,cc=delta_stats(rr,out[:,0],out[:,j],dict(trial=t,analysis_set=name,calendar=k));sens+=ss;senscounts.append(dict(**cc,change_type=kind))
   senscontrast+=model(rr,out,dict(trial=t,analysis_set=name))
   for e in strict:
    if e not in rr:excluded.append(dict(trial=t,scenario=name,participant_id=e['participant_id'],eye=e['eye'],group=e['group'],max_gap=e['elig']['max_gap'],annual_day=e['elig']['annual_day']))
  # Phase shift: internal phase+m*j, retain0/48/52; common sample across all phases.
  for k,phases in [(8,[0,4]),(12,[0,4,8])]:
   cals={ph:sorted({0,48,52}|{w for w in range(ph,48,k) if 0<w<48}) for ph in phases}
   common=[e for e in strict if all(set(c)<=set(e['assigned']) for c in cals.values())]
   for ph,cal in cals.items():
    rr=[e for e in strict if set(cal)<=set(e['assigned'])];phasecounts.append(dict(trial=t,calendar=k,phase=ph,target_weeks=','.join(map(str,cal)),eligible_eyes=len(rr),common_eyes=len(common)))
    yy4=np.array([outcome(e,4) for e in common]);yy=np.array([outcome(e,cal) for e in common]);meta=dict(trial=t,analysis_set=f'phase{k}_{ph}',calendar=k)
    ss,cc=delta_stats(common,yy4,yy,meta);sens+=ss;senscounts.append(dict(**cc,change_type='calendar_phase'))
    phase_model=model(common,np.column_stack([yy4,yy]),dict(trial=t,analysis_set=f'phase{k}_{ph}'))
    if k==12:
     for c in phase_model:c['metric']=c['metric'].replace('_8','_12')
    senscontrast.extend(phase_model)
    base=np.array([outcome(e,cals[0]) for e in common]);change=np.abs(yy-yy4)-np.abs(base-yy4)
    phasediff.append(statrow(dict(trial=t,calendar=k,phase=ph,eyes=len(common)),'MAE_change_vs_phase0',np.mean(change),boot(common).mean(change)))
  # Frozen wide sets: no per-calendar complete-case change.
  for k,expected in [(8,465 if t=='T' else 187),(12,486 if t=='T' else 203)]:
   rr=[e for e in allrows if e['elig'][f'wide_{k}']=='True'];check(t+f' wide{k}',expected,len(rr));y4=np.array([outcome(e,4) for e in rr]);yk=np.array([outcome(e,k) for e in rr]);ss,cc=delta_stats(rr,y4,yk,dict(trial=t,analysis_set='wide',calendar=k));sens+=ss;senscounts.append(dict(**cc,change_type='population'))
   # For a two-calendar wide comparison, label the sparse model by its actual calendar.
   co=model(rr,np.column_stack([y4,yk]),dict(trial=t,analysis_set=f'wide{k}'))
   if k==12:
    for c in co:c['metric']=c['metric'].replace('_8','_12')
   senscontrast+=co
   strictkeys={(e['participant_id'],e['eye']) for e in strict};widekeys={(e['participant_id'],e['eye']) for e in rr}
   overlap.append(dict(trial=t,calendar=k,strict=len(strict),wide=len(rr),intersection=len(strictkeys&widekeys),strict_only=len(strictkeys-widekeys),wide_only=len(widekeys-strictkeys)))
   hp,hi,hc=heldout(rr,k);ht,hm=holdout_summary(rr,k,hc,hi,'wide');holdtables+=ht;holdmetrics+=hm
  common=[e for e in allrows if all(e['elig'][f'wide_{k}']=='True' for k in [8,12])];check(t+' widecommon',439 if t=='T' else 175,len(common))
  for aset,rr in [('wide_common',common),('strict_common',strict)]:
   yy=np.array([[outcome(e,k) for k in [4,8,12]] for e in rr]);bb=boot(rr)
   for j,k in enumerate([8,12],1):
    ss,cc=delta_stats(rr,yy[:,0],yy[:,j],dict(trial=t,analysis_set=aset,calendar=k));sens+=ss;senscounts.append(dict(**cc,change_type='population_common'))
   for metric,d in [('abs_error12_minus8',abs(yy[:,2]-yy[:,0])-abs(yy[:,1]-yy[:,0])),('signed_error12_minus8',yy[:,2]-yy[:,1])]:phasediff.append(statrow(dict(trial=t,calendar=12,phase=aset,eyes=len(rr)),metric,np.mean(d),bb.mean(d)))
   if aset=='wide_common':senscontrast+=model(rr,yy,dict(trial=t,analysis_set=aset))
 # Additional deterministic bootstrap unit check for bilateral opposite-arm enrollment.
 example=[dict(trial='AC',participant_id='1',eye='OD',group='AC-A'),dict(trial='AC',participant_id='1',eye='OS',group='AC-B'),dict(trial='AC',participant_id='2',eye='OD',group='AC-A')]
 bb=boot(example);check('bilateral opposite arms share multiplicity',True,bool(np.array_equal(bb.w[:,0],bb.w[:,1])))
 for name,rows in [('main_results.tsv',summary),('direction_counts.tsv',signs),('treatment_contrasts.tsv',contrasts),('per_eye_strict.tsv',per),('numeric_tail_identity.tsv',qa),('holdout_points_strict.tsv',holdpoints),('holdout_intervals_strict.tsv',holdintervals),('holdout_2x2.tsv',holdtables),('holdout_validation.tsv',holdmetrics),('sensitivity_results.tsv',sens),('sensitivity_counts.tsv',senscounts),('sensitivity_contrasts.tsv',senscontrast),('phase_eligibility.tsv',phasecounts),('paired_calendar_comparison.tsv',phasediff),('sample_selection.tsv',selection),('set_overlap.tsv',overlap),('sensitivity_excluded_eyes.tsv',excluded),('validation_checks.tsv',CHECKS),('bootstrap_log.tsv',BOOTLOG)]:write(name,rows)
 (OUT/'environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,bootstrap_replicates=B,seed=SEED,config_sha256=hashlib.sha256((P/'RUN_CONFIG.json').read_bytes()).hexdigest(),real_TIR_authorized=True,manuscript_authorized=False),indent=2))
 print('Phase1 complete',OUT,'checks',len(CHECKS),'tail identities',len(qa),'populations',len(BOOTLOG))
if __name__=='__main__':main()
