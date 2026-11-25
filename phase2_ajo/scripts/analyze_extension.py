"""POST_PHASE1_EXPLORATORY_EXTENSION; immutable historical inputs; no new cohort."""
from pathlib import Path
import csv,json,sys,hashlib,collections,math,datetime
import numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent
sys.path.insert(0,str(R/'phase1/scripts'))
import analyze as A
C=json.loads((P/'RUN_CONFIG_EXTENSION.json').read_text());B=C['bootstrap']['B'];TOL=1e-8
A.B=B;A.SEED=20261124
O=P/'results';O.mkdir(exist_ok=True)
CHECKS=[];BOOT=[]
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write(name,rr):
 if not rr:return
 keys=list(dict.fromkeys(k for r in rr for k in r))
 with (P/name).open('w',newline='') as f:
  w=csv.DictWriter(f,keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def check(label,ok,detail=''):
 CHECKS.append(dict(check=label,passed=bool(ok),detail=str(detail)));assert ok,(label,detail)
def q(y,w,p=.9):
 order=np.argsort(y,kind='stable');ys=np.asarray(y)[order];ww=w[:,order];cs=np.cumsum(ww,axis=1);den=ww.sum(axis=1)
 h=(den-1)*p;l=np.floor(h);u=np.ceil(h)
 li=np.argmax(cs>l[:,None],axis=1);ui=np.argmax(cs>u[:,None],axis=1)
 ans=ys[li]+(h-l)*(ys[ui]-ys[li]);ans[den<=0]=np.nan;return ans

def mean(y,w):
 den=w.sum(axis=1);return np.divide(w@np.asarray(y,float),den,out=np.full(len(w),np.nan),where=den>0)
def sd(y,w):
 y=np.asarray(y,float);den=w.sum(axis=1);mu=mean(y,w);ss=w@(y*y)-den*mu*mu
 return np.sqrt(np.divide(np.maximum(ss,0),den-1,out=np.full(len(w),np.nan),where=den>1))
def row(meta,metric,v,bs,unit,ci_allowed=True,**extra):
 z=np.asarray(bs,float);z=z[np.isfinite(z)]
 if len(z) and ci_allowed:lo,hi=np.quantile(z,[.025,.975],method='linear')
 else:lo=hi=math.nan
 out=dict(label='POST_PHASE1_EXPLORATORY_EXTENSION',**meta,metric=metric,unit=unit,estimate=float(v),ci_low=float(lo),ci_high=float(hi),bootstrap_valid=len(z),bootstrap_failed=B-len(z),ci_status='descriptive_unadjusted' if ci_allowed and len(z) else 'point_only_small_stratum' if not ci_allowed else 'not_estimable',**extra)
 return out

def curve(e,k):
 a=e['assigned'];pts=[(a[w]['day'],a[w]['va']) for w in A.CALENDARS[k]];x=np.array([v[0] for v in pts]);y=np.array([v[1] for v in pts]);check_tail=x[-1]<364
 if check_tail:
  assert 357<=x[-1]<=363;x=np.append(x,364);y=np.append(y,y[-1])
 xx=np.r_[0,x[(x>0)&(x<364)],364];yy=np.interp(xx,x,y)
 assert len(np.unique(xx))==len(xx) and xx[0]==0 and xx[-1]==364
 return xx,yy

def tir_xy(x,y,h):
 total=0.
 for a,b,v,z in zip(x,x[1:],y,y[1:]):
  if v>=h and z>=h:total+=b-a
  elif v<h and z<h:continue
  elif v!=z:
   cross=a+(h-v)*(b-a)/(z-v);total+=(b-cross) if z>v else (cross-a)
 return float(total)

def load():
 buck=collections.defaultdict(list)
 for r in read(R/'phase1/normalized/visits.tsv'):
  for k in ['week','day','va','raw_actual_day']:r[k]=float(r[k]) if r[k] else None
  buck[r['trial'],r['participant_id'],r['eye']].append(r)
 elig={(r['trial'],r['participant_id'],r['eye']):r for r in read(R/'phase1/eligibility/eye_eligibility_real.tsv')}
 rows=[]
 for r in read(R/'phase1/normalized/eyes.tsv'):
  key=r['trial'],r['participant_id'],r['eye'];r['assigned'],_=A.assign_planned(buck[key]);r['elig']=elig[key]
  if r['elig']['strict_8']=='True':rows.append(r)
 return sorted(rows,key=lambda r:(r['trial'],int(r['participant_id']),r['eye']))

def fit(rows,yy,t):
 ix=[i for i,e in enumerate(rows) if t=='AC' or e['group'] in ['Aflibercept','Bevacizumab']]
 rr=[rows[i] for i in ix];y=np.asarray(yy)[ix];bb=A.boot(rr)
 X=np.column_stack([np.ones(len(rr)),[r['group'] in ['Aflibercept','AC-A'] for r in rr],[float(r['baseline_va'])-69 for r in rr]])
 if t=='AC':X=np.column_stack([X,[r['bilateral']=='True' for r in rr]])
 mat=np.einsum('bi,ij,ik->bjk',bb.w,X,X,optimize=True);rhs=np.einsum('bi,ij,ik->bjk',bb.w,X,y,optimize=True)
 good=np.linalg.matrix_rank(mat)==X.shape[1];bs=np.full((B,y.shape[1]),np.nan);bs[good]=np.linalg.solve(mat[good],rhs[good])[:,1,:]
 point=np.linalg.lstsq(X,y,rcond=None)[0][1,:]
 BOOT.append(dict(trial=t,population='treatment',eyes=len(rr),people=len(bb.ids),population_sha256=bb.key,replicates=B,rank_failures=int((~good).sum())))
 return rr,point,bs

def summary(rows,yy,t,endpoint,unit):
 bb=A.boot(rows);w=bb.w;ret=[];den=np.std(yy[:,0],ddof=1);bden=sd(yy[:,0],w)
 for j,k in enumerate([8,12],1):
  d=yy[:,j]-yy[:,0];ad=abs(d);base=dict(trial=t,endpoint=endpoint,calendar=k,eyes=len(rows),people=len(bb.ids))
  vals=[('reference_mean',np.mean(yy[:,0]),mean(yy[:,0],w)),('reference_sd',den,bden),('sparse_mean',np.mean(yy[:,j]),mean(yy[:,j],w)),('sparse_sd',np.std(yy[:,j],ddof=1),sd(yy[:,j],w)),('mean_signed',np.mean(d),mean(d,w)),('mean_abs',np.mean(ad),mean(ad,w)),('median_abs',np.median(ad),q(ad,w,.5)),('p90_abs',np.quantile(ad,.9),q(ad,w,.9))]
  for name,v,bs in vals:ret.append(row(base,name,v,bs,unit))
  val=np.mean(ad)/den if den>0 else math.nan;bs=np.divide(mean(ad,w),bden,out=np.full(B,np.nan),where=bden>0)
  ret.append(row(base,'standardized_discrepancy',val,bs,'reference_SD'))
 return ret

def main():
 assert C['label']=='POST_PHASE1_EXPLORATORY_EXTENSION'
 historical=read(R/'phase1/results/per_eye_strict.tsv');old={(r['trial'],r['participant_id'],r['eye']):r for r in historical}
 rows=load();check('T401',sum(e['trial']=='T' for e in rows)==401);check('AC157',sum(e['trial']=='AC' for e in rows)==157)
 # Hand-checkable synthetic examples before real endpoint calculations.
 check('AUC constant change0',abs(np.trapezoid(np.array([69,69])-69,[0,364])/364)<TOL)
 check('AUC irregular trapezoids',abs(np.trapezoid([0,10,0],[0,100,364])/364-5)<TOL)
 check('TIR 72,66,72',abs(tir_xy([0,28,56],[72,66,72],69)-28)<TOL)
 check('TIR constant69',tir_xy([0,364],[69,69],69)==364)
 endpoint=[];contrasts=[];cross=[];multi=[];thresholdclass=[];subgroups=[];clinical=[];held=[];holdthreshold=[];pereyes=[];trajectories=[];example=[]
 points=read(R/'phase1/results/holdout_points_strict.tsv')
 for t in ['T','AC']:
  rr=[e for e in rows if e['trial']==t];bb=A.boot(rr);w=bb.w;n=len(rr)
  check(t+' people',len(bb.ids)==(401 if t=='T' else 139));BOOT.append(dict(trial=t,population='strict_all',eyes=n,people=len(bb.ids),population_sha256=bb.key,replicates=B,rank_failures=0))
  ytir=np.array([[float(old[t,e['participant_id'],e['eye']][f'T{k}_days']) for k in [4,8,12]] for e in rr]);auc=[];ann=[];threshold_arrays={h:[] for h in C['thresholds']};cl=[]
  for i,e in enumerate(rr):
   baseline=float(e['baseline_va']);cv={k:curve(e,k) for k in [4,8,12]}
   au=[np.trapezoid(y-baseline,x)/364 for x,y in cv.values()];an=e['assigned'][52]['va']-baseline
   auc.append(au);ann.append([an]*3)
   for k,j in zip([4,8,12],range(3)):check(f'{t} TIR69 reproduce {i}:{k}',abs(tir_xy(*cv[k],69)-ytir[i,j])<TOL)
   observed=[e['assigned'][a]['va'] for a in A.CALENDARS[4] if e['assigned'][a]['day']<=364]
   state=np.array(observed)>=69;down=int(np.sum(state[:-1]&~state[1:]));up=int(np.sum(~state[:-1]&state[1:]));total=down+up
   x,y=cv[4];allabove=bool(min(y)>=69);allbelow=bool(max(y)<69)
   cls='all_reference_knots_attain' if allabove else 'all_reference_knots_below' if allbelow else f'{total}_observed_conversion'+('s' if total==2 else '') if total in [1,2] else '3plus_observed_conversions' if total>=3 else 'boundary_only_change_or_no_observed_conversion'
   cl.append(cls)
   trajectories.append(dict(trial=t,participant_id=e['participant_id'],eye=e['eye'],observed_measurements=len(observed),down=down,up=up,total=total,all_reference_knots_attain=allabove,all_reference_knots_below=allbelow,trajectory_class=cls))
   pe=dict(trial=t,participant_id=e['participant_id'],eye=e['eye'],group=e['group'],baseline_va=baseline,trajectory_class=cls,down=down,up=up,total_conversions=total)
   for j,k in enumerate([4,8,12]):pe.update({f'TIR{k}_days':ytir[i,j],f'AUC{k}_letters':au[j],f'anniversary{k}_letters':an})
   pereyes.append(pe)
   for h in C['thresholds']:
    vals=[tir_xy(*cv[k],h) for k in [4,8,12]];threshold_arrays[h].append(vals)
  auc=np.array(auc);ann=np.array(ann)
  for name,y,unit in [('TIR_days',ytir,'days'),('TIR_fraction',ytir/364,'proportion_of_year'),('time_averaged_VA_change',auc,'letters'),('anniversary_VA_change',ann,'letters')]:
   endpoint+=summary(rr,y,t,name,unit)
   modelrows,pt,bs=fit(rr,y,t)
   for j,k in enumerate([4,8,12]):
    meta=dict(trial=t,endpoint=name,calendar=k,eyes=len(modelrows),people=len({e['participant_id'] for e in modelrows}),direction='aflibercept_minus_bevacizumab' if t=='T' else 'aflibercept_monotherapy_minus_bevacizumab_first')
    contrasts.append(row(meta,'adjusted_contrast',pt[j],bs[:,j],unit))
    if j:contrasts.append(row(meta,'paired_contrast_change',pt[j]-pt[0],bs[:,j]-bs[:,0],unit))
  for j,k in enumerate([8,12],1):
   sdT=np.std(ytir[:,0],ddof=1);sdA=np.std(auc[:,0],ddof=1);dT=abs(ytir[:,j]-ytir[:,0]);dA=abs(auc[:,j]-auc[:,0]);bT=mean(dT,w)/sd(ytir[:,0],w);bA=mean(dA,w)/sd(auc[:,0],w)
   cross.append(row(dict(trial=t,calendar=k,eyes=n,people=len(bb.ids)),'standardized_TIR_minus_AUC',dT.mean()/sdT-dA.mean()/sdA,bT-bA,'reference_SD',TIR_reference_sd=sdT,AUC_reference_sd=sdA,TIR_standardized=dT.mean()/sdT,AUC_standardized=dA.mean()/sdA))
   for limit in C['discrepancy_days']:
    indicator=dT>=limit;clinical.append(row(dict(trial=t,calendar=k,eyes=n,people=len(bb.ids),threshold_days=limit),'absolute_TIR_difference_ge',indicator.mean(),mean(indicator,w),'proportion',numerator=int(indicator.sum()),denominator=n))
  for h,raw in threshold_arrays.items():
   yy=np.array(raw);ptrows,pt,bs=fit(rr,yy,t)
   fullabove=[];fullbelow=[];obschange=[]
   for e in rr:
    x,y=curve(e,4);fullabove.append(min(y)>=h);fullbelow.append(max(y)<h)
    ob=np.array([e['assigned'][a]['va'] for a in A.CALENDARS[4] if e['assigned'][a]['day']<=364])>=h
    obschange.append(bool(np.any(ob[1:]!=ob[:-1])))
   thresholdclass.append(dict(trial=t,threshold=h,eyes=n,people=len(bb.ids),reconstructed_all_attain=int(sum(fullabove)),reconstructed_all_below=int(sum(fullbelow)),observed_any_state_conversion=int(sum(obschange)),classification_basis=C['threshold_floor_ceiling']))
   for j,k in enumerate([8,12],1):
    meta=dict(trial=t,threshold=h,calendar=k,eyes=n,people=len(bb.ids));d=yy[:,j]-yy[:,0];ad=abs(d)
    for name,v,rep,unit in [('reference_mean',yy[:,0].mean(),mean(yy[:,0],w),'days'),('sparse_mean',yy[:,j].mean(),mean(yy[:,j],w),'days'),('mean_abs',ad.mean(),mean(ad,w),'days'),('mean_signed',d.mean(),mean(d,w),'days'),('p90_abs',np.quantile(ad,.9),q(ad,w),'days'),('zero_fraction',np.mean(ad<=TOL),mean(ad<=TOL,w),'proportion')]:
     multi.append(row(meta,name,v,rep,unit))
    multi.append(row({**meta,'model_eyes':len(ptrows)},'treatment_contrast_change',(pt[j]-pt[0])/7,(bs[:,j]-bs[:,0])/7,'weeks'))
  # Baseline and observed-trajectory explanations; resample original participant clusters.
  base=np.array([float(e['baseline_va']) for e in rr]);cl=np.array(cl)
  masks=[('baseline_VA','below69',base<69),('baseline_VA','atleast69',base>=69)]+[('reference_trajectory',v,cl==v) for v in C['trajectory_classes']]
  for kind,level,mask in masks:
   ix=np.where(mask)[0];eyes=len(ix);people=len({rr[i]['participant_id'] for i in ix});allow=eyes>=20 and people>=20
   for j,k in enumerate([8,12],1):
    meta=dict(trial=t,stratifier=kind,stratum=level,calendar=k,eyes=eyes,people=people)
    if not eyes:
     subgroups.append(row(meta,'empty_stratum',math.nan,np.full(B,np.nan),'none',False));continue
    ww=w[:,ix];d=ytir[ix,j]-ytir[ix,0];ad=abs(d)
    vals=[('mean_abs',ad.mean(),mean(ad,ww),'days'),('mean_signed',d.mean(),mean(d,ww),'days'),('p90_abs',np.quantile(ad,.9),q(ad,ww),'days')]
    vals += [(f'abs_ge{lim}',np.mean(ad>=lim),mean(ad>=lim,ww),'proportion') for lim in [14,28]]
    for name,v,rep,unit in vals:subgroups.append(row(meta,name,v,rep,unit,allow,numerator=int(sum(ad>=int(name[6:]))) if name.startswith('abs_ge') else '',denominator=eyes))
  # Transparent complex examples, not a claim of representativeness.
  candidates=[r for r in pereyes if r['trial']==t];chosen=sorted(candidates,key=lambda r:(-r['total_conversions'],-abs(r['TIR8_days']-r['TIR4_days']),int(r['participant_id']),r['eye']))[0]
  ex=next(e for e in rr if e['participant_id']==chosen['participant_id'] and e['eye']==chosen['eye'])
  example.append(dict(example=f'Example {t}',**chosen))
  write(f'results/example_{t}_points.tsv',[dict(example=f'Example {t}',week=week,actual_day=ex['assigned'][week]['day'],letters=ex['assigned'][week]['va'],retained8=week in A.CALENDARS[8],retained12=week in A.CALENDARS[12]) for week in A.CALENDARS[4]])
  # Frozen point set, continuous error and every threshold on the same observations.
  idx={(e['participant_id'],e['eye']):i for i,e in enumerate(rr)}
  for k in [8,12]:
   hp=[p for p in points if p['trial']==t and p['calendar']==str(k)];check(f'{t}{k} hidden count',len(hp)==(n*6 if k==8 else n*8))
   eyeix=[idx[p['participant_id'],p['eye']] for p in hp];wp=w[:,eyeix];obs=np.array([float(p['actual_va']) for p in hp]);pred=np.array([float(p['predicted_va']) for p in hp]);d=pred-obs;ad=abs(d)
   meta=dict(trial=t,calendar=k,eyes=n,people=len(bb.ids),points=len(hp),weighting='pooled_points')
   for name,v,rep,unit in [('mean_signed',d.mean(),mean(d,wp),'letters'),('mean_abs',ad.mean(),mean(ad,wp),'letters'),('p90_abs',np.quantile(ad,.9),q(ad,wp),'letters')]+[(f'abs_ge{lim}',np.mean(ad>=lim),mean(ad>=lim,wp),'proportion') for lim in [5,10]]:
    held.append(row(meta,name,v,rep,unit,numerator=int(sum(ad>=int(name[6:]))) if name.startswith('abs_ge') else '',denominator=len(hp)))
   # all eyes have equal number of hidden points: pooled mean/proportion equals eye-equal mean/proportion.
   for h in C['thresholds']:
    mismatch=(obs>=h)!=(pred>=h);low=obs<h;high=obs>=h
    for metric,num,den in [('overall_disagreement',mismatch,np.ones(len(hp),bool)),('missed_below_threshold',(pred>=h)&low,low),('missed_attainment',(pred<h)&high,high)]:
     nd=int(sum(den));nn=int(sum(num));bwden=wp@den;rep=np.divide(wp@num,bwden,out=np.full(B,np.nan),where=bwden>0)
     holdthreshold.append(row(dict(trial=t,calendar=k,threshold=h,eyes=n,people=len(bb.ids),points=len(hp)) ,metric,nn/nd if nd else math.nan,rep,'proportion',numerator=nn,denominator=nd))
  # release large point-weight arrays before the next trial
  del wp
 # Historical estimates and limits must reproduce, not just round similarly.
 oldmain=read(R/'phase1/results/main_results.tsv');oldco=read(R/'phase1/results/treatment_contrasts.tsv')
 for r in endpoint:
  if r['endpoint']=='TIR_days' and r['metric'] in ['mean_abs','mean_signed','p90_abs']:
   src=next(x for x in oldmain if x['trial']==r['trial'] and x['analysis_set']=='strict' and int(x['calendar'])==r['calendar'] and x['metric']==r['metric'] and x['unit']=='days')
   for col in ['estimate','ci_low','ci_high']:check(f'historical {r["trial"]}{r["calendar"]}{r["metric"]} {col}',abs(r[col]-float(src[col]))<TOL)
 for r in contrasts:
  if r['endpoint']=='TIR_days':
   metric=('contrast_' if r['metric']=='adjusted_contrast' else 'contrast_change_')+str(r['calendar']);src=next(x for x in oldco if x['trial']==r['trial'] and x['analysis_set']=='strict' and x['metric']==metric and x['unit']=='weeks')
   for col in ['estimate','ci_low','ci_high']:check(f'historical treatment {r["trial"]}{metric} {col}',abs(r[col]/7-float(src[col]))<TOL)
 check('Negative control exact zero',all(r['estimate']==0 for r in endpoint if r['endpoint']=='anniversary_VA_change' and r['metric'] in ['mean_signed','mean_abs','median_abs','p90_abs']))
 for name,rs in [('ENDPOINT_COMPARISON.tsv',endpoint),('results/ENDPOINT_TREATMENT_CONTRASTS.tsv',contrasts),('results/CROSS_ENDPOINT_COMPARISON.tsv',cross),('MULTI_THRESHOLD_RESULTS.tsv',multi),('results/THRESHOLD_FLOOR_CEILING.tsv',thresholdclass),('CLINICAL_SUBGROUP_RESULTS.tsv',subgroups),('results/CLINICAL_DISCREPANCY_SCALES.tsv',clinical),('CONTINUOUS_HOLDOUT_VALIDATION.tsv',held),('results/MULTI_THRESHOLD_HOLDOUT.tsv',holdthreshold),('results/PER_EYE_EXTENSION.tsv',pereyes),('results/TRAJECTORY_COUNTS.tsv',trajectories),('qa/EXAMPLE_SELECTION_KEY_LOCAL_ONLY.tsv',example),('qa/NUMERIC_CHECKS.tsv',CHECKS),('qa/BOOTSTRAP_LOG.tsv',BOOT)]:write(name,rs)
 outputs=endpoint+contrasts+cross+multi+clinical+held
 failures=[dict(trial=r['trial'],metric=r['metric'],endpoint=r.get('endpoint',''),calendar=r.get('calendar',''),failed=r['bootstrap_failed']) for r in outputs if r['bootstrap_failed']]
 (P/'qa/ANALYSIS_STATUS.json').write_text(json.dumps(dict(completed=datetime.datetime.now().isoformat(),checks=len(CHECKS),all_checks_passed=all(c['passed'] for c in CHECKS),full_population_bootstrap_failures=failures,bootstrap_replicates=B,config_sha256=hashlib.sha256((P/'RUN_CONFIG_EXTENSION.json').read_bytes()).hexdigest()),indent=2))
 print('Extension numerical analysis complete; checks',len(CHECKS),'full-population failed intervals',len(failures),flush=True)
if __name__=='__main__':main()
