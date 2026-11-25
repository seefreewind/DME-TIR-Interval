"""Frozen AJO V2 exploratory analyses. Historical inputs are read-only."""
from pathlib import Path
import csv,json,sys,hashlib,math,collections,time
import numpy as np
R=Path(__file__).resolve().parents[2];P=R/'phase2_ajo_revision';O=P/'results'
sys.path.insert(0,str(R/'phase1/scripts'));import analyze as A
C=json.loads((P/'RUN_CONFIG_REVISION.json').read_text());B=C['bootstrap']['B'];A.B=B;A.SEED=C['bootstrap']['seed']
def rd(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def wr(name,rr):
 if not rr:return
 keys=list(dict.fromkeys(k for r in rr for k in r))
 with (P/name).open('w') as f:
  w=csv.DictWriter(f,keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def fq(x,w,p):
 ix=np.argsort(x,kind='stable');v=x[ix];cum=np.cumsum(w[ix]);n=cum[-1]
 if n<=0:return math.nan
 h=(n-1)*p;l=np.searchsorted(cum,math.floor(h),side='right');u=np.searchsorted(cum,math.ceil(h),side='right');return v[l]+(h-math.floor(h))*(v[u]-v[l])
def stat(meta,metric,v,bs,unit='days',**extra):
 bs=np.asarray(bs,float);good=np.isfinite(bs);lo,hi=np.quantile(bs[good],[.025,.975]) if good.any() else [math.nan]*2
 return dict(identity='POST_HOC_AJO_V2_EXTENSION',**meta,metric=metric,estimate=float(v),ci_low=float(lo),ci_high=float(hi),unit=unit,bootstrap_valid=int(good.sum()),bootstrap_failed=int((~good).sum()),**extra)
def quantrows(trial,what,x):return [dict(trial=trial,distribution=what,quantile=q,value=float(np.quantile(x,q))) for q in [0,.01,.05,.25,.5,.75,.95,.99,1]]
NUM=['age','baseline_VA','diabetes_duration_approx','HbA1c'];CAT=['sex','group','diabetes_type','prior_DME_treatment']
REF={'sex':'F','group':None,'diabetes_type':'Type 2','prior_DME_treatment':'No'}
def model_setup(rows,t):
 vals={name:np.array([float(r[name]) for r in rows]) for name in NUM};cats={k:np.array([r[k] if r[k] and (k!='prior_DME_treatment' or r[k] in ['Yes','No']) else 'Missing' for r in rows]) for k in CAT}
 levels={k:sorted(set(v)) for k,v in cats.items()};refs={**REF,'group':'Bevacizumab' if t=='T' else 'AC-B'}
 first=np.array([i for i,r in enumerate(rows) if r['participant_id'] not in {z['participant_id'] for z in rows[:i]}]);mis={k:~np.isfinite(v) for k,v in vals.items()}
 def design(m):
  cols=[np.ones(len(rows))];names=['Intercept'];imps={}
  for k in NUM:
   ix=first if k!='baseline_VA' else np.arange(len(rows));ix=ix[np.isfinite(vals[k][ix])];median=fq(vals[k][ix],m[ix],.5)
   if not np.isfinite(median):raise ValueError('No finite covariate '+k)
   imps[k]=median;x=np.where(mis[k],median,vals[k]);scale=10 if k!='HbA1c' else 1;cols.append((x-median)/scale);names.append(k)
   if mis[k].any():cols.append(mis[k].astype(float));names.append(k+'_missing')
  for k in CAT:
   for lev in levels[k]:
    if lev!=refs[k]:cols.append((cats[k]==lev).astype(float));names.append(k+'='+lev)
  if t=='AC':cols.append(np.array([r['bilateral']=='True' for r in rows],float));names.append('original_bilateral')
  return np.column_stack(cols),names,imps
 return design

def logistic(X,y,m):
 beta=np.zeros(X.shape[1]);q=np.sum(m*y)/sum(m);beta[0]=np.log(q/(1-q));pen=np.r_[0,np.full(X.shape[1]-1,1e-6)]
 def obj(b):
  z=X@b;return np.sum(m*(y*z-np.logaddexp(0,z)))-.5*np.sum(pen*b*b)
 for it in range(100):
  z=X@beta;p=1/(1+np.exp(-np.clip(z,-700,700)));grad=X.T@(m*(y-p))-pen*beta;H=X.T@((m*p*(1-p))[:,None]*X)+np.diag(pen)
  step=np.linalg.solve(H,grad)
  if np.max(abs(step))<1e-9 and np.max(abs(grad))<1e-7:return beta,p,True,it+1
  fac=1.;old=obj(beta)
  while fac>2**-25 and obj(beta+fac*step)<old-1e-10:fac*=.5
  beta+=fac*step
 return beta,1/(1+np.exp(-np.clip(X@beta,-700,700))),False,100

def fit_weights(rows,design,y,m):
 X,names,imp=design(m);be,p,ok,it=logistic(X,y,m);sel=y.astype(bool);q=np.sum(m*y)/sum(m);raw=q/np.maximum(p[sel],1e-300);cut=np.array([fq(raw,m[sel],v) for v in [.01,.99]]);w=np.clip(raw,*cut)
 return dict(X=X,names=names,imputation=imp,beta=be,p=p,converged=ok,iterations=it,raw=raw,trim=w,cutoffs=cut)

def selection(base,per):
 out=[];diag=[];balance=[];prob=[];bootlog=[];members=[];coef=[];imps=[];status=[]
 for t in ['T','AC']:
  rows=[r for r in base if r['trial']==t];y=np.array([r['included']=='True' for r in rows],float);sel=y==1;n=int(y.sum());design=model_setup(rows,t);m=np.ones(len(rows));f=fit_weights(rows,design,y,m);bb=A.boot(rows)
  lookup={(r['participant_id'],r['eye']):r for r in per if r['trial']==t};pr=[lookup[r['participant_id'],r['eye']] for r in rows if r['included']=='True'];d=np.array([[float(r[f'TIR{k}_days'])-float(r['TIR4_days']) for k in [8,12]] for r in pr]);ad=abs(d)
  point={k:np.average(ad,weights=f[k],axis=0) for k in ['raw','trim']};un=ad.mean(axis=0);rep={k:np.full((B,2),np.nan) for k in ['raw','trim','unweighted','raw_delta','trim_delta']}
  for b,m in enumerate(bb.w):
   err=''
   try:
    z=fit_weights(rows,design,y,m)
    if not z['converged']:raise ValueError('fit_not_converged')
    sm=m[sel];rep['unweighted'][b]=np.average(ad,weights=sm,axis=0)
    for k in ['raw','trim']:rep[k][b]=np.average(ad,weights=sm*z[k],axis=0);rep[k+'_delta'][b]=rep[k][b]-rep['unweighted'][b]
    assert np.all(np.isfinite(rep['trim'][b]))
   except Exception as e:err=type(e).__name__+':'+str(e)
   bootlog.append(dict(trial=t,replicate=b+1,success=not bool(err),error=err,population_sha256=bb.key,seed=20261124,**{name+str(k):rep[name][b,j] for name in rep for j,k in enumerate([8,12])}))
  for k,v in f['imputation'].items():imps.append(dict(trial=t,variable=k,median=v))
  for k,v in zip(f['names'],f['beta']):coef.append(dict(trial=t,term=k,coefficient=v))
  # Balance uses original eye population SD, including missing indicators.
  before=[];after=[]
  for j,name in enumerate(f['names'][1:],1):
   x=f['X'][:,j];mu=x.mean();sd=np.std(x,ddof=0)
   if sd<=1e-12:continue
   a=(x[sel].mean()-mu)/sd;c=(np.average(x[sel],weights=f['trim'])-mu)/sd;before.append(abs(a));after.append(abs(c));balance.append(dict(trial=t,term=name,original_mean=mu,selected_mean=x[sel].mean(),weighted_selected_mean=np.average(x[sel],weights=f['trim']),target_SD=sd,SMD_unweighted=a,SMD_weighted=c))
  ess=f['trim'].sum()**2/np.sum(f['trim']**2);valid=sum(r['success'] for r in bootlog if r['trial']==t)
  gates=dict(point_converged=f['converged'],positivity_min_all_ge05=bool(min(f['p'])>=.05),raw_max_median_le10=bool(max(f['raw'])/np.median(f['raw'])<=10),trim_ESS_gehalf=bool(ess>=.5*n),weighted_max_abs_SMD_le10=bool(max(after)<=.1),mean_abs_SMD_improved=bool(np.mean(after)<np.mean(before)+1e-8),bootstrap_all5000=valid==B)
  qualified=all(gates.values());status.append(dict(trial=t,status='QUALIFIED' if qualified else 'NOT_QUALIFIED',original_eyes=len(rows),selected_eyes=n,original_people=len(bb.ids),min_probability_all=min(f['p']),raw_weight_max_median=max(f['raw'])/np.median(f['raw']),ESS=ess,ESS_fraction=ess/n,max_abs_SMD_before=max(before),max_abs_SMD_after=max(after),mean_abs_SMD_before=np.mean(before),mean_abs_SMD_after=np.mean(after),bootstrap_valid=valid,bootstrap_failed=B-valid,failed_gates=';'.join(k for k,v in gates.items() if not v)))
  diag+=quantrows(t,'probability_all',f['p'])+quantrows(t,'probability_selected',f['p'][sel])+quantrows(t,'weight_stabilized_raw',f['raw'])+quantrows(t,'weight_stabilized_truncated',f['trim'])
  diag += [dict(trial=t,distribution='weight_truncation_cutoff',quantile=q,value=v) for q,v in zip([.01,.99],f['cutoffs'])]
  si=0
  for i,r in enumerate(rows):
   members.append(dict(trial=t,participant_id=r['participant_id'],eye=r['eye'],group=r['group'],original_bilateral=r['bilateral'],selected=bool(y[i]),probability=f['p'][i],raw_weight=f['raw'][si] if y[i] else '',truncated_weight=f['trim'][si] if y[i] else '',model_complete=True))
   if y[i]:si+=1
  for j,k in enumerate([8,12]):
   for model in ['unweighted','raw','trim']:
    v=un[j] if model=='unweighted' else point[model][j];out.append(stat(dict(trial=t,calendar=k,weighting=model,selected_eyes=n,selected_people=len({r['participant_id'] for r in pr}),qualification='not_applicable' if model=='unweighted' else ('QUALIFIED' if qualified else 'NOT_QUALIFIED')),'mean_abs',v,rep[model][:,j]))
    if model!='unweighted':out.append(stat(dict(trial=t,calendar=k,weighting=model,selected_eyes=n,selected_people=len({r['participant_id'] for r in pr}),qualification='QUALIFIED' if qualified else 'NOT_QUALIFIED'),'change_from_unweighted',v-un[j],rep[model+'_delta'][:,j]))
 for name,rr in [('SELECTION_WEIGHTING_RESULTS',out),('SELECTION_DIAGNOSTICS',diag),('SELECTION_BALANCE',balance),('SELECTION_QUALIFICATION',status),('SELECTION_EYE_RECORDS_LOCAL_ONLY',members),('SELECTION_COEFFICIENTS',coef),('SELECTION_IMPUTATION',imps),('SELECTION_BOOTSTRAP_LOG',bootlog)]:wr('results/'+name+'.tsv',rr)
 return status

def common(per):
 results=[];distribution=[];elig=[];arr={};bsarr={};bins={};bsbins={};binrows=[];std=[];comp=[];qual=[];stdboot=[]
 for t in ['T','AC']:
  rr=[r for r in per if r['trial']==t];bb=A.boot(rr);va=np.array([float(r['baseline_va']) for r in rr]);mask=(va>=24)&(va<=68);ix=np.where(mask)[0];v=va[ix];ww=bb.w[:,ix];d=np.array([[float(r[f'TIR{k}_days'])-float(r['TIR4_days']) for k in [8,12]] for r in rr])[ix];ad=abs(d);meta=dict(trial=t,range_letters='24-68',eyes=len(ix),people=len({rr[i]['participant_id'] for i in ix}))
  for i,r in enumerate(rr):elig.append(dict(trial=t,participant_id=r['participant_id'],eye=r['eye'],baseline_VA=va[i],common_range=mask[i],bin='24-53' if 24<=va[i]<=53 else '54-68' if 54<=va[i]<=68 else 'outside'))
  distribution.append(dict(**meta,minimum=min(v),q25=np.quantile(v,.25),median=np.median(v),q75=np.quantile(v,.75),maximum=max(v),mean=v.mean(),SD=np.std(v,ddof=1)))
  zmean=lambda z,w: (w@z)/w.sum(axis=1)
  for j,k in enumerate([8,12]):
   z=ad[:,j];cols=[('mean_abs',z.mean(),zmean(z,ww),'days'),('mean_signed',d[:,j].mean(),zmean(d[:,j],ww),'days'),('p90_abs',np.quantile(z,.9),np.array([fq(z,w,.9) for w in ww]),'days')]
   for h in [14,28]:cols.append((f'abs_ge{h}',np.mean(z>=h),zmean(z>=h,ww),'proportion'))
   for met,val,reps,unit in cols:results.append(stat(dict(**meta,calendar=k),met,val,reps,unit,numerator=int(sum(z>=int(met[6:]))) if met.startswith('abs_ge') else '',denominator=len(ix)))
  bins[t]=[];bsbins[t]=[];arr[t]=[];bsarr[t]=[]
  for low,high in C['common_support']['bins_inclusive']:
   ii=(v>=low)&(v<=high);n=int(ii.sum());pp=len({rr[ix[i]]['participant_id'] for i in np.where(ii)[0]});qualified=n>=20 and pp>=20;qual.append(dict(trial=t,bin=f'{low}-{high}',eyes=n,people=pp,qualified=qualified));bins[t].append(n/len(v));bsbins[t].append(ww[:,ii].sum(axis=1)/ww.sum(axis=1));vals=np.column_stack([ad[ii].mean(axis=0),d[ii].mean(axis=0)]).ravel();reps=np.column_stack([zmean(ad[ii,j],ww[:,ii]) for j in [0,1]]+[zmean(d[ii,j],ww[:,ii]) for j in [0,1]])
   # Explicit order: MAE8, MAE12, signed8, signed12.
   vals=np.r_[ad[ii].mean(axis=0),d[ii].mean(axis=0)];arr[t].append(vals);bsarr[t].append(reps)
   for pos,(metric,k) in enumerate([('mean_abs',8),('mean_abs',12),('mean_signed',8),('mean_signed',12)]):binrows.append(stat(dict(trial=t,bin=f'{low}-{high}',eyes=n,people=pp,calendar=k),metric,vals[pos],reps[:,pos]))
 good=all(r['qualified'] for r in qual);target=.5*(np.array(bins['T'])+np.array(bins['AC']));bstarget=.5*(np.array(bsbins['T'])+np.array(bsbins['AC']))
 wr('results/COMMON_STANDARDIZATION_QUALIFICATION.tsv',qual)
 wr('results/COMMON_STANDARDIZATION_TARGET.tsv',[dict(bin=f'{lo}-{hi}',target_proportion=target[i],T_proportion=bins['T'][i],AC_proportion=bins['AC'][i],status='QUALIFIED' if good else 'NOT_QUALIFIED') for i,(lo,hi) in enumerate(C['common_support']['bins_inclusive'])])
 if good:
  est={};reps={}
  for t in ['T','AC']:
   est[t]=target@np.array(arr[t]);reps[t]=np.einsum('jb,jbk->bk',bstarget,np.array(bsarr[t]))
   stdboot += [dict(trial=t,replicate=b+1,mean_abs8=reps[t][b,0],mean_abs12=reps[t][b,1],mean_signed8=reps[t][b,2],mean_signed12=reps[t][b,3],target_low=bstarget[0,b],target_high=bstarget[1,b]) for b in range(B)]
   for pos,(metric,k) in enumerate([('mean_abs',8),('mean_abs',12),('mean_signed',8),('mean_signed',12)]):std.append(stat(dict(trial=t,calendar=k,qualification='QUALIFIED'),metric,est[t][pos],reps[t][:,pos]))
  for pos,(metric,k) in enumerate([('mean_abs',8),('mean_abs',12),('mean_signed',8),('mean_signed',12)]):comp.append(stat(dict(trial='AC_minus_T',calendar=k,qualification='QUALIFIED'),metric,est['AC'][pos]-est['T'][pos],reps['AC'][:,pos]-reps['T'][:,pos]))
 for nm,z in [('COMMON_SUPPORT_RESULTS',results),('COMMON_BASELINE_DISTRIBUTION',distribution),('COMMON_EYE_QUALIFICATION_LOCAL_ONLY',elig),('COMMON_BIN_RESULTS',binrows),('COMMON_STANDARDIZED_RESULTS',std),('COMMON_STANDARDIZED_CONTRAST',comp)]:wr('results/'+nm+'.tsv',z)
 wr('results/COMMON_BOOTSTRAP_ESTIMATES.tsv',stdboot)
 return good

def examples(per):
 keys=[];pts=[]
 sys.path.insert(0,str(R/'phase2_ajo/scripts'));import analyze_extension as E
 rows=E.load()
 for t in ['T','AC']:
  rr=[r for r in per if r['trial']==t];md8=np.median([abs(float(r['TIR8_days'])-float(r['TIR4_days'])) for r in rr]);md12=np.median([abs(float(r['TIR12_days'])-float(r['TIR4_days'])) for r in rr])
  high=sorted(rr,key=lambda r:(-int(r['total_conversions']),-abs(float(r['TIR8_days'])-float(r['TIR4_days'])),int(r['participant_id']),r['eye']))[0]
  central=sorted(rr,key=lambda r:(abs(abs(float(r['TIR8_days'])-float(r['TIR4_days']))-md8),abs(abs(float(r['TIR12_days'])-float(r['TIR4_days']))-md12),int(r['participant_id']),r['eye']))[0]
  for kind,r in [('fluctuant',high),('near_median',central)]:
   label=f'{t} '+('fluctuant' if kind=='fluctuant' else 'near median');keys.append(dict(example=label,kind=kind,median_absD8=md8,median_absD12=md12,**r));e=next(e for e in rows if e['trial']==t and e['participant_id']==r['participant_id'] and e['eye']==r['eye'])
   for week in A.CALENDARS[4]:pts.append(dict(example=label,trial=t,kind=kind,week=week,actual_day=e['assigned'][week]['day'],letters=e['assigned'][week]['va'],retained8=week in A.CALENDARS[8],retained12=week in A.CALENDARS[12],TIR4_days=r['TIR4_days'],TIR8_days=r['TIR8_days'],TIR12_days=r['TIR12_days']))
 wr('qa/EXAMPLE_KEYS_LOCAL_ONLY.tsv',keys);wr('results/EXAMPLE_POINTS_DEIDENTIFIED.tsv',pts)

def main():
 start=time.time();baseline=rd(P/'results/BASELINE_RECORDS_LOCAL_ONLY.tsv');historical=rd(R/'phase2_ajo/results/BASELINE_RECORDS_LOCAL_ONLY.tsv');assert baseline==historical
 eyes={(r['trial'],r['participant_id'],r['eye']):r for r in rd(R/'phase1/normalized/eyes.tsv')}
 for r in baseline:r['bilateral']=eyes[r['trial'],r['participant_id'],r['eye']]['bilateral']
 baseline.sort(key=lambda r:(r['trial'],int(r['participant_id']),r['eye']));per=rd(R/'phase2_ajo/results/PER_EYE_EXTENSION.tsv')
 s=selection(baseline,per);g=common(per);examples(per)
 (P/'qa/ANALYSIS_STATUS.json').write_text(json.dumps(dict(elapsed_seconds=time.time()-start,selection_qualifications=s,common_standardization_qualified=g,config_sha256=hashlib.sha256((P/'RUN_CONFIG_REVISION.json').read_bytes()).hexdigest()),indent=2));print(json.dumps(s,indent=2));print('Common standardization:',g,'seconds',time.time()-start)
if __name__=='__main__':main()
