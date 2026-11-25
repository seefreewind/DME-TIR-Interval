from pathlib import Path
import sys,csv,json,hashlib,zipfile,numpy as np
R=Path(__file__).resolve().parents[2];P=R/'phase2_ajo_revision';sys.path.insert(0,str(P/'scripts'));import analyze_revision as E
checks=[]
def ck(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)));assert ok,(name,detail)
for r in E.rd(P/'sources/BASELINE_INPUT_VERIFICATION.tsv'):
 with zipfile.ZipFile(R/'data/raw'/r['trial']/('Protocol_'+r['trial']+'.zip')) as z:ck('Raw baseline member '+r['member'],hashlib.sha256(z.read(r['member'])).hexdigest()==r['sha256']==hashlib.sha256(Path(r['file']).read_bytes()).hexdigest())
ck('Baseline rebuilt unchanged',E.rd(P/'results/BASELINE_RECORDS_LOCAL_ONLY.tsv')==E.rd(R/'phase2_ajo/results/BASELINE_RECORDS_LOCAL_ONLY.tsv'))
base=E.rd(P/'results/BASELINE_RECORDS_LOCAL_ONLY.tsv');eyes={(r['trial'],r['participant_id'],r['eye']):r for r in E.rd(R/'phase1/normalized/eyes.tsv')};per=E.rd(R/'phase2_ajo/results/PER_EYE_EXTENSION.tsv');pe={(r['trial'],r['participant_id'],r['eye']):r for r in per};members=E.rd(P/'results/SELECTION_EYE_RECORDS_LOCAL_ONLY.tsv');boot=E.rd(P/'results/SELECTION_BOOTSTRAP_LOG.tsv');sres=E.rd(P/'results/SELECTION_WEIGHTING_RESULTS.tsv')
for r in base:r['bilateral']=eyes[r['trial'],r['participant_id'],r['eye']]['bilateral']
base.sort(key=lambda r:(r['trial'],int(r['participant_id']),r['eye']))
for t in ['T','AC']:
 rr=[r for r in base if r['trial']==t];y=np.array([r['included']=='True' for r in rr],float);sel=y==1;design=E.model_setup(rr,t);fit=E.fit_weights(rr,design,y,np.ones(len(rr)));bb=E.A.boot(rr)
 pp=[pe[t,r['participant_id'],r['eye']] for r in rr if r['included']=='True'];d=np.array([[abs(float(r[f'TIR{k}_days'])-float(r['TIR4_days'])) for k in [8,12]] for r in pp])
 ck(t+' original count',len(rr)==(660 if t=='T' else 312));ck(t+' strict count',int(y.sum())==(401 if t=='T' else 157))
 for rep in [0,23,4999]:
  # Literal cluster expansion independently validates frequency-weight refitting and medians.
  ix=np.repeat(np.arange(len(rr)),bb.w[rep].astype(int));expanded=[rr[i] for i in ix]
  # Unique IDs remain duplicated by resampling; medians must use participant multiplicities, not dedup across draws.
  m=bb.w[rep];f=E.fit_weights(rr,design,y,m);X,_,_=design(m);exX=X[ix];exy=y[ix];beta,prob,ok,it=E.logistic(exX,exy,np.ones(len(ix)))
  ck(t+' literal bootstrap coefficients '+str(rep),np.max(abs(beta-f['beta']))<1e-5)
  ads=d[np.repeat(np.arange(len(d)),m[sel].astype(int))];raw=(exy.mean())/prob[exy==1];lo,hi=np.quantile(raw,[.01,.99]);literal=np.average(ads,weights=np.clip(raw,lo,hi),axis=0)
  row=next(r for r in boot if r['trial']==t and int(r['replicate'])==rep+1)
  for j,k in enumerate([8,12]):ck(t+' literal weighted MAE '+str(rep)+str(k),abs(literal[j]-float(row['trim'+str(k)]))<1e-6)
 for r in [z for z in sres if z['trial']==t]:
  key=('trim_delta' if r['weighting']=='trim' else 'raw_delta') if r['metric']=='change_from_unweighted' else r['weighting'];vals=np.array([float(z[key+r['calendar']]) for z in boot if z['trial']==t]);lims=np.quantile(vals,[.025,.975]);ck('Weighting CI '+t+r['calendar']+key,np.max(abs(lims-np.array([float(r['ci_low']),float(r['ci_high'])])))<1e-10)
 old=[r for r in E.rd(R/'phase2_ajo/ENDPOINT_COMPARISON.tsv') if r['trial']==t and r['endpoint']=='TIR_days' and r['metric']=='mean_abs']
 for r in old:ck('Historical MAE intact '+t+r['calendar'],abs(float(r['estimate'])-np.mean(d[:,0 if r['calendar']=='8' else 1]))<1e-10)
common=E.rd(P/'results/COMMON_STANDARDIZATION_QUALIFICATION.tsv');ck('All common bins qualified',all(r['qualified']=='True' for r in common))
std=E.rd(P/'results/COMMON_STANDARDIZED_RESULTS.tsv');bs=E.rd(P/'results/COMMON_BOOTSTRAP_ESTIMATES.tsv')
for r in std:
 vals=np.array([float(z[r['metric']+r['calendar']]) for z in bs if z['trial']==r['trial']]);ck('Standardization CI '+r['trial']+r['calendar']+r['metric'],np.max(abs(np.quantile(vals,[.025,.975])-np.array([float(r['ci_low']),float(r['ci_high'])])))<1e-10)
critical=[]
for name,met,tn,dn in [('CLINICAL_DISCREPANCY_SCALES','absolute_TIR_difference_ge',[106,58],[401,157]),('CONTINUOUS_HOLDOUT_VALIDATION','abs_ge5',[481,196],[2406,942])]:
 for i,t in enumerate(['T','AC']):
  rs=E.rd(R/'phase2_ajo'/('results/'+name+'.tsv' if name=='CLINICAL_DISCREPANCY_SCALES' else name+'.tsv'));r=next(r for r in rs if r['trial']==t and r['calendar']=='8' and r['metric']==met and (name!='CLINICAL_DISCREPANCY_SCALES' or r['threshold_days']=='14'));ck('Core count '+name+t,int(r['numerator'])==tn[i] and int(r['denominator'])==dn[i]);critical.append(r)
E.wr('qa/CORE_RESULTS_AUDIT.tsv',critical);E.wr('qa/INDEPENDENT_CHECKS.tsv',checks);(P/'qa/AUDIT_STATUS.json').write_text(json.dumps(dict(passed=True,checks=len(checks),historical_input_files='checked by public runner final-table digests'),indent=2));print(len(checks),'checks passed')
