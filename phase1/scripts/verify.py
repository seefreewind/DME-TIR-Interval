"""Independent breakpoint integration, bootstrap multiplicity and output cross-checks."""
import pathlib,csv,sys,json,collections,hashlib
import numpy as np
import analyze as A
P=pathlib.Path(__file__).resolve().parents[1];O=P/'verification';O.mkdir(exist_ok=True)
def rd(p):return A.rd(p)
def write(name,rr):
 with (O/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
# Alternative integration: split each observed interval at threshold crossing,
# evaluate midpoint on each subsegment, add the lengths of attaining pieces.
def independent(pts,end=364,method='linear',inclusive=True):
 tt=np.array([p[0] for p in pts],float);yy=np.array([p[1] for p in pts],float)
 if tt[-1]<end:
  assert end-tt[-1]<=7
  tt=np.append(tt,end);yy=np.append(yy,yy[-1])
 cuts={0.,float(end)}|{float(t) for t in tt if 0<t<end}
 if method=='linear':
  for a,b,x,y in zip(tt,tt[1:],yy,yy[1:]):
   if (x-69)*(y-69)<0:
    cross=a+(69-x)*(b-a)/(y-x)
    if 0<cross<end:cuts.add(float(cross))
 cs=sorted(cuts);total=0.
 for a,b in zip(cs,cs[1:]):
  mid=(a+b)/2
  val=np.interp(mid,tt,yy) if method=='linear' else yy[np.searchsorted(tt,mid,side='right')-1]
  total+=(b-a)*((val>=69) if inclusive else (val>69))
 return total
bucket=collections.defaultdict(list)
for r in rd(P/'normalized/visits.tsv'):
 for c in ['week','day','va','raw_actual_day']:r[c]=float(r[c]) if r[c] else None
 bucket[r['trial'],r['participant_id'],r['eye']].append(r)
eyes=rd(P/'normalized/eyes.tsv');elig={(r['trial'],r['participant_id'],r['eye']):r for r in rd(P/'eligibility/eye_eligibility_real.tsv')}
rows=[]
for e in eyes:
 key=e['trial'],e['participant_id'],e['eye'];e['assigned'],_=A.assign_planned(bucket[key]);e['elig']=elig[key]
 if e['elig']['strict_8']=='True':rows.append(e)
rows.sort(key=lambda r:(r['trial'],int(r['participant_id']),r['eye']))
per={(r['trial'],r['participant_id'],r['eye']):r for r in rd(P/'results/per_eye_strict.tsv')}
checks=[]
for e in rows:
 for k in [4,8,12]:
  pts=[(e['assigned'][w]['day'],e['assigned'][w]['va']) for w in A.CALENDARS[k]]
  expected=independent(pts);actual=float(per[e['trial'],e['participant_id'],e['eye']][f'T{k}_days']);err=abs(expected-actual)
  checks.append(dict(trial=e['trial'],participant_id=e['participant_id'],eye=e['eye'],calendar=k,expected_days=expected,observed_days=actual,absolute_error=err,passed=err<1e-8))
  assert err<1e-8
  for method,inc in [('locf',True),('linear',False)]:assert abs(independent(pts,method=method,inclusive=inc)-A.integral(pts,method=method,inclusive=inc))<1e-8
write('independent_integration.tsv',checks)
# Bootstrap empirical quantiles and weighted ANCOVA vs literal participant duplication.
unit=[]
for t in ['T','AC']:
 rr=[r for r in rows if r['trial']==t];bb=A.boot(rr);y=np.array([float(per[r['trial'],r['participant_id'],r['eye']]['D8_days']) for r in rr]);mean=bb.mean(y);q=bb.quant(abs(y),.9)
 for j in [0,1,2,10,4999]:
  literal=np.repeat(np.arange(len(rr)),bb.w[j].astype(int));assert abs(mean[j]-y[literal].mean())<1e-10;assert abs(q[j]-np.quantile(abs(y[literal]),.9,method='linear'))<1e-10
  unit.append(dict(check='literal_cluster_repeat_mean_and_quantile',trial=t,replicate=j,passed=True))
 # Two-arm T or AC original-strategy baseline + original bilateral covariate.
 rr=[r for r in rr if t=='AC' or r['group']!='Ranibizumab'];bb=A.boot(rr);x=np.column_stack([np.ones(len(rr)),[r['group'] in ['Aflibercept','AC-A'] for r in rr],[float(r['baseline_va'])-69 for r in rr]])
 if t=='AC':x=np.column_stack([x,[r['bilateral']=='True' for r in rr]])
 y=np.array([[float(per[r['trial'],r['participant_id'],r['eye']][f'T{k}_days']) for k in [4,8,12]] for r in rr])
 for j in [0,1,10,4999]:
  literal=np.repeat(np.arange(len(rr)),bb.w[j].astype(int));literal_co=np.linalg.lstsq(x[literal],y[literal],rcond=None)[0];weighted=np.linalg.solve(x.T@(bb.w[j,:,None]*x),x.T@(bb.w[j,:,None]*y));assert np.max(abs(weighted-literal_co))<1e-8
  unit.append(dict(check='literal_cluster_repeat_ANCOVA',trial=t,replicate=j,passed=True))
write('bootstrap_literal_checks.tsv',unit)
# Wide per-eye outputs, eligibility sensitivity grid, and qualified masked score counts.
wide=[];grid=[]
for e in eyes:
 for k in [8,12]:
  if e['elig'][f'wide_{k}']=='True':
   a=A.outcome(e,4);b=A.outcome(e,k);wide.append(dict(trial=e['trial'],participant_id=e['participant_id'],eye=e['eye'],group=e['group'],calendar=k,T4_days=a,T_sparse_days=b,D_days=b-a,T4_weeks=a/7,T_sparse_weeks=b/7,D_weeks=(b-a)/7))
for t in ['T','AC']:
 for k in [8,12]:
  for m in [12,13,14]:
   for gap in [56,63,84]:
    rr=[e for e in eyes if e['trial']==t and int(e['elig']['n_valid'])>=m and float(e['elig']['max_gap'])<=gap and set(A.CALENDARS[k])<=set(e['assigned']) and e['elig']['boundary'] in ['bracket_or_exact','tail_le7']]
    grid.append(dict(trial=t,calendar=k,min_points=m,max_gap_days=gap,eyes=len(rr),people=len({r['participant_id'] for r in rr})))
write('per_eye_wide.tsv',wide);write('qualification_grid.tsv',grid)
summary=dict(independent_real_integrations=len(checks),max_integration_error_days=max(r['absolute_error'] for r in checks),literal_bootstrap_checks=len(unit),all_passed=True)
(O/'verification_summary.json').write_text(json.dumps(summary,indent=2));print(summary)
