from pathlib import Path
import sys,csv,hashlib,json,collections,math,datetime
import numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent

def rd(p,sep='\t'):return list(csv.DictReader(Path(p).open(),delimiter=sep))
def wr(p,rr):
 keys=list(dict.fromkeys(k for r in rr for k in r))
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def number(x):
 try:return float(x)
 except:return float('nan')
def valid(x):return np.isfinite(x)
def checked(x):return str(x).strip().lower() in ['true','1','yes']
def unique(rr,key):
 d={}
 for r in rr:
  k=key(r);assert k not in d,(k,'ambiguous baseline form');d[k]=r
 return d
# Written before reading baseline values; mapping does not depend on group summaries.
maprows=[]
for t,h,l,a,b in [('T','THxBaseline','TLab','AgeAsOfEnrollDt','Baseline'),('AC','ACBaselineOcuMedHx','ACLab','AgeAsofEnrollDt','Randomization')]:
 fields=[('age','PtRoster',a,'participant','years','Numeric exact age only; censored strings missing'),('sex',h,'Gender' if t=='T' else 'Sex','participant','category','M/F as reported'),('diabetes_type',h,'DiabetesType','participant','category','Type 1/Type 2/Uncertain as reported'),('diabetes_duration_approx',h,'DiabAge + PtRoster.'+a,'participant','years','Age minus age at diagnosis; numeric and 0<=diagnosis<=age; no substitution for censored ages'),('HbA1c',l,'HgbA1c; HbA1cLabVal','participant','percent','Central numeric result unless HgbA1cNA checked; otherwise local 4..19 unless HbA1cNotDone checked; source retained'),('baseline_VA','phase1/normalized/eyes.tsv','baseline_va','eye','E-ETDRS letters','Frozen valid baseline unchanged'),('prior_DME_treatment',h,'DMETrtHxOD / DMETrtHxOS','eye','Yes/No','Use exact eye-specific Yes/No; other missing'),('CST','RCOCTGradingResult','OCTCentSubCalc; OCTCentSubCalcStatus; OCTMachine','eye','micrometers native by device','Baseline visit only; positive <=2000 and not Ungradable; use only one unique (machine,value) pair across all gradable records; discordant duplicate pairs missing; no cross-device conversion or pooling')]
 for v,tb,col,unit,measure,rule in fields:maprows.append(dict(trial=t,variable=v,table=tb,columns=col,join_key='PtID'+(' + OCTScanEye' if v=='CST' else ' + study eye suffix' if v=='prior_DME_treatment' else ''),baseline_visit=b,analysis_unit=unit,measurement_unit=measure,rule=rule,dictionary='TReadMe.txt' if t=='T' else 'DataGlossary.txt'))
mp=P/'FIELD_MAPPING_BASELINE.tsv'
if mp.exists():assert rd(mp)==maprows
else:wr(mp,maprows)
(P/'sources/BASELINE_MAPPING_FREEZE.json').exists() or (P/'sources/BASELINE_MAPPING_FREEZE.json').write_text(json.dumps(dict(timestamp=datetime.datetime.now().astimezone().isoformat(),mapping_sha256=hashlib.sha256(mp.read_bytes()).hexdigest(),note='Rules fixed before first baseline descriptive summaries; duplicate handling is conservative, not trial adjudication.'),indent=2))
manifest=rd(R/'phase1/input_member_verification.tsv');expected={r['fresh_file']:r for r in manifest}
inputs=[];records=[];ambiguities=[]
eyes=rd(R/'phase1/normalized/eyes.tsv');elig={(r['trial'],r['participant_id'],r['eye']):r for r in rd(R/'phase1/eligibility/eye_eligibility_real.tsv')}
for t,sub,h,l,agecol,visit in [('T','Data Tables - Text files','THxBaseline','TLab','AgeAsOfEnrollDt','Baseline'),('AC','Data Tables','ACBaselineOcuMedHx','ACLab','AgeAsofEnrollDt','Randomization')]:
 d=R/'phase1/extracted'/t/sub
 data={}
 for name in ['PtRoster',h,l,'RCOCTGradingResult']:
  p=d/(name+'.txt');sha=hashlib.sha256(p.read_bytes()).hexdigest();er=expected[str(p)];assert sha==er['sha256'];inputs.append(dict(trial=t,file=str(p),sha256=sha,member=er['member'],zip_sha256=er['zip_sha256'],matches_verified_zip_member=True));data[name]=rd(p,'|')
 roster=unique(data['PtRoster'],lambda r:r['PtID']);hist=unique([r for r in data[h] if r['Visit']==visit],lambda r:r['PtID']);lab=unique([r for r in data[l] if r['Visit']==visit],lambda r:r['PtID']);octs=collections.defaultdict(list)
 for r in data['RCOCTGradingResult']:
  if r['Visit']==visit:octs[r['PtID'],r['OCTScanEye']].append(r)
 for e in [x for x in eyes if x['trial']==t]:
  p=e['participant_id'];hh=hist[p];ll=lab[p];age=number(roster[p][agecol]);diag=number(hh['DiabAge']);dur=age-diag if valid(age) and valid(diag) and 0<=diag<=age else math.nan
  cen=number(ll.get('HgbA1c'));local=number(ll.get('HbA1cLabVal'));a1c=cen if valid(cen) and not checked(ll.get('HgbA1cNA')) else local if valid(local) and 4<=local<=19 and not checked(ll.get('HbA1cNotDone')) else math.nan
  source='central' if valid(a1c) and valid(cen) and a1c==cen else 'local' if valid(a1c) else 'missing'
  oo=octs[p,e['eye']];pairs={(r['OCTMachine'],number(r['OCTCentSubCalc'])) for r in oo if r['OCTCentSubCalcStatus']!='Ungradable' and valid(number(r['OCTCentSubCalc'])) and 0<number(r['OCTCentSubCalc'])<=2000 and r['OCTMachine']}
  machine,cst=next(iter(pairs)) if len(pairs)==1 else ('',math.nan)
  ambiguities.append(dict(trial=t,participant_id=p,eye=e['eye'],baseline_oct_rows=len(oo),unique_gradable_pairs=len(pairs),status='accepted_consensus' if len(pairs)==1 else 'discordant_duplicates' if len(pairs)>1 else 'no_gradable_CST'))
  records.append(dict(trial=t,participant_id=p,eye=e['eye'],group=e['group'],included=elig[t,p,e['eye']]['strict_8']=='True',age=age,sex=hh.get('Gender',hh.get('Sex','')),diabetes_type=hh['DiabetesType'],diabetes_duration_approx=dur,HbA1c=a1c,HbA1c_source=source,baseline_VA=number(e['baseline_va']),prior_DME_treatment=hh.get('DMETrtHx'+e['eye'],''),CST=cst,OCT_device=machine))
wr(P/'sources/BASELINE_INPUT_VERIFICATION.tsv',inputs);wr(P/'results/BASELINE_RECORDS_LOCAL_ONLY.tsv',records);wr(P/'qa/BASELINE_OCT_AMBIGUITY.tsv',ambiguities)
output=[];membership=[]
for t in ['T','AC']:
 rr=[r for r in records if r['trial']==t];by=collections.defaultdict(list)
 for r in rr:by[r['participant_id']].append(r)
 persons=[]
 for pid,z in by.items():
  q=z[0].copy();q['included']=any(x['included'] for x in z);persons.append(q)
  membership.append(dict(trial=t,participant_id=pid,included_eye_count=sum(x['included'] for x in z),excluded_eye_count=sum(not x['included'] for x in z),participant_group='included_at_least_one_eye' if q['included'] else 'excluded_all_eyes'))
 for var,unit,cat in [('age','participant',False),('sex','participant',True),('diabetes_type','participant',True),('diabetes_duration_approx','participant',False),('HbA1c','participant',False),('HbA1c_source','participant',True),('baseline_VA','eye',False),('prior_DME_treatment','eye',True),('OCT_device','eye',True)]+[('CST|'+dev,'eye',False) for dev in sorted({r['OCT_device'] for r in rr if r['OCT_device']})]:
  src=persons if unit=='participant' else rr
  if var.startswith('CST|'):dev=var.split('|')[1];src=[r for r in src if r['OCT_device']==dev];field='CST'
  else:field=var
  subsets={'original':src,'included':[r for r in src if r['included']],'excluded':[r for r in src if not r['included']]}
  levels=sorted({r[field] for r in src if r[field] and (field!='prior_DME_treatment' or r[field] in ['Yes','No'])}) if cat else ['continuous']
  for level in levels:
   stats={}
   for pop,ss in subsets.items():
    vals=[r[field] for r in ss if (r[field] and (field!='prior_DME_treatment' or r[field] in ['Yes','No']))] if cat else [r[field] for r in ss if valid(r[field])]
    num=sum(v==level for v in vals) if cat else len(vals);mean=num/len(vals) if cat and vals else np.mean(vals) if vals else math.nan;sd=math.sqrt(mean*(1-mean)) if cat and vals else np.std(vals,ddof=1) if len(vals)>1 else math.nan
    stats[pop]=dict(n_total=len(ss),n_valid=len(vals),n_missing=len(ss)-len(vals),n_category=num if cat else '',estimate=float(mean),sd=float(sd))
   den=math.sqrt((stats['included']['sd']**2+stats['excluded']['sd']**2)/2);smd=(stats['included']['estimate']-stats['excluded']['estimate'])/den if valid(den) and den>0 else math.nan
   for pop,st in stats.items():output.append(dict(label='POST_PHASE1_EXPLORATORY_EXTENSION',trial=t,variable=var,level=level,analysis_unit=unit,population=pop,**st,SMD_included_minus_excluded=smd,participant_selection='At least one strict eye versus all eyes excluded; disjoint participants' if unit=='participant' else 'Eye inclusion; no inferential baseline P values',CST_missing_note='Device-specific denominator excludes unassigned/ambiguous CST; see overall OCT_device missing' if var.startswith('CST|') else ''))
wr(P/'BASELINE_SELECTION_AUDIT.tsv',output);wr(P/'results/BASELINE_PARTICIPANT_MEMBERSHIP.tsv',membership)
print('Baseline audit complete',len(output),'rows',collections.Counter((r['trial'],r['status']) for r in ambiguities))
