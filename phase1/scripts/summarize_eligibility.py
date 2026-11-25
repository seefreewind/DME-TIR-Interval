"""Structural counts, boundary and rule losses only. Never calls tir or fits outcomes."""
import pathlib,csv,sys,collections,json,math,argparse
RUN=pathlib.Path(__file__).resolve().parents[1];ROOT=next(q for q in [RUN,*RUN.parents] if (q/'data/raw').is_dir())
sys.path.insert(0,str(ROOT/'phase0/scripts'))
from phase0_core import assign_planned,eligibility,boundary,CALENDARS
from audit_normalized import num,write

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--normalized',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 rd=lambda n:list(csv.DictReader((a.normalized/n).open(),delimiter='\t'))
 eyes=rd('eyes.tsv');vv=rd('visits.tsv');bucket=collections.defaultdict(list)
 for r in vv:
  for c in ['week','day','va','raw_actual_day']:r[c]=num(r[c])
  bucket[r['trial'],r['participant_id'],r['eye']].append(r)
 details=[];ruleflags=[];masks=[];holdouts=[]
 for e in eyes:
  key=e['trial'],e['participant_id'],e['eye'];raw=bucket[key];ass,au=assign_planned(raw)
  assert not any(x[0] in ['nonmonotone_actual_dates','conflicting_duplicate'] for x in au)
  valid=sorted(ass);days=[ass[w]['day'] for w in valid];actual_days=[ass[w]['raw_actual_day'] for w in valid];gap=max((v-u for u,v in zip(actual_days,actual_days[1:])),default=0)
  d={c:e[c] for c in ['trial','participant_id','eye','group']};d.update(linked=bool(raw),baseline=0 in ass,annual=52 in ass,complete14=len(ass)==14,n_valid=len(ass),baseline_actual_day=num(e['baseline_actual_day']),annual_day=ass.get(52,{}).get('day'),max_gap=gap,boundary=boundary(ass))
  for k in [8,12]:
   retained=[ass[w]['day'] for w in CALENDARS[k] if w in ass]
   held=[w for w in valid if w not in CALENDARS[k] and 0<=ass[w]['day']<=364 and retained and min(retained)<ass[w]['day']<max(retained)]
   for strict in [True,False]:
    name=('strict' if strict else 'wide')+'_'+str(k)
    rules=dict(linked=bool(raw),baseline=0 in ass,actual_dates_order=all(x is not None for x in days) and all(y>x for x,y in zip(days,days[1:])),point_count=len(ass)>=(14 if strict else 12),sparse_targets=set(CALENDARS[k])<=set(ass),boundary=boundary(ass) in ['bracket_or_exact','tail_le7'],holdout=len(held)>0)
    if not strict:rules['max_gap63']=gap<=63
    okay=all(rules.values());core_ok,reason,n=eligibility(ass,k,strict=strict);assert core_ok==okay
    d[name]=okay;d[name+'_points']=n if okay else 0;d[name+'_no_tail']=eligibility(ass,k,strict=strict,allow_tail=False)[0]
    for rule,passed in rules.items():ruleflags.append(dict(**{c:d[c] for c in ['trial','participant_id','eye','group']},analysis_set=name,rule=rule,passed=passed,only_remove_recovers=(not passed and all(v for kk,v in rules.items() if kk!=rule))))
   for w in CALENDARS[4]:masks.append(dict(**{c:d[c] for c in ['trial','participant_id','eye','group']},calendar=k,week=w,status=('original_missing_or_invalid' if w not in ass else 'retained' if w in CALENDARS[k] else 'artificial_hidden'),measurement_id=ass.get(w,{}).get('measurement_id','')))
   for w in held:holdouts.append(dict(**{c:d[c] for c in ['trial','participant_id','eye','group']},calendar=k,week=w,actual_day=ass[w]['day'],measurement_id=ass[w]['measurement_id'],strict_eligible=d['strict_'+str(k)],wide_eligible=d['wide_'+str(k)]))
  details.append(d)
 flow=[];boundaryrows=[];loss=[];gaps=[];sens=[]
 quant=lambda xs,q:sorted(xs)[int((len(xs)-1)*q)]+((len(xs)-1)*q%1)*(sorted(xs)[min(int((len(xs)-1)*q)+1,len(xs)-1)]-sorted(xs)[int((len(xs)-1)*q)]) if xs else None
 for t,g in sorted({(r['trial'],r['group']) for r in details}|{(r['trial'],'ALL') for r in details}):
  rr=[r for r in details if r['trial']==t and (g=='ALL' or r['group']==g)]
  metrics={'randomized':rr,'linked':[r for r in rr if r['linked']],'valid_baseline':[r for r in rr if r['baseline']],'annual_EVA':[r for r in rr if r['annual']],'complete14':[r for r in rr if r['complete14']]}
  for typ in ['strict','wide']:
   for k in [8,12]:
    name=typ+'_'+str(k);rows=[r for r in rr if r[name]];metrics[name]=rows;metrics[name+'_holdout']=rows
    flow.append(dict(trial=t,group=g,stage=name+'_holdout_points',people='',eyes='',points=sum(r[name+'_points'] for r in rr)))
    metrics[name+'_no_tail']=[r for r in rr if r[name+'_no_tail']]
   metrics[typ+'_common3']=[r for r in rr if r[typ+'_8'] and r[typ+'_12']]
  for name,rows in metrics.items():flow.append(dict(trial=t,group=g,stage=name,people=len({r['participant_id'] for r in rows}),eyes=len(rows),points=''))
  predicates={'before364':lambda d:d<364,'equal364':lambda d:d==364,'after364':lambda d:d>364,'357_363':lambda d:357<=d<=363,'308_356':lambda d:308<=d<=356,'365_420':lambda d:365<=d<=420,'outside308_420':lambda d:not 308<=d<=420}
  for label,fn in predicates.items():
   rows=[r for r in rr if r['annual_day'] is not None and fn(r['annual_day'])];boundaryrows.append(dict(trial=t,group=g,metric='annual_'+label,people=len({r['participant_id'] for r in rows}),eyes=len(rows)))
  for typ in ['strict','wide']:
   for k in [8,12]:
    name=typ+'_'+str(k);rows=[r for r in rr if r[name] and not r[name+'_no_tail']]
    boundaryrows.append(dict(trial=t,group=g,metric=name+'_lost_if_no_tail',people=len({r['participant_id'] for r in rows}),eyes=len(rows)))
    additional=[r for r in rr if r['wide_'+str(k)] and not r['strict_'+str(k)]];strict_only=[r for r in rr if r['strict_'+str(k)] and not r['wide_'+str(k)]]
    if typ=='wide':
     for label,sel in [('additional_wide_eyes',additional),('strict_only_eyes',strict_only)]:sens.append(dict(trial=t,group=g,calendar=k,metric=label,eyes=len(sel)))
    zero=[r for r in rr if r[name] and r['baseline_actual_day']==0];sens.append(dict(trial=t,group=g,calendar=k,metric=name+'_actual_day0_baseline_only',eyes=len(zero)))
  ff=[r for r in ruleflags if r['trial']==t and (g=='ALL' or r['group']==g)]
  for aset,rule in sorted({(r['analysis_set'],r['rule']) for r in ff}):
   ss=[r for r in ff if r['analysis_set']==aset and r['rule']==rule];failed=[r for r in ss if not r['passed']];recovered=[r for r in ss if r['only_remove_recovers']]
   loss.append(dict(trial=t,group=g,analysis_set=aset,rule=rule,standalone_failed_eyes=len(failed),standalone_failed_people=len({r['participant_id'] for r in failed}),only_remove_recovers_eyes=len(recovered),only_remove_recovers_people=len({r['participant_id'] for r in recovered})))
  xx=[r['max_gap'] for r in rr if r['n_valid']>=2]
  for label,val in [('n',len(xx)),('min',min(xx)),('median',quant(xx,.5)),('p90',quant(xx,.9)),('max',max(xx)),('n_gt63',sum(v>63 for v in xx))]:gaps.append(dict(trial=t,group=g,statistic=label,value=val))
 for name,rows in [('paired_flow_real.tsv',flow),('boundary_counts.tsv',boundaryrows),('rule_losses.tsv',loss),('max_gap_distribution.tsv',gaps),('baseline_and_wide_sensitivity.tsv',sens),('eye_eligibility_real.tsv',details),('rule_pass_flags.tsv',ruleflags),('calendar_masks.tsv',masks),('holdout_inventory.tsv',holdouts)]:write(a.out/name,rows)
 print('Structural summaries generated:',a.out)
if __name__=='__main__':main()
