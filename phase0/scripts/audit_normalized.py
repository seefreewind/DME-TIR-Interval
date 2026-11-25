"""Audit counts ONLY. Input is not an official DRCR schema.
Requires eyes.tsv (trial,participant_id,eye,group,randomized,baseline_va)
and visits.tsv (trial,participant_id,eye,week,day,va,measurement_id).
Do not use until mapping, visit codes and dates have been verified against raw dictionary.
No TIR, classification accuracy, treatment effect or P value is calculated here.
"""
import argparse,csv,json,pathlib,collections,statistics
from phase0_core import assign_planned,eligibility,boundary,CALENDARS,finite

def write(p,rows):
    if not rows:return
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def num(v):
    try:return float(v)
    except (TypeError,ValueError):return None

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--eyes',type=pathlib.Path,required=True);p.add_argument('--visits',type=pathlib.Path,required=True)
    p.add_argument('--mapping-evidence',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True)
    a=p.parse_args()
    if not a.mapping_evidence.is_file():p.error('actual mapping evidence required')
    if a.out.exists():p.error('use a fresh output directory to preserve prior files')
    eyes=list(csv.DictReader(a.eyes.open(),delimiter='\t'));visits=list(csv.DictReader(a.visits.open(),delimiter='\t'))
    key=lambda r:tuple(r[k] for k in ('trial','participant_id','eye'))
    if len({key(r) for r in eyes})!=len(eyes):raise ValueError('duplicate eye roster')
    bucket=collections.defaultdict(list)
    for r in visits:
        r['week']=num(r['week']);r['day']=num(r['day']);r['va']=num(r['va']);bucket[key(r)].append(r)
    details=[];issues=[];gaps=[];flags=[]
    for e in eyes:
        if e['randomized']!='1':continue
        rows=bucket[key(e)]; assigned,au=assign_planned(rows)
        for reason,detail in au:issues.append(dict(trial=e['trial'],participant_id=e['participant_id'],eye=e['eye'],reason=reason,detail=str(detail)))
        base=finite(num(e['baseline_va'])) and 0<=num(e['baseline_va'])<=100 and 0 in assigned
        if base and num(e['baseline_va'])!=assigned[0]['va']:base=False;issues.append(dict(trial=e['trial'],participant_id=e['participant_id'],eye=e['eye'],reason='baseline_mismatch',detail='quarantine'))
        d=dict(trial=e['trial'],participant_id=e['participant_id'],eye=e['eye'],group=e['group'],linked=bool(rows),baseline=base,n_valid=len(assigned),annual=52 in assigned,boundary=boundary(assigned),annual_day=assigned.get(52,{}).get('day'))
        days=sorted(r['day'] for r in assigned.values() if finite(r.get('day')))
        gap=max((v-u for u,v in zip(days,days[1:])),default=None)
        gaps.append(dict(trial=e['trial'],group=e['group'],participant_id=e['participant_id'],eye=e['eye'],max_actual_gap_days=gap))
        for k in (8,12):
            for strict in (True,False):
                label=f'{"strict" if strict else "wide"}_{k}'
                ok,reason,n=eligibility(assigned,k,strict=strict)
                d[label]=ok and base;d[label+'_holdout_points']=n if ok and base else 0
                d[label+'_no_tail']=eligibility(assigned,k,strict=strict,allow_tail=False)[0] and base
                if not ok:issues.append(dict(trial=e['trial'],participant_id=e['participant_id'],eye=e['eye'],reason=label+':'+reason,detail='sequential_first_failure'))
            for minimum in (12,13,14):
                for limit in (56,63,84):
                    flags.append(dict(trial=e['trial'],group=e['group'],eye_key='|'.join(key(e)),calendar=k,min_points=minimum,max_gap=limit,eligible=base and eligibility(assigned,k,strict=False,min_points=minimum,max_gap=limit)[0]))
        details.append(d)
    flow=[]
    for trial,group in sorted({(r['trial'],r['group']) for r in details}|{(r['trial'],'ALL') for r in details}):
        rr=[r for r in details if r['trial']==trial and (group=='ALL' or r['group']==group)]
        d=dict(trial=trial,group=group,randomized_people=len({r['participant_id'] for r in rr}),randomized_eyes=len(rr),linked_eyes=sum(r['linked'] for r in rr),baseline_eyes=sum(r['baseline'] for r in rr),annual_eyes=sum(r['annual'] for r in rr),complete14=sum(r['n_valid']==14 and r['baseline'] for r in rr))
        for typ in ('strict','wide'):
            for k in (8,12):
                label=f'{typ}_{k}';d[label]=sum(r[label] for r in rr);d[label+'_holdout_eyes']=sum(r[label] and r[label+'_holdout_points']>0 for r in rr);d[label+'_holdout_points']=sum(r[label+'_holdout_points'] for r in rr);d[label+'_no_tail']=sum(r[label+'_no_tail'] for r in rr)
            d[typ+'_common3']=sum(r[typ+'_8'] and r[typ+'_12'] for r in rr)
        for name,fn in [('before364',lambda x:x<364),('on364',lambda x:x==364),('after364',lambda x:x>364),('357to363',lambda x:357<=x<=363)]:
            dd=[r for r in rr if finite(r['annual_day']) and fn(r['annual_day'])];d['annual_'+name+'_eyes']=len(dd);d['annual_'+name+'_people']=len({r['participant_id'] for r in dd})
        flow.append(d)
    a.out.mkdir(parents=True)
    for name,rows in [('paired_flow',flow),('eye_eligibility',details),('exclusions',issues),('max_gaps',gaps),('tolerance_grid',flags)]:write(a.out/(name+'.tsv'),rows)
    print('Counts written; no real-eye TIR/effect computed.')
if __name__=='__main__':main()
