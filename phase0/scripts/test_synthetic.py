import csv,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).parent))
from phase0_core import *
OUT=pathlib.Path(__file__).resolve().parents[1]/'tests'; OUT.mkdir(exist_ok=True)
results=[]
def test(name,expected,fn):
    try:
        actual=fn(); good=math.isclose(expected,actual,abs_tol=1e-9) if isinstance(expected,(int,float)) else actual==expected
    except Exception as e:actual=repr(e);good=False
    results.append(dict(test=name,expected=expected,actual=actual,passed=good))
def T(p,**kw):return tir(p,provenance='synthetic',**kw)
test('constant_69_ge_linear_52weeks',52,lambda:T([(0,69),(364,69)]))
test('constant_69_gt_linear_zero',0,lambda:T([(0,69),(364,69)],inclusive=False))
test('constant_69_ge_locf',52,lambda:T([(0,69),(364,69)],method='locf'))
test('72_66_72_full',4,lambda:T([(0,72),(28,66),(56,72)],end=56))
test('72_66_72_hidden_week4',8,lambda:T([(0,72),(56,72)],end=56))
test('reverse_hidden_peak_full',4,lambda:T([(0,66),(28,72),(56,66)],end=56))
test('reverse_hidden_peak_sparse',0,lambda:T([(0,66),(56,66)],end=56))
test('irregular_dates_crossings_10days',10/7,lambda:T([(0,72),(10,66),(20,72)],end=20))
test('364_clip_using_late_annual',50,lambda:T([(0,72),(336,72),(392,60)]))
test('tail_from357_7days',52,lambda:T([(0,70),(336,70),(357,70)],allow_tail=True))
def reject(p,**kw):
    try:T(p,**kw);return False
    except ValueError:return True
test('no_tail_reject357',True,lambda:reject([(0,70),(336,70),(357,70)]))
test('tail8days_rejected',True,lambda:reject([(0,70),(356,70)],allow_tail=True))
def recs(missing=(),annual=364):
    return [dict(week=w,day=annual if w==52 else w*7,va=70,measurement_id=str(w)) for w in CALENDARS[4] if w not in missing]
def aligned():
    a,_=assign_planned(list(reversed(recs((4,)))))
    m=mask_calendar(a,8)
    return [r['week'] for r in m if r['status']=='retained']
test('missing4_does_not_shift8_calendar',CALENDARS[8],aligned)
def maskcounts():
    a,_=assign_planned(recs((4,)));m=mask_calendar(a,8)
    return [sum(r['status']==s for r in m) for s in ('retained','original_missing','artificial_hidden')]
test('missing_vs_hidden_separate',[8,1,5],maskcounts)
a,_=assign_planned(recs());test('strict8_6holdout',(True,'eligible',6),lambda:eligibility(a,8))
test('strict12_8holdout',(True,'eligible',8),lambda:eligibility(a,12))
b,_=assign_planned(recs((4,12)));test('wide12of14_gap56_8eligible',(True,'eligible',4),lambda:eligibility(b,8,strict=False))
c,_=assign_planned(recs((4,8)));test('wide_gap84_rejected',(False,'gap',0),lambda:eligibility(c,12,strict=False))
d,_=assign_planned(recs((8,)));test('missing_sparse_target',(False,'sparse_target_missing',0),lambda:eligibility(d,8,strict=False))
e,_=assign_planned(recs(annual=357));test('boundary357','tail_le7',lambda:boundary(e))
f,_=assign_planned(recs(annual=356));test('boundary356','unsupported_end',lambda:boundary(f))
g,_=assign_planned(recs(annual=370));test('boundary370','bracket_or_exact',lambda:boundary(g))
h,_=assign_planned(recs()+[dict(week=4,day=28,va=65)]);test('conflicting_duplicate_quarantined',False,lambda:4 in h)
i,au=assign_planned(recs()+[recs()[1]]);test('exact_duplicate_collapsed',14,lambda:len(i))
rows=[dict(participant_id='P1',eye='R',group='A',calendars=[4,8,12]),dict(participant_id='P1',eye='L',group='B',calendars=[4,8,12]),dict(participant_id='P2',eye='R',group='B',calendars=[4,8,12])]
def clustercheck():
    for seed in range(30):
        draw=cluster_draw(rows,seed)
        for cluster in set(r['bootstrap_cluster'] for r in draw):
            block=[r for r in draw if r['bootstrap_cluster']==cluster]
            if block[0]['participant_id']=='P1' and {(r['eye'],r['group'],tuple(r['calendars'])) for r in block}!={('R','A',(4,8,12)),('L','B',(4,8,12))}:return False
    return True
test('bilateral_crossgroup_cluster_intact',True,clustercheck)
def guard():
    try:tir([(0,70),(364,70)],provenance='real');return False
    except PermissionError:return True
test('real_outcome_guard',True,guard)
with (OUT/'synthetic_tests.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['test','expected','actual','passed'],delimiter='\t');w.writeheader();w.writerows(results)
(OUT/'synthetic_tests.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
print(f'{sum(r["passed"] for r in results)}/{len(results)} synthetic tests passed')
assert all(r['passed'] for r in results)
