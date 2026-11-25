"""Phase 0: reusable calendar/eligibility logic; TIR entry point accepts synthetic only.
No official-file adapter exists until actual dictionaries are verified.
Python >=3.9, standard library only.
"""
import math, random
CALENDARS={4:list(range(0,53,4)),8:[0,8,16,24,32,40,48,52],12:[0,12,24,36,48,52]}

def finite(v):
    return isinstance(v,(int,float)) and math.isfinite(v)

def assign_planned(records):
    """Already dictionary-verified planned weeks only. Never infer from row position.
    Identical duplicates collapse; conflicting duplicate values/dates quarantine week.
    Missing day is allowed for a separate nominal-only analysis.
    """
    buckets={}; audit=[]
    for r in records:
        w=r['week']
        if w not in CALENDARS[4]:
            audit.append(('outside_calendar',r)); continue
        if not finite(r.get('va')) or not 0<=r['va']<=100:
            audit.append(('missing_or_invalid_va',r)); continue
        buckets.setdefault(w,[]).append(r)
    assigned={}
    for w,rows in buckets.items():
        vals={(r.get('day'),r['va']) for r in rows}
        if len(vals)!=1:
            audit.append(('conflicting_duplicate',w)); continue
        assigned[w]=rows[0]
        if len(rows)>1:audit.append(('identical_duplicates_collapsed',len(rows)-1))
    dated=[(w,r['day']) for w,r in sorted(assigned.items()) if finite(r.get('day'))]
    if any(b[1]<=a[1] for a,b in zip(dated,dated[1:])):
        audit.append(('nonmonotone_actual_dates',dated))
        return {},audit
    ids=[r.get('measurement_id') for r in assigned.values() if r.get('measurement_id') is not None]
    if len(ids)!=len(set(ids)):raise ValueError('measurement assigned more than once')
    return assigned,audit

def mask_calendar(assigned,k):
    out=[]
    for w in CALENDARS[4]:
        r=assigned.get(w)
        status='original_missing' if r is None else ('retained' if w in CALENDARS[k] else 'artificial_hidden')
        out.append({'week':w,'status':status,'record':r})
    return out

def boundary(assigned,allow_tail=True):
    if not all(w in assigned and finite(assigned[w].get('day')) for w in (0,48,52)):
        return 'missing_boundary_or_date'
    d0,d48,d52=(assigned[w]['day'] for w in (0,48,52))
    if d0!=0 or not 0<d48<d52:return 'invalid_order_or_baseline'
    if d48<=364<=d52:return 'bracket_or_exact'
    if allow_tail and 357<=d52<364:return 'tail_le7'
    return 'unsupported_end'

def eligibility(assigned,k,strict=True,min_points=12,max_gap=63,allow_tail=True):
    valid=sorted(assigned)
    if 0 not in assigned:return False,'missing_baseline',0
    if not all(finite(assigned[w].get('day')) for w in valid):return False,'nominal_only',0
    if len(valid)<(14 if strict else min_points):return False,'point_count',0
    days=[assigned[w]['day'] for w in valid]
    if any(b<=a for a,b in zip(days,days[1:])):return False,'date_order',0
    if not strict and max(b-a for a,b in zip(days,days[1:]))>max_gap:return False,'gap',0
    if not set(CALENDARS[k])<=set(valid):return False,'sparse_target_missing',0
    if boundary(assigned,allow_tail) not in ('bracket_or_exact','tail_le7'):return False,'boundary',0
    kept=[assigned[w]['day'] for w in CALENDARS[k]]
    held=[assigned[w] for w in valid if w not in CALENDARS[k] and 0<=assigned[w]['day']<=364 and min(kept)<assigned[w]['day']<max(kept)]
    if not held:return False,'no_holdout',0
    return True,'eligible',len(held)

def cluster_draw(rows,seed=20261124):
    ids=sorted({r['participant_id'] for r in rows}); rng=random.Random(seed); out=[]
    for j,pid in enumerate(rng.choices(ids,k=len(ids))):
        for r in rows:
            if r['participant_id']==pid:out.append(dict(r,bootstrap_cluster=j))
    return out

def tir(points,end=364,method='linear',inclusive=True,allow_tail=False,*,provenance):
    if provenance!='synthetic':raise PermissionError('PHASE1_AUTHORIZATION=FALSE: real-eye outcome computation disabled')
    pts=sorted(points)
    if not pts or pts[0][0]!=0 or any(not finite(t) or not finite(y) for t,y in pts):raise ValueError('invalid start/data')
    if any(b[0]<=a[0] for a,b in zip(pts,pts[1:])):raise ValueError('duplicate dates')
    if pts[-1][0]<end:
        if not allow_tail or end-pts[-1][0]>7:raise ValueError('unsupported end')
        pts.append((end,pts[-1][1]))
    ok=lambda y:y>=69 if inclusive else y>69
    total=0.
    for (t0,y0),(t1,y1) in zip(pts,pts[1:]):
        if t0>=end:break
        stop=min(t1,end); dt=stop-t0
        if method=='locf':total+=dt*ok(y0);continue
        if method!='linear':raise ValueError('unknown interpolation')
        yend=y0+(y1-y0)*(stop-t0)/(t1-t0)
        if y0==yend:total+=dt*ok(y0)
        elif min(y0,yend)>=69:total+=dt
        elif max(y0,yend)<=69:pass
        else:
            cross=(69-y0)/(yend-y0)
            total+=dt*((1-cross) if yend>y0 else cross)
    return total/7
