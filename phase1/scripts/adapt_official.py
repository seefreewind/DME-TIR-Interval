"""Official T/AC -> normalized eligibility inputs. No derived visual outcomes.
Method: EVA-ETDRS only; preserve raw day, official protocol baseline [-8,0] anchored at0.
"""
import csv,json,pathlib,datetime,collections,argparse,sys
ROOT=pathlib.Path(__file__).resolve().parents[2];RUN=pathlib.Path(__file__).resolve().parents[1]
def read(p):
 with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='|'))
def write(p,rr):
 if not rr:return
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def date(s):return datetime.datetime.strptime(s.split(' ')[0],'%m/%d/%Y').date()
def number(s):
 try:return float(s)
 except (TypeError,ValueError):return None
def unique(rows,field):
 assert len(rows)==len({r[field] for r in rows}),(field,'nonunique')
 return {r[field]:r for r in rows}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 eyes=[];visits=[];checks={};filters=[];mapping=[];schemas=[]
 for t,folder,base in [('T','Data Tables - Text files','Baseline'),('AC','Data Tables','Randomization')]:
  d=RUN/'extracted'/t/folder
  if not d.exists():checks[t]={'status':'BLOCKED_BY_ACCESS'};continue
  for p in sorted(d.glob('*.txt')):
   if p.name.startswith('._'):continue
   rr=read(p);schemas.append(dict(trial=t,table=p.name,rows=len(rr),columns='|'.join(rr[0]) if rr else ''))
  roster=unique(read(d/'PtRoster.txt'),'PtID'); rand=unique(read(d/(t+'Randomization.txt')),'PtID')
  assert set(roster)==set(rand)
  va=read(d/(t+'VAReVATest.txt'));vi=read(d/(t+'VisitInfo.txt'))
  unique(va,'RecID');unique(vi,'ParentLoginVisitID')
  vimap={(r['PtID'],r['ParentLoginVisitID']):r for r in vi};assert len(vimap)==len(vi)
  badlinks=[r['RecID'] for r in va if (r['PtID'],r['ParentLoginVisitID']) not in vimap];assert not badlinks
  for r in va:assert r['Visit']==vimap[r['PtID'],r['ParentLoginVisitID']]['Visit']
  bases=unique([r for r in va if r['Visit']==base],'PtID');assert set(bases)==set(roster)
  groupmap=unique(read(d/'TTrtGroupUnmasked.txt'),'PtID') if t=='T' else None
  origin={};trialeyes=[]
  for pid,r in roster.items():
   rr=rand[pid];rv=vimap[pid,rr['ParentLoginVisitID']];assert rv['Visit']==base
   assert bases[pid]['ParentLoginVisitID']==rr['ParentLoginVisitID']
   if t=='T':
    origin[pid]=int(rv['VisitDtDaysFromEnroll']);assert origin[pid]==0
    assert r['StudyEye']==rr['RandEye']==groupmap[pid]['RandEye'];ee=[rr['RandEye']]
   else:
    origin[pid]=date(r['RandDt']);assert date(r['EnrollDt'])==origin[pid]==date(r['VisitSchedStartDt'])==date(rv['VisitDt'])
    assert r['StudyEye']==rr['ACRandEye'];ee=['OD','OS'] if r['StudyEye']=='OU' else [r['StudyEye']]
   for eye in ee:
    group=groupmap[pid]['TTrtGroup'] if t=='T' else r['TrtGroup'+eye]
    assert group in (['Aflibercept','Bevacizumab','Ranibizumab'] if t=='T' else ['AC-A','AC-B'])
    b=bases[pid];rawday=int(b['VAExamDtDaysFromEnroll'])-origin[pid] if t=='T' else (date(b['VAExamDt'])-origin[pid]).days
    value=number(b['Etdrs'+eye]);valid=(b['VAMethod'+eye]=='EVA-ETDRS' and b['VACompleted'] in ['Both eyes','Right eye only' if eye=='OD' else 'Left eye only'] and value is not None and 0<=value<=100 and -8<=rawday<=0)
    e=dict(trial=t,participant_id=pid,eye=eye,group=group,randomized=1,baseline_va=value if valid else '',baseline_actual_day=rawday,bilateral=len(ee)==2)
    eyes.append(e);trialeyes.append(e)
  for r in va:
   s=r['Visit'];week=0 if s==base else (int(s.split()[0]) if s.lower() in [f'{w} week' for w in range(4,53,4)] else None)
   if week is None:continue
   pid=r['PtID'];rawday=(number(r['VAExamDtDaysFromEnroll'])-origin[pid] if number(r['VAExamDtDaysFromEnroll']) is not None else None) if t=='T' else ((date(r['VAExamDt'])-origin[pid]).days if r['VAExamDt'] else None)
   for e in (e for e in trialeyes if e['participant_id']==pid):
    eye=e['eye'];v=number(r['Etdrs'+eye]);method=r['VAMethod'+eye]
    valid=v is not None and 0<=v<=100 and method=='EVA-ETDRS' and r['VACompleted'] in ['Both eyes','Right eye only' if eye=='OD' else 'Left eye only']
    if week==0:valid=valid and rawday is not None and -8<=rawday<=0
    else:valid=valid and rawday is not None and rawday>0
    reason='valid' if valid else ('non_EVA_ETDRS' if method!='EVA-ETDRS' else 'missing_invalid_or_uncompleted')
    filters.append(dict(trial=t,participant_id=pid,eye=eye,week=week,method=method,valid=valid,reason=reason))
    visitday=number(vimap[pid,r['ParentLoginVisitID']]['VisitDtDaysFromEnroll']) if t=='T' else (date(vimap[pid,r['ParentLoginVisitID']]['VisitDt'])-origin[pid]).days
    visits.append(dict(trial=t,participant_id=pid,eye=eye,week=week,day=(0 if week==0 and valid else rawday),va=(v if valid else ''),measurement_id=t+':'+r['RecID']+':'+eye,raw_actual_day=rawday,raw_table=t+'VAReVATest.txt',raw_recid=r['RecID'],visit_id=r['ParentLoginVisitID'],raw_visit=s,method=method,visit_actual_day=visitday,visit_missed=vimap[pid,r['ParentLoginVisitID']].get('VisitMiss') or '',baseline_anchored=week==0 and valid and rawday!=0))
  counts=collections.Counter(e['group'] for e in trialeyes)
  if t=='T':assert counts=={'Aflibercept':224,'Bevacizumab':218,'Ranibizumab':218}
  else:
   assert counts=={'AC-A':158,'AC-B':154}
   bi=[pid for pid,r in roster.items() if r['StudyEye']=='OU'];assert len(bi)==42
   assert all({r['TrtGroupOD'],r['TrtGroupOS']}=={'AC-A','AC-B'} for pid,r in roster.items() if pid in bi)
  checks[t]=dict(people=len(roster),eyes=len(trialeyes),groups=dict(counts),randomization_to_visit_links='ALL_MATCH',VA_to_visit_links='ALL_MATCH',baseline_offsets=dict(collections.Counter(e['baseline_actual_day'] for e in trialeyes)),time_origin=('randomization linked baseline VisitDtDaysFromEnroll=0' if t=='T' else 'date(VAExamDt)-date(PtRoster.RandDt); common deidentification offset cancels'),bilateral_people=(42 if t=='AC' else 0))
  entries=[('participant_id','PtRoster.txt','PtID','PtID','试验内受试者唯一ID'),('eye',t+'Randomization.txt','RandEye' if t=='T' else 'ACRandEye','PtID','与PtRoster.StudyEye逐人一致；OU展开两眼'),('random_group','TTrtGroupUnmasked.txt' if t=='T' else 'PtRoster.txt','TTrtGroup' if t=='T' else 'TrtGroupOD/TrtGroupOS','PtID+eye','保留随机分组，禁止使用当前药物/换药重分类'),('planned_visit',t+'VAReVATest.txt','Visit','PtID+ParentLoginVisitID','Baseline/Randomization→0；4 week…52 week；不按行号'),('letters',t+'VAReVATest.txt','EtdrsOD/EtdrsOS','PtID+ParentLoginVisitID+eye','按研究眼选择0–100数值'),('method',t+'VAReVATest.txt','VAMethodOD/VAMethodOS;VACompleted','PtID+ParentLoginVisitID+eye','仅EVA-ETDRS且该眼已完成；ETDRS Chart不混入主计数'),('exam_time',t+'VAReVATest.txt','VAExamDtDaysFromEnroll' if t=='T' else 'VAExamDt','PtID+ParentLoginVisitID','actual detection day, not visit/injection day'),('randomization_origin',t+'VisitInfo.txt' if t=='T' else 'PtRoster.txt','VisitDtDaysFromEnroll' if t=='T' else 'RandDt','PtID+randomization.ParentLoginVisitID' if t=='T' else 'PtID','T随机化表链接基线VisitInfo=0；AC减随机化日历日期'),('missing',t+'VisitInfo.txt','VisitMiss/VisitMissReason','PtID+ParentLoginVisitID','原始缺访单独保留；不转为人为隐藏'),('duplicate',t+'VAReVATest.txt','RecID;ParentLoginVisitID','PtID+Visit+eye','计划周重复检查；实测ID唯一'),('baseline',t+'VAReVATest.txt','Visit;EtdrsOD/EtdrsOS;VAExamDtDaysFromEnroll' if t=='T' else 'Visit;EtdrsOD/EtdrsOS;VAExamDt','PtID+eye','协议允许随机化前8天基线；原始负日保留，日0锚点显式记录')]
  for concept,table,field,key,rule in entries:mapping.append(dict(trial=t,concept=concept,raw_table=table,raw_columns=field,join_key=key,time_unit=('days_from_enrollment' if t=='T' else 'deidentified_calendar_date_or_datetime') if concept in ['exam_time','randomization_origin','baseline'] else 'not_applicable',normalized_time_unit='days_from_randomization' if concept in ['exam_time','randomization_origin','baseline'] else 'not_applicable',rule=rule,dictionary='TReadMe.docx' if t=='T' else 'DataGlossary.rtf'))
 write(a.out/'eyes.tsv',eyes);write(a.out/'visits.tsv',visits);write(a.out/'record_validity.tsv',filters);write(a.out/'field_mapping_real.tsv',mapping);write(a.out/'raw_table_schema.tsv',schemas)
 (a.out/'mapping_evidence.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2));print(json.dumps(checks,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
