"""Code-only reproduction; all individual-level outputs remain in a private local run."""
from pathlib import Path
import argparse,subprocess,sys,shutil,json,hashlib,os
R=Path(__file__).resolve().parent

def run(root,script,*args):
 print('Running',script,flush=True)
 subprocess.run([sys.executable,str(root/script),*map(str,args)],cwd=root,check=True)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--synthetic-only',action='store_true');p.add_argument('--formal-only',action='store_true');p.add_argument('--out',type=Path)
 a=p.parse_args()
 if a.synthetic_only:
  run(R,'phase0/scripts/test_synthetic.py');run(R,'phase0/scripts/test_normalized_cli.py');return
 if a.out is None:p.error('--out is required for real-data reproduction')
 out=a.out.resolve()
 if not out.is_relative_to(R) or out==R:p.error('--out must be a new directory inside this checkout')
 out.mkdir(parents=True,exist_ok=False)
 for d in ['phase0','phase0_real_2026-11-24','phase1','phase2_ajo','phase2_ajo_revision','provenance']:
  # Explicit source/config allowlist; never copy generated individual-level outputs.
  for src in (R/d).rglob('*'):
   rel=src.relative_to(R)
   if src.is_file() and not any(x.startswith('._') or x=='__pycache__' for x in rel.parts) and (src.suffix=='.py' or src.name.startswith('RUN_CONFIG') or src.name=='requirements.txt' or src.parent.name=='provenance' or src.name in ['ANALYSIS_RULES_FROZEN.md','PROTOCOL_ADDENDUM_2026-11-24.md','BASELINE_IMPLEMENTATION_NOTE_2026-11-24.md'] or str(rel).startswith('phase0/tests/normalized_fixture/') and 'audit_run' not in str(rel)):
    target=out/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,target)
 for t in ['T','AC']:
  source=R/'data/raw'/t/f'Protocol_{t}.zip'
  if not source.is_file():raise FileNotFoundError(f'Obtain official package separately: {source}')
  target=out/'data/raw'/t/source.name;target.parent.mkdir(parents=True,exist_ok=True);target.symlink_to(source)
 for d in ['phase2_ajo','phase2_ajo_revision']:
  for sub in ['results','qa','sources']:(out/d/sub).mkdir(parents=True,exist_ok=True)
 run(out,'phase0/scripts/test_synthetic.py');run(out,'phase0/scripts/test_normalized_cli.py')
 run(out,'phase1/scripts/prepare.py');run(out,'phase1/scripts/analyze.py');run(out,'phase1/scripts/verify.py');run(out,'phase1/scripts/supplement.py')
 if not a.formal_only:
  run(out,'phase2_ajo/scripts/baseline_audit.py');run(out,'phase2_ajo/scripts/analyze_extension.py');run(out,'phase2_ajo_revision/scripts/rebuild_baseline.py');run(out,'phase2_ajo_revision/scripts/analyze_revision.py');run(out,'phase2_ajo_revision/scripts/audit_revision.py')
 run(R,'plot_results.py','--run',out)
 expected=json.loads((R/'provenance/final_table_sha256.json').read_text());checks=[]
 for rel,h in expected.items():
  f=out/rel
  if f.is_file():checks.append(dict(file=rel,match=hashlib.sha256(f.read_bytes()).hexdigest()==h))
 status=dict(generated_tables_checked=len(checks),byte_identical=sum(x['match'] for x in checks),all_generated_match=all(x['match'] for x in checks),checks=checks,not_generated=[k for k in expected if not (out/k).is_file()],python=sys.version)
 (out/'release_check.json').write_text(json.dumps(status,indent=2)+'\n')
 print(json.dumps({k:v for k,v in status.items() if k not in ['checks','not_generated','python']},indent=2))
 if not status['all_generated_match']:raise AssertionError('Frozen-table digest mismatch; inspect release_check.json')
if __name__=='__main__':main()
