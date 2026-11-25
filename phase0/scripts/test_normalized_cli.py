"""Re-run normalized synthetic CLI integration; preserves all earlier audit runs."""
from pathlib import Path
import subprocess,sys,csv,json
p=Path(__file__).resolve().parents[1];d=p/'tests/normalized_fixture'
i=1
while (d/f'audit_run_v{i}').exists():i+=1
out=d/f'audit_run_v{i}'
subprocess.run([sys.executable,str(p/'scripts/audit_normalized.py'),'--eyes',str(d/'eyes.tsv'),'--visits',str(d/'visits.tsv'),'--mapping-evidence',str(d/'mapping_SYNTHETIC_ONLY.json'),'--out',str(out)],check=True)
r=next(r for r in csv.DictReader((out/'paired_flow.tsv').open(),delimiter='\t') if r['group']=='ALL')
expected={'randomized_people':'2','randomized_eyes':'3','strict_8':'3','strict_12':'3','strict_8_holdout_points':'18','strict_12_holdout_points':'24','strict_8_no_tail':'2','annual_357to363_eyes':'1'}
checks=[dict(check=k,expected=v,actual=r[k],passed=r[k]==v) for k,v in expected.items()]
(out/'integration_checks.json').write_text(json.dumps(checks,indent=2))
assert all(x['passed'] for x in checks)
print('Synthetic normalized CLI integration: 8/8 passed; '+str(out))
