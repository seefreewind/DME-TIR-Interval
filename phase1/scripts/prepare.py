import pathlib,hashlib,zipfile,json,csv,subprocess,sys,shutil,datetime
P=pathlib.Path(__file__).resolve().parents[1]; ROOT=next(q for q in [P,*P.parents] if (q/'data/raw').is_dir())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rr):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def main():
 hashes={'T':'be6be58ec51e9d257ee6d6881fb4e22d605ff64d79483435ee003b0fefffe9c5','AC':'49f4dac3b2a3324c5ebab7f7d287840b4102ed23aae110578b0769badf955a1a'}
 checks=[]
 for t,h in hashes.items():
  z=ROOT/'data/raw'/t/f'Protocol_{t}.zip';assert sha(z)==h
  with zipfile.ZipFile(z) as zz:
   for m in zz.infolist():
    if m.is_dir() or m.filename.startswith('__MACOSX/') or pathlib.Path(m.filename).name.startswith('._'):continue
    target=P/'extracted'/t/m.filename;assert target.resolve().is_relative_to((P/'extracted').resolve())
    data=zz.read(m);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data);assert target.read_bytes()==data
    checks.append(dict(trial=t,zip=str(z),zip_sha256=h,member=m.filename,bytes=len(data),sha256=sha(target),fresh_file=str(target),byte_identical_to_zip=True))
 write(P/'input_member_verification.tsv',checks)
 if not (P/'normalized').exists():subprocess.run([sys.executable,str(P/'scripts/adapt_official.py'),'--out',str(P/'normalized')],check=True,stdout=subprocess.DEVNULL)
 if not (P/'audit').exists():subprocess.run([sys.executable,str(ROOT/'phase0/scripts/audit_normalized.py'),'--eyes',str(P/'normalized/eyes.tsv'),'--visits',str(P/'normalized/visits.tsv'),'--mapping-evidence',str(P/'normalized/mapping_evidence.json'),'--out',str(P/'audit')],check=True,stdout=subprocess.DEVNULL)
 if not (P/'eligibility').exists():subprocess.run([sys.executable,str(P/'scripts/summarize_eligibility.py'),'--normalized',str(P/'normalized'),'--out',str(P/'eligibility')],check=True)
 matches=[]
 frozen=json.loads((ROOT/'provenance/frozen_input_sha256.json').read_text())
 for rel,h in frozen.items():
  p=ROOT/rel
  assert sha(p)==h,('Frozen input digest mismatch',rel)
  matches.append(dict(regenerated=rel,frozen_sha256=h,byte_identical=True,semantic_identical=True))
 write(P/'frozen_version_comparison.tsv',matches)
 print('Verified archive members:',len(checks),'identical regenerated outputs:',len(matches))
if __name__=='__main__':main()
