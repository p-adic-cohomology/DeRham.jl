"""Separate native-kernel changes from Julia wrapper cost for the small threefold."""
import json,os,statistics,subprocess,sys,tempfile,time
from pathlib import Path
root=Path(__file__).resolve().parents[2];workshop=Path(sys.argv[1])
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=workshop,text=True).strip()=='a5f802454f425ed6eae2aea197a6b3fee65ddba8'
case=next(r for r in json.loads((root/'test/e2e/examples/native_hypersurfaces.json').read_text())['cases'] if r['id']=='native_dwork_n5_p13')
model=case['variety']['model']
data='13 5 35\n'+'\n'.join(' '.join(map(str,(*e,c))) for e,c in zip(model['monomials'],model['coeffs']))+'\n'
with tempfile.TemporaryDirectory() as temp:
 original=Path(temp)/'original'
 subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17','-DNV=5',str(workshop/'dwork/native/dense_frobenius.cpp'),'-o',str(original)],check=True)
 paths={'workshop':original,'shipping':root/'deps/usr/bin/dense_dwork_5'}
 samples={name:[] for name in paths};diagnostics={};expected=None
 for repetition in range(8):
  for name in (list(paths) if repetition%2==0 else list(reversed(paths))):
   start=time.perf_counter();result=subprocess.run([str(paths[name])],input=data,text=True,capture_output=True,check=True);elapsed=time.perf_counter()-start
   if expected is None:expected=result.stdout
   assert result.stdout==expected,'native matrices differ'
   if repetition:samples[name].append(elapsed)
   diagnostics[name]=result.stderr
 for name in paths:
  print(name,'matrix_median_seconds='+str(statistics.median(samples[name])),'samples='+str(samples[name]))
  print(diagnostics[name])
