"""Matched serial benchmark of the pinned workshop direct native route, no dispatcher."""
import json,os,statistics,sys,time
from pathlib import Path
root,executable,case_id=Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3]
sys.path.insert(0,str(root/'dwork'))
import dense_native
from projective_zeta import primitive_polynomial,_assemble
dense_native.kernel=lambda n:executable
rows=json.loads((Path(__file__).resolve().parents[2]/'test/e2e/examples/native_hypersurfaces.json').read_text())['cases']
case=next(c for c in rows if c['id']==case_id);n=case['variety']['dim']+2;p=case['results'][0]['p']
s=case['variety']['model']['monomials'];c=case['variety']['model']['coeffs']
expected=list(map(int,case['results'][0]['Lpoly']['coeffs_asc']));zeta=_assemble(expected,p,0,n-2)
def solve():
 for N in (range(3,6 if p<=47 else 5) if n==4 else [5]):
  _,A,_=dense_native.dense_frobenius(p,s,c,N,measure_rss=False)
  poly=primitive_polynomial(A,p,N,n-2)
  if poly is not None:
   result=_assemble(poly,p,N,n-2)
   assert poly==expected and result.numerator==zeta.numerator and result.denominator==zeta.denominator
   return
 raise ArithmeticError('ambiguous lift')
start=time.perf_counter();solve();first=time.perf_counter()-start
print('first_solve_seconds='+str(first),flush=True)
if os.environ.get('NATIVE_COLD')=='1':sys.exit(0)
times=[]
for _ in range(5):
 start=time.perf_counter();solve();times.append(time.perf_counter()-start)
 print('warm_repetition_seconds='+str(times[-1]),flush=True)
print('median_seconds='+str(statistics.median(times)))
print('full_zeta_equal=true mode=workshop case='+case_id)
