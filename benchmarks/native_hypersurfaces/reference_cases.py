"""Generate a fixed-seed native-reference cohort for the Julia public-API check.
Usage: sage -python reference_cases.py WORKSHOP OUTPUT [COUNT]
Builds the pinned, unchanged workshop kernel once per shape. Each output row
contains coefficients and the full primitive polynomial, or an explicit refusal.
"""
import itertools,os,random,subprocess,sys,tempfile
from pathlib import Path
root=Path(sys.argv[1]);output=Path(sys.argv[2]);count=int(sys.argv[3]) if len(sys.argv)>3 else 100
sys.path.insert(0,str(root/'dwork'))
import dense_native
from projective_zeta import primitive_polynomial
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()=='a5f802454f425ed6eae2aea197a6b3fee65ddba8'
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
with tempfile.TemporaryDirectory() as directory,output.open('w') as out:
 for n in (4,5,6):
  executable=Path(directory)/f'original_{n}'
  subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',f'-DNV={n}',str(root/'dwork/native/dense_frobenius.cpp'),'-o',str(executable)],check=True)
  dense_native.kernel=lambda dimension,exe=executable:exe
  p=13 if n==4 else 11;d=4 if n==4 else 3
  support=[e for e in itertools.product(range(d+1),repeat=n) if sum(e)==d]
  good=0;attempt=0
  while good<count and attempt<2*count:
   seed=2026112300+10000*n+attempt;attempt+=1
   rng=random.Random(seed);c=[rng.randrange(1,p) for _ in support]
   status='ambiguous';poly=None
   try:
    for N in (range(3,6) if n==4 else [5]):
     _,A,_=dense_native.dense_frobenius(p,support,c,N,measure_rss=False)
     poly=primitive_polynomial(A,p,N,n-2)
     if poly is not None:status='ok';good+=1;break
   except ValueError:status='singular'
   out.write('\t'.join([str(n),str(p),str(seed),status,','.join(map(str,c)),','.join(map(str,poly or []))])+'\n');out.flush()
  assert good==count,(n,good,attempt)
  print(f'n={n} smooth={good} attempts={attempt}',flush=True)
