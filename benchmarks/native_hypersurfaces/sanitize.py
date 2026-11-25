"""ASan/UBSan on every shape's largest accepted prime, plus refusals.
python3 benchmarks/native_hypersurfaces/sanitize.py
"""
import itertools,json,math,os,random,subprocess,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
cases=json.loads((root/'test/e2e/examples/native_hypersurfaces.json').read_text())['cases']
def payload(p,N,s,c):return f'{p} {N} {len(s)}\n'+'\n'.join(' '.join(map(str,(*e,a))) for e,a in zip(s,c))+'\n'
with tempfile.TemporaryDirectory() as temp:
 for n,p,N in [(4,251,4),(5,13,5),(5,23,5),(6,167,5),(6,13,5)]:
  exe=Path(temp)/f'sanitize_{n}'
  subprocess.run([os.environ.get('CXX','c++'),'-O1','-g','-std=c++17',f'-DNV={n}','-fsanitize=address,undefined','-fno-sanitize-recover=all',str(root/'deps/native/dense_frobenius.cpp'),'-o',str(exe)],check=True)
  case=next(c for c in cases if c['id']==f'native_dwork_n{n}_p{11 if n==6 and p==13 else p}')
  s=case['variety']['model']['monomials'];c=case['variety']['model']['coeffs']
  if n==6 and p==13:
   rng=random.Random(113);c=[rng.randrange(1,p) for _ in s]
  result=subprocess.run([str(exe)],input=payload(p,N,s,c),capture_output=True,text=True,timeout=900)
  assert result.returncode==0,result.stderr
  rank={4:21,5:10,6:22}[n];lines=result.stdout.splitlines()
  assert len(lines)==2*rank+1 and int(lines[0])==rank
  assert all(0<=int(x)<p**N for line in lines[rank+1:] for x in line.split())
  print(f'PASS sanitizer n={n} p={p} N={N}',flush=True)
  for q,digits in [(7,N),(257 if n==4 else 29 if n==5 else 173,N),(p,6)]:
   result=subprocess.run([str(exe)],input=payload(q,digits,s,[1]*len(s)),capture_output=True,text=True)
   assert result.returncode==2,result.stderr
   print(f'PASS refused n={n} p={q} N={digits}',flush=True)
  d=4 if n==4 else 3
  singular=[math.factorial(d)//math.prod(math.factorial(x) for x in e)%11 for e in s]
  result=subprocess.run([str(exe)],input=payload(11,N,s,singular),capture_output=True,text=True)
  assert result.returncode==3,result.stderr
  print(f'PASS native singular refusal n={n}',flush=True)
