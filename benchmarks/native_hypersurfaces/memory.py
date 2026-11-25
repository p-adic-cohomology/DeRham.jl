"""Peak RSS of the shipping native child alone (one matrix, no Julia host).
python3 benchmarks/native_hypersurfaces/memory.py
"""
import json,os,platform,random,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[2]
rows=json.loads((root/'test/e2e/examples/native_hypersurfaces.json').read_text())['cases']
out=root/'benchmarks/native_hypersurfaces/results'
for n,p in [(5,13),(6,11),(6,17),(4,151),(6,13)]:
 case=next(r for r in rows if r['id']==f'native_dwork_n{n}_p{11 if n==6 and p==13 else p}')
 model=case['variety']['model'];support=model['monomials'];coefficients=model['coeffs']
 if n==6 and p==13:
  rng=random.Random(113);coefficients=[rng.randrange(1,p) for _ in support]
 N=4 if n==4 else 5
 data=f'{p} {N} {len(support)}\n'+'\n'.join(' '.join(map(str,(*e,c))) for e,c in zip(support,coefficients))+'\n'
 binary=root/'deps/usr/bin'/('dense_dwork' if n==4 else f'dense_dwork_{n}')
 run=subprocess.run(['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v',str(binary)],input=data,text=True,capture_output=True,env=dict(os.environ,OMP_NUM_THREADS='1'),check=True)
 rank={4:21,5:10,6:22}[n]
 assert len(run.stdout.splitlines())==2*rank+1
 (out/f'child-memory-n{n}-p{p}.txt').write_text(f'n={n} p={p} N={N} one_matrix=true\n'+run.stderr)
 print(f'PASS native child memory n={n} p={p}',flush=True)
