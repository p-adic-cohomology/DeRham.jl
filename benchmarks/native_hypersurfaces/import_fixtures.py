"""Import reproducible full-polynomial fixtures from workshop a5f8024.
Usage: python3 benchmarks/native_hypersurfaces/import_fixtures.py WORKSHOP
No workshop or Python is required by the installed calculator.
"""
import ast,itertools,json,random,re,subprocess,sys
from pathlib import Path
work=Path(sys.argv[1]);root=Path(__file__).resolve().parents[2]
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=work,text=True).strip()=='a5f802454f425ed6eae2aea197a6b3fee65ddba8'
def support(n,d):
 return [list(e) for e in itertools.product(range(d+1),repeat=n) if sum(e)==d]
def rows(file):return [json.loads(x) for x in (work/file).read_text().splitlines()]
fixtures=[]
def add(n,d,p,seed,P,source):
 s=support(n,d);rng=random.Random(seed);c=[rng.randrange(1,p) for _ in s]
 fixtures.append(dict(id=f'native_dwork_n{n}_p{p}',n=n,d=d,p=p,support=s,coefficients=c,expected=P,source=source,seed=seed))
text=(work/'experiments/check_threefold_transport.txt').read_text()
for p,seed,poly in re.findall(r'p =\s*(\d+).*?seed (\d+).*?\n\s*P\(T\) = (\[[^\n]+\])',text):
 if int(p) in (13,19,23):add(5,3,int(p),int(seed),ast.literal_eval(poly),'experiments/check_threefold_transport.txt')
for p in (17,23):
 row=next(r for r in rows('experiments/check_fourfold_window_panels.jsonl') if r['p']==p and not r.get('jacobi'))
 add(6,3,p,p,row['polynomial'],'experiments/check_fourfold_window_panels.jsonl')
for p in (89,167):
 row=next(r for r in rows('experiments/check_fourfold_divided_fibres.jsonl') if r['prime']==p and 'polynomial' in r)
 add(6,3,p,row.get('replacement_seed') or p,row['polynomial'],'experiments/check_fourfold_divided_fibres.jsonl')
row=next(r for r in rows('experiments/check_dense_fourfold_socle.jsonl') if r['prime']==11 and r['terms']==56)
add(6,3,11,125,row['polynomial'],'experiments/check_dense_fourfold_socle.jsonl')
row=next(r for r in rows('experiments/check_window_prime_range.jsonl') if r['p']==151)
add(4,4,151,row['seed'],row['polynomial'],'experiments/check_window_prime_range.jsonl')
# The generic generator samples the mixed support even when requesting every term.
rng=random.Random(row['seed']);s=support(4,4);rng.sample([e for e in s if max(e)<4],31)
fixtures[-1]['coefficients']=[rng.randrange(1,151) for _ in s]
# For the extended boundary, expand the archived sparse K3 in its recorded frame.
p=251;row=next(r for r in rows('experiments/check_window_prime_range.jsonl') if r['p']==p and r['d']==4 and r['case']=='same surface in new coordinates')
rng=random.Random(8129+p)
def invertible(A):
 A=[r[:] for r in A]
 for j in range(4):
  i=next((i for i in range(j,4) if A[i][j]%p),None)
  if i is None:return False
  A[j],A[i]=A[i],A[j];v=pow(A[j][j],-1,p);A[j]=[x*v%p for x in A[j]]
  for i in range(j+1,4):
   v=A[i][j];A[i]=[(x-v*y)%p for x,y in zip(A[i],A[j])]
 return True
while True:
 frame=[[rng.randrange(p) for _ in range(4)] for _ in range(4)]
 if invertible(frame):break
def mul(a,b):
 out={}
 for e,c in a.items():
  for f,d in b.items():
   t=tuple(x+y for x,y in zip(e,f));out[t]=(out.get(t,0)+c*d)%p
 return {e:c for e,c in out.items() if c}
linear=[{tuple(int(i==j) for i in range(4)):frame[j][k] for j in range(4)} for k in range(4)]
f={};product={(0,0,0,0):2}
for line in linear:
 power={(0,0,0,0):1}
 for _ in range(4):power=mul(power,line)
 for e,c in power.items():f[e]=(f.get(e,0)+c)%p
 product=mul(product,line)
for e,c in product.items():f[e]=(f.get(e,0)+c)%p
add(4,4,p,8129+p,row['polynomial'],'experiments/check_window_prime_range.jsonl; changed_coordinates')
fixtures[-1]['coefficients']=[f[tuple(e)] for e in fixtures[-1]['support']]
assert all(all(r['coefficients']) for r in fixtures)
# Keep one canonical JSON fixture source in the E2E schema.
cases=[]
for r in fixtures:
 cases.append(dict(id=r['id'],variety=dict(coeff_domain={'kind':'integer'},dim=r['n']-2,
  non_middle_factors={'kind':'projective_lefschetz','middle_factor_content':'primitive'},
  model=dict(kind='projective_hypersurface',vars=[f'x{i}' for i in range(r['n'])],monomials=r['support'],coeffs=r['coefficients'])),
  results=[dict(p=r['p'],Lpoly={'coeffs_asc':list(map(str,r['expected']))})],
  notes=f"Workshop a5f802454f425ed6eae2aea197a6b3fee65ddba8; {r['source']}; seed {r['seed']}"))
(root/'test/e2e/examples/native_hypersurfaces.json').write_text(json.dumps({'schema_version':'3','cases':cases},indent=2)+'\n')
# Julia fixture data avoids adding a JSON runtime dependency.
out=['const NATIVE_HYPERSURFACE_FIXTURES = [']
for r in fixtures:
 out.append('(id='+json.dumps(r['id'])+f', n={r["n"]}, d={r["d"]}, p={r["p"]}, support='+repr(r['support'])+', coefficients='+repr(r['coefficients'])+', expected=parse.(BigInt,'+json.dumps(list(map(str,r['expected'])))+')),')
out.extend([']','function native_fixture(case)',' R, xs = polynomial_ring(GF(case.p), ["x$i" for i in 1:case.n])',' return sum(R(c)*prod(xs.^e) for (e,c) in zip(case.support,case.coefficients))','end'])
(root/'benchmarks/native_hypersurfaces/fixtures.jl').write_text('\n'.join(out)+'\n')
print('Imported',len(fixtures),'archived fixtures')
