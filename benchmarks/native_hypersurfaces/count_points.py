"""Independent small-prime point counts for the new odd/even-dimensional fixtures."""
import json,os,subprocess,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
cases=json.loads((root/'test/e2e/examples/native_hypersurfaces.json').read_text())['cases']
with tempfile.TemporaryDirectory() as temp:
 exe=Path(temp)/'count'
 subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',str(root/'benchmarks/native_hypersurfaces/count_points.cpp'),'-o',str(exe)],check=True)
 for name in ['native_dwork_n5_p13','native_dwork_n5_p19','native_dwork_n6_p11','native_dwork_n6_p17']:
  case=next(c for c in cases if c['id']==name);model=case['variety']['model'];dim=case['variety']['dim'];n=dim+2;p=case['results'][0]['p'];s=model['monomials'];c=model['coeffs']
  data=f'{p} {n} 3 {len(s)}\n'+'\n'.join(' '.join(map(str,(*e,a))) for e,a in zip(s,c))+'\n'
  observed=int(subprocess.check_output([str(exe)],input=data,text=True))
  expected=sum(p**i for i in range(dim+1))-(-1)**dim*int(case['results'][0]['Lpoly']['coeffs_asc'][1])
  assert observed==expected,(name,observed,expected)
  print(name,'direct_points='+str(observed),'matches=true',flush=True)
