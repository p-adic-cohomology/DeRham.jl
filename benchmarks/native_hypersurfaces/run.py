"""Serial exact-output benchmark; compile outside timers; retain failures/timeouts.
python3 benchmarks/native_hypersurfaces/run.py WORKSHOP --python 'sage -python'
"""
import argparse,json,os,platform,shlex,signal,subprocess,tempfile,time
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('workshop',type=Path)
parser.add_argument('--python',default='python3')
parser.add_argument('--julia',default='julia +1.13.2')
parser.add_argument('--cases',nargs='+',default=['native_dwork_n5_p13','native_dwork_n6_p11','native_dwork_n6_p17','native_dwork_n4_p151'])
parser.add_argument('--modes',nargs='+',choices=['workshop','native','upstream-default','upstream-depthfirst'],default=['workshop','native'])
parser.add_argument('--timeout',type=int,default=600)
parser.add_argument('--cold',action='store_true')
parser.add_argument('--label',default='paired')
args=parser.parse_args();root=Path(__file__).resolve().parents[2];output=root/'benchmarks/native_hypersurfaces/results';output.mkdir(exist_ok=True)
env=dict(os.environ,JULIA_NUM_THREADS='1',JULIA_NUM_PRECOMPILE_TASKS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',DWORK_POWER_METHOD='auto',NATIVE_COLD='1' if args.cold else '0')
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.workshop,text=True).strip()=='a5f802454f425ed6eae2aea197a6b3fee65ddba8'
records=[]
with tempfile.TemporaryDirectory() as temp:
 if 'workshop' in args.modes:
  for n in (4,5,6):
   subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17',f'-DNV={n}',str(args.workshop/'dwork/native/dense_frobenius.cpp'),'-o',str(Path(temp)/f'original_{n}')],check=True)
 for i,case in enumerate(args.cases):
  # Alternate the order of the paired methods across fixtures.
  for mode in (args.modes if i%2==0 else list(reversed(args.modes))):
   n=int(case.split('_n')[1].split('_')[0])
   command=shlex.split(args.python)+['benchmarks/native_hypersurfaces/workshop.py',str(args.workshop),str(Path(temp)/f'original_{n}'),case] if mode=='workshop' else shlex.split(args.julia)+['--project=.','benchmarks/native_hypersurfaces/benchmark.jl',mode,case]
   load_start=os.getloadavg();start=time.perf_counter();name=f'{args.label}-{case}-{mode}'
   with (output/f'{name}.txt').open('w') as log:
    process=subprocess.Popen(['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v']+command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    try:status=process.wait(timeout=args.timeout)
    except subprocess.TimeoutExpired:
     os.killpg(process.pid,signal.SIGTERM)
     try:process.wait(timeout=10)
     except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
     status='timeout'
   record=dict(case=case,mode=mode,status=status,cold=args.cold,wall_seconds=time.perf_counter()-start,load_average_start=load_start,load_average_end=os.getloadavg())
   records.append(record);print(json.dumps(record),flush=True)
 (output/f'{args.label}.json').write_text(json.dumps(records,indent=2)+'\n')
raise SystemExit(0 if all(r['status']==0 for r in records) else 1)
