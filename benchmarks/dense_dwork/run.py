"""Serial benchmark driver: python3 benchmarks/dense_dwork/run.py WORKSHOP
Use --python 'sage -python' when cypari2 is provided by Sage.
"""
import argparse, json, os, platform, shlex, signal, subprocess, tempfile, time
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('workshop',type=Path)
parser.add_argument('--python',default='python3')
parser.add_argument('--julia',default='julia')
parser.add_argument('--modes',nargs='+',default=['workshop','native','upstream'])
parser.add_argument('--primes',type=int,nargs='+',default=[37,61,67,101])
parser.add_argument('--timeout',type=int,default=600)
parser.add_argument('--cold',action='store_true',help='one solve in a fresh process; includes startup')
parser.add_argument('--output',type=Path,default=Path('benchmarks/dense_dwork/results'))
args=parser.parse_args()
repo=Path(__file__).resolve().parents[2]
args.output.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,JULIA_NUM_THREADS='1',JULIA_NUM_PRECOMPILE_TASKS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',DENSE_BENCH_COLD='1' if args.cold else '0')
# Compilation is explicitly outside every timer; the pinned source is unchanged.
with tempfile.TemporaryDirectory() as temp:
    executable=Path(temp)/'original'
    if 'workshop' in args.modes:
        subprocess.run([os.environ.get('CXX','c++'),'-O3','-std=c++17','-DNV=4',
                        str(args.workshop/'dwork/native/dense_frobenius.cpp'),'-o',str(executable)],check=True)
    records=[]
    for mode in args.modes:
        for p in args.primes:
            command=(shlex.split(args.python)+['benchmarks/dense_dwork/workshop.py',str(args.workshop),str(executable),str(p)] if mode=='workshop' else
                     shlex.split(args.julia)+['--project=.','benchmarks/dense_dwork/benchmark.jl',mode,str(p)])
            prefix=['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v']
            label=f'{mode}-{p}'+('-cold' if args.cold else '')
            start=time.perf_counter()
            with (args.output/f'{label}.txt').open('w') as log:
                process=subprocess.Popen(prefix+command,cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                try:
                    status=process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGTERM)
                    try: process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid,signal.SIGKILL);process.wait()
                    status='timeout'
            record=dict(mode=mode,p=p,status=status,process_wall_seconds=time.perf_counter()-start,cold=args.cold)
            records.append(record)
            print(json.dumps(record),flush=True)
    (args.output/('cold.json' if args.cold else 'runs.json')).write_text(json.dumps(records,indent=2)+'\n')
    raise SystemExit(0 if all(record['status'] == 0 for record in records) else 1)
