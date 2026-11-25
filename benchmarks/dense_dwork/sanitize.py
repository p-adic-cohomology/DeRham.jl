"""Native boundary and undefined-behavior check; no workshop dependency.
python3 benchmarks/dense_dwork/sanitize.py
"""
import os, subprocess, tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
support=[(a,b,c,4-a-b-c) for a in range(5) for b in range(5-a) for c in range(5-a-b)]
def data(p,N):
    coefficients=[(i**3+3*i+1)%(p-1)+1 for i in range(1,36)]
    return f'{p} {N} 35\n'+''.join(' '.join(map(str,(*e,c)))+'\n' for e,c in zip(support,coefficients))
with tempfile.TemporaryDirectory() as temp:
    binary=Path(temp)/'dense-sanitized'
    subprocess.run([os.environ.get('CXX','c++'),'-O1','-g','-std=c++17',
                    '-fsanitize=address,undefined','-fno-sanitize-recover=all',
                    str(root/'deps/native/dense_frobenius.cpp'),'-o',str(binary)],check=True)
    for p,N in [(11,5),(13,5),(47,5),(101,4),(47,6),(53,5),(257,4)]:
        result=subprocess.run([str(binary)],input=data(p,N),text=True,capture_output=True)
        expected=0 if (p,N) in [(11,5),(13,5),(47,5),(101,4)] else 2
        assert result.returncode==expected,(p,N,result.returncode,result.stderr)
        assert 'runtime error:' not in result.stderr and 'ERROR: AddressSanitizer' not in result.stderr
        print(f'p={p} N={N} exit={result.returncode} PASS',flush=True)
    bad=data(11,4).replace('0 0 0 4','-1 0 0 5',1)
    assert subprocess.run([str(binary)],input=bad,text=True,capture_output=True).returncode==2
    print('malformed exponent rejection PASS')
