"""Original workshop native solver, same input and exact projective output.
Run with the pinned workshop and a prebuilt NV=4 executable as arguments.
The kernel is unchanged; overriding kernel() only excludes build discovery.
On macOS use Sage's Python (cypari2), and /usr/bin/time -l externally for RSS.
"""
import json, os, statistics, sys, time
from pathlib import Path
root, executable, prime = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
sys.path.insert(0, str(root / 'dwork'))
from dense_native import dense_frobenius
import dense_native
from projective_zeta import primitive_polynomial, _assemble
# Both binaries use the same compiler/serial flags; all batching is unmodified.
dense_native.kernel = lambda n: executable
cases = json.loads((Path(__file__).resolve().parents[2] / 'test/e2e/examples/dense_dwork.json').read_text())['cases']
case = next(c for c in cases if c['results'][0]['p'] == prime)
support = case['variety']['model']['monomials']
coefficients = case['variety']['model']['coeffs']
expected = list(map(int, case['results'][0]['Lpoly']['coeffs_asc']))
expected_zeta = _assemble(expected, prime, 0, 2)
def solve():
    for N in range(3, 6 if prime <= 47 else 5):
        _, A, _ = dense_frobenius(prime, support, coefficients, N, measure_rss=False)
        polynomial = primitive_polynomial(A, prime, N, 2)
        if polynomial is not None:
            zeta = _assemble(polynomial, prime, N, 2)
            assert polynomial == expected
            assert zeta.numerator == expected_zeta.numerator and zeta.denominator == expected_zeta.denominator
            return
    raise ArithmeticError('ambiguous lift')
start=time.perf_counter(); solve(); first=time.perf_counter()-start
if os.environ.get("DENSE_BENCH_COLD") == "1":
    print(json.dumps(dict(p=prime,first_solve_seconds=first)))
    sys.exit(0)
times=[]
for _ in range(3):
    start=time.perf_counter(); solve(); times.append(time.perf_counter()-start)
print(json.dumps(dict(p=prime, first_solve_seconds=first, warm_seconds=times,
                     median_seconds=statistics.median(times), full_projective_zeta_equal=True)))
