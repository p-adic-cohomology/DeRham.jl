# Reproduce the dense Dwork checks

Source: MathAgents/andrew-huchala at
`708c1efbf00b9675ca3857f92382aaf200063551`; upstream baseline
`d9633cf1847fed1f2ccf7a48ab855d4ce3ca4bde`.
`fixtures.jl` contains the four archived equations and all 22 expected primitive
coefficients, plus the independently checked small inputs. The E2E JSON stores
the same data in the upstream catalogue format. The archived equations use
Python Random seeds 20261160, 20261184, 20261190, 20261224 in lexicographic
monomial order, exactly as the workshop's `scoping_q1.py`.

## Build and verify

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build("DeRham")'
julia --project=. benchmarks/dense_dwork/verify.jl
python3 benchmarks/dense_dwork/sanitize.py
julia +1.13.2 --project=. -e 'using Pkg; Pkg.test(test_args=["--workflow=presubmit"])'
# Also run each new record by name (11, 13, 37, 61, 67, 101):
julia +1.13.2 --project=. -e 'using Pkg; Pkg.test(test_args=["--name=dense_dwork_p37"])'
```

`verify.jl` exercises the public Julia methods, all 22 archived coefficients,
ambient factors, the projective frontend, direct projective point counts at
11 and 13, independent upstream coefficients on those inputs, singular and
unsupported equations, the maximum native precision and rejected next
precision. `sanitize.py` compiles a temporary ASan/UBSan build and checks the
arithmetic boundaries, with no workshop dependency. Runtime solves never
invoke a compiler.

## Timing boundaries

```sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 benchmarks/dense_dwork/run.py /path/to/pinned/workshop \
  --python 'sage -python' --julia 'julia +1.13.2' \
  --primes 11 13 37 61 67 101 --timeout 600
python3 benchmarks/dense_dwork/run.py /path/to/pinned/workshop \
  --python 'sage -python' --julia 'julia +1.13.2' \
  --primes 101 --cold --timeout 600
julia +1.13.2 --project=. benchmarks/dense_dwork/profile.jl
```

A normal Python interpreter with cypari2 also works (`--python python3`). The
workshop binary is compiled **unchanged**, once outside all timers, with the
same C++17 compiler, `-O3`, four variables and serial execution as the port.
Its OpenMP pragmas are inactive in this matched serial build. `workshop.py`
calls the original `dense_frobenius`, `primitive_polynomial` and `_assemble`;
only `kernel()` is replaced with the prebuilt executable path to exclude
compiler discovery. This measures the original native solver and exact
reconstruction, **not** the workshop's Linux-only automatic dispatcher or
its optional coordinate-frame search. Native Jacobian smoothness checks remain
included. Julia additionally validates with Oscar.

The Julia timer calls **`DeRham.dense_dwork_zeta(f)`**, including validation,
serialization, native process launches, every precision attempt, parsing,
integer reconstruction and ambient factors. Upstream uses its current default
solver with `use_threads=false, use_gpu=false`. Every timed call checks the same
exact full projective rational function; the original Python path also checks
each primitive coefficient. Fixture construction and expected-result assembly
are outside solve timers. Warm values are medians of three complete repeated
solves after a first call; no answer cache is introduced.

The cold command measures one solve in a **fresh process**, including imports,
Julia JIT and fixture construction, after installation/precompilation. It does
not flush OS filesystem caches. `/usr/bin/time` records peak process RSS
(bytes on macOS, KiB on Linux); this is not the sum of simultaneously resident
parent/child processes. Native children are much smaller than the Julia host.
Do not run other benchmarks or package precompilation concurrently.

## Results

Machine: Apple M3 Max, 16 CPU cores, 48 GiB RAM, arm64, macOS 26.6.2.
Final timings use Julia 1.13.2, Oscar 1.8.2, Apple Clang 21.0.0 (`-O3
-std=c++17`), CPython 3.14.3 from SageMath 10.9, cypari2 2.2.4 and PARI
2.17.1. All three solvers use one CPU thread. No GPU is used.

The four archived cases use the Harvey window backend in both native solvers.
The small p=11,13 cases make setup and wrapper overhead visible; both solvers
choose their half-power coefficient backend there. The profile reports
validation, native-process/matrix time, exact reconstruction and assembly
separately.

Warm medians of three complete solves; all outputs checked equal:

| p | Original workshop (s) | Public Julia API (s) | Upstream solver (s) |
| --- | ---: | ---: | ---: |
| 11 | 0.037329 | 0.035548 | 46.480284 |
| 13 | 0.048017 | 0.044821 | 50.847685 |
| 37 | 0.446078 | 0.440611 | 29.456245 |
| 61 | 1.081772 | 1.052626 | 41.580703 |
| 67 | 1.184845 | 1.167399 | 45.205936 |
| 101 | 2.063286 | 2.058480 | 71.615613 |

The public Julia warm medians are within 7% of the workshop's on these six
inputs (slightly faster in this run). Initially forcing the window backend at
p=11 cost about 135 ms versus about 34 ms for the workshop. Profiling and
retaining its bounded half-power backend at p=11,13 removed that regression.
The final p=11 stage medians are: validation 1.404 ms, native 31.622 ms, reconstruction 2.803 ms, assembly 0.072 ms.
Stage medians are diagnostic and need not sum to the independently timed
public-call median. Raw stage measurements are in `results/profile-11.txt`.

Fresh-process measurements after installation/precompilation:

| Solver, p=101 | Fresh-process wall (s) | First solve within process (s) | Peak process RSS (GiB) |
| --- | ---: | ---: | ---: |
| Original workshop | 3.733 | 3.229 | 0.047 |
| Public Julia API | 50.750 | 24.305 | 2.769 |
| Upstream solver | 100.264 | 77.731 | 2.537 |

**Cold startup is substantially more expensive through Julia.** Warm solve
parity does not imply startup or memory parity with the Python/native driver.
Peak process RSS across the six repeated-run processes was
workshop 0.046–0.050 GiB; public Julia 2.461–2.899 GiB; upstream 2.557–2.700 GiB. These are process maxima reported by
`/usr/bin/time`, not a simultaneous process-tree sum or solve-only allocations.

The upstream p=101 warm repetitions ranged from 66.073 to 86.183 seconds;
its median should not be read as a tightly estimated speed ratio. Runs were
serial in solver order, not interleaved, and thermal state was not controlled.

Every final benchmark process exited successfully. `results/runs.json` and
`results/cold.json` contain exit status and wall time; adjacent text files
contain every repetition and OS resource measurements. No automatic backend
selection rule is proposed from this sample.

## Validation status

- **Passed:** installation/instantiation and native build; the public-entry-point verification script (37 assertions); all 22 primitive coefficients on each of the four archived dense quartics at p=37,61,67,101; full ambient factors and `ProjectiveHypersurface` dispatch.
- **Passed:** independent upstream full-zeta comparisons on all six inputs, in the repeated benchmark calls. Direct projective counts at p=11 and p=13 are 123 and 183, and all primitive coefficients agree with upstream there.
- **Passed:** singular input, unsupported fields/support/primes, ambiguous/inconsistent reconstruction guards, and the maximum native precision paths. ASan/UBSan passed at (p,N)=(11,5),(13,5),(47,5),(101,4); the next unsupported precision/prime pairs and malformed exponents were rejected.
- **Passed:** Julia 1.13.2 presubmit, including 167 new short reconstruction/arithmetic checks, 315 E2E assertions, 23 public-API checks and 18 quality-gate checks. All six new exact-name E2E invocations passed; p=11,13 were rerun after the small-prime optimization, and the final presubmit reran all six.
- **Skipped:** three GPU catalogue records because this machine has no CUDA: `gpu_p7_depthfirst__dim_1_deg_4_random_1`, `gpu_s012__dim_2_deg_3_fermat_1`, `gpu_s012__dim_2_deg_3_fermat_2`. The separate `gpu_primer` skip in the log is an intentional harness test. No local Linux/Windows or documentation-site build was run.
- **Failed:** none in the final validation/benchmark suite. Preliminary 180-second upstream timing runs hit the time limit; the recorded final runs use 600 seconds and all completed.
