# Native dense hypersurface integration

This experiment extends the original native K3 integration to cubic threefolds
and fourfolds, and widens the K3 prime range. It uses **dense Dwork–Cartier with
Harvey window coefficient extraction**, with scalar serial C++ kernels and a
Julia/Oscar public interface. All coefficient requests stay in one native
matrix call. It does not select itself through `zeta_function(f)`.

## Reproducible sources and supported inputs

- Upstream base: `d9633cf1847fed1f2ccf7a48ab855d4ce3ca4bde`.
- Original K3 integration: `61d0ecc64fc03dfe80b71734ab7452513c973c9f` (PR #99).
- Workshop reference: `MathAgents/andrew-huchala` at
  `a5f802454f425ed6eae2aea197a6b3fee65ddba8`.
- Source attribution and protocol: [`deps/native/PROVENANCE.md`](../../deps/native/PROVENANCE.md).
- Shipping proof and complete word bounds:
  [`native_hypersurfaces.md`](../../docs/src/experimental/native_hypersurfaces.md).

Every monomial of the stated degree must occur with nonzero coefficient.
Inputs must be geometrically smooth over a prime field.

| Input | Prime range | Primitive degree | Native precision |
| --- | --- | ---: | --- |
| Quartic in four variables | 11–251 | 21 | N=3,4; N=5 only for p≤47 |
| Cubic in five variables | 11–23 | 10 | N=5 |
| Cubic in six variables | 11–167 | 22 | N=5 (internally α/p modulo p⁴) |

The table specifies accepted arithmetic, not a promise that every input has a
unique lift at the precision limit. Ambiguous lifts fail explicitly. Singular,
partial-support, extension-field and out-of-range inputs are refused.

`DeRham.dense_dwork_coefficients(f)` returns descending integer primitive
coefficients. `DeRham.dense_dwork_zeta(f)` returns the full exact projective zeta
function, including every ambient factor. `ProjectiveHypersurface` is accepted.
The [documentation](../../docs/src/experimental/native_hypersurfaces.md#runnable-example)
contains a standalone cubic-threefold example.

## Build and test commands

From the DeRham checkout, with Julia 1.13.2 and a C++17 GCC/Clang compiler:

```sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
julia +1.13.2 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build("DeRham")'
julia +1.13.2 --project=. benchmarks/native_hypersurfaces/verify_archives.jl
julia +1.13.2 --project=. benchmarks/native_hypersurfaces/verify_guards.jl
python3 benchmarks/native_hypersurfaces/count_points.py
python3 benchmarks/native_hypersurfaces/sanitize.py
julia +1.13.2 --project=. benchmarks/dense_dwork/verify.jl
python3 benchmarks/dense_dwork/sanitize.py
julia +1.13.2 --project=. -e 'using Pkg; Pkg.test(test_args=["--workflow=presubmit"])'
julia +1.13.2 --project=. -e 'using Pkg; Pkg.test(test_args=["--name=native_dwork_n6_p167"])'
```

The build happens once; installed solves require no compiler, Python or workshop
checkout. Native binaries are not committed. Build failure leaves the rest of
DeRham usable and this opt-in backend unavailable. JLL packaging is future work.

`import_fixtures.py WORKSHOP` imports ten archived full polynomials into the E2E
v3 catalogue and generates `fixtures.jl`. The workshop is required only to
regenerate reference data and benchmark the original implementation. Four
archived K3s at p=37,61,67,101 remain in the original K3 tests. New archives
cover threefolds at p=13,19,23; fourfolds at p=11,17,23,89,167; and K3s at
p=151,251. The p=251 record is the archived sparse quartic after a change of coordinates over F₂₅₁,
not the separate generic p=251 row.
The JSON notes give the exact archive source for each equation.

For the fixed-seed cohort, use a Python environment with `cypari2` (these runs
used Sage's Python):

```sh
sage -python benchmarks/native_hypersurfaces/reference_cases.py WORKSHOP /tmp/native-cohort.tsv 100
julia +1.13.2 --project=. benchmarks/native_hypersurfaces/check_cohort.jl /tmp/native-cohort.tsv
```

The reference script compiles the pinned workshop kernel unchanged, once per
shape. It saves equations, complete primitive coefficients and refusals.
Native-to-native agreement is a regression check, not an independent proof.
The projective point-count program and upstream solver provide separate checks.

## Validation results

Correctness checks (see the storage-change recheck below):

- 300/300 smooth fixed-seed comparisons through the public Julia coefficient API:
  100 quartic K3s at p=13, 100 cubic threefolds at p=11, and 100 cubic fourfolds
  at p=11. All primitive coefficients matched the unchanged workshop source.
  All 26 singular reference cases were refused (9 K3, 11 threefold, 6 fourfold).
  Cohort TSV SHA-256:
  `4aca5dc7a449697498677d8c1e774209e3868a866d48d69909ec0864ec28a97b`.
- Ten archived equations: all 11/22/23 primitive coefficients matched, and each
  full zeta function matched an independently assembled numerator/denominator.
- 21 public guard/frontend checks passed, including singular, partial-support,
  wrong degree, extension-field, unsupported prime and native precision cases.
- 624 short reconstruction/range checks passed (167 original plus 457 new).
- Direct projective point counts matched: threefolds F₁₃=2341 and F₁₉=7259;
  fourfolds F₁₁=16160 and F₁₇=88486. These check the first primitive trace and
  the dimension-dependent sign independently of the matrix code.
- Documenter built the complete documentation successfully; deployment was
  skipped as expected for a local build.

- ASan/UBSan passed for `(n,p,N)=(4,251,4),(5,13,5),(5,23,5),(6,167,5),(6,13,5)`,
  plus unsupported-prime/precision and native singularity refusals.
- All 37 original K3 checks passed, including the archived p=37,61,67,101
  polynomials, small-prime direct counts and independent upstream controlled
  reduction at p=11,13, and maximum-precision guards.
- The exact documented Julia example executed successfully.

After selecting the smaller threefold storage, all 100 threefold cohort cases
and 11 singular refusals were repeated successfully, as were the 21 public
interface checks. The subset checksum is
`9600432ba152997b3c72f76bc9f9f76e848694581cdc1e90d7fed9deca6bfd4f`.
The two other shapes retain 16-bit storage and the same arithmetic.
All E2E records are exercised again on the final build below.

Final package checks passed:

- `Pkg.test(test_args=["--workflow=presubmit"])`: 321 E2E assertions,
  23 public-API checks and 18 quality-gate checks, plus the conventional tests.
- All 16 exact-name invocations passed: the ten new records and six re-timed
  original K3 records. Run `Pkg.test(test_args=["--name=NAME"])` for each name
  in `test/e2e/records/{native_hypersurfaces,dense_dwork}.toml` to reproduce.
- The first exact-name call ran serially in a fresh Julia process. Its 18
  available priming records totaled 78.511 seconds, below the 300-second
  budget; two GPU priming records were unavailable. Remaining exact-name
  correctness checks used at most three independent processes. All reported
  benchmark medians and catalogue runtime metadata were measured serially.
- No final native validation or package-test failures. The initial upstream
  benchmark configuration error and subsequent timeout/skip statuses are
  reported separately below; they are not successful equality comparisons.
- Presubmit skipped `gpu_p7_depthfirst__dim_1_deg_4_random_1`,
  `gpu_s012__dim_2_deg_3_fermat_1`, and `gpu_s012__dim_2_deg_3_fermat_2`
  because CUDA was unavailable. The `gpu_primer` skip is an intentional
  harness test. The full all-runtime catalogue, Linux and Windows were not run.

Package-test logs in `results/presubmit.txt` and `results/focused.txt` summarize
the completed runs, omitting dependency inventories and unrelated runtime
rows. `results/validation.json` records final statuses. Benchmark logs retain
stdout/stderr and timing data; carriage returns and trailing spaces are removed.


## Matched performance

Machine: Apple M3 Max, 16 CPU cores, 48 GiB RAM, macOS 26.6.2 arm64.
Software: Julia 1.13.2, Oscar 1.8.2, Apple Clang 21.0.0; workshop reference
Python 3.14.3, cypari2 2.2.4, PARI 2.17.1 (SageMath 10.9).
The upstream comparisons explicitly pass `S=[0,1,2]` for cubics, as in its
existing cubic catalogue; an unconfigured call fails before computation with
`Length of S must be <= 3`. Both native sources are compiled with `-O3 -std=c++17 -DNV=n`, without OpenMP,
SIMD or architecture-specific flags. Julia, BLAS and OpenMP thread counts are
one. Each warm row is the median of five complete exact-zeta calls, after one
first solve; every answer is compared with the archived full zeta function.
Julia timings include Oscar smoothness, process startup, all requested
precisions, arbitrary-precision reconstruction and rational-function assembly.
Compilation is outside all timings. Julia GC runs before, not inside, each
warm timer. Original Python/native timings use the corresponding complete
matrix/reconstruction/assembly route, not just the C++ kernel.

This comparison deliberately calls the workshop's **direct native** route;
it is not a comparison with its Picard–Fuchs or automatic dispatcher. Those
methods have different setup and certification boundaries and are not shipped.

```sh
python3 benchmarks/native_hypersurfaces/run.py WORKSHOP --python 'sage -python'
python3 benchmarks/native_hypersurfaces/run.py WORKSHOP --python 'sage -python' \
  --modes upstream-default --cases native_dwork_n5_p13 --timeout 360 --label upstream-threefold
python3 benchmarks/native_hypersurfaces/run.py WORKSHOP --python 'sage -python' \
  --modes upstream-default upstream-depthfirst \
  --cases native_dwork_n6_p11 --timeout 360 --label upstream-fourfold
python3 benchmarks/native_hypersurfaces/run.py WORKSHOP --python 'sage -python' \
  --modes workshop native --cases native_dwork_n6_p11 \
  --cold --timeout 360 --label cold
```

The driver alternates paired method order across equations. Raw stdout/stderr,
`/usr/bin/time -l` memory data, process wall times, exit/timeout status and host
load averages are retained under `results/`. Cold means a fresh interpreter
process with existing package caches, not a clean Julia installation. The first
solve timing starts after imports; the cold process wall time includes imports.
Peak RSS from macOS `time` includes the largest child maximum, not a sum of
simultaneously live host/child RSS; it is not a whole-process-tree peak.

The first small-input comparison exposed a storage regression: threefold
half-power tables used 32-bit words even at p=11,13, where 16 bits suffice.
The profile attributed about 60 ms to native coefficient extraction; Julia's
smoothness and reconstruction took about 1.3 ms. The shipping kernel now
selects 16-bit tables when `Q<=65536`, preserving 32 bits for larger threefold
primes. Preliminary logs are marked `pre-width-*`; the main tables below use
the corrected implementation.

| Equation / backend | Workshop warm (s) | Julia warm (s) | Workshop peak RSS (MiB) | Julia peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| Cubic threefold, p=13; half-power / 16-bit | 0.338 | 0.348 | 45.1 | 3094.9 |
| Cubic fourfold, p=11; half-power / grouped residues | 4.117 | 4.026 | 92.6 | 3167.7 |
| Cubic fourfold, p=17; shared scalar windows | 3.462 | 3.458 | 46.5 | 3117.0 |
| Quartic K3, p=151; scalar windows | 3.517 | 3.570 | 51.4 | 3125.7 |

The final warm medians differ from the workshop by at most 3% on these four
inputs. Julia/Oscar process lifetime peak RSS is about 3.0–3.1 GiB, much larger
than the Python/native route. This is not a claim of cold-start or memory parity.
The native child alone peaked at 16.7 MiB (threefold p=13), 91.9 MiB (fourfold
p=11), 8.6 MiB (fourfold p=17), and 4.3 MiB (K3 p=151). The separate fourfold
p=13 storage-boundary probe peaked at 237.6 MiB. Reproduce these child-only
measurements with `python3 benchmarks/native_hypersurfaces/memory.py`.
All measurements were on a desktop with normal background activity, not an
isolated benchmark machine; driver JSON retains the load averages.

The upstream default threefold solve at p=13 did not return before the
360-second **fresh-process** deadline. Its default dispatch is `:depthfirst`
for this dimension, so the redundant explicit-depthfirst attempt was canceled
and recorded as skipped. For the fourfold at p=11, both `:default` and
`:depthfirst` reached the same deadline without completing their first exact
solve. Therefore there is **no independent full-polynomial equality result or
warm median for these new dense cubics**. The completed independent upstream
comparisons are the original K3 checks at p=11,13; the new cubic first traces
also have the independent point-count checks above. No speedup factor is
inferred from the timeouts.

| Fresh process, fourfold p=11 | Total wall time (s) | First solve after imports (s) | Peak RSS (MiB) |
| --- | ---: | ---: | ---: |
| workshop | 4.69 | 4.336 | 92.9 |
| native | 38.20 | 6.634 | 3131.7 |

The upstream cold solve was not repeated: each already-timed-out comparison
started in a fresh process and had not completed its first solve. Its cold
result remains censored at 360 seconds, with no final RSS measurement. Package
loading and JIT are substantial on the Julia route; warm parity does not imply
cold parity. No automatic selection rule is introduced.


## Limitations and review focus

- The all-coordinate-strata trace comparison with ordinary projective primitive
  cohomology is the central mathematical review point. The implementation is
  opt-in and experimental; tests do not replace that argument.
- Conservative complete-splitting cutoffs are used throughout. No sharper
  empirical reduced-slab precision or unproved Hodge-column shortcut is used.
- Exact bounds cover all accepted native widths, including delayed sums and
  half-power tables. The largest accepted primes are exercised under ASan and
  UBSan. Sanitizers supplement the bounds, rather than prove them.
- Prime fields and full support only; no extension fields, sparse dispatch,
  quartic threefolds, double-plane K3s, toric reduction or transport.
- The recurrence advances linearly in its exponent. No square-root-in-p runtime
  or controlled-reduction characterization is claimed.
- Portable serial GCC/Clang build; local execution was on arm64 macOS. Linux,
  Windows and GPU execution were not tested here. MSVC without 128-bit integer
  support is unsupported.
- Files under `benchmarks/dense_dwork/results/` are historical PR #99 results;
  the current validation and measurements are under this directory.
