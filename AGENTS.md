# AGENTS.md

Before adding or changing an E2E test, read `test/e2e/README.md`.

The E2E README documents executor registration, expected-value sources,
notes-tag rules, records, and runtime measurement.

## Invariant testing

- Keep unit tests short and use precomputed invariant data.
- Put longer examples computed from the defining polynomial in the E2E
  catalogue, with an accurate measured runtime.
- Each selected E2E record runs once, timed, and its normalized result is
  checked against the expected value. The fixed priming set runs before `ci`
  and `full` workflows and before `--name=` / `--tag=` selections; `presubmit`
  runs the priming records for its own configurations. Keep the set below five
  minutes and cover every configuration whose cheapest record is below
  `ci_seconds`; see `test/e2e/README.md` for the configuration key and how to
  add an entry.

## Test tiers

- `presubmit` has a budget of one minute after priming. It runs the unit test
  files marked `:presubmit` in `UNIT_TEST_TIERS` (`test/runtests.jl`) and the
  E2E records below `presubmit_seconds`. `ci` runs every unit test file, the
  quality gates and the records below `ci_seconds`; `full` runs everything.
- New unit test files go in `UNIT_TEST_TIERS`. Mark a file `:presubmit` only if
  it takes a few seconds on the shared host under its usual load.
- Declared E2E runtimes are measured on the shared host under its usual load.
- Besides presubmit, run the tests for what you changed:
  `Pkg.test(test_args=["--unit=<file>", ...])` for the unit test files that
  cover the changed code (for example `--unit=reduction/AKR.jl` for
  `src/reduction/policies/AKR.jl`), and `--name=` / `--tag=` for the E2E
  records that exercise it. Report which targeted tests you ran.
- When a change affects execution speed, re-time affected E2E records in the
  same change.
- Derive an expected value from stored `Lpoly` data whenever that data
  determines the invariant. Use the `derham-invariants` notes tag only for
  values that `Lpoly` does not determine; see the E2E README for its format.

## Experimental and proven algorithms

- Proven algorithms have a proof of correctness. Experimental algorithms have
  no proof of correctness; this includes algorithms whose precision choices
  are not certified.
- Experimental algorithms are documented in the "Experimental" section of
  the documentation, not alongside the proven invariants.
- Symbols for experimental algorithms are public but not exported: declare
  them with `public` in `src/DeRham.jl`, never with `export`.
- For `zeta_function` and its variants (`zeta_coefficients`, `L_polynomial`,
  and other coefficient-based entry points), experimental algorithms may take
  `series_precision` and `algorithm_precision` inputs. Proven algorithms
  calculate their own precision and do not take those inputs.

## Source layout

`src/` is one module, `DeRham`. Directories group files by component; they are
not submodules. The zeta pipeline runs graded ring → Frobenius → pole reduction
→ final reduction → invariants, and the precision code chooses the working
precision for each stage.

- `src/graded_ring/`: linear algebra in graded rings. Monomial bases of graded
  pieces (`GradedPieces.jl`, `GradedExpCache.jl`) and conversion between
  polynomials and coefficient vectors. `maps/` has one file per map between
  graded pieces (cohomology basis, controlled map, classical standard-reduction
  map, smoothness test); each builds the map's matrix and solves its linear
  algebra problem (pseudo-inverse, Hensel lift). Shared solvers live in
  `LinearAlgebra.jl`.
- `src/frobenius/`: applies Frobenius to the cohomology basis, producing
  truncated sums of terms with poles.
- `src/finalreduction/`: reduces the output of pole reduction to the cohomology
  basis and assembles the Frobenius matrix.
- `src/recurrence/`: polynomial recurrences and memory management. Partially
  evaluated polynomial matrices (`AbstractPEP` and its eager, lazy and cached
  implementations), their evaluation along a reduction step, the choice of how
  the R_uv matrices are stored, and the CPU/GPU/Karatsuba linear-algebra
  backends (`backend/`, the `my_*` wrappers).
- `src/reduction/`: pole (controlled) reduction. The reduction context, the
  shared steps (direction choice, Costa-data terms), construction of the R_uv
  matrices, the `reducetransform` dispatcher, and the interchangeable policies
  in `policies/` (P-chunk, var-by-var, depth-first, AKR). The default policy is
  chosen by a deliberate heuristic.
- `src/precision/`: precision calculations. p-adic and series precision bounds,
  the default series truncation, and the precision data the pipeline consumes.
- `src/invariants/`: public wrappers. `zeta_coefficients` runs the pipeline;
  the other invariants (L-polynomial, point counts, Newton and Hodge polygons,
  Hasse–Witt matrix, a-number, Picard rank, heights, ...) are computed from it
  or from its coefficients.

Files at the top of `src/` (options, utilities, L-polynomial reconstruction,
slope polygons, example polynomials, the K3/𝔽₃ wrapper) do not belong to one
component yet. Unused historical source files have been removed.
`NemoDerivedWrappers.jl` (BSD-2-Clause)
and `SeriesPrecisionDefaults.jl` (GPL) carry their own license headers; do not
merge other code into them.
