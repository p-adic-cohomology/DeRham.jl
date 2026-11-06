# AGENTS.md

Before adding or changing an E2E test, read `test/e2e/README.md`.

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
