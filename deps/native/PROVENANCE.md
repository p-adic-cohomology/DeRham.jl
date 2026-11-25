# Provenance and build boundary

The original K3 integration below is extended by the current
[broader integration](#broader-direct-integration). Its complete accepted ranges
and arithmetic bounds are in `docs/src/experimental/native_hypersurfaces.md`.

The copyright holder authorized contribution under **MIT OR GPL-3.0-or-later**
(the repository license). Source: Andrew Huchala's workshop,
[MathAgents/andrew-huchala at 708c1efbf00b9675ca3857f92382aaf200063551](https://github.com/MathAgents/andrew-huchala/tree/708c1efbf00b9675ca3857f92382aaf200063551).

- `dense_frobenius.cpp`: adapted from `dwork/native/dense_frobenius.cpp`.
  Retains Jacobian reduction, Teichmüller lifting, grouped Cartier expansion
  and batched windows. Retains the faster half-power path only at p=11,13,
  where complete sums fit in uint64 by explicit bounds. Removes the unused
  fewnomial and top-chart paths, unverified basis fallback, environment dispatch
  and OpenMP pragmas.
  Fixes four variables and validates the supported input/precision range.
- `harvey_window.hpp`: the workshop's batched recurrence, with attribution.
- `src/invariants/DenseDwork.jl`: Julia interface and reconstruction adapted
  from `dwork/dense_native.py` and `dwork/projective_zeta.py`.
- E2E fixtures: `experiments/dense_word_range.txt` and the seeded equations
  in `experiments/scoping_q1.py`, replayed by `benchmark_window_snapshot.py`.
- Proof starting point: `results/projective_soundness.html`; the shipping
  cutoff and word bounds are documented separately in
  `docs/src/experimental/dense_dwork.md`.

This is dense Dwork–Cartier with Harvey window coefficient extraction.
Harvey's recurrence is a mathematical algorithm, not a claim of square-root
runtime for this implementation. No code from PR #97 is required.

## Process protocol, version 1 (original K3 instance)

`deps/build.jl` compiles C++17 once to `deps/usr/bin/dense_dwork` (or `.exe`).
It uses `CXX` or `c++`, `-O3 -std=c++17`, without architecture-specific flags
or an OpenMP runtime. GCC/Clang's 128-bit integer extension is required;
MSVC alone is unsupported. There are no external native libraries.
The installed binary can be run without a compiler, Python, or the workshop.
If compilation is unavailable, existing DeRham functions still work; this
opt-in function reports the missing backend. A future JLL artifact can replace
the build step without changing the protocol. Binaries are not committed.

One process computes one complete matrix (441 residues), retaining all
coefficient requests in batched C++ windows. Processes isolate global state;
concurrent Julia calls cannot share it. No untrusted shell interpolation occurs.
Input on stdin is `p N 35`, then 35 lines of four exponents and a nonzero
canonical coefficient modulo p. Terms must be distinct quartic monomials.
Output is `21`, then 21 lines `k nu1 nu2 nu3 nu4`, then 21 matrix rows.
Row i gives the coordinates of alpha applied to basis vector i, modulo p^N.
Stderr contains diagnostic timings. Exit 0 means success, 3 a failed geometric
smoothness rank check, 2 unsupported input/arithmetic. Other process failures
are errors. The Julia caller checks dimensions and canonical residues.

## Broader direct integration

Revision `a5f802454f425ed6eae2aea197a6b3fee65ddba8` supplies the shared fourfold
readouts, small-prime residue grouping, socle certificate, and divided-matrix implementation in
`dwork/native/dense_frobenius.cpp`. The scalar window's guarded modulus is
extended to 52 bits with the complete audit in
`docs/src/experimental/native_hypersurfaces.md`. No external native library,
SIMD instruction requirement, toric implementation, or transport implementation
is added. Native code remains under MIT OR GPL-3.0-or-later with the copyright
holder's authorization. The Julia weight-dependent reconstruction is adapted
from the same revision of `dwork/projective_zeta.py`.

Protocol v1 is also built as `dense_dwork_5` and `dense_dwork_6`, for full cubic
threefolds (35 terms, rank 10) and fourfolds (56 terms, rank 22). The quartic
executable still uses `dense_dwork` (35 terms, rank 21). There are r basis lines
and r matrix rows after the rank header. Fourfold output residues fit 64 bits
and are parsed into arbitrary-precision integers. The optional build handles
each executable separately; unavailable backends fail explicitly.

The broader branch retains generic half-power convolution for cubic threefolds
and the p=11,13 fourfold residue packing from the pinned revision. The complete
sum bounds and word types are documented in the current audit. All SIMD and
OpenMP code is excluded from these portable builds.
