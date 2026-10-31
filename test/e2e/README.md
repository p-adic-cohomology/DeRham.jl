# E2E test catalogue

A declarative, metadata-driven catalogue for DeRham's end-to-end
zeta-coefficient and invariant tests. Example data (the mathematical variety
plus stored `Lpoly` results and optional invariant notes) is separate from
execution records (which invariant executor to run, with which parameters,
tags, and runtime). `Pkg.test()` runs the catalogue as part of the conventional
test suite; see `test/runtests.jl`.

## Add local JSON example data

Add a file under `test/e2e/examples/` (any subdirectory — placement has no
execution semantics) containing one or more grouped cases in the
[zeta_test_suite v3 schema](https://github.com/jjgarzella/zeta_test_suite)
(`schema.json` there is the normative definition):

```json
{
  "schema_version": "3",
  "cases": [
    {
      "id": "my_new_example",
      "variety": {
        "coeff_domain": { "kind": "integer" },
        "dim": 1,
        "non_middle_factors": {
          "kind": "projective_lefschetz",
          "middle_factor_content": "full"
        },
        "model": {
          "kind": "plane_curve",
          "pretty": "x1^3 + x2^3 + x3^3 = 0",
          "vars": ["x1", "x2", "x3"],
          "monomials": [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
          "coeffs": [1, 1, 1]
        }
      },
      "results": [
        { "p": 7, "Lpoly": { "coeffs_asc": [1, 1, 7] } }
      ],
      "notes": "why this case exists"
    }
  ]
}
```

The `zeta_coefficients` executor supports `model.kind` of `plane_curve` (3
variables) or `projective_hypersurface` (>= 2 variables), and
`coeff_domain.kind == "integer"`. DeRham returns the **primitive** middle
factor. Set `middle_factor_content` to `"primitive"` for surfaces. Odd-
dimensional hypersurfaces may use `"primitive"` or `"full"` because those
factors coincide. `Lpoly.coeffs_asc` is **ascending** (constant term first);
DeRham's own `zeta_coefficients` return value is **descending** (leading
coefficient first) — the executor reverses the stored coefficients and rejects
an accidentally-reversed `coeffs_asc` (its first entry must be `1`).

Example ids must be nonempty and globally unique across every local file and
the pinned artifact.

## Add a TOML execution record

Add (or append to) a `[[test]]` array-of-tables entry under
`test/e2e/records/` (any subdirectory):

```toml
[[test]]
name = "my_new_example_default"
type = "zeta_coefficients"
tags = []
runtime = 0.1
example = "my_new_example"
parameters = { algorithm = "depthfirst", fastevaluation = true }
```

- `name` must be nonempty and globally unique across every record file.
- `runtime` is the measured warm execution time in seconds, used only for
  `--workflow=` membership. Use Julia 1.13.2; update affected records in the
  same change when an implementation change affects their speed.
- `example` must match an example `id`. If that example has more than one
  prime result, also set `p = <prime>` to disambiguate.
- `tags` must be a subset of `config.toml`'s `[tags] allowed` list (currently
  `gpu`, `multithreading`, `uint`, `bigint`, `karatsuba`). A `gpu` tag (or
  `parameters.use_gpu = true`) marks the record as requiring a functional
  CUDA device.
- For `Lpoly`-derived records, common computation options are passed through as
  keyword arguments to the selected DeRham invariant. `point_counts` and
  `k3_picard_rank` also require a positive integer `k`; the notes-based
  `a_number` and `artin_invariant` executors accept no parameters.

### Supported record types

| `type` | Expected value and input | Record-specific parameters |
| --- | --- | --- |
| `zeta_coefficients` | Stored `Lpoly`; plane curves or projective hypersurfaces | None required |
| `newton_polygon` | Stored `Lpoly` and the example polynomial | None required |
| `point_counts` | Stored `Lpoly` and the example dimension | `k`: positive extension degree |
| `picard_rank_bound` | Stored `Lpoly`; projective surface | None required |
| `k3_geometric_picard_rank` | Stored `Lpoly`; quartic K3 surface | None required |
| `k3_picard_rank` | Stored `Lpoly`; quartic K3 surface | `k`: positive extension degree |
| `k3_picard_realization_degree` | Stored `Lpoly`; quartic K3 surface | None required |
| `artin_mazur_height` | Stored `Lpoly`; quartic K3 surface | None required |
| `a_number` | `derham-invariants` notes value; plane curve | None |
| `artin_invariant` | `derham-invariants` notes value; supersingular quartic K3 over an odd prime | None; note value must be in `1:10` |

All `Lpoly`-derived record types accept the shared DeRham computation options
recognized by the catalogue: `S`, `verbose`, `changef`, `givefrobmat`,
`algorithm`, `termorder`, `vars_reversed`, `fastevaluation`,
`always_use_bigints`, `use_gpu`, and `use_threads`. `algorithm` accepts
`default`, `depthfirst`, `varbyvar`, or `pchunk`; `S` must be an array and
`termorder` must be a string.

## Reference an example from the pinned artifact

Artifact examples (from the `zeta_test_suite` Julia artifact, resolved lazily
— `using DeRham` alone never downloads it) load into the same catalogue as
local examples, keyed by their own `id`. Reference one from a TOML record
exactly as you would a local example; there is no separate syntax.

## Run tests by name, tag, or workflow

`Pkg.test` passes `test_args` straight through to `ARGS`. Exactly one
selector may be given; a bare `Pkg.test()` defaults to `--workflow=ci`.
An invocation with only a `--name=` or `--tag=` selector runs the selected E2E
record(s) without running the conventional testsets, E2E harness tests, or
quality gates. Bare `Pkg.test()` and `--workflow=` invocations still run those
tests and gates.

```julia
Pkg.test("DeRham")                                    # --workflow=ci (default)
Pkg.test("DeRham"; test_args = ["--workflow=presubmit"])
Pkg.test("DeRham"; test_args = ["--workflow=full"])
Pkg.test("DeRham"; test_args = ["--name=my_new_example_default"])
Pkg.test("DeRham"; test_args = ["--tag=bigint"])
```

`--workflow=` selects by each record's `runtime` against the configured
`presubmit_seconds` / `ci_seconds` thresholds (`presubmit ⊆ ci ⊆ full`).
`--name=` and `--tag=` ignore runtime thresholds; a `--tag=` selection (like
a workflow) silently skips cases whose capability (e.g. `gpu`) is
unavailable, while an unavailable `--name=` selection fails.

Each selected record runs twice: a warm-up execution followed by a timed
execution. The harness normalizes both results, checks each against the
expected value, and requires the two results to match. The `[e2e]` case line
reports `first_runtime` and `actual_runtime`; the runtime summary compares the
declared warm runtime with the timed second execution.


## Adding an executor

Add one module under `test/e2e/executors/` and register it with one `include`
and one entry in `EXECUTOR_REGISTRY` in `E2E.jl`. Keep the module's functions
focused on one invariant:

- `validate_parameters(name, parameters)` rejects unsupported record
  parameters before execution.
- `expected(example, record)` returns the expected value in the form that
  executor compares.
- `execute(input, parameters)` computes the invariant from the example's
  defining polynomial.
- `normalize(raw)` and `compare(expected, actual)` make result representation
  and equality explicit.
- `INVARIANT_NAME` identifies the invariant for notes-tag validation.
- `EXPECTED_SOURCE` is `:lpoly` when `expected` derives its value from the
  stored `Lpoly`, or `:notes` when it reads a value from the notes tag.

`Common.jl` supplies prime/result resolution, polynomial construction,
parameter conversion, capability checks, notes-tag parsing, and the stored
descending-order `Lpoly` accessor with its orientation check. Use these shared
helpers instead of duplicating catalogue behavior.

Derive expected values from the stored `Lpoly` whenever that polynomial
determines the invariant, applying DeRham's data-form function. Do not store a
second expected value for a derivable invariant. Put only values not determined
by `Lpoly` in the optional notes line:

    derham-invariants: {"5": {"a_number": 1, "artin_invariant": 2}}

The JSON object is keyed by prime strings, then invariant names. Free prose may
appear on other notes lines. Catalogue validation rejects malformed tags,
primes without a case result, unregistered invariant names, and values for
invariants derived from `Lpoly`. Do not add provenance fields or new example
file kinds.

Keep short tests on precomputed data in `test/invariants/`. Put longer
examples computed from the defining polynomial in E2E records, with the warm
measured runtime. Run `Pkg.test(test_args=["--name=<record>"])` for each record
added or changed, plus the prescribed workflow test tier. The end-of-run
runtime table flags cases whose measured time is more than twice or less than
half their declared runtime; this is informational and does not change a test
result.
