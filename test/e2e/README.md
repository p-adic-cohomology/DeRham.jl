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
- `parameters` are passed straight through as keyword arguments to
  `DeRham.zeta_coefficients`.

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
