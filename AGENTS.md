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
  does not run it. Keep the set below five minutes and cover every configuration
  whose cheapest record is below `ci_seconds`; see `test/e2e/README.md` for the
  configuration key and how to add an entry.
- When a change affects execution speed, re-time affected E2E records in the
  same change.
- Derive an expected value from stored `Lpoly` data whenever that data
  determines the invariant. Use the `derham-invariants` notes tag only for
  values that `Lpoly` does not determine; see the E2E README for its format.
