# AGENTS.md

Before adding or changing an E2E test, read `test/e2e/README.md`.

The E2E README documents executor registration, expected-value sources,
notes-tag rules, records, and runtime measurement.

## Invariant testing

- Keep unit tests short and use precomputed invariant data.
- Put longer examples computed from the defining polynomial in the E2E
  catalogue, with an accurate warm runtime.
- Each E2E record runs twice: a warm-up followed by a timed execution. Both
  normalized results are checked against the expected value and against each
  other; runtime summaries use the timed second execution.
- When a change affects execution speed, re-time affected E2E records in the
  same change.
- Derive an expected value from stored `Lpoly` data whenever that data
  determines the invariant. Use the `derham-invariants` notes tag only for
  values that `Lpoly` does not determine; see the E2E README for its format.
