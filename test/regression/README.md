# Regression fixtures

This directory holds the numerical-equivalence safety net for the
architectural refactor. The integrator is deterministic at full Float64
precision, so the invariant checked here is **zero drift**: re-running the
captured problem must reproduce the stored trajectory bit-for-bit (tolerance
1e-12 on every `t`, `q`, `p` component).

## Files

- [`capture.jl`](capture.jl) — one-shot script that builds four representative
  problems against the current public API and serialises each to a `.jld2`
  file under `fixtures/`. Re-run only when intentionally replacing the
  reference trajectories (commit the new fixtures in the same PR that changes
  the integrator).
- [`validate.jl`](validate.jl) — loads every fixture, rebuilds the problem
  using the **current** public API, solves it, and compares against the stored
  trajectory. Any drift beyond 1e-12 is a hard failure. Included automatically
  in `Pkg.test()` via the `regression fixtures` testset in `runtests.jl`.
- [`fixtures/*.jld2`](fixtures/) — committed reference trajectories produced
  by `capture.jl`.

## Fixture coverage

| Fixture                    | Covers                                                    |
|----------------------------|-----------------------------------------------------------|
| `twobody_ellipse`          | Unregularized symmetric-projection integrator, finite c   |
| `threebody_mixed`          | Unregularized 3-body with mixed-sign charges              |
| `close_approach_lifted`    | `:lifted_pair` Levi-Civita regularization fires (2D)      |
| `zollner_offmatch`         | `:adaptive_cartesian` + Zöllner κ ≠ 1 in substep path     |

## Regenerating

```bash
julia test/regression/capture.jl
```

Requires `JLD2` in the active environment (the package's default dev env is
fine; JLD2 is also listed in `[extras]` so the `Pkg.test()` path sees it).

## Validating

```bash
julia --project=. -e 'using Pkg; Pkg.test()'   # includes the regression testset
julia test/regression/validate.jl              # standalone
```

## When a refactor breaks equivalence

`validate.jl` reads only the fixture setup dict and calls the current public
API in `rebuild_problem`. When the API shape changes (e.g. Phase 1 of the
refactor renames `HamiltonianProblem` → `HamiltonianProblem`), update
`rebuild_problem` to use the new constructors while keeping the fixtures
themselves untouched. If the refactor genuinely changes numerical behaviour,
document the reason and regenerate the fixtures in the same PR.
