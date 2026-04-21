# κ Storage Refactor — Future-Work Prompt

A self-contained prompt for a future session. Hand the contents below to a Claude Code agent to pick this up.

---

## Task

Refactor κ (per-pair Zöllner coupling) out of the global `params` vector so that it lives with the term that owns it (a `NamedTerm`-carried field or a per-term context), rather than being sliced out of `params` via `_pair_index`. Downstream code (statistics, plots, tests) should ask the term for its κ; the compiled EOMs should no longer require κ inside the global `params`.

This is the cascading half of what the `_research/NextSteps.md` "Follow-ups on deliberate refactor deferrals" section deferred. The `kinetic_term` / `coulomb_term` half has already been shipped.

## Why this is non-trivial

The current layout is:

```
params = [m₁…mₙ, q₁…qₙ, c, κ₁₂…κ_{N-1,N}]
length = 2N + 1 + N·(N-1)/2
```

It is defined at [src/types.jl:505-510](../src/types.jl#L505-L510) and referenced from the compiled-EOM signature `dq_dt_compiled(out, q, p, t, params)` / `dp_dt_compiled(out, q, p, t, params)` compiled at [src/hamiltonian_system.jl:135-142](../src/hamiltonian_system.jl#L135-L142). Those EOMs are on the integration hot path, called from:

- [src/solve.jl:4-52](../src/solve.jl#L4-L52) (`strang_splitting_flow!`)
- [src/solve.jl:55-95](../src/solve.jl#L55-L95) (`compute_constraint_residual!`)
- [src/solve.jl:97-243](../src/solve.jl#L97-L243) (`_projected_cartesian_step!`)
- Every regularized-step path in `src/solve.jl` and [src/integrators/regularized.jl](../src/integrators/regularized.jl)

`_pair_index` is defined at [src/hamiltonian_system.jl:86-89](../src/hamiltonian_system.jl#L86-L89) and called from:

| File | Line(s) | Purpose |
|------|---------|---------|
| `src/hamiltonian_system.jl` | 69, 71 | inside `weber_term` / `zollner_term` symbolic builders |
| `src/hamiltonian/builders/weber.jl` | 69, 104 | symbolic + `_weber_pair_decomposition` runtime slice |
| `src/hamiltonian/builders/zollner.jl` | 71, 102 | symbolic + `_zollner_pair_decomposition` runtime slice |
| `src/statistics/forces.jl` | 252 | `compute_pair_force_timeseries` reads κ via `kappas(prob)[_pair_index(...)]` |
| `src/statistics/energy.jl` | 276 | `compute_energy_timeseries` reads κ the same way |
| `ext/WeberElectrodynamicsMakieExt.jl` | 101 | animation dashboard |
| `test/test_hamiltonian_system.jl` | 153-175 | ~19 test cases that construct params manually |
| `docs/src/internals.md` | 62-64 | documented layout |
| `CLAUDE.md` | 58-59 | "Direct calls to `dq_dt_compiled`/`dp_dt_compiled` **must** include κ entries" |

Changing the compiled-EOM signature is a **breaking change** cascading through all of the above plus ~100 test cases that construct `HamiltonianProblem` / raw params.

## Suggested workflow

**Phase 1 — Re-exploration (single Explore agent).** Call sites drift. Before planning, confirm the table above is still accurate:

> Search for every use of `_pair_index` and every direct indexing into `params` at offset `≥ 2N+1`. Return a file/line table. Also verify the compiled-EOM signature has not already been changed.

**Phase 2 — Choose the end-state (single Plan agent).** The two realistic designs:

- **A. Extra positional arg:** `dq_dt_compiled(out, q, p, t, params, kappas)`. Minimal semantic change; `Symbolics.build_function` accepts multiple vector args. Downstream callers pass `kappas(prob)` explicitly.
- **B. Context struct:** `dq_dt_compiled(out, q, p, t, ctx::EOMContext)` where `ctx` bundles masses/charges/c/kappas (and anything future terms need). Cleaner extension point for time-dependent or parameter-promoted terms; bigger ripple.

Both remove κ from the flat `params` vector. Plan agent should pick one and lay out the migration order (build_function signature → solve.jl hot path → regularization → statistics → Makie ext → tests → docs).

**Phase 3 — Mechanical migration (foreground, incremental).** Do not parallelize; every step breaks the package until consistent. After each file, run the single-file test invocation from `CLAUDE.md` to catch regressions early.

**Phase 4 — Term-owned κ surface.** Once κ is out of `params`, give `NamedTerm` (or a small sibling) an optional `kappas::Vector{Float64}` slot. Rewrite `kappas(prob)` to delegate to the `:zollner` term. Update `_weber_pair_decomposition` and `_zollner_pair_decomposition` to close over that κ rather than slicing from `params`. This is the payoff step — downstream code now says `get_term(sys, :zollner).kappas[pair_idx]` (or `kappas(term, i, j)`) instead of `params[2N+1+_pair_index(i, j, N)]`.

## Acceptance criteria

- `grep _pair_index src/ ext/` returns only the definition in `hamiltonian_system.jl` and internal call sites inside the κ-owning term — no references in `solve.jl`, `statistics/`, or the Makie ext.
- `kappas(prob)` still works on `HamiltonianProblem` (construct-time accessor).
- Full test suite passes: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- `CLAUDE.md` §"Params vector layout" is rewritten to describe the new layout and the new EOM signature. `docs/src/internals.md` matches.
- `RELEASENOTES.md` carries a `### Breaking changes` section describing the EOM signature change and the upgrade path.

## What not to do

- Don't add a deprecation shim that keeps the old `params`-with-κ layout working alongside the new one. Pre-1.0 Julia semver treats minor bumps as breaking — a clean break is cheaper than a two-layout codebase.
- Don't bundle this with unrelated refactors (term-owned mass/charge, time-dependent Hamiltonian plumbing). Those are separate items in `_research/NextSteps.md` and each cascades differently.
- Don't skip the Plan-agent step. Choosing between designs A and B sets the scope of the migration; picking wrong makes Phase 3 painful.
