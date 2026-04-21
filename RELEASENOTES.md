<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->

### Breaking changes

This release completes the layered **System → Problem → Algorithm → Callbacks** architectural refactor. The public API is reshaped end-to-end; no deprecation shims are provided.

- **Renamed top-level types.** `WeberSystem` → `HamiltonianSystem`, `WeberProblem` → `HamiltonianProblem`, `WeberSolution` → `HamiltonianSolution`, `WeberIntegrator` → `HamiltonianIntegrator`, `WeberAlgorithm` → `HamiltonianAlgorithm`.
- **`HamiltonianSystem` is now built from terms.** `HamiltonianSystem(n, dims)` still works (defaults to Weber + Zöllner). A new generic constructor `HamiltonianSystem(H, q, p; params, t)` and term builders `weber_term(...)` and `zollner_term(...)` let users assemble custom Hamiltonians. Each system exposes `system.terms::Vector{NamedTerm}` with per-term `pair_decomposition` closures driving statistics.
- **`HamiltonianProblem` fields dropped.** `masses`, `charges`, `c`, `kappas`, `regularization`, `zollner` are no longer fields on the problem. Construct with kwargs as before; access via exported accessors: `masses(prob)`, `charges(prob)`, `speed_of_light(prob)`, `kappas(prob)`, `regularization(prob)`, `zollner(prob)`, `params(prob)`, `n_particles(sys)`, `dims(sys)`.
- **Regularization is now an algorithm wrapper, not a problem option.** Replace `regularization=RegularizationOptions(...)` on `HamiltonianProblem` with `RegularizedIntegrator(SymmetricProjectionIntegrator(); r_on_factor=..., r_off_factor=..., backend=..., ...)` passed to `solve`.
- **Collision bounce is now a callback.** Replace `regularization=RegularizationOptions(collision_bounce_radius=r)` with `solve(prob, alg; callbacks=CollisionBounce(r))`.
- **Compiled EOM signature includes `t`.** `dq_dt_compiled(out, q, p, t, params)` and `dp_dt_compiled(out, q, p, t, params)` — direct callers must thread the time argument.
- **κ (Zöllner per-pair coupling) removed from `params`.** The parameter vector is now `[m₁…mₙ, q₁…qₙ, c]` (length `2N+1`); κ values live on `HamiltonianProblem.kappas` and are accessed via `kappas(prob)` (same accessor as before, now a direct field read). Compiled EOMs take κ as a separate positional argument: `dq_dt_compiled(out, q, p, t, params, kappas)`, `dp_dt_compiled(out, q, p, t, params, kappas)`, and `hamiltonian_compiled(q, p, t, params, kappas)`. New per-pair accessor `kappa(prob, i, j)` replaces manual `_pair_index` arithmetic in downstream code. Code that spliced κ onto `params` literals, sliced `params[(2N+2):end]`, or called the compiled EOMs directly must be updated.

### Added

- `RegularizedIntegrator{BaseAlg}` algorithm wrapper preserving symplectic projection via the base algorithm's `projection_kernel`.
- `CollisionBounce(radius)` pre-step callback. `CallbackSet(...)` composition.
- Accessor API on `HamiltonianSystem` and `HamiltonianProblem` to insulate downstream code from struct layout.
- `NamedTerm{S}` with per-term compiled EOMs and optional `pair_decomposition(i, j, q, p, params, kappas)` for pair-wise statistics.
- `kappa(prob, i, j)` per-pair accessor, co-located with `_pair_index` in `hamiltonian_system.jl`.
- Symbolic builders: `weber_term`, `zollner_term`.
- Regression fixtures (`test/regression/fixtures/*.jld2`) with 1e-12 numerical-equivalence validation across every refactor phase.

### Migration guide

```julia
# Before
prob = WeberProblem(sys, (0.0, T), q0, p0;
    masses=ms, charges=qs, c=1.0,
    regularization=RegularizationOptions(enabled=true, backend=:lifted_pair,
                                          collision_bounce_radius=0.02))
sol = solve(prob, SymmetricProjectionIntegrator())

# After
prob = HamiltonianProblem(sys, (0.0, T), q0, p0; masses=ms, charges=qs, c=1.0)
alg  = RegularizedIntegrator(SymmetricProjectionIntegrator(); backend=:lifted_pair)
sol  = solve(prob, alg; callbacks=CollisionBounce(0.02))

# Field access
ms  = masses(prob)          # was: prob.masses
κs  = kappas(prob)          # was: prob.kappas (field access; slicing params[tail] no longer works)
κij = kappa(prob, 1, 2)     # per-pair — preferred over kappas(prob)[_pair_index(...)]
n   = n_particles(sys)      # was: sys.n_particles

# Direct compiled-EOM calls now take kappas as a separate positional arg
sys.dp_dt_compiled(out, q, p, t, params(prob), kappas(prob))
```
