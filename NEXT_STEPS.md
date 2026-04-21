# Next Steps

Forward-looking notes after the System → Problem → Algorithm → Callbacks refactor
landed in v0.5.0. Captures intentional deferrals, debts surfaced during the
refactor, and unblocked work items.

These are candidates, not commitments — prioritise as needs arise.

## Follow-ups on deliberate refactor deferrals

- **Finish the Phase-2 interface extraction.** `SymmetricProjectionIntegrator`
  and its buffers still live in [src/types.jl](src/types.jl). A clean
  [src/integrators/symmetric_projection.jl](src/integrators/symmetric_projection.jl)
  with `step!(cache, prob, alg, t, dt)` dispatch and a `projection_kernel` trait
  only pays off when a second base algorithm lands (ImplicitMidpoint,
  GaussLegendre). Defer until then.
- **Ship `kinetic_term` and `coulomb_term` builders.** Small, orthogonal, and
  unlocks non-Weber custom-Hamiltonian demos (`H = kinetic_term(...) +
  coulomb_term(...)`). Good first-issue size.
- **Move κ out of the `params` vector to term-owned storage.** Current layout
  `[masses…, charges…, c, κ₁₂…]` is insulated by `kappas(prob)` but still
  couples downstream code via `_pair_index`. Term-owned κ would let us delete
  that coupling; requires coordinated changes to the compiled-EOM signature.
- **Split [src/callbacks.jl](src/callbacks.jl) into a directory** when a second
  callback lands (e.g. `EnergyObserver`, `RegularizationDiagnosticsObserver`).

## Unblocked by the refactor (previously out of scope)

- **Lagrange-multiplier / holonomic constraints.** The layered architecture
  was explicitly designed to prepare for these.
- **Time-dependent Hamiltonians.** `t` is threaded through every compiled-EOM
  signature but has zero real users today. A minimal time-dependent term would
  exercise the plumbing.
- **Custom user Hamiltonians beyond Weber/Zöllner.** `HamiltonianSystem(H, q,
  p; params, t)` is the intended entry point; needs a worked example and a
  regression test with non-Weber physics.
- **Promoting physical parameters (mass, charge) to generalized coordinates.**
  One of the motivating use cases for the refactor; still not done.

## Research-sandbox hygiene

Four `_research/` files still construct `params` manually with
`[masses; charges; [c]; kappas]` and call `hamiltonian_compiled` directly.
They work today because the layout is unchanged, but they'd break the moment
κ moves term-owned:

- [_research/homology/03_contact_reeb_numerics/reeb_2body.jl](_research/homology/03_contact_reeb_numerics/reeb_2body.jl)
- [_research/homology/03_contact_reeb_numerics/star_center_search.jl](_research/homology/03_contact_reeb_numerics/star_center_search.jl)
- [_research/FourBodyTwoPlusTwoMinus/02_symmetry_reduction/verify_symmetries.jl](_research/FourBodyTwoPlusTwoMinus/02_symmetry_reduction/verify_symmetries.jl)
- [_research/FourBodyTwoPlusTwoMinus/11_contact_reeb/star_shaped_check.jl](_research/FourBodyTwoPlusTwoMinus/11_contact_reeb/star_shaped_check.jl)

Migrate them to `params(prob)` and the accessor API.
