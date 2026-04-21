# Next Steps

Forward-looking notes after the System → Problem → Algorithm → Callbacks refactor
landed in v0.5.0. Captures intentional deferrals, debts surfaced during the
refactor, candidates for promotion from the research sandbox, and gaps in the
theory/code mapping.

These are candidates, not commitments — prioritise as needs arise.

## Follow-ups on deliberate refactor deferrals

- **Finish the Phase-2 interface extraction.** `SymmetricProjectionIntegrator`
  and its buffers still live in [src/types.jl](../src/types.jl). A clean
  [src/integrators/symmetric_projection.jl](../src/integrators/symmetric_projection.jl)
  with `step!(cache, prob, alg, t, dt)` dispatch and a `projection_kernel` trait
  only pays off when a second base algorithm lands (ImplicitMidpoint,
  GaussLegendre). Defer until then.
- **Move κ out of the `params` vector to term-owned storage.** Pulled out
  into a standalone work prompt: [KappaStorageRefactor.md](KappaStorageRefactor.md).
  Breaking change to the compiled-EOM signature, not justified until another
  caller needs it.
- **Split [src/callbacks.jl](../src/callbacks.jl) into a directory** when a second
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

## Theory/code gaps worth closing

The code implements more than [theory/](../theory/) specifies. Candidates for new
spec documents:

1. **Parameter vector specification.** Layout `[m₁…m_N, q₁…q_N, c, κ₁₂…]`,
   the `_pair_index` formula, and construction protocol. Today only in
   CLAUDE.md + [docs/src/internals.md](../docs/src/internals.md); no spec-level
   reference for users writing custom ICs.
2. **Symmetric-projection defaults & diagnostics.** Document that
   `relaxation=0.25` in [src/types.jl:34](../src/types.jl#L34) is a default, not
   a derived constant — `SemiExplicitIntegrator.md:144,162` derives 1/4 but
   the code generalises it to a tunable parameter. Also tolerance/max-iter
   choices and residual semantics.
3. **Regularization backend selection guide.** `:lifted_pair` (2D only) vs
   `:adaptive_cartesian` (all dims), Weber-velocity-term limitation, the 3D
   KS "diagnostic-only" status, interaction with collision bounce. Currently
   scattered across [docs/src/regularization.md](../docs/src/regularization.md)
   and [theory/RegularizedIntegrationDesign.md](../theory/RegularizedIntegrationDesign.md).
4. **Zöllner predictions & measurement.** Extend
   [theory/ZollnerElectrogravitationalTheory.md](../theory/ZollnerElectrogravitationalTheory.md)
   with testable numerical predictions and how to infer `a` from observed
   orbital behaviour.
5. **Collision-singularity classification.** Head-on (ℓ=0, regularisable) vs
   spiralling (ℓ≠0, not regularisable). Partially covered in
   [theory/NonZeroRadialVelocityBoundICs.md](../theory/NonZeroRadialVelocityBoundICs.md)
   §2.2; deserves its own note.
