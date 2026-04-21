# Next Steps

Forward-looking notes after the System → Problem → Algorithm → Callbacks refactor
landed in v0.5.0. Captures intentional deferrals, debts surfaced during the
refactor, candidates for promotion from the research sandbox, and gaps in the
theory/code mapping.

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

## `_research/` promotion candidates for `theory/`

From [_research/investigations/](_research/investigations/):

| Source | Verdict | Proposed `theory/` name |
|--------|---------|-------------------------|
| [AngularMomentumRegularization.md](_research/investigations/AngularMomentumRegularization.md) | Ready — topological-obstruction theorem, seven failed approaches catalogued | `AngularMomentumNonRegularizability.md` |
| [TetheringImpossibility.md](_research/investigations/TetheringImpossibility.md) | Ready — non-existence theorem + numerical validation | `ExternalChargeStabilizationImpossibility.md` |
| [CriticalRadiusAndLikeChargeAttraction.md](_research/investigations/CriticalRadiusAndLikeChargeAttraction.md) | Borderline — foundational, cited by four other investigations; move or mirror | `CriticalRadiusAndSubCriticalAttraction.md` |
| [TransformedWeberHamiltonians.md](_research/investigations/TransformedWeberHamiltonians.md) | Needs cleanup — catalogue style, add synthesis section before promotion | `WeberHamiltonianTransformations.md` |

Also from [_research/investigations/](_research/investigations/) (studies rather than theorems):

| Source | Verdict | Proposed `theory/` name |
|--------|---------|-------------------------|
| [CollisionBounceRegularization.md](_research/investigations/CollisionBounceRegularization.md) | Ready — already cited from [src/callbacks.jl:42](src/callbacks.jl#L42); validated convergence | `CollisionBounceRegularization.md` |
| [ThreeBodyBoundStates.md](_research/investigations/ThreeBodyBoundStates.md) | Ready — ±± vs ±∓ dichotomy, Case B is canonical demo | `ThreeBodyBoundStates.md` |
| [FourPositiveChargeCrossInvestigation.md](_research/investigations/FourPositiveChargeCrossInvestigation.md) | Ready after trim — lead with collinear success; cross result to a paragraph | `FourBodySubcriticalBoundStates.md` |
| [ThreePositiveChargeInvestigation.md](_research/investigations/ThreePositiveChargeInvestigation.md) | Ready after trim — lead with the impossibility result; scans → appendix | `ThreeBodyCollinearAnalysis.md` |
| [RustImplementation.md](_research/investigations/RustImplementation.md) | Keep in `_research/` — speculative design, not yet implemented | — |
| [HypergeometricStructure.md](_research/investigations/HypergeometricStructure.md) | Symbolic-algebra study — no finalised theorems yet | — |

`_research/homology/` and the sub-critical exploration notes in
`_research/investigations/SubCriticalWeberExploration.md` have no
finalised theorems yet.

## Theory/code gaps worth closing

The code implements more than [theory/](theory/) specifies. Candidates for new
spec documents:

1. **Parameter vector specification.** Layout `[m₁…m_N, q₁…q_N, c, κ₁₂…]`,
   the `_pair_index` formula, and construction protocol. Today only in
   CLAUDE.md + [docs/src/internals.md](docs/src/internals.md); no spec-level
   reference for users writing custom ICs.
2. **Symmetric-projection defaults & diagnostics.** Document that
   `relaxation=0.25` in [src/types.jl:34](src/types.jl#L34) is a default, not
   a derived constant — `SemiExplicitIntegrator.md:144,162` derives 1/4 but
   the code generalises it to a tunable parameter. Also tolerance/max-iter
   choices and residual semantics.
3. **Regularization backend selection guide.** `:lifted_pair` (2D only) vs
   `:adaptive_cartesian` (all dims), Weber-velocity-term limitation, the 3D
   KS "diagnostic-only" status, interaction with collision bounce. Currently
   scattered across [docs/src/regularization.md](docs/src/regularization.md)
   and [theory/RegularizedIntegrationDesign.md](theory/RegularizedIntegrationDesign.md).
4. **Zöllner predictions & measurement.** Extend
   [theory/ZollnerElectrogravitationalTheory.md](theory/ZollnerElectrogravitationalTheory.md)
   with testable numerical predictions and how to infer `a` from observed
   orbital behaviour.
5. **Collision-singularity classification.** Head-on (ℓ=0, regularisable) vs
   spiralling (ℓ≠0, not regularisable). Partially covered in
   [theory/NonZeroRadialVelocityBoundICs.md](theory/NonZeroRadialVelocityBoundICs.md)
   §2.2; deserves its own note, ideally alongside the promoted
   `AngularMomentumNonRegularizability.md` above.

## Cleared in this pass

- Accessor functions (`masses`, `charges`, `speed_of_light`, `kappas`,
  `params`, `dims`, `degrees_of_freedom`) added to
  [docs/src/api/problem.md](docs/src/api/problem.md) and
  [docs/src/api/system.md](docs/src/api/system.md).
- [CLAUDE.md](CLAUDE.md) corrected: `regularization(prob)` and `zollner(prob)`
  were listed as accessors but do not exist.
- `benchmarks/` removed. It was a one-off Phase 0 gate artifact from commit
  `91922a0`; no CI or test consumer. Re-add with CI integration when real
  benchmarking work starts.
