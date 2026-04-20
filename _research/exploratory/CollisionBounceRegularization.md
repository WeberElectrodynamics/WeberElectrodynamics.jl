# Collision Bounce & Sub-Critical Like-Charge Oscillation

## Physics: Sub-Critical Regime (r₀ < ρ)

- **Critical radius**: ρ = q₁q₂/(μc²). For m1=m2=1, q1=q2=1, c=4: ρ = 0.125.
- Below ρ, effective inertial mass μ_eff = μ(1 - ρ/r) is **negative** → like charges attract.
- Particles oscillate radially between r₀ and r=0 with ṙ² → 2c² as r→0.
- Energy E = k/r₀ > 0 is positive, yet the system is permanently bound.
- Period estimate: T ≈ 2√2·r₀/c.

## Regularizability Constraint (Frauenfelder & Weber 2024)

- **ℓ = 0 (head-on)**: collision at finite speed √2·c, C⁰-continuable. Regularizable.
- **ℓ ≠ 0 (spiralling)**: reaches r=0 at **infinite speed** in finite time. NOT regularizable by any smooth coordinate-time transform. No periodic orbits exist inside ρ.
- For point charges in our system, ℓ=0 is the physically correct case.

## Why Levi-Civita (lifted_pair) Fails for Weber

The LC backend transforms coordinates (q→u) but evaluates forces in **Cartesian space** then transforms back. This does NOT regularize the velocity-dependent Weber force terms, which still diverge as 1/r². Result: energy error drifts and accumulates (non-symplectic behavior). The LC approach works for Coulomb/Kepler (position-only forces) but not for Weber.

## Why Adaptive Cartesian Also Struggles

The implicit midpoint fixed-point iteration diverges near r=0 due to the 1/r² force singularity. The iteration simply doesn't converge when the force changes too rapidly across the step.

## Solution: Collision Bounce (Pre-Step Reflection)

### Concept
When pair separation r < bounce_radius, reflect the relative coordinate through the origin: q_rel → -q_rel (momenta unchanged). This analytically continues the C⁰ trajectory without integrating through the singularity.

### Implementation
- `CollisionBounce(r)` callback; pass via
  `solve(prob, alg; callbacks=CollisionBounce(r))` (omit callback to disable)
- `_apply_collision_bounces!` called at START of each `step!` in solve.jl
- Works for any dimension via `_current_pair_r` and `_reflect_pair!` helpers
- Energy is exactly preserved by the reflection (only positions change, not momenta)

### Why Symplectic + Bounce Works Best
The unregularized Strang splitting integrator preserves a modified Hamiltonian (symplectic), giving **bounded** (non-accumulating) energy error. Combined with the bounce, this gives constant ~0.01% error over 100+ periods. The LC backend's error drifts because the coordinate transform breaks symplecticity.

## Validated Parameters

| Parameter | Run 1 (r₀=0.05) | Run 2 (r₀=0.10) |
|---|---|---|
| dt | 1e-4 (fixed, absolute) | 1e-4 (fixed, absolute) |
| bounce_r | 0.02 | 0.02 |
| Algorithm | SymmetricProjectionIntegrator | SymmetricProjectionIntegrator |
| Energy error (100 periods) | 0.014% | 0.005% |
| max |ṙ|/c | 1.20 | 1.38 |

### Parameter Selection Lessons
- **dt must be absolute** (e.g. 1e-4), NOT scaled to period (T/40). Period-scaled dt fails because dt becomes too large for midpoint convergence near the bounce radius.
- **bounce_r = 0.02** is the sweet spot: large enough that the midpoint iteration converges comfortably, small enough that max|ṙ|/c approaches √2 ≈ 1.414.
- bounce_r = 0.01 works with dt=T/2000 for r₀=0.05 but is fragile for larger r₀.
- bounce_r = 0.03 gives even better energy (0.001%) but truncates the trajectory further from r=0.

## Second-Order Convergence

Energy error scales as O(dt²) — halving dt quarters the error:

| dt | Energy error (%) |
|---|---|
| T/250 | 0.51 |
| T/500 | 0.18 |
| T/1000 | 0.033 |
| T/2000 | 0.005 |

This confirms the integrator is second-order (Strang splitting).

## Test Coverage
- Test "Sub-critical like-charge oscillation (collision bounce)" in test_regularization.jl
- Uses r₀=0.05, bounce_r=0.01, dt=T/2000, 10 periods
- Asserts: Success, all finite, r ≤ r₀, energy < 1%, ≥10 oscillation minima

## Files Modified
- `src/types.jl` — `collision_bounce_radius` field + constructor kwarg
- `src/solve.jl` — `_current_pair_r`, `_reflect_pair!`, `_apply_collision_bounces!`, pre-step bounce in `step!`
- `test/test_regularization.jl` — sub-critical oscillation test
- `examples/critical_radius_reference.ipynb` — Runs 1 & 2 use bounce, Run 2 changed from ℓ≠0 to ℓ=0

## References

- Frauenfelder, U., Weber, J. "A mathematical description of the Weber nucleus as a classical and quantum mechanical system." *Anal. Math. Phys.* **14**:31 (2024). [DOI: 10.1007/s13324-024-00891-5](https://doi.org/10.1007/s13324-024-00891-5)
- See also: [../theory/CriticalRadiusAndLikeChargeAttraction.md](../theory/CriticalRadiusAndLikeChargeAttraction.md), [../theory/Regularization.md](../theory/Regularization.md)
