# Can Three Positive Charges Form a Bound State?

## Question

Can three positive charges, all below the critical radius, mutually attract and
form a stable oscillating triple? This extends the two-body sub-critical bound
pair (§5 of [CriticalRadiusAndLikeChargeAttraction.md](../theory/CriticalRadiusAndLikeChargeAttraction.md))
to three bodies.

## Short Answer

**No.** While three positive charges do attract each other below the critical
radius, a stable bound triple is impossible in practice. The system requires
perfect collinear symmetry (angular momentum ℓ=0 for every pair), and any
transverse perturbation — even as small as 10⁻⁸ — triggers exponentially
growing instability that destroys the system. In 1D (forced collinearity), the
system can be marginally stabilized by a large collision bounce radius, but with
poor energy conservation (~7%) that does not converge under dt refinement.

## Theoretical Analysis

### The ℓ=0 Constraint

Sub-critical bound oscillation (the "molecular state") requires **purely radial
motion** (ℓ=0) for each like-charge pair. This is established in
[CriticalRadiusAndLikeChargeAttraction.md](../theory/CriticalRadiusAndLikeChargeAttraction.md), §6:

> For ℓ ≠ 0, trajectories inside ρ spiral into the origin at infinite speed in
> finite time. These collisions are not regularizable.

For three particles A, B, C, the condition ℓ=0 for *every* pair constrains each
particle to move along the line connecting it to every other particle. Particle A
must simultaneously move along line AB and line AC. Unless these lines are
collinear, A cannot satisfy both constraints. Therefore **all three particles
must be collinear**.

### No Long-Range Binding

Unlike the ++- planetary atom (which has Coulomb attraction at large distances),
the +++ system has no attractive force above the critical radius. If any pair
separation exceeds ρ, that pair repels. There is no restoring mechanism to
re-confine escaped particles.

### Zöllner Parameter Has No Effect

The Zöllner coupling factor κ_ij = 1+a applies only to **unlike-sign** charge
pairs. For three positive charges, all pairs are like-sign and have κ=1
regardless of the mismatch parameter a. Zöllner cannot help stabilize a +++
system.

### Regularization Cannot Help

Neither Levi-Civita nor adaptive Cartesian regularization can regularize Weber's
velocity-dependent force (see [CollisionBounceRegularization.md](CollisionBounceRegularization.md)).
More fundamentally, the ℓ≠0 non-regularizability is a property of the
differential equation itself, not the numerical method: no coordinate transform
removes the singularity when angular momentum is nonzero (Frauenfelder & Weber
2024, Theorem 2.1).

## Numerical Experiments

Parameters throughout: c = 4, ρ = q₁q₂/(μc²). Integration uses the
unregularized symplectic integrator with collision bounce.

### 1D Experiments (Forced Collinearity)

**Setup**: `HamiltonianSystem(3, 1)`, three particles at positions [-s, 0, s], zero
momenta. Tests whether the radial dynamics alone can sustain a bound state.

#### Equal mass/charge baseline

m₁ = m₂ = m₃ = 1, q₁ = q₂ = q₃ = 1, s = 0.05.
All pairs sub-critical (inner: 0.05 < ρ=0.125, outer: 0.10 < ρ=0.125).

| bounce_r | dt | tspan | retcode | E_err max (%) | r range | Crossed ρ? |
|---|---|---|---|---|---|---|
| 0.02 | 1e-5 | (0,5) | Failure | 3962 | [0.001, 0.182] | YES |
| 0.03 | 1e-5 | (0,5) | Success | 91.5 | [0.014, 0.147] | YES |
| 0.035 | 1e-5 | (0,5) | Success | 23.0 | [0.026, 0.109] | no |
| 0.04 | 1e-5 | (0,5) | Success | 15.6 | [0.033, 0.103] | no |
| **0.045** | **1e-5** | **(0,20)** | **Success** | **7.4** | **[0.041, 0.100]** | **no** |

The best case (bounce_r=0.045) survived 2,000,000 steps over t=20 with all pairs
staying in [0.041, 0.100]. Energy oscillates rather than drifts monotonically:

| Time | Energy error (%) |
|---|---|
| t=1 | 2.5 |
| t=5 | 0.7 |
| t=10 | 0.1 |
| t=15 | 2.6 |
| t=20 | 5.9 |

However, **dt convergence is absent**: halving dt from 5e-5 to 5e-6 only
improves the max error from 16.3% to 13.9%. The error is dominated by the
collision bounce mechanism (each bounce is a non-symplectic perturbation), not by
discretization. This means the ~7% error floor cannot be reduced by refining the
time step.

#### Mass asymmetry (heavy center)

masses=[1, 10, 1], charges=[1, 1, 1], s=0.05.
ρ_heavy-light = 0.069, ρ_light-light = 0.125. All pairs sub-critical.

| bounce_r | dt | tspan | retcode | E_err max (%) | Outcome |
|---|---|---|---|---|---|
| 0.02 | 1e-5 | (0,5) | Failure | 1082 | Cascade bounces |

**Worse than equal mass.** The asymmetric reduced masses create mismatched
critical radii. The heavy-light pairs have smaller ρ and oscillate differently
from the light-light pair, leading to faster desynchronization.

#### Charge asymmetry (super-critical outer pair)

charges=[0.5, 3.0, 0.5], masses=[1, 1, 1], s=0.05.
ρ_center-outer = 0.188 (sub-critical), ρ_outer-outer = 0.031 (outer pair at
2s=0.10 is **super-critical** and repulsive).

| bounce_r | dt | tspan | retcode | E_err max (%) | Outcome |
|---|---|---|---|---|---|
| 0.02 | 1e-5 | (0,5) | Failure | 2840 | Center-outer pairs cross ρ |

**Also fails.** The strong center charge (q=3) creates large forces on the outer
particles, which get pushed hard against each other. Despite the outer-outer pair
being repulsive, the force imbalance drives the center-outer pairs past their
critical radius.

Combined mass+charge asymmetry (M_center=10, q_center=3, m_outer=1, q_outer=0.5)
also failed (15955% energy error at t=0.8).

### 2D Transverse Stability Tests

**Setup**: `HamiltonianSystem(3, 2)`, collinear initial positions along x-axis with
bounce_r=0.045 (the best 1D case). Small perturbation added in y-direction.

#### Position perturbation (middle particle nudged in y)

| ε_y | retcode | t_final | E_err max (%) | y_max |
|---|---|---|---|---|
| 0 (baseline) | Success | 20.0 | 7.2 | 0 |
| 1e-8 | Failure | 2.36 | 2522 | 2.6e-2 |
| 1e-6 | Failure | 0.073 | 5548 | 1.4e-2 |
| 1e-4 | Failure | 0.043 | 7316 | 1.4e-2 |
| 1e-3 | Failure | 0.028 | 4811 | 1.2e-2 |
| 1e-2 | Failure | 0.014 | 7839 | 1.4e-2 |

#### Momentum perturbation (middle particle kicked in y)

| ε_p | retcode | t_final | E_err max (%) | y_max |
|---|---|---|---|---|
| 1e-6 | Failure | 2.00 | 6165 | 1.4e-2 |
| 1e-4 | Failure | 0.071 | 6967 | 1.2e-2 |
| 1e-3 | Failure | 0.056 | 5251 | 1.2e-2 |
| 1e-2 | Failure | 0.042 | 8599 | 1.4e-2 |
| 0.1 | Failure | 0.027 | 8843 | 1.2e-2 |

**Every perturbation, no matter how small, destroys the system.** The transverse
displacement grows exponentially: at ε_y=1e-8, it amplifies by a factor of
~10⁶ (from 1e-8 to 2.6e-2) before the integrator fails. This confirms
that the collinear +++ configuration is **linearly unstable** to transverse
perturbations.

The instability mechanism: any transverse displacement gives pairs nonzero
angular momentum (ℓ≠0). Below the critical radius, ℓ≠0 trajectories spiral
inward at infinite speed in finite time. The integrator cannot resolve this
singularity, and energy conservation is immediately destroyed.

## Comparison with the Two-Body Case

The two-body sub-critical pair works because:
1. In 1D, it is a single pair with no cascade bounces
2. The bounce energy error converges as O(dt²) — reaching 0.005% at dt=T/2000
3. In 2D, there is no transverse instability to worry about (a single pair
   oscillating head-on has no mechanism for angular momentum generation)

The three-body case adds two fatal complications:
1. **Cascade bounces**: when pair (1,2) bounces, particle 1 moves toward
   particle 3, triggering a bounce of pair (1,3), and so on. Each bounce is
   non-symplectic, and the cascade accumulates error that cannot be reduced by
   dt refinement.
2. **Transverse instability**: any departure from perfect collinearity gives
   some pair ℓ≠0, triggering the non-regularizable inward spiral.

## Conclusion

Three positive charges below the critical radius do mutually attract (via the
negative effective mass mechanism), but they **cannot form a stable bound
triple**. The combination of the ℓ=0 requirement (forcing collinearity), the
absence of long-range binding, cascade collision bounces, and violent transverse
instability make this configuration fundamentally non-viable.

This contrasts with the ++- planetary atom ([ThreeBodyBoundStates.md](ThreeBodyBoundStates.md)),
which succeeds because:
- The orbiter binding is Coulomb-like (ℓ≠0 allowed at super-critical distances)
- Timescale separation (T_nuc ≪ T_orb) prevents 3-body close encounters
- Only one pair (the nucleus) requires collision bounce

## References

- [CriticalRadiusAndLikeChargeAttraction.md](../theory/CriticalRadiusAndLikeChargeAttraction.md) — critical radius theory
- [ThreeBodyBoundStates.md](ThreeBodyBoundStates.md) — original +++ investigation (Case A) and ++- success (Case B)
- [CollisionBounceRegularization.md](CollisionBounceRegularization.md) — bounce mechanism and limitations
- [FourPositiveChargeCrossInvestigation.md](FourPositiveChargeCrossInvestigation.md) — 4-body investigation: cross configuration (fails) and collinear chain (**stable bound state found**)
- Frauenfelder, U., Weber, J. "A mathematical description of the Weber nucleus." *Anal. Math. Phys.* **14**:31 (2024)
