# Three-Body Bound States in Weber Electrodynamics

## Overview

We systematically investigated whether compact, bound 3-body states exist in
Weber electrodynamics. Two configurations were studied:
- **Case A (+++)**: Three positive charges, all pairs sub-critical
- **Case B (++-)**: Two positive charges forming a nucleus + one negative orbiter

**Main result**: Case A does not produce stable bound states. Case B — the
"planetary atom" — produces **robust bound orbits** that persist for 200+ time
units (55+ orbiter periods). The Zöllner mismatch parameter `a` can circularize
eccentric orbits.

Parameters throughout: m₁ = m₂ = m₃ = 1, |q| = 1, c = 4, ρ = 0.125.

## Case A: Three Positive Charges (+++)

### Equilateral Triangle (s < ρ, p = 0)

All 3 pairs sub-critical. Particles released from rest at equilateral vertices
with side length s. COM at origin.

| s | retcode | E_err (%) | max r | Steps |
|---|---|---|---|---|
| 0.05 | Failure | 111% | 0.050 | 37 |
| 0.10 | Failure | 304% | 0.100 | 169 |

**Failure mode**: The integrator diverges after ~37 steps (t ≈ 0.004). All 3
particles collapse toward the COM simultaneously. The 3-body near-collision
creates forces that the midpoint iteration cannot handle. Varying dt (1e-4 to
1e-5) and bounce_r (0.02 to 0.04) does not help — smaller dt actually produces
**worse** energy errors (up to 1372%).

Adding rigid-body rotation (ω = 0.1 to 10.0) does not stabilize the system.

### Collinear Configuration

Particles at (-s, 0), (0, 0), (+s, 0). Pairs (1,2) and (2,3) at distance s,
pair (1,3) at distance 2s.

| s | 2s vs ρ | retcode (t=2) | E_err (%) | Bound? |
|---|---|---|---|---|
| 0.05 | 0.10 < ρ | Failure | 251% | — |
| 0.07 | 0.14 > ρ | Failure | 139% | — |
| 0.10 | 0.20 > ρ | Success | 0.002% | Yes (short-lived) |
| 0.11 | 0.22 > ρ | Success | 0.001% | Yes (short-lived) |

At s = 0.10–0.11, the inner pairs oscillate sub-critically while the outer pair
is super-critical. This looks like a "linear molecule" with the middle particle
shared. However, **long-duration runs fail**: s=0.10 at t=5.0 produces 81%
energy error. The symmetry eventually breaks, causing one pair to reach very
small r and triggering cascade bounces.

### Case A Conclusion

Three-body sub-critical dynamics are fundamentally unstable for our integrator.
The simultaneous multi-pair close encounters create conditions where:
1. The midpoint iteration diverges
2. Cascade collision bounces (bouncing pair A pushes into pair B) accumulate error
3. Even symmetry-protected configurations eventually break symmetry numerically

**No stable +++ bound state was found.**

## Case B: Planetary Atom (++-)

### Configuration

- **Nucleus**: Particles 1, 2 (charges +1, +1) at ±r_nuc/2 on x-axis, zero momenta
- **Orbiter**: Particle 3 (charge -1) at (0, R) with tangential velocity v = η·v_circ
- v_circ = √(Q_eff / (μ_orb · R)), Q_eff = 2, μ_orb = 2m·m₃/(2m+m₃) = 2/3
- All positions/momenta adjusted for COM frame and zero total momentum
- Integration: dt = 1e-4, bounce_r = 0.02, regularization disabled

### Eta Scan (R = 1.0, a = 0.1, t = 50)

| η | retcode | E_err (%) | Orbiter range | Nucleus max | Status |
|---|---|---|---|---|---|
| 0.3 | Failure | 7.7 | [0.07, 1.00] | 0.1252 | Nucleus broken |
| 0.4 | Failure | 740 | [0.02, 3.71] | 0.0528 | Failure |
| 0.5 | Success | 3.4 | [0.29, 1.01] | 0.0511 | **Bound** |
| 0.6 | Success | 11.1 | [0.52, 1.01] | 0.0506 | **Bound** |
| 0.7 | Success | 1.5 | [0.93, 1.04] | 0.0504 | **Bound** |
| 0.8 | Success | 1.1 | [1.00, 1.80] | 0.0503 | **Bound** |
| 0.9 | Success | 1.0 | [1.00, 4.67] | 0.0503 | **Bound** |
| 1.0 | Failure | 1.0 | [1.00, 29.4] | 0.0503 | Escaped |
| ≥1.1 | Failure | ~1% | [1.00, 30+] | 0.050 | Escaped |

**Bound range: η ∈ [0.5, 0.9]**. Below η = 0.5, the orbiter falls too close
to the nucleus and disrupts it. Above η = 1.0, the orbiter escapes.

The nearly-circular orbit at **η = 0.7** has orbiter range [0.93, 1.04] — only
±5% variation from R.

### Zöllner Parameter Scan (R = 1.0, t = 100)

**At η = 0.7 (tight orbit):**

| a | E_err (%) | Orbiter range | Effect |
|---|---|---|---|
| 0.00 | 1.19 | [0.995, 1.178] | Baseline |
| 0.05 | 1.32 | [0.993, 1.073] | Tighter |
| 0.10 | 1.46 | [0.926, 1.039] | Tighter |
| 0.20 | 1.78 | [0.779, 1.025] | Eccentric |
| 0.50 | 5.69 | [0.522, 1.240] | Very eccentric |
| 1.00 | 6.72 | [0.315, 2.400] | Too eccentric |

**At η = 0.8 (wider orbit):**

| a | E_err (%) | Orbiter range | Effect |
|---|---|---|---|
| 0.00 | 0.99 | [0.999, 2.455] | Wide ellipse |
| 0.05 | 1.06 | [0.999, 2.076] | Tighter |
| 0.10 | 1.14 | [0.999, 1.798] | Moderate |
| 0.20 | 1.33 | [0.998, 1.421] | Compact |
| 0.50 | 2.12 | [0.849, 1.034] | **Near-circular** |
| 1.00 | 20.3 | [0.504, 1.074] | Eccentric, poor energy |

**Key finding**: Zöllner circularizes the orbit. At η = 0.8, a = 0.5 transforms
a wide ellipse [1.0, 2.5] into a near-circular orbit [0.85, 1.03]. This is the
Zöllner gravitational residual acting as additional binding.

### Long-Duration Stability (η = 0.8, a = 0.1, t = 200)

- **retcode: Success** (2,000,001 steps)
- Nucleus intact: r₁₂ ∈ [0.017, 0.051] (well below ρ = 0.125)
- Orbiter bounded: r₁₃ ∈ [0.69, 2.67]
- Energy error: 5.3% (drifting but bounded over 200 time units)
- 24,865 periapsis passages counted

The planetary atom survives for **200 time units** (~55 orbiter periods). The
energy drift is due to the collision bounce at the nucleus — each bounce
introduces a small non-symplectic perturbation that accumulates.

### Timescale Separation

The nucleus oscillation period T_nuc ≈ 2√2·r_nuc/c = 0.035 is 100× faster
than the orbiter orbital period T_orb ≈ 2πR/v_circ ≈ 3.6. This natural
timescale separation keeps the orbiter far from the nucleus during most of its
orbit, preventing the 3-body close-encounter instability that plagues Case A.

## Physical Interpretation

The planetary atom model (Case B) realizes the configuration described in
§10 of CriticalRadiusAndLikeChargeAttraction.md: a sub-critical nucleus of two
like charges, orbited by an unlike charge under combined electrodynamic and
(optionally) Zöllner gravitational forces.

**What makes it work:**
1. The nucleus binding is purely Weber (sub-critical molecular oscillation)
2. The orbiter binding is Coulomb-like (attractive force from net charge 2q)
3. Scale separation (R >> r_nuc) prevents the orbiter from disrupting the nucleus
4. Zöllner coupling (a > 0) strengthens the orbiter binding without affecting the
   nucleus (like-charge pairs have κ = 1 regardless of a)

**Why Case A fails:**
1. No natural hierarchy of length scales — all pairs at similar separations
2. No attractive pairs to provide long-range binding
3. 3-body simultaneous close encounters overwhelm the integrator
4. Even with collision bounce, cascade bounces destroy energy conservation

## Notebook Configurations

For the demonstration notebook, the following configurations are recommended:

| Run | η | a | R | Duration | Shows |
|---|---|---|---|---|---|
| 1 | 0.7 | 0.0 | 1.0 | 100 | Pure Weber planetary atom (near-circular) |
| 2 | 0.8 | 0.0 | 1.0 | 100 | Weber with wider orbit |
| 3 | 0.8 | 0.1 | 1.0 | 100 | Zöllner tightens the orbit |
| 4 | 0.8 | 0.5 | 1.0 | 100 | Zöllner circularization (a=0 vs a=0.5) |

All with r_nuc = 0.05, dt = 1e-4, bounce_r = 0.02.

## References

- Weber, Sixth Memoir (1871), §§9.8–9.17
- Frauenfelder & Weber, *Anal. Math. Phys.* **14**:31 (2024)
- See [CriticalRadiusAndLikeChargeAttraction.md](../theory/CriticalRadiusAndLikeChargeAttraction.md) §10
- See [CollisionBounceRegularization.md](CollisionBounceRegularization.md) for 2-body bounce validation
- See [ThreePositiveChargeInvestigation.md](ThreePositiveChargeInvestigation.md) for extended +++ investigation with mass/charge asymmetry and transverse stability analysis
- See [FourPositiveChargeCrossInvestigation.md](FourPositiveChargeCrossInvestigation.md) for 4-body investigation: cross (fails) and collinear chain (**stable +++ bound state found**)
