# 09 Three-Body Weber Bound Orbit Search

Agent 09/10 -- first-ever systematic 3-body Weber electrodynamics survey.

## Summary

179 total runs across 5 categories. 27 runs classified as "good bound" (retcode=Success,
all pair distances < 10x initial, energy drift < 1%).

| Category   | Runs | Success | Bound | Good bound (E_drift<1%) |
|------------|------|---------|-------|-------------------------|
| 2+/1-      |   69 |      17 |    11 |                       1 |
| 1+/2-      |   28 |      12 |     9 |                       1 |
| +++        |   24 |       2 |     2 |                       0 |
| asym       |   28 |      15 |    11 |                      10 |
| perturb    |   30 |      23 |    15 |                      15 |

**Key finding**: asymmetric (helium-like) configurations dominate the good-bound count.
The equal-mass symmetric configurations (equilateral, collinear, isosceles) are overwhelmingly
unstable. Same-sign (all-repulsive) configurations produced zero good-bound orbits.

## Bug fix

The original script had a type error: `angle_sep::Float64 = pi` failed because Julia's `pi` is
`Irrational{:pi}`, not `Float64`. Fixed by wrapping with `Float64(pi)` in both
`ic_generators_3body.jl` and `three_body_survey.jl`.

## Results by category

### Category 1: 2+/1- (equal masses, [+1,+1,-1])

69 runs. Only 1 good-bound orbit (run 69: planetary atom, R=2.0, eta=0.9, c=4.0, E_drift=0.87%).

- **Equilateral rotating** (24 runs): 0 good-bound. Most fail the integrator (close encounters).
  c=1 is particularly bad; c=4 sometimes completes but particles escape (max_r >> initial).
  Low eta (0.1) at c=4 produces unbound escape orbits. Higher eta at c=1 causes integrator failure.
- **Collinear +/-/+** (9 runs): 0 good-bound. All fail. The Euler-like configuration is unstable.
- **Planetary atom** (36 runs): 11 bound, 1 good-bound. The planetary atom geometry (two
  like-charge particles close together, one opposite-charge orbiter far away) is the only
  2+/1- geometry that produces bound orbits. However, energy drift is typically 1-23%,
  with only the highest eta (0.9) and largest orbit (R=2.0) achieving sub-1% drift.
  **All bound planetary orbits require c=4** -- every c=1 run fails immediately (12 steps only,
  likely due to the nucleus pair at separation 0.05 being far below the critical distance).

### Category 2: 1+/2- (equal masses, [-1,-1,+1])

28 runs. By charge symmetry, equal-mass [-1,-1,+1] and [+1,+1,-1] produce identical physics
(just relabeled). The equilateral results are identical to Category 1 (as expected). Planetary
atom results mirror Category 1 exactly (same E_drift values). 1 good-bound (run 97: R=2.0,
eta=0.9, c=4.0, E_drift=0.87%).

**Note on symmetric vs non-symmetric**: Categories 1 and 2 use equal masses [1,1,1], making
the charge-sign flip an exact symmetry. The results confirm this -- every run pair produces
identical retcodes and energy drifts. There is no physical distinction between 2+/1- and 1+/2-
at equal mass.

### Category 3: Same-sign +++ (all repulsive)

24 runs. 0 good-bound.

- **Rotating equilateral** (18 runs): All fail. Sub-critical distances (0.05--0.10) with
  bounce radius 0.02 cannot stabilize the 3-body repulsive system.
- **Collinear breathing** (4 runs): 2 nominally "bound" (collin_breath at d=0.05, c=1 and c=4)
  but with catastrophic energy drift (109--152%). These are numerical artifacts, not real bound
  states.
- **Super-critical rotating** (2 runs): Both fail quickly.

**Conclusion**: Three same-sign charges show no evidence of bound states, even at sub-critical
Weber distances. This contrasts with the 2-body case where sub-critical same-sign pairs do form
bound orbits -- the third body destabilizes the system.

### Category 4: Asymmetric (helium-like, unequal charges/masses)

28 runs. 10 good-bound -- **the most productive category by far**.

Best configurations:

| Run | Config                        | Charges   | Masses    | c   | eta  | R   | E_drift% | max_r  |
|-----|-------------------------------|-----------|-----------|-----|------|-----|----------|--------|
| 132 | He_q2m10_R2.0_api_eta0.5      | [+2,-1,-1]| [10,1,1]  | 1.0 | 0.50 | 2.0 | 0.00003  | 4.000  |
| 136 | He_q2m10_R2.0_api_eta0.75     | [+2,-1,-1]| [10,1,1]  | 1.0 | 0.75 | 2.0 | 0.00000  | 4.000  |
| 133 | He_q2m10_R2.0_a2pi3_eta0.5    | [+2,-1,-1]| [10,1,1]  | 1.0 | 0.50 | 2.0 | 0.020    | 5.193  |
| 137 | He_q2m10_R2.0_a2pi3_eta0.75   | [+2,-1,-1]| [10,1,1]  | 1.0 | 0.75 | 2.0 | 0.008    | 9.829  |
| 146 | He_c4_R0.5_eta0.7             | [+2,-1,-1]| [10,1,1]  | 4.0 | 0.70 | 0.5 | 0.009    | 5.568  |
| 147 | He_c4_R1.0_eta0.7             | [+2,-1,-1]| [10,1,1]  | 4.0 | 0.70 | 1.0 | 0.001    | 14.695 |
| 148 | He_c4_R0.5_eta0.8             | [+2,-1,-1]| [10,1,1]  | 4.0 | 0.80 | 0.5 | 0.003    | 8.915  |
| 149 | He_c4_R1.0_eta0.8             | [+2,-1,-1]| [10,1,1]  | 4.0 | 0.80 | 1.0 | 0.001    | 11.861 |
| 128 | He_q2m10_R2.0_api_eta0.25     | [+2,-1,-1]| [10,1,1]  | 1.0 | 0.25 | 2.0 | 0.220    | 26.409 |
| 143 | inv_He_q-2_R2.0_eta0.75       | [-2,+1,+1]| [1,1,1]   | 1.0 | 0.75 | 2.0 | 0.003    | 29.238 |

**Runs 132 and 136 are the standout results**: energy drift effectively zero (< 0.001%) over
t=50 with dt=1e-3. These are heavy-nucleus helium-like configurations with two light "electrons"
in anti-podal orbits (angle_sep=pi), orbit radius R=2.0. The heavy nucleus (m=10) stays nearly
stationary, providing a stable potential well.

**Non-symmetric highlights**:
- Run 133 (2pi/3 angle separation instead of pi): still bound with E_drift=0.02%, showing the
  orbit survives breaking the electron-electron symmetry.
- Run 137 (2pi/3, eta=0.75): bound with E_drift=0.008% but max_r expands to ~10.
- Run 143 (inverted helium [-2,+1,+1], all equal mass): bound but with max_r=29.2, suggesting
  a loosely bound state. This is the only non-heavy-nucleus bound config.

**Inverted helium** ([-2,+1,+1], equal mass): only 1 of 6 runs bound (run 143). Much less
stable than the heavy-nucleus version -- the equal-mass nucleus cannot anchor the system.

### Category 5: Perturbation tests

30 runs testing stability of the 12 best bound orbits under random perturbations
(delta = 0.01, 0.05, 0.1 applied to both q and p).

**Strongly stable** (survive all perturbation levels):
- He_q2m10_R2.0_api_eta0.75 (runs 150--152): all 3 perturbation levels remain bound with
  E_drift < 0.003%. Extremely robust.
- inv_He_q-2_R2.0_eta0.75 (runs 165--167): all 3 levels bound, E_drift < 0.031%.
- He_q2m10_R2.0_a2pi3_eta0.75 (runs 168--170): all 3 levels bound, E_drift < 0.031%.
- He_c4_R1.0_eta0.8 (runs 156--158): bound at delta=0.01 and 0.05, marginally unbound at
  delta=0.1 (max_r=23.6 just exceeds 10x threshold).

**Fragile** (break under perturbation):
- He_q2m10_R2.0_api_eta0.5 (runs 153--155): bound at delta=0.01 but fails at delta=0.05.
  The eta=0.5 configuration is near the stability boundary.
- He_q2m10_R2.0_api_eta0.25 (runs 177--179): fails at all perturbation levels. Marginal
  orbit that was already borderline (max_r=26.4 in the unperturbed case).

## Key findings

1. **Asymmetric mass is essential for 3-body Weber bound states.** The heavy-nucleus
   helium-like configuration (m_nucleus=10, m_electron=1) is by far the most successful.
   Equal-mass configurations are overwhelmingly unstable.

2. **Unlike charges remain necessary.** Consistent with Agent 05's 2-body finding (804 bound
   orbits, all unlike-charge), the 3-body search finds bound states only when unlike charges
   are present. Same-sign triple (+++) produced zero good-bound orbits.

3. **Higher c stabilizes.** c=4 runs succeed far more often than c=1. The Weber correction
   (velocity-dependent force) scales as 1/c^2, so larger c means the system is closer to
   pure Coulomb, which has well-known stable orbits. The most impressive results (runs
   132, 136) are at c=1, but they require large R=2.0 and heavy nucleus.

4. **Anti-podal electron orbits (angle_sep=pi) are most stable**, but non-symmetric
   angle_sep=2pi/3 also works at higher eta, just with larger orbital excursions.

5. **Planetary atom geometry works but with poor energy conservation** in the 2+/1- equal-mass
   case. The tiny nucleus separation (0.05) at c=1 is below the critical distance, causing
   immediate integrator failure. At c=4, orbits survive but drift 1--23%.

6. **Perturbation robustness**: The eta=0.75, R=2.0 helium-like orbits are genuinely stable
   attractors, surviving 10% perturbations. The eta=0.5 orbits are marginally stable.

## Comparison: symmetric vs non-symmetric initial conditions

| Property            | Symmetric (equilateral, equal mass)     | Non-symmetric (helium-like, unequal mass) |
|---------------------|-----------------------------------------|-------------------------------------------|
| Good-bound count    | 0 out of 40 equilateral runs            | 10 out of 28 asymmetric runs              |
| Best E_drift        | n/a (none bound)                        | 0.00003% (run 132)                        |
| Perturbation robust | n/a                                     | Yes, up to delta=0.1 for best orbits      |
| Mechanism           | --                                      | Heavy nucleus creates stable potential well|

The symmetric configurations suffer from a fundamental instability: all three particles have
comparable inertia and respond to perturbations on similar timescales, leading to resonant
energy transfer and escape. The asymmetric helium-like case suppresses this by anchoring the
heavy nucleus, turning the problem into an effective 2-body problem (two electrons orbiting
a nearly fixed center) with a perturbative electron-electron interaction.

## Files

- `ic_generators_3body.jl` -- initial condition generators (6 geometries)
- `three_body_survey.jl` -- survey driver script (179 runs)
- `survey_results.csv` -- full results table
- `NOTES.md` -- this file
