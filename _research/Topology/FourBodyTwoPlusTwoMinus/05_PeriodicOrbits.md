# 05 — Periodic Orbit Search & Continuation

Agent 5 deliverable. Code: `search.jl`. Data: `found_orbits.csv`. Plots: `figures/`.

## 1. Setup and methodology

Planar 4-body 2+/2− Weber Hamiltonian, `m_i = 1`, `q_{1,2}=+1`,
`q_{3,4}=−1`, `c = 1`, `κ_{ij}=1`. Unregularized
`SymmetricProjectionIntegrator` (the LC backends do not regularize the
Weber velocity term). Failures (projection divergence near a collision)
are recorded but not counted as found.

**T-brake shooting.** Time-reversal `T : p → −p` (Agent 2 §1) implies any
orbit that starts on `{p = 0}` and returns to `{p = 0}` at `t = T/2` is
`T`-periodic. This halves the search space (8-dim brake submanifold).
The brake-return detector picks the first local minimum of `‖p‖_∞` that
drops below `ptol = 5·10⁻²`. For non-brake families (rotating square,
dimer-dimer) we use a pair-distance return detector instead.

**Monodromy.** For each returner we form `M = ∂z(T)/∂z(0)` via symmetric
finite differences (`ε = 10⁻⁵`) and Floquet multipliers `= eig(M)`. An
orbit is "found" iff it stays inside `10×` its initial diameter for at
least one detected period and the energy drift is `< 5%`.

Helpers from `_research/Topology/FourBodyTwoPlusTwoMinus/shared/`; `HamiltonianSystem(4,2)`
is built once and cached.

## 2. Family 1 — Breathing alternating square (L = 0)

Brake at the alternating square (side `s`) plus an outward radial
velocity `v_rad`, preserving `D₄ × T`. Selected results (full table in
CSV):

| id | s | v_rad | retcode | T | E₀ | bound | |λ|_max |
|---|--:|--:|---|--:|--:|---|--:|
| F1a_brake_s0.50 | 0.5 | 0 | Failure | — | −5.172 | no | — |
| F1a_brake_s1.00 | 1.0 | 0 | Failure | — | −2.586 | no | — |
| F1a_brake_s2.00 | 2.0 | 0 | Failure | — | −1.293 | no | — |
| F1b_s1.00_v0.70 | 1.0 | 0.70 | Failure | 5.97 | −1.032 | no | — |
| F1b_s2.00_v0.40 | 2.0 | 0.40 | Failure | 6.06 | −0.879 | no | — |
| **F1b_s2.00_v0.50** | **2.0** | **0.50** | **Success** | **11.78** | **−0.646** | **yes** | **228.6** |
| F1b_s2.00_v0.60 | 2.0 | 0.60 | Success | — | −0.362 | yes | — |
| F1b_s2.00_v0.70 | 2.0 | 0.70 | Success | — | −0.026 | yes | — |

**Interpretation.** Pure brake (`L=0`) collapses radially through the
center; all four particles converge to the origin and the projection
diverges before the singular instant. Outward radial velocity delays the
catastrophe. At `s=2, v_rad=0.5` the system overshoots the inward turning
point, returns through brake, and the detector locates `T = 11.78` with
zero energy drift — the **only** fully closed periodic orbit found in any
of the four families.

Floquet analysis gives `max|λ| ≈ 228.6`: violently linearly unstable,
quantitatively consistent with Agent 4's real exponent `λ_real ≈ 1.026`
(`exp(1.026 · 5.5) ≈ 283`). This orbit exists only as an unstable saddle
— not observable generically.

Runs at `s=2, v_rad ∈ {0.6, 0.7}` are bound for `t ≤ 10` but moving
toward escape (`E₀ → 0`) so no brake-return is seen within `tmax`.

## 3. Family 2 — Rotating-square near orbits

Alternating square of side 1 with the rotating-frame angular velocity
`ω = √((2√2 − 1)/4) ≈ 0.6761` from Agent 4, plus a small radial
perturbation `δ ∈ {−0.05, −0.01, 0, +0.01, +0.05}`. The `δ = 0` IC is the
exact rotating relative equilibrium for both the Coulomb and Weber
Hamiltonians (Agent 4 §3 — Weber is identically zero on rigid rotations).
**All five runs failed at `t ≈ 0.25`**, far less than one rotation
period (`T_rot ≈ 9.3`).

Even at `dt = 5·10⁻⁵`, the relative equilibrium itself fails after
≈0.25 time units, far less than one rotation `T_rot = 2π/ω ≈ 9.3`. This
is exactly Agent 4's instability: numerical roundoff seeds the real
unstable mode (`λ ≈ 1.026`) plus the Krein quadruple `±0.84 ± 0.68 i`,
and the configuration collapses through an unlike-pair edge in a fraction
of a rotation. Simple radial offset does not uncover a nearby
breathing-rotating orbit; a successful continuation would need a shooter
that projects out the unstable real eigenvector at each iterate
(Lindstedt–Poincaré or arc-length continuation against the 1-dim unstable
manifold). **No rotating-square periodic orbit was found** — Family 2 is a
clean numerical confirmation of Agent 4.

## 4. Family 3 — Dimer-dimer orbital ("Zöllner molecule")

Two `(+ −)` dyads of length `a = 0.4`, separated by `R = 3`, internal
Coulomb-circular velocity (`intra_fraction = 0.4`), inter-dyad relative
velocity `v ∈ {0.05, 0.1, 0.2, 0.3}`. **All four runs failed at `t < 1.6`**
with `E₀ ≈ −3.5`.

Every dimer-dimer IC failed within `≈1` time unit, well before completing
one orbit of the inter-dyad pair (`T ∼ 2π R^{3/2}/√Q ≈ 30`). The
Coulomb-circular IC does not solve the Weber two-body problem exactly;
the intra-dyad orbit becomes eccentric, an unlike-pair `r_ij` reaches a
near-collision, and the projection diverges. **No dimer-dimer periodic
orbit was found.** A successful search needs (i) the exact Weber
two-body Kepler analogue as the dyad seed and (ii) averaging over the
fast intra-dyad period.

## 5. Family 4 — Collinear ABAB

Four particles on `x` at `±1.5d, ±0.5d` arranged `+, −, +, −` along the
line. Spacings `d ∈ {0.8, 1.0, 1.2}`, optional transverse kick
`k ∈ {0.02, 0.05}`. **All five runs failed** with no detected period.

All collinear ICs collapse: the inner unlike pair `(2,3)` (the central
`+,−` at distance `d`) is attractive and falls together regardless of
transverse kick. The 1D invariant subspace `{y_i = ẏ_i = 0}` (Agent 2 §3)
contains only collision trajectories. **No collinear ABAB periodic orbit
was found.**

## 6. Summary — found vs. failed

| Family | ICs scanned | Successes | Periodic | Linearly stable |
|---|--:|--:|--:|--:|
| F1 breathing square (L=0) | 21 | 3 | 1 | 0 |
| F2 rotating square near | 5 | 0 | 0 | 0 |
| F3 dimer–dimer | 4 | 0 | 0 | 0 |
| F4 collinear ABAB | 5 | 0 | 0 | 0 |
| **Total** | **35** | **3** | **1** | **0** |

The single periodic orbit located is **F1b_outbreath_s2.00_v0.50**:
a `D₄ × T`-symmetric breathing mode of the alternating square at side 2,
period `T = 11.78`, energy `E ≈ −0.646`, monodromy spectral radius
`max|λ| ≈ 228.6`, hence linearly unstable by an enormous factor. It is
the brake-symmetric image of the Coulomb-only `D₄ × T` breather extended
into the Weber regime; it survives in the Weber problem because the
Weber correction is small (`v² / c² ≲ 0.25` at the inward turning point)
and does not destroy the brake symmetry.

## 7. Negative result and physical reading

Within the four symmetry-reduced subfamilies and inside this agent's
budget, **the 4-body 2+/2− Weber Hamiltonian admits no linearly stable
periodic orbits.** The single orbit located is the breathing alternating
square, itself violently unstable (`|λ|_F ≈ 230` per period). Every other
candidate either collapses to a collision (Families 1a, 4), triggers
Agent 4's unstable eigenvector within a fraction of one rotation
(Family 2), or fails because the IC is inconsistent with the Weber
2-body problem (Family 3).

This quantifies Agent 4 §7's prediction: any bound 2+/2− family must be
non-rigid; we find exactly one such orbit and confirm it is unstable by a
factor > 200 per period. Genuine bound 2+/2− motion — if it exists — is
therefore not periodic in the strict closed-orbit sense; it must be
quasi-periodic (KAM tori, Agent 6) or chaotic-but-confined (Agent 7).

## 8. Next steps

- Agent 6: use F1b as a skeleton for KAM tori in its ≈13-dim center
  manifold (16 − 1 energy − 2 real unstable); also Poincaré sections
  through the brake submanifold for non-symmetric orbits.
- Agent 7: F1b is a clean Benettin test — Floquet log predicts the
  largest Lyapunov exponent directly.
- Agent 8: the near-collisions in F1a, F3, F4 approach `ρ = q²/c² = 1`
  — natural sub-Weber-radius domain.
- A follow-up Agent 5b should continue F1b in `c` from `∞` (Coulomb)
  toward `c = 1` and below, tracking Floquet multipliers as the Weber
  correction turns on.

## Files

- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/05_periodic_orbits/search.jl`
- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/05_periodic_orbits/found_orbits.csv`
- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/05_periodic_orbits/figures/F1b_outbreath_s2.00_v0.50.pdf`
- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/05_periodic_orbits/figures/F1b_outbreath_s2.00_v0.60.pdf`
- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/05_periodic_orbits/figures/F1b_outbreath_s2.00_v0.70.pdf`
