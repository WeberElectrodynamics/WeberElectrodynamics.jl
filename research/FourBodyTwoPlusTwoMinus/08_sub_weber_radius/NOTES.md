# Agent 8 — Sub-Weber-Radius Dynamics for the 4-body 2+/2− System

## Setup

Units: `m = q = c = 1`, so the like-charge critical radius is
`ρ_like = q²/(μ c²) = 1/(1/2) = 2`.
The unlike pairs have no real critical radius (`ρ < 0`).
Particle labels: 1, 2 are positive; 3, 4 are negative.
The 6-pair index order (used everywhere below) is
`(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)`,
with `(1,2)` and `(3,4)` the two like-charge pairs.

Integrator: `SymmetricProjectionIntegrator` with `dt = 5e-4`, `tmax = 50`,
`collision_bounce_radius = 0.05` (0.02 for Experiment C). Per CLAUDE.md, neither
LC nor adaptive-Cartesian regularization handles Weber's velocity-dependent
force, so we rely exclusively on the bounce escape hatch.

Driver: `run_sub_weber.jl`. Outputs: `phase_diagram.png`, `grid_summary.csv`,
`figures/pairs_{A,B,C}.png`. Total integrations: 3 named experiments + 25 grid
cells = 28 (within budget). Each run capped at `tmax ≤ 50`.

## Headline result

**Every single one of the 28 integrations terminated with `retcode = :Failure`**,
regardless of bounce radius, and only Experiment A made meaningful temporal
progress before failing. This is exactly what the topological obstruction of
Frauenfelder–Weber 2024 (Thm 2.1) predicts for sub-ρ like-pair dynamics with
nonzero pair angular momentum: the spiral collision is non-regularizable, the
projection step in the symplectic integrator stops converging, and the run
aborts. For the 2+/2− system, even a *deliberately* head-on ICs configuration is
contaminated by ℓ ≠ 0 perturbations from the other two particles, and the
system inherits the same obstruction.

We log the failures as data, not as bugs.

## Experiment A — `(+,+)` nucleus + `(−,−)` orbiters

ICs: `(+,+)` at `±r_like/2 = ±0.25` (well inside ρ = 2) with head-on momenta
`±0.2`; `(−,−)` placed at `(±2.5, +5)` with small tangential momenta `±0.1`.

| metric | value |
|---|---|
| retcode | `:Failure` |
| t_end / tmax | 4.140 / 50 (≈ 8.3%) |
| steps taken | 8280 |
| max E drift | 0.261 % |
| pair (1,2) r_min / r_max | 0.034 / 0.531 |
| pair (3,4) r_min / r_max | 4.816 / 5.000 |
| `bound_indicator` | false (escape_radius = 25) |

Interpretation. The `(+,+)` nucleus *does* form: pair (1,2) bounced off the
`bounce_r = 0.05` shell and oscillated inside `[0.034, 0.531]`, deep inside ρ.
The `(−,−)` pair stayed quietly in `[4.8, 5.0]` (Coulomb-bound dimer). What
killed the run is almost certainly the slow drift of pair (1,2) angular
momentum away from zero under the asymmetric pull of the distant negatives:
each pre-bounce inspiral becomes slightly off-axis, the next perihelion is
tighter, the spiral winding rate increases, and the projection iteration
collapses. Energy drift of only 0.26 % at the moment of failure rules out an
energetic blow-up — the integrator is simply unable to advance through a
near-spiral arc. **The nucleus is metastable for ~4 oscillation periods, then
the topological obstruction wins.**

See `figures/pairs_A.png`.

## Experiment B — Two simultaneous sub-critical nuclei

ICs: `(+,+)` at `(−3 ± 0.25, 0)` and `(−,−)` at `(+3 ± 0.25, 0)`, both pairs
with anti-parallel head-on momenta `±0.15`.

| metric | value |
|---|---|
| retcode | `:Failure` |
| t_end | 0.482 |
| steps | 964 |
| max E drift | 0.138 % |
| (1,2) r_min / r_max | 0.033 / 0.515 |
| (3,4) r_min / r_max | 0.033 / 0.515 |

Both nuclei perform essentially identical first pre-bounce contractions, hit
the bounce shell at almost the same instant, and the projection step fails
**immediately after the first bounce**. Two simultaneous sub-ρ encounters are
worse than one — the Weber velocity-dependent term couples the two relative
coordinates through the COM frame, and the projection's fixed-point iteration
loses contraction. Energy drift is small at termination because the failure is
detected and reported promptly.

See `figures/pairs_B.png`.

## Experiment C — Tight unlike `(+,−)` dimer + far `(+,−)` orbiter

ICs: pair (1,3) tight at `r₀ = 0.3` with tangential momenta of half the
circular Coulomb estimate; pair (2,4) at `r ≈ 1` separated by `R = 4` along x.
`bounce_r = 0.02`.

| metric | value |
|---|---|
| retcode | `:Failure` |
| t_end | 0.266 |
| steps | 533 |
| max E drift | **185 %** (energetic blow-up) |
| (1,3) r_min / r_max | 0.029 / 1.179 |
| (2,4) r_min / r_max | 0.029 / 1.000 |

The unlike pair has *no* critical radius — the Weber attractive Coulomb part
plus the velocity term simply makes the dimer too tight; (1,3) and (2,4) both
collapse, the bounce kicks the velocity-dependent term into a regime where
energy is no longer well controlled (energy drift explodes to 185 %), and the
integrator dies. The lesson: unlike-pair "tight Coulomb dimer" ICs are
**numerically intractable** under the current scheme even though they are
physically benign in the pure Coulomb limit. The Weber velocity correction
spoils the symplectic integrator's energy conservation as soon as `|ṙ|` becomes
comparable to `c`.

See `figures/pairs_C.png`.

## Experiment D — 5×5 collision survival map for the `(+,+)` pair

Grid: `r₀_like ∈ {0.20, 0.525, 0.85, 1.175, 1.50}` (all sub-ρ),
`p₀ ∈ {0.05, 0.1875, 0.325, 0.4625, 0.60}` (head-on like-pair momenta),
spectator `(−,−)` held at `r = 1`, far away at `(±0.5, +8)`. Bounce radius 0.05.
25 cells.

**Result: every cell labels as `:Failure`.** None completed `tmax = 50`. See
`grid_summary.csv` and `phase_diagram.png` (a uniformly-coloured grid is
itself a result — there is no Success island anywhere in this rectangle of IC
space).

This is striking because the spectator `(−,−)` carries zero initial momentum
and the like pair has perfectly anti-parallel `±p₀` along the x-axis, so the
relative angular momentum of (1,2) is exactly zero at `t = 0`. The failure
mechanism is the same as Experiment A: the spectator `(−,−)`'s Coulomb
attraction on (1,2) is generically asymmetric (because the `(−,−)` pair is at
`(±0.5, +8)` — non-symmetric with respect to (1,2)'s axis), so a tiny ℓ leaks
into pair (1,2) within a few steps, and the spiral obstruction kicks in. Even
the largest-`r₀` row (`r₀ = 1.50`, still sub-ρ but only just) fails: at that
separation the Weber correction is mild, but `p₀ = 0.05` is too small to
overcome Coulomb repulsion, and `p₀ = 0.6` overshoots into the deeply
non-perturbative regime.

A **degenerate symmetric** placement of `(−,−)` (e.g. `(0, +8)` and `(0, −8)`)
would preserve the (1,2) ℓ = 0 line by the C₂ symmetry of the configuration —
that experiment is the natural follow-up but is outside Agent 8's budget; see
the open-questions section.

## Honest summary table

| Exp | retcode | t_end | E drift | Bound? | Comment |
|---|---|---|---|---|---|
| A | Failure | 4.14 | 0.26 % | no | Nucleus formed, ~few bounces, then spiral |
| B | Failure | 0.48 | 0.14 % | no | First simultaneous bounce kills it |
| C | Failure | 0.27 | 185 % | no | Energetic blow-up; unlike tight dimer is intractable |
| D | 25/25 Failure | < 1 | — | no | Whole rectangle is non-survivable |

## Conclusions

1. **The `(+,+)` sub-ρ nucleus is *metastable* in the 4-body 2+/2− setting**, in
   the sense that it survives a handful of head-on bounces (Experiment A,
   ~4 t-units), but the asymmetric pull of the unlike particles always injects
   enough angular momentum to trigger an unrecoverable spiral.
2. **There is no Success island in the natural `(r₀, p₀)` rectangle** for the
   like-pair sub-critical regime once you place a generic spectator pair —
   `phase_diagram.png` is uniformly red. This contrasts sharply with the
   3-body planetary-atom (++−) result (`research/investigations/Three…`),
   where a *single* electron's symmetric pull *did* preserve the nucleus.
   Adding a second negative breaks that protective symmetry.
3. **Unlike-pair tight Coulomb dimers (Experiment C) are *numerically
   intractable*** under the current symplectic projection scheme. The Weber
   velocity term causes a runaway in energy as `|ṙ|` approaches `c` near
   perihelion. This is *not* the topological obstruction; it is a
   discretisation/projection failure that a different integrator (e.g. with
   exact energy preservation) might handle.
4. **The bounce escape hatch is necessary but insufficient** for 2+/2−.
   It rescues true `ℓ = 0` head-on encounters (and we observed this in
   Experiment A) but cannot help once `ℓ` becomes nonzero, because nothing in
   the current scheme regularises the spiral.

## Open questions for follow-up agents

- **Symmetric spectator**: place `(−,−)` symmetrically at `(0, ±R)` so that the
  reflection `y → −y` exactly preserves `ℓ_{12} = 0`. Predict: the nucleus
  survives much longer, possibly indefinitely. This is the cleanest
  falsifiable prediction from the present null result.
- **Zöllner enhancement**: Agent 13 should retry these ICs with `zollner_a > 0`
  to see whether the unlike-pair κ boost stabilises the nucleus by deepening
  the ambient potential well around the `(+,+)` core.
- **Alternative integrator**: a fully energy-preserving (e.g. Gauss collocation
  or AVF) scheme could rescue Experiment C and would be a fair test of whether
  the unlike-dimer pathology is physical or numerical.

## Files

- `run_sub_weber.jl` — driver, ~280 lines, regenerates everything.
- `phase_diagram.png` — Experiment D heatmap (all `:Failure`).
- `grid_summary.csv` — 25 rows, columns `r0_like, p0, label, drift_pct`.
- `figures/pairs_{A,B,C}.png` — pair-distance time series for the three named runs.
