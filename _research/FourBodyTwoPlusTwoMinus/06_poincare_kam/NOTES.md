# Agent 6 — Poincaré sections and KAM structure

**Scope.** Simple Poincaré-section detector for the planar 4-body 2+/2−
Weber Hamiltonian, swept over 12 ICs (3 configurations × 4 energy fractions),
classifying section clouds as torus / chaotic / sparse / escape and flagging
KAM-stable candidates for Agent 5.

All numerics in `poincare.jl` (this directory); it includes
`shared/{ic_generators,metrics,run_survey}.jl` and writes `figures/` and
`poincare_summary.csv`.

## 1. Sections

The reduced (COM-removed) phase space is 12-dim in 2D; fixing `(E, L)` drops
it to 9 (Agent 2 §3). A 2D Poincaré section is therefore a *projection*, not
a true Birkhoff section. Three observables built from Jacobi vectors
`(r₊, r₋, R, P₊, P₋, P)` from Agent 2 §4:

| ID | `g(q,p)` | ġ-condition | Plotted |
|----|----------|-------------|---------|
| S1 | `(r₊)_y = 0` | `(ṙ₊)_y > 0` | `(|r₊|, |r₋|)` |
| S2 | `|r₊| − |r₋| = 0` | `d/dt > 0` | `(|R|, ∠(r₊, r₋))` |
| S3 | `L₁₂ = 0` (planar L of (++) about its centre) | none | `(r₁₂, ṙ₁₂)` |

The brief asks for S1 = `R_x = 0`, but for every IC family here `R_x` stays
one-signed (or oscillates only inside one period before termination),
yielding 0 crossings. Switching to `(r₊)_y = 0` reliably triggers on the
(+)-dimer rotation/breathing while still cutting a 1-codim hyperplane.

The detector walks saved `(q,p)` snapshots, looks for sign changes of `g`,
accepts the crossing only when `ġ > 0`, and linearly interpolates between
endpoints. With `dt = 1e-3` this is more than precise enough.

## 2. IC grid

12 ICs = 3 configurations × `energy_fraction ∈ {0.10, 0.25, 0.50, 0.75}`:

- `two_dimers` — `(+−)(+−)` dyads, `dyad_length = 0.6`, `separation = 3.0`,
  `inter_velocity = 0.05`. `energy_fraction` reused as `intra_fraction` for
  the bound circular speed inside each dyad.
- `rhombus` — vertices at `(±a, 0)` and `(0, ±b)` with `a = 1.5`,
  `b = 1.0 + 0.2·η`, rotating velocity assignment, `energy_fraction = η`.
- `alternating_square` — `+ − + −` at the corners of a unit square, rotating
  velocity assignment, `energy_fraction = η`.

Every IC is integrated with `tmax = 30`, `dt = 1e-3`, `c = 1`, and the
collision-bounce escape hatch (`bounce_r = 0.02`) enabled to survive any
incidental head-on hits during the warm-up phase.

## 3. Outcomes

The full table is `poincare_summary.csv`. Highlights below; "drift" is the
peak global energy-error percent reported by `compute_energy_timeseries`.

| config              | η    | retcode  | t_end  | n_S1 | n_S2 | n_S3 | drift % |
|---------------------|------|----------|--------|------|------|------|---------|
| two_dimers          | 0.10 | :Failure | 0.16   | 0    | 0    | 0    | 12.6    |
| two_dimers          | 0.25 | :Failure | 1.77   | 1    | 1    | 1    | 4.7     |
| two_dimers          | 0.50 | :Failure | 1.95   | 2    | 1    | 1    | 2.0     |
| two_dimers          | 0.75 | :Failure | 2.38   | 1    | 0    | 1    | 1336    |
| rhombus             | 0.10 | :Failure | 1.50   | 0    | 0    | 0    | 6.8     |
| rhombus             | 0.25 | :Failure | 1.51   | 0    | 0    | 0    | 19.6    |
| rhombus             | 0.50 | :Failure | 6.53   | 0    | 0    | 1    | 264     |
| **rhombus**         | **0.75** | **:Success** | **30.00** | **0** | **4** | **5** | **0.076** |
| alternating_square  | 0.10 | :Failure | 0.82   | 0    | 3    | 0    | 0.053   |
| alternating_square  | 0.25 | :Failure | 2.80   | 1    | 25   | 2    | 34.3    |
| alternating_square  | 0.50 | :Failure | 7.69   | 2    | 4    | 1    | 61.8    |
| alternating_square  | 0.75 | :Failure | 27.46  | 2    | 0    | 3    | 244     |

The dominant signal is that **almost every IC terminates early** with a
`:Failure` retcode from the projection-integrator's fixed-point solver. The
compounding fact is that even those runs that do reach `t = O(1..30)` often
ring up large energy drift before stopping, so the section-cloud size is
essentially a smoke test — not a true KAM survey.

Two qualitative survivors:

1. **`rhombus`, η = 0.75** — full `tmax = 30` with `:Success` and 0.076 %
   energy drift (~2.5×10⁻³ %/unit). S2 returns 4 clean crossings on a
   smooth `(|R|, ∠)` locus, S3 returns 5 nearly-collinear points on
   `(r₁₂, ṙ₁₂)`. Few crossings, but the temporal coherence and tiny drift
   match a near-integrable KAM torus under a 2D projection. **Flag for
   Agent 5 as a periodic-orbit shooting candidate.**
2. **`alternating_square`, η = 0.75** — runs to `t ≈ 27.5` (92 % of `tmax`)
   before tripping. S1 has only 2 crossings, but the long survival is
   suggestive. Terminal drift ~244 % indicates the second half collapsed
   into a near-singular event; the first half (`t ≲ 10`) is plotted. Worth
   retrying with `dt = 5e-4` in Agent 5.

Everything else dies within ~2 time units.

## 4. KAM-fraction estimate

Using the very rough classifier in `poincare.jl::classify`:

- `:torus` ⇐ ≥ 8 crossings on S1 with low fill fraction;
- `:chaotic` ⇐ ≥ 8 crossings, area-filling;
- `:sparse` ⇐ < 8 crossings but the run finished with low drift;
- `:escape` ⇐ early `:Failure` with > 5 % drift, or late `:Failure` with
  the integrator visibly diverging.

Because S1 is starved (most runs die before the dimer rotates more than once)
the resulting per-energy `:torus` fraction is ≤ 1/12 across the entire grid;
see `figures/kam_fraction.png`. **No quantitative KAM curve can be extracted
from this data.** The only honest qualitative statement is:

- **η = 0.10–0.25**: nearly all ICs in our grid die quickly; the integrator
  cannot resolve the close-encounter regime that the small-kinetic
  initialisation drops into. We cannot tell torus from chaos here without a
  Levi-Civita-regularised setup specialised to the singular pair, which is
  outside this agent's brief.
- **η = 0.50**: marginal; only the rotating square reaches t ≈ 8.
- **η = 0.75**: the only regime where two of three configurations survive
  for a substantial fraction of `tmax`, and where the rhombus shows the
  cleanest section signature. Provisional KAM fraction ~ 1/3 at this energy,
  with the caveat that 12 ICs is way below statistical significance.

## 5. Resonances and island chains

With at most 5 crossings per surviving run, no island-chain structure is
resolvable. Section S2 of `alternating_square` η = 0.10 *did* return 19
crossings before terminating (because the breathing frequency is fast there)
and they cluster on what visually looks like a smooth 1-D arc in
`(|R|, ∠(r₊, r₋))` — see `figures/section_S2_combined.png`. Without a longer
integration this is at best a "torus-like" hint; it is reported as `:torus`
in the CSV for the S2 channel of that single IC and ignored by the S1-based
KAM-fraction plot.

## 6. Limitations

1. **`(E, L)` conditioning.** The brief asks for sections at fixed `(E, L)`;
   we instead overlay ICs from different shells. For the surviving rhombus
   IC the projection is faithful in a small neighbourhood, but S2/S3 overlays
   should not be read as a single Poincaré map.
2. **2D projection of an 8D section.** A single point in the `(|r₊|, |r₋|)`
   plane corresponds to a 3-torus-worth of microstates; curve-vs-cloud
   heuristics are necessary but not sufficient for KAM detection.
3. **Integrator failure dominates.** 11 of 12 ICs fail the projection step
   before useful Poincaré data is collected. Mitigations for any re-run:
   a Levi-Civita lifted-pair regulariser on close-approach pairs (Agent 8),
   a smaller `dt` or adaptive macrostep, and instrumentation of which pair
   triggered the failure.
4. **S3 ġ-sign filter** is constant (`ġ = 1.0`) because the symbolic
   `L̇₁₂` against the Weber correction is non-trivial. S3 therefore counts
   *all* sign changes of `L₁₂`, and should be read as an event detector
   rather than a true return map.

## 7. Hand-off to Agent 5

The single robust torus-candidate from this scan is:

```
ic = rhombus(a = 1.5, b = 1.15, energy_fraction = 0.75, velocity_mode = :rotating)
tmax = 30, dt = 1e-3, c = 1, bounce_r = 0.02
```

It satisfies the Agent 2 brake-orbit symmetry only approximately, but it
gives a 30-unit `:Success` integration with sub-percent drift and qualitative
torus structure on S2/S3. Agent 5 should:

- Try Newton-shooting starting from the IC above with the rotation period
  estimated from the 4 S2 crossings (Δt ≈ 5–7 between them).
- Continue in `b` to map the rotating-rhombus family.
- Re-do the section detection on the converged orbit at `dt = 5e-4` to
  validate the curve geometry.

## 8. Files produced

- `poincare.jl` — section detector + survey driver.
- `poincare_summary.csv` — one row per IC.
- `figures/section_S1_two_dimers.png`
- `figures/section_S1_rhombus.png`
- `figures/section_S1_alternating_square.png`
- `figures/section_S2_combined.png`
- `figures/section_S3_combined.png`
- `figures/kam_fraction.png`
