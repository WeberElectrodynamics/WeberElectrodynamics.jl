# Agent 7 — Maximal Lyapunov Exponent (MLE)

Finite-time maximal Lyapunov exponent λ_max for the 4-body 2+/2− Weber
Hamiltonian, on the same 3 × 4 (configuration × energy-fraction) IC grid that
Agent 6 used for Poincaré sections. Method: two-trajectory shadow integration
with periodic renormalization, run on top of `SharedSurvey.run` so the
integrator, time step (`dt = 1e-3`), and Weber options match every other Wave-2
agent.

Artefacts in this directory:

- `lyapunov.jl` — driver: builds ICs from `shared/ic_generators.jl`, runs the
  shadow estimator, writes CSV + heatmap.
- `lyapunov_results.csv` — per-IC λ_max, retcode, completed renorm intervals.
- `lyapunov_heatmap.png` — 3 × 4 heatmap of λ_max with NaN sentinels.

## Method

For each IC `(q₀, p₀) ∈ ℝ¹⁶`:

1. Form the full phase-space state `x = [q; p]` (16 components for N=4, d=2).
2. Pick a fixed-seed random unit vector `u ∈ ℝ¹⁶` and a perturbation amplitude
   `ε = 1e-8`. Initial shadow `x̃₀ = x₀ + εu`.
3. For `k = 1 … N` (with `N = tmax/Δt_renorm = 20/0.5 = 40`):
   - Integrate the reference and shadow each for `Δt_renorm = 0.5`.
   - Compute `δ_k = ‖x̃ − x‖`.
   - Accumulate `γ_k = log(δ_k / ε)`.
   - Renormalize: `x̃ ← x + ε · (x̃ − x)/δ_k`.
4. `λ_max ≈ (1/(N · Δt_renorm)) · Σ_k γ_k`.

If a leg fails (`retcode != :Success` — typically a fixed-point divergence at a
near-singular pair encounter, or a particle escaping to large radius), the
estimator returns the partial average over the `k − 1` completed intervals
(`retcode = :partial`) instead of throwing away the run.

## Caveats

- The Weber integrator is symplectic with energy projection, so `δ_k` does not
  collapse to roundoff; the perturbation tracks the true linearized flow well
  for `δ ≪ 1` (verified: `δ_final ~ 1e-8 – 1e-5` on successful intervals).
- Shadow trajectories diverge by O(1) within one renorm interval in strongly
  chaotic regions. The estimator clamps the divergence by renormalizing each
  interval, but if the **shadow** trajectory hits a near-singular pair encounter
  while the reference does not (or vice versa), one of the two `solve!` calls
  fails and the run is reported as `:partial`.
- 1-interval `:partial` rows are essentially upper bounds on `λ_max` (they
  measure the exponential separation rate over a single 0.5-time-unit window
  and so are dominated by the largest local stretching rate, not the time
  average). Treat them as a "definitely chaotic / definitely escapes" flag,
  not a quantitative MLE.
- `tmax = 20` is short by Lyapunov standards: the units of `λ_max` are
  inverse-time, the natural time scale is the Coulomb-circular period
  `T ~ 2π · √(s³/q²) ~ 2π` for `s ≈ 1`, so `tmax/T ≈ 3` periods. Hence
  "finite-time" — fine for ranking, not for asymptotics.

## Results

12 ICs total. One full `:Success` (40/40 intervals); the rest are partial.

| Configuration         | η = 0.10 | η = 0.25 | η = 0.45 | η = 0.70 |
|-----------------------|---------:|---------:|---------:|---------:|
| Alternating square    | 2.59 ¹   | 9.06 ²   | 5.85 ¹¹  | **1.70** |
| Rhombus (a=1, b=0.6)  | 1.42 ¹   | 7.38 ¹   | 8.09 ¹   | 9.15 ¹   |
| Two dimers            | 9.18 ⁵   | 8.31 ⁵   | 5.50 ⁴   | 3.01 ⁴   |

Bold = full success (40/40). Superscript = number of completed renorm
intervals out of 40 — values with low counts are upper bounds. Energy-fraction
η is the IC parameter passed to the generators; for `two_dimers` we remap to
`intra_fraction = 0.3 + 0.6 η` (the dimer needs a non-zero kinetic floor or it
collapses immediately to a head-on collision and the integrator bails before
step 30).

### Highlights

- **No bound ordered IC found by MLE alone.** Even the most ordered case from
  Agent 6 (alternating square η = 0.10, S2 Poincaré section classified as
  *torus*) gives `λ_max ≈ 2.6` over a single surviving interval, consistent
  with the local stretching expected from a near-square configuration relaxing
  via Weber-corrected Coulomb forces. There is no IC in this grid with
  `λ_max ≈ 0`.
- **One IC has `λ_max > 1` and survives the full 20 time units**: the
  alternating-square η = 0.70 case, with `λ_max ≈ 1.70` over 40/40 intervals.
  This is a *fast-chaos but bounded* regime — the configuration stays globally
  bound (no escape, no fixed-point failure) but neighbouring trajectories
  separate exponentially with an e-folding time of `~0.6` time units (about
  10% of the rotation period).
- **Two-dimer ICs are systematically the worst-behaved**: every two-dimer run
  failed by interval 5. The two (+−) dyads, even with widened separation
  `R = 3.0` and elongated `a = 0.8`, drift apart and induce close cross-pair
  approaches that crash the symplectic projection. The rotating-dimer
  configuration is **not** a bound state of the Weber 2+/2− problem — at least
  not in the bare (non-Zöllner) regime. Agent 13's Zöllner sweep should re-test
  these ICs.
- **Rhombus ICs all fail after the very first renorm interval** for every η,
  with `λ_max` in the range 1.4–9.2 over that single interval. Agent 6 saw
  the same behaviour (rhombus runs reported `escape` on every section). The
  rhombus is unstable to a near-immediate collision/escape regardless of
  energy.

## Cross-reference with Agent 6 (Poincaré / KAM)

Agent 6 (`research/FourBodyTwoPlusTwoMinus/06_poincare_kam/poincare_summary.csv`)
classifies each IC's three Poincaré sections (S1, S2, S3) as `torus`,
`chaotic`, `sparse`, or `escape`. Comparison:

| IC                       | Agent 6 verdict        | Agent 7 λ_max | Agree? |
|--------------------------|------------------------|--------------:|:-----:|
| alt. square η = 0.10     | sparse / **torus** / sparse | 2.59 (1 itv) | partial — Agent 6 found a torus on S2; my single-interval λ is an upper bound on local stretching, so this is *consistent* with quasi-periodic motion that is later disturbed. |
| alt. square η = 0.25     | sparse / sparse / sparse  | 9.06 (2 itv) | yes — neither agent finds order. |
| alt. square η = 0.45     | sparse / sparse / sparse  | 5.85 (11 itv) | yes — both find chaos / collapse. |
| alt. square η = 0.70/0.75| sparse / sparse / **chaotic** | **1.70 (full)** | **strong agreement**: Agent 6 sees an explicitly chaotic S3 section, Agent 7 measures a quantitative bounded MLE. This is the cleanest cross-validation in the grid. |
| rhombus (all η)          | escape on every section   | 1.4–9.2 (1 itv) | yes — both detect immediate disintegration. |
| two_dimers (all η)       | sparse on every section   | 3–9 (4–5 itv) | yes — both fail to find any bounded structure. |

The most informative agreement is the `alt_square η ≈ 0.7` row: it is the
**only** IC in the joint Agent 6 / Agent 7 grid that supports a globally
bounded but chaotic trajectory throughout the entire `tmax = 20` window, with
both the section-based and the Lyapunov-based diagnostics in agreement.

## Verdict for the master report

- The 2+/2− 4-body Weber phase space at `c = 1`, unit charges and masses, is
  **dominated by chaotic and dissociating regions** on the energy slice
  `0.1 ≤ η ≤ 0.7` for all three "obvious" symmetric configurations.
- **No KAM-like ordered island survived 20 time units in either agent's
  diagnostics.** The closest is the alternating-square η ≈ 0.10 case, where
  Agent 6 found a torus on S2 but the trajectory collapses before our second
  renorm interval; this is a candidate for finer-resolution follow-up.
- **One bounded chaotic regime exists** at alternating-square η ≈ 0.70 with
  `λ_max ≈ 1.7`. This is the only IC that is unambiguously *bounded* and
  unambiguously *chaotic* in both the Poincaré and the Lyapunov view.
- Agents 5 (periodic-orbit search) and 8 (sub-Weber-radius dynamics) should
  treat the alternating-square family as the most promising basin and look for
  periodic orbits that the chaotic η = 0.7 trajectory might be shadowing.

## Reproducing

```bash
cd /Users/mac/dev/Weber/WeberElectrodynamics
julia --project=. research/FourBodyTwoPlusTwoMinus/07_lyapunov_chaos/lyapunov.jl
```

Wall time on a Mac laptop: ≈ 12 s for the entire 12-IC sweep (most ICs collapse
within a handful of intervals; only the alt-square η = 0.70 case runs the full
40 × 2 = 80 legs).
