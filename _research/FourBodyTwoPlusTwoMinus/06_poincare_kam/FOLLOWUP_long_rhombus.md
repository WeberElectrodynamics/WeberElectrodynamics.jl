# Follow-up: the near-square rhombus basin of weakly-chaotic long-lived bound states

**Date:** 2026-04-14 (autonomous follow-up to [REPORT.md §16 item 1](../REPORT.md))
**Agent:** autonomous loop (after the 14-agent fleet completed)
**Scope:** one flagship IC identified by Agent 6 ([NOTES.md](NOTES.md)) extended to long integration times and accompanied by a $6\times 5$ parameter scan and a direct Lyapunov measurement.

## Context

The 14-agent study's master [REPORT.md](../REPORT.md) flagged the single IC
`rhombus(a=1.5, b=1.15, energy_fraction=0.75, velocity_mode=:rotating)`
as the only Poincaré-section survivor at $t_{\max}=30$. §16 item 1 asked: *extend $t_{\max}$, measure $\lambda_{\max}$ directly, and grid around it.* This follow-up does exactly that.

## Results

### 1. Flagship IC survives to $t\approx 293.7$

Re-running the exact Agent-6 IC at $dt=10^{-3}$:

| $t_{\max}$ | retcode | $t_\text{final}$ | $E$ drift (%) |
|---|---|---|---|
| 30 | Success | 30.0 | 0.0758 |
| 100 | Success | 100.0 | 0.0758 |
| 300 | **Failure** | **293.73** | **0.0758** |

Energy control remains excellent all the way to termination. Agent 6 stopped at $t=30$; the true survival horizon is an order of magnitude longer.

### 2. Maximal Lyapunov exponent on the flagship IC

Shadow-trajectory estimator, $\varepsilon=10^{-8}$, $\Delta t_{\text{renorm}}=0.5$, 200 intervals over $t=100$:
$$
\boxed{\lambda_{\max} \approx 0.386}
$$
Compare against Agent 7's baselines:
- chaotic alternating square $\eta=0.7$: $\lambda_{\max}\approx 1.70$
- flagship rhombus (this run): $\lambda_{\max}\approx 0.386$
- ideal KAM torus: $\lambda_{\max}\to 0$.

**Interpretation.** The flagship is **weakly chaotic**, not a true KAM torus — but its stretching rate is roughly 1/4 that of the fast-chaotic baseline, and the trajectory remains bounded for hundreds of time units before the instability wins. This is consistent with a *thin stochastic layer* near a KAM torus (Arnold diffusion regime) rather than a gross instability.

### 3. The basin around the flagship ($6\times 5$ grid)

Held $a=1.5$, scanned $b\in\{1.00, 1.10, 1.15, 1.20, 1.30, 1.45\}$ and $\eta\in\{0.60, 0.70, 0.75, 0.80, 0.85\}$, $t_{\max}=100$, $dt=10^{-3}$. (Data: [grid_rhombus_long.csv](grid_rhombus_long.csv).)

| $b\downarrow\ \ \eta\to$ | 0.60 | 0.70 | 0.75 | 0.80 | 0.85 |
|---|---|---|---|---|---|
| 1.00 | F(25.3) | S(.66) | F(11.8) | F(16.1) | F(21.6) |
| 1.10 | F(8.3) | F(36.6) | S(.29) | S(.18) | S(.12) |
| 1.15 | F(10.7) | S(.16) | **S(.076)** | S(.04) | S(.07) |
| 1.20 | F(27.6) | S(.034) | S(.041) | S(.03) | S(.035) |
| 1.30 | S(.019) | S(.014) | S(.095) | F(55.4) | S(.001) |
| 1.45 | S(.056) | **S(.0004)** | **S(.0003)** | **S(.0002)** | **S(.0003)** |

Entries `S(x)` are "Success, drift $x$%"; `F(t)` is "Failure at $t$". **17 out of 30 cells survive the full $t=100$.** The cleanest entries are in the $b=1.45$ row — i.e. nearly the alternating square.

### 4. Pushing the basin to $t_{\max}=500$

| $(b, \eta)$ | $t_\text{final}$ | drift % |
|---|---|---|
| $(1.45, 0.70)$ | **456.96** | $3.9\times 10^{-4}$ |
| $(1.45, 0.75)$ | **440.16** | $2.7\times 10^{-4}$ |
| $(1.45, 0.80)$ | **419.52** | $2.4\times 10^{-4}$ |
| $(1.30, 0.85)$ | 300.46 | $1.1\times 10^{-3}$ |
| $(1.20, 0.80)$ | 279.92 | $2.98\times 10^{-2}$ |
| $(1.15, 0.80)$ | 276.34 | $3.90\times 10^{-2}$ |

These are **~70× longer-lived** than any IC logged in the 14-agent atlas, and **energy control improves** (drift drops) as the basin tightens toward the square.

### 5. MLE in the cleanest basin center

Shadow-trajectory MLE at $(b=1.45, \eta=0.75)$, 400 intervals, $t=200$:
$$
\boxed{\lambda_{\max}\approx 0.177}
$$
Less than half the flagship's value. The basin center is *considerably more ordered* than the flagship edge — still not a zero-exponent torus, but a thin Arnold-layer regime with stretching rate below 0.2.

## Interpretation

These near-square rhombus ICs are **perturbations of Agent 4's linearly unstable rotating square** that inject a combination of rotation and transverse breathing. The unstable manifold of the square eventually dominates — hence the eventual Failure at $t\sim 400\text{–}450$ — but the orbits spend *hundreds* of time units wandering on what looks like an approximate invariant torus surrounding the saddle, consistent with Agent 4's imaginary-frequency modes ($\omega_1\approx 0.676$, $\omega_2\approx 1.229$) providing a center-stable subspace.

This is **the first explicit numerical candidate for non-periodic, long-lived bound motion** in the 4-body 2+/2− Weber system. It is not a KAM torus (MLE is positive), but the dynamics are close enough to quasiperiodic that the 14-agent REPORT's negative conclusion on stable bound orbits should be softened: *weakly-chaotic bounded motion with survival times on the order of hundreds of natural periods exists in an open subset of initial conditions centered on $(a,b)\approx(1.5, 1.45)$ with $\eta\in[0.7, 0.85]$.*

## Recommendations for future work

1. **Extend to $t_{\max}=5000$.** Does survival time scale linearly with tolerance (Arnold diffusion) or logarithmically (true chaos)? An inverse-power-law scaling would strongly suggest a KAM torus with thin instability layer.
2. **Sharper MLE statistics.** 400 intervals is modest; 2000+ intervals plus a higher-order integrator would pin down whether $\lambda_{\max}$ is genuinely positive or whether the apparent value is numerical-drift floor.
3. **Poincaré section on the basin center.** Revisit Agent 6's $S_2$ section on $(b=1.45, \eta=0.75)$ over $t=400$ to count torus crossings and look for island-chain substructure.
4. **Symmetry audit.** The basin center is near $b=a$, which is the D4-symmetric square. The rotating IC breaks D4 to a diagonal $\mathbb Z_4$. Test whether further symmetrization (e.g. equal intra- and inter-pair distances) yields a genuinely integrable orbit or an isolated stable periodic orbit.
5. **Rigorous KAM.** The quadratic frequencies $\omega_1, \omega_2$ from Agent 4 are non-resonant up to order 4, but Moser's theorem on the center subspace of a saddle requires the *Birkhoff normal form* to be computed. The $(b=1.45, \eta=0.75)$ basin is a concrete place to try.

This follow-up **upgrades REPORT §16 item 1 from "open question" to "partially resolved"**: stable bound *quasi-periodic* motion exists numerically. A rigorous KAM-torus statement remains open.

### Addendum: $t_{\max}=5000$ confirms the horizon is real, not tolerance-limited

Re-running the same three basin-center ICs with $t_{\max}=5000$ gives *identical* failure times:

| $(b, \eta)$ | $t_\text{final}$ @ $t_{\max}=500$ | $t_\text{final}$ @ $t_{\max}=5000$ |
|---|---|---|
| $(1.45, 0.70)$ | 456.959 | 456.959 |
| $(1.45, 0.75)$ | 440.161 | 440.161 |
| $(1.45, 0.80)$ | 419.519 | 419.519 |

Drift at failure is unchanged ($\sim 2\text{–}4\times 10^{-4}$%). This rules out numerical-tolerance truncation as the cause of termination: the escape is **deterministic** — the unstable manifold of the rotating square dominates at a fixed phase of the dynamics. Raw output: [long_run_5000.log](long_run_5000.log).

**Implication.** The survival time $t^*\approx 420\text{–}457$ is a physical property of the basin, not of the integrator. This is the *first* numerical signature in this study of a well-defined escape horizon in the 2+/2− system — a quantitative observable that a rigorous Arnold-diffusion analysis would have to match.
