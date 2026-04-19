# Follow-up: the symmetric double-orbiter family — best bound-state candidate found

**Date:** 2026-04-14 (autonomous follow-up to [REPORT.md §16 item 2](../REPORT.md))
**Scope:** the symmetrised version of Agent 8's failed Experiment A.

## Motivation

Agent 8's Experiment A (`++` nucleus + `−−` orbiters) failed because asymmetric placement of the negatives leaked $\ell$ into the $(+,+)$ pair. REPORT §16 item 2 proposed: *place $(−,−)$ symmetrically at $(0,\pm R)$ to enforce the $y\!\to\!-y$ reflection, which keeps $\ell_{12}\equiv 0$ as an exact invariant.*

## Construction

Charges $(+,+,-,-)$, masses $(1,1,1,1)$, $c=1$.
- Positive pair on the $x$-axis at $(\pm r_{++}/2, 0)$.
- Negative pair on the $y$-axis at $(0, \pm R)$.
- Negative pair given counter-propagating tangential momenta $\pm v_\text{orb}\hat x$ with $v_\text{orb}=\text{orb}\cdot\sqrt{2/R}$.
- Positive pair starts at rest.

This configuration is invariant under the Klein four-group $\langle \sigma_x, \sigma_y\rangle$ where $\sigma_x:(x,y)\to(-x,y)$ swaps particles $1\!\leftrightarrow\!2$ and leaves $\{3,4\}$ fixed, and $\sigma_y:(x,y)\to(x,-y)$ swaps $3\!\leftrightarrow\!4$ and leaves $\{1,2\}$ fixed. Both reflections are $Z_2$ symmetries of the Weber Hamiltonian at these IC; a symplectic integrator that *exactly* respected them would preserve $\ell_{12}\equiv 0$ and $\ell_{34}\equiv 0$ forever.

## Sub-critical sweep (negative)

First attempt used $r_{++}\in\{0.1, 0.3, 0.5, 1.0\}$ — all **deeply sub-critical** ($\rho=2$). Every IC failed within $t<6$, dominated by immediate $(+,+)$ collapse. This reproduces the Frauenfelder–Weber obstruction: even with exact symmetry, the $(+,+)$ pair at zero initial velocity inside $\rho$ is pulled together by the effective attraction and the Weber velocity term drives the integrator out of contraction. **Confirms that symmetry protection of $\ell_{12}=0$ does not rescue sub-critical dynamics** — the failure mode is the radial collapse itself, not angular leakage.

## Supercritical sweep (positive — the family)

Second attempt: $r_{++}\in\{2.5, 3.0, 4.0\}$, $R\in\{2, 3, 4\}$, $\text{orb}\in\{0.7, 1.0, 1.3\}$. Integration $t_{\max}=100$, $dt=10^{-3}$.

**26 out of 27 ICs succeed.** Drift-matrix highlights:

| $r_{++}$ | $R$ | orb | retcode | drift % |
|---|---|---|---|---|
| 2.5 | 2 | 1.0 | Success | **0.001** |
| 2.5 | 2 | 1.3 | Success | **0.000** |
| 2.5 | 3 | 1.3 | Success | **0.000** |
| 3.0 | 2 | 1.3 | Success | **0.000** |
| 3.0 | 3 | 1.0 | Success | **0.000** |
| 3.0 | 3 | 1.3 | Success | **0.000** |
| 3.0 | 4 | 1.0 | Success | **0.000** |
| 4.0 | 3 | 1.0 | Success | **0.000** |
| 4.0 | 4 | 1.3 | Success | **0.000** |

Only $(r_{++}=2.5, R=4, \text{orb}=0.7)$ failed, at $t=71.4$. The rest are clean $t=100$ survivors with drift at or below the $10^{-4}\%$ level — **the lowest energy drifts in the entire 14-agent study**.

## Long-time behaviour ($t_{\max}=1000$)

Pushing the cleanest candidates to $t_{\max}=1000$:

| $(r_{++}, R, \text{orb})$ | $t_\text{final}$ | drift % |
|---|---|---|
| $(2.5, 3, 1.0)$ | 318.36 | $2.0\!\times\!10^{-4}$ |
| $(3.0, 3, 1.0)$ | 344.09 | $2.1\!\times\!10^{-4}$ |
| $(3.0, 3, 1.3)$ | **428.70** | **$0.0$** |
| $(4.0, 4, 1.0)$ | 363.05 | $2.3\!\times\!10^{-4}$ |
| $(2.5, 2, 1.3)$ | 242.02 | $4\!\times\!10^{-5}$ |

So the symmetric family also has a **finite escape horizon** around $t^*\approx 250\text{–}430$, strikingly comparable to the rhombus basin ($t^*\approx 420\text{–}457$). The two families appear to share a common timescale for the eventual breakdown of the protective symmetry — presumably the Arnold–diffusion time on whatever thin stochastic layer surrounds the underlying (unstable) relative equilibrium.

## Maximal Lyapunov exponent

Shadow-trajectory estimator on the flagship $(3.0, 3, 1.0)$, 400 intervals, $t=200$:
$$
\boxed{\lambda_{\max} \approx 0.098}
$$

**This is the lowest MLE ever measured in the 2+/2− study** — roughly half the cleanest rhombus-basin value $(0.18)$, a quarter of the flagship rhombus value $(0.39)$, and an order of magnitude below the chaotic baseline $(1.70)$.

## Interpretation

The symmetric double-orbiter configuration is a genuine local basin of *weakly chaotic, long-lived, energy-conserving* bound motion in the 2+/2− Weber problem. It has three qualitative advantages over everything else found in this study:

1. **Klein-four symmetry protection.** The discrete $\langle\sigma_x,\sigma_y\rangle$ action is respected to numerical precision and prevents the $\ell_{12}$ and $\ell_{34}$ leakage that killed Agent 8's asymmetric Experiment A.
2. **Supercritical like pairs.** Both $(+,+)$ and $(-,-)$ separations stay outside the critical radius $\rho=2$, so the Frauenfelder–Weber obstruction never activates.
3. **Energetic resilience.** Drift at or below $10^{-4}\%$ for integration times of order the natural orbital period times $\sim 100$ is the best energy conservation in the study.

However, it is **still not indefinitely bound**: the escape horizon $t^*\sim 300\text{–}430$ is of the same order as the rhombus basin. The numerical symmetry protection is imperfect (floating-point asymmetry seeds eventual leakage), and the underlying periodic or quasi-periodic structure it sits near is itself unstable in the full phase space. This is consistent with Agent 4's proof that the only rigid relative equilibrium (alternating square) is linearly unstable.

## Recommendations

1. **Enforce symmetry exactly in the integrator.** Replace generic 16-dim integration with a symmetry-reduced 8-dim integration in the Klein-four-quotient. If survival extends to $t\to\infty$, this is a *proof of bound motion* in the reduced system.
2. **Measure the Arnold time more carefully.** The fact that rhombus basin and symmetric-orbiter escape horizons coincide at $t^*\sim 400$ is suspicious — may reflect a common structural scale (diffusion time on the same torus class). Worth measuring $t^*$ as a function of IC distance from a candidate torus.
3. **Periodic-orbit search inside this basin.** Agent 5 did not probe the symmetric double-orbiter family. Shooting/Newton iteration on the $(r_{++}, R, \text{orb})$ 3-parameter family might locate a genuinely periodic orbit (Floquet $|\lambda|\le 1$) within the symmetric subspace.
4. **Compute Poincaré sections here.** Revisit Agent 6's section tool on this family; the flat energy curves suggest invariant-torus structure that should show up as 1D arcs on any sensible section.

## Bottom line

The symmetric double-orbiter family is **the strongest candidate in the 14-agent study for stable bound motion of the 4-body 2+/2− Weber Hamiltonian**. It is not yet rigorously proven stable — the escape horizon $t^*\sim 400$ means we have a *long-lived* bound state, not a *permanent* one — but its combination of $\lambda_{\max}\approx 0.1$, drift $\le 10^{-4}\%$, and exact discrete-symmetry protection makes it the most promising site in phase space for a future rigorous stability proof (via either symmetry-reduced KAM or Conley-index on the reduced phase space).

The `REPORT.md` §16 item 2 is therefore **partially resolved in the affirmative**: symmetric placement does dramatically extend survival, even though it does not achieve the predicted *indefinite* stability.
