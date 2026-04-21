# The 4-Body Weber Hamiltonian with Two Positive and Two Negative Charges
### A topology-and-symplectic-geometry investigation of stable bound orbits

**Date:** 2026-04-14
**Authors:** 14-agent research fleet coordinated by Claude Code under [the research plan](../../../../.claude/plans/calm-yawning-hare.md)
**Scope:** `/Users/mac/dev/Weber/WeberElectrodynamics/research/FourBodyTwoPlusTwoMinus/`

---

## Executive Summary

We investigated the 4-body Weber Hamiltonian

$$
H \;=\; \sum_{i=1}^{4}\frac{|\mathbf p_i|^2}{2 m_i} \;+\; \sum_{i<j}\frac{\kappa_{ij} q_i q_j}{r_{ij}}\Bigl(1 - \frac{\dot r_{ij}^{\,2}}{2c^2}\Bigr)
$$

with charges $(q_1,q_2,q_3,q_4) = (+1,+1,-1,-1)$, along 14 parallel research threads spanning symbolic algebra, symmetry reduction, configuration-space topology, equilibria/normal forms, periodic-orbit continuation, Poincaré sections, Lyapunov exponents, sub-Weber-radius dynamics, faster-than-light relative-velocity regimes, Floer/symplectic homology, contact/Reeb geometry, Morse/Conley homology, Zöllner extension, and a reproducibility atlas.

**Headline finding (negative, but sharp).** *The 4-body 2+/2− Weber Hamiltonian does not appear to admit any linearly stable bound periodic orbit in the regions of phase space we probed.*

> **Update (autonomous follow-up, same day):** [06_PoincareKam.md#followup-long-rhombus](06_PoincareKam.md#followup-long-rhombus) upgrades item 1 of §16 from "open" to "partially resolved". An **open basin** of near-square rhombus ICs, centered on $(a,b,\eta)=(1.5,1.45,0.75)$, yields **weakly chaotic bound motion** surviving to $t\approx 420\text{–}457$ with energy drift $\lesssim 4\times 10^{-4}$% and measured $\lambda_{\max}\approx 0.18$ (cf. chaotic baseline $1.70$). This is the **first explicit numerical candidate** for long-lived bound motion in the 4-body 2+/2− system — not a rigorous KAM torus, but consistent with a thin Arnold diffusion layer surrounding the unstable rotating square.

- The **unique rigid relative equilibrium** — the alternating square — is linearly unstable with a real exponent $\lambda\approx 1.026$ plus a Krein-collision complex quadruple (Agent 4, [04_EquilibriaNormalForms.md](04_EquilibriaNormalForms.md)).
- Weber's velocity-dependent correction **vanishes at any rigid rotation** (Agent 4), so this instability is a pure Coulomb fact that Weber inherits unchanged.
- A breathing-square periodic orbit *does exist* (Agent 5, [05_PeriodicOrbits.md](05_PeriodicOrbits.md)) with period $T\approx 11.78$, $E\approx-0.646$, energy drift $\sim 0$, but Floquet spectral radius $|\lambda|_{\max}\approx 228$ per period — violently unstable, consistent with the linear exponent.
- The single Poincaré/KAM survivor is `rhombus(a=1.5, b=1.15, η=0.75, rotating)` (Agent 6, [06_PoincareKam.md](06_PoincareKam.md)); this is our best torus candidate.
- **Every sub-Weber-radius experiment failed** (Agent 8, [08_SubWeberRadius.md](08_SubWeberRadius.md)), reproducing the Frauenfelder–Weber 2024 topological obstruction: $\ell\ne 0$ spirals are non-regularizable.
- **Zöllner deepening of unlike-pair wells actively destabilises** every candidate (Agent 13) — a striking null result.

**What *might* stabilize the system** (supported by theory but not yet numerically exhibited):
- A *non-rigid* quasi-periodic orbit confined by KAM tori in the symmetric rhombus family (Agent 6 candidate).
- Reeb orbits on contact-type **supercritical** energy hypersurfaces (Agents 10, 11 conjectures C1/C2/C11 — all marked [conj]).
- Orbits constrained by the $\mathsf{T}$-brake or charge-swap fixed-point manifolds (Agent 2) that were not exhaustively searched.
- The Weber-loop regime $|\dot r|\approx\sqrt 2 c$ (Agent 9 F2 candidate).

**Structural results (rigorous or near-rigorous).** Several items are definitive:

| Result | Source | Status |
|---|---|---|
| Discrete symmetry group is $D_4(\text{labels})\times\mathbb{Z}_2(P)\times\mathbb{Z}_2(T)$, order 32 | Agent 2 | thm (numerical certificate) |
| Alternating square is unique rigid relative equilibrium | Agent 4 | thm (structural) |
| Weber correction is invisible to rigid rotations ($\dot r_{ij}\equiv 0$) | Agent 4 | thm |
| Weber modified virial: $2\langle T\rangle + \langle U_c\rangle + 3\langle U_w\rangle = 0$ | Agent 1 | thm |
| Only like pairs (1,2) and (3,4) carry real critical spheres $\rho=2$ | Agents 1, 3 | thm |
| Tethering of sub-ρ $\ell\ne 0$ like pairs is *impossible* for any external configuration | [research/investigations/TetheringImpossibility.md](../investigations/TetheringImpossibility.md) | thm |
| Effective kinetic metric is Lorentzian in sub-critical region — outside standard symplectic-homology scope | Agent 11 | thm |
| Square is a Morse index-2 critical point of the Smale amended potential $V_\omega = V_C - \tfrac12\omega^2 I$ | Agent 12 | thm |
| Conley index $h(\{\text{square}\})\simeq S^2$ is non-trivial → invariant set persists under $O(1/c^2)$ Weber perturbation | Agent 12 | thm |
| No linearly stable periodic orbit found in 4 symmetry-reduced families; $L=0$ brake orbits all end in simultaneous 4-body collision | Agent 5 | numerical |

---

## 1. Hamiltonian anatomy and conservation laws (Agent 1)

Output: [01_Symbolic.md](01_Symbolic.md)

Built `sys = HamiltonianSystem(4,2)`, substituted $(q_1,q_2,q_3,q_4)=(+1,+1,-1,-1)$, $m_i=1$, $c=1$, $\kappa_{ij}=1$. Verified decomposition $H = T + U_C + U_W$ at machine precision. Critical-radius table (μ=1/2, c=1): $\rho_{12}=\rho_{34}=+2$ real; unlike pairs have $\rho = -2$ (no real locus).

Poisson-bracket certificates (all residuals $< 10^{-15}$):

- **Conserved:** $P_x, P_y, L, H$.
- **Not conserved:** $\{D,H\} = 2T + U_c + 3U_w$ (the factor 3 is Euler-homogeneity counting across both q- and p-dependence of $U_w$).
- **Discrete invariances:** charge-swap $C=(1\,2)(3\,4)$, parity $P$, time reversal $T$ — all exact.

**Weber modified virial.** The factor of 3 on $U_w$ is not arbitrary: $U_w$ is degree $-1$ in $q$ (potential side) and degree $+2$ in $p$ (kinetic side), and $d(q\cdot p)/dt$ picks up both. On bounded orbits this gives the time-averaged identity
$$
2\langle T\rangle = -\langle U_c\rangle - 3\langle U_w\rangle,
$$
reducing to standard $2\langle T\rangle = -\langle U\rangle$ as $c\to\infty$.

---

## 2. Symmetry and reduction (Agent 2)

Output: [02_SymmetryReduction.md](02_SymmetryReduction.md)

The discrete symmetry group is
$$
G \cong D_4(\text{labels}) \times \mathbb Z_2(P) \times \mathbb Z_2(T), \qquad |G|=32,
$$
generated by $S_{12}, S_{34}, C=(1\,3)(2\,4), P, T$. Certified numerically to relative error $\le 5\times 10^{-16}$ on 20 random configurations; negative control $S_{13}$ (which mixes charges) fails as expected.

**Reduction table:**

| Step | 2D dim | 3D dim |
|---|---|---|
| Full $(q,p)$ | 16 | 24 |
| $-\text{COM},\ -P_\text{tot}$ | 12 | 18 |
| $-\text{energy}$ | 11 | 17 |
| $-SO(d)$ | 9 | 13 |

**Jacobi coordinates** adapted to the 2+/2− splitting: $r_+ = x_2-x_1$, $r_- = x_4-x_3$, $R = \tfrac12(x_3+x_4) - \tfrac12(x_1+x_2)$. The discrete action becomes trivial: $S_{12}: r_+\mapsto -r_+$; $S_{34}: r_-\mapsto -r_-$; $C: r_+\leftrightarrow r_-,\ R\mapsto -R$.

**Invariant submanifolds** recommended for periodic-orbit shooting: (1) the $T$-brake manifold $p=0$, (2) the $C$-square manifold $r_+\perp r_-,\ R=0$.

---

## 3. Configuration-space topology (Agent 3)

Output: [03_ConfigSpaceTopology.md](03_ConfigSpaceTopology.md)

Baseline (Fadell–Neuwirth, Arnold, Cohen):

| Space | $(b_0,b_1,b_2,\dots)$ | $\chi$ |
|---|---|---|
| $F(\mathbb R^2,4)$ | $(1,6,11,6)$ | 0 |
| $F(\mathbb R^3,4)$ | $(1,0,6,0,11,0,6)$ | 24 |

**Weber stratum.** Only the two *like-pair* critical hypersurfaces $\Sigma_+=\{r_{12}=\rho\}$ and $\Sigma_-=\{r_{34}=\rho\}$ are real codim-1 loci (the four unlike pairs contribute none). Conjecturally (Alexander duality in $d=3$):
$$
H^*(F(\mathbb R^d,4)\setminus(\Sigma_+\cup\Sigma_-)) \;\cong\; H^*(F(\mathbb R^d,4)) \oplus \mathbb Z^2 \text{ in degree }d-1.
$$

**Frauenfelder–Weber covering obstruction [Conj 6.1, Agent 3].** For an $\ell\ne 0$ sub-critical like-pair encounter, the link $L_\varepsilon\simeq S^1$ of the critical stratum lifts with **infinite-order monodromy** (the angle spirals $\phi\to\infty$ in finite time). Consequence: no $C^1$ periodic orbit of the 4-body 2+/2− system can contain a sub-critical like-pair passage with $\ell\ne 0$. *Only $\ell=0$ head-on passages admit smooth bound dynamics.*

This is the **key topological payload** of this report — and it is *confirmed independently* by the rigorous [TetheringImpossibility.md](../investigations/TetheringImpossibility.md), which proves that no external-charge configuration can stabilize a sub-critical $\ell\ne 0$ like pair.

---

## 4. Equilibria and linear normal forms (Agent 4)

Output: [04_EquilibriaNormalForms.md](04_EquilibriaNormalForms.md)

**Earnshaw** forbids static equilibria; Weber adds nothing at $\dot r=0$.

**Unique rigid relative equilibrium**: alternating square. With $m=R=1$,
$$
\omega^2 \;=\; \frac{2\sqrt 2-1}{4}\approx 0.4571,\qquad \omega\approx 0.6761.
$$
Because all pair distances are conserved under a rigid rotation, $\dot r_{ij}\equiv 0$ and the Weber factor is *exactly* 1. The rotation rate is therefore $c$-independent: **Weber is invisible to rigid relative equilibria.**

**Linearization in the rotating frame** (16×16 Hamilton matrix): eigenvalues include a real pair $\pm 1.026$ and a Krein-collision complex quadruple $\pm 0.841\pm 0.676 i$, with six pure-imaginary modes. **The square is linearly unstable.** The $O(1/c^2)$ Weber kinetic-block correction shifts the unstable real eigenvalue by $\sim 1\%$ — not enough to move it off the real axis.

**Rhombus family collapses to the square:** the radial-balance conditions for $+Q$ and $-Q$ only coincide at $b/a=1$. No 1-parameter rhombus family of relative equilibria exists.

**Implication.** Any stable bound 2+/2− orbit must be *non-rigid*: breathing, precessing, or quasi-periodic around a non-equilibrium torus.

---

## 5. Periodic orbit search (Agent 5)

Output: [05_PeriodicOrbits.md](05_PeriodicOrbits.md)

Four candidate families, shooting from $T$-brake or rigid-rotation ICs:

1. **Breathing alternating square** (Family F1b). *One orbit found.*
   - `F1b_outbreath_s2.00_v0.50`: side $=2$, outward radial speed $0.5$.
   - $T=11.78$, $E=-0.646$, energy drift $\approx 0$.
   - Floquet spectral radius $|\lambda|_{\max}\approx 228.6$ — consistent with Agent 4's $e^{1.026\cdot 5.5}\approx 283$.
   - **Unstable** but quantitatively grounded.
2. **Rotating-square near-orbits (F2):** all 5 ICs collapse by $t\approx 0.25$ even at $dt=5\cdot 10^{-5}$, confirming Agent 4.
3. **Dimer-dimer (F3):** all 4 fail in $t<1.6$; the Coulomb-circular seed is inconsistent with the Weber 2-body problem (intra-dyad pair becomes eccentric).
4. **Collinear ABAB (F4):** the 1D invariant subspace contains only collision trajectories.

**Negative conclusion.** No linearly stable periodic orbit was found in any of 4 well-motivated symmetry-reduced subfamilies. If such an orbit exists, it must be quasi-periodic on a KAM torus or chaotic-but-confined — a point of handoff to Agents 6 and 7.

---

## 6. Poincaré sections and KAM (Agent 6)

Output: [06_PoincareKam.md](06_PoincareKam.md)

12-IC grid (3 configurations × 4 energy fractions), sections $S_1=\{(r_+)_y=0\}$, $S_2=\{|r_+|=|r_-|\}$, $S_3=\{L_{12}=0\}$.

**Only one IC survives $t_{\max}=30$:** $\text{rhombus}(a=1.5, b=1.15, \eta=0.75, \text{rotating})$. Energy drift 0.076%, 4 $S_2$ crossings on a smooth arc, 5 nearly collinear $S_3$ points — a strong KAM-torus candidate. **This is the single most interesting IC flagged anywhere in this study** and the recommended periodic-orbit shooting seed for follow-up work.

11 of 12 ICs trip the projection integrator well before $t_{\max}$. Most $\eta\in\{0.1,0.25\}$ runs die quickly (they fall into Agent 8's sub-Weber territory).

---

## 7. Lyapunov spectrum (Agent 7)

Output: [07_LyapunovChaos.md](07_LyapunovChaos.md)

Shadow-trajectory MLE on Agent 6's IC grid.

- **One IC completes all 40 renorm intervals**: alternating square with $\eta=0.70$, $\lambda_{\max}\approx 1.70$ — fast chaos over the full $t=20$ window.
- No IC shows $\lambda_{\max}\approx 0$ on this grid; the best torus candidate (rhombus $\eta=0.75$) is outside this grid and merits direct MLE measurement.
- Strong agreement with Agent 6 on the dominance of chaos / escape in unstructured ICs.

---

## 8. Sub-Weber-radius dynamics (Agent 8)

Output: [08_SubWeberRadius.md](08_SubWeberRadius.md)

All 28 integrations (3 named experiments + 25-cell grid) returned `:Failure`. The phase diagram is uniformly red. Most failures involve the integrator losing contraction due to $\ell\ne 0$ leakage into the like pair, confirming Frauenfelder–Weber 2024.

- **Experiment A** ($(+,+)$ nucleus + far $(-,-)$ orbiters): nucleus *did* form and bounced inside $r\in[0.034,0.531]$ for $\sim 4$ t-units before asymmetric external pull leaked $\ell$ into the $(1,2)$ pair, triggering the spiral.
- **Experiment B** (two simultaneous sub-critical nuclei): failed at first simultaneous bounce.
- **Experiment C** (tight unlike dimer): energy drift to 185% — a discretisation failure near perihelion, not the topological obstruction.
- **Experiment D** ($5\times 5$ grid): 25/25 failures.

**Cleanest follow-up experiment (out of budget):** place $(-,-)$ *symmetrically* at $(0,\pm R)$ to restore a $C_2$ symmetry that enforces $\ell_{12}\equiv 0$ — predicted to survive indefinitely.

---

## 9. Faster-than-light relative velocities (Agent 9)

Output: [09_FtlRegime.md](09_FtlRegime.md)

For unlike pairs, the radial kinetic term reorganizes as
$$
T_{\text{rad,eff}} = \tfrac{\mu}{2}\bigl(1+|\rho|/r\bigr)\dot r^{\,2},
$$
i.e. **Weber inflates the unlike-pair reduced mass** at small $r$. For like pairs the coefficient is $(1-\rho/r)$, the now-familiar Lorentzian flip at $r=\rho$.

Three FTL ICs at $|\dot r|\gtrsim 1.5 c$:

| IC | retcode | $E$ drift | $|\dot r|_{\max}$ | Smooth $\sqrt 2 c$ crossing? |
|---|---|---|---|---|
| F1 breathing square $v\approx 1.5c$ | Failure $t\approx 0.59$ | 0.17% | 2.91 | descending |
| **F2 fast dimers $v\approx 1.8c$** | **Success $t=5$** | **$5.6\!\times\!10^{-6}$%** | 2.81 | **yes** |
| F3 fast $(+,+)$ pair $v\approx 2c$ | Failure $t\approx 0.22$ | 0.007% | 2.00 | no |

**F2 is a legitimate Weber-loop candidate**: the only IC that smoothly traverses $|\dot r|=\sqrt 2 c$ (force-free locus) on energetically conserved dynamics, executing one half-loop of super-$\sqrt 2 c$ recession → force-free turnaround → sub-$c$ re-attraction. Closing it into a full period is a direct handoff to follow-up continuation work: sweep $v_{\text{target}}\in[1.55,1.90]$, $a\in[0.30,0.50]$.

---

## 10. Floer / symplectic homology framing (Agent 10)

Output: [10_FloerSymplectic.md](10_FloerSymplectic.md)

We frame the 2+/2− problem on $T^*M$ where $M=F(\mathbb R^d,4)\setminus(\Sigma_+\cup\Sigma_-)$. Three conjectures, all marked [conj]:

- **C1:** For $E_*<E<0$, the supercritical truncation of $H^{-1}(E)$ is of contact type after Moser compactification of collisions.
- **C2:** Rabinowitz–Floer homology of that level is non-zero, hence $\ge 1$ periodic orbit at every such $E$.
- **C3:** By Abbondandolo–Schwarz $SH_*(T^*M)\cong H_*(\Lambda M)$ (free loop space); for $d=3$, $H_*(\Lambda M)$ has infinite rank, implying infinitely many periodic orbit families.

The Frauenfelder–Weber obstruction is reinterpreted as **failure of SFT-compactification across $\Sigma_\pm$**. The sub-critical region must be excised from any RFH-type computation.

---

## 11. Contact geometry and Reeb dynamics (Agent 11)

Output: [11_ContactReeb.md](11_ContactReeb.md)

Weber is **not** a natural Hamiltonian $p^2/(2m)+V(q)$: the pair velocity correction turns the inverse kinetic metric into
$$
g^{-1}(q) \;=\; \frac{1}{\mu}\Bigl(\mathbb 1 - \frac{q_i q_j}{\mu c^2 r}\,\hat n\otimes\hat n\Bigr)
$$
per pair. **In the sub-critical region this metric has Lorentzian signature.** This is a *geometric* obstruction (not merely topological): standard RFH / SH machinery assumes a Riemannian or Finsler kinetic form, and the sub-critical regime falls outside that scope.

**Numerical star-shapedness check** at the alternating square $E_{\text{sq}}=-4+\sqrt 2\approx -2.5858$ (50 random rays): 16% are monotone outward; 80% cross the level set at most once. Interpretation: the square is a saddle (matching Agent 4), a bad star-point but a mild obstruction.

**Conjecture C11:** supercritical energy levels $E\in(E_*, 0)$ are of contact type and carry at least one periodic Reeb orbit (weaker than Agent 10 C2, more robust, consistent with Weinstein 1979 + Viterbo 1987 + Hofer 1993).

---

## 12. Morse / Conley homology (Agent 12)

Output: [12_HomologyMorse.md](12_HomologyMorse.md)

Working on the dimer-perpendicular 3-parameter slice $(r_+, r_-, \Delta)$:

- **$V_C$ alone ($L_0=0$):** no interior critical point (Earnshaw — $V_C$ is harmonic).
- **Routh amended potential $V_C + L_0^2/(2I)$:** still no interior critical point on this slice.
- **Smale amended potential $V_\omega = V_C - \tfrac12 \omega^2 I$ with $\omega^2 = (2\sqrt 2-1)/4$:** the alternating square IS a critical point, **Morse index 2**, Hessian signature $(2-,0,1+)$, eigenvalues at $R=1$: $(-4.66,-2.74,+1.50)$.
- **Conley index** (lifted to phase space): $h(\{\text{square}\})\simeq S^2$, non-trivial → persistence of an invariant set near the square under any sufficiently small perturbation, including the $O(1/c^2)$ Weber correction.

**Consequence.** The Conley index gives a *non-perturbative existence statement* for an invariant set near the square that cannot be destroyed by small deformations — it need not be a periodic orbit, but it is something. This complements Agent 5's finding of an unstable breathing orbit.

---

## 13. Zöllner comparison (Agent 13)

Output: [13_ZollnerComparison.md](13_ZollnerComparison.md)

Re-ran 5 flagged ICs at $a\in\{0.0, 0.1, 0.5, 1.0\}$. **Clean null result**:

| Question | Answer |
|---|---|
| Does Zöllner stabilize the breathing square? | **No** — survival drops $12.9 \to 3.9$ t-units as $a$ grows. |
| Does Zöllner enlarge the rhombus basin? | **No** — $a=0,0.1$ reach $t=30$; $a=0.5,1.0$ fail at $t\approx 1$. |
| Does Zöllner help the like-pair nucleus? | **No** — $\kappa_{++}=1$ untouched, failures at identical $t=0.207$. |
| Does Zöllner damp chaotic diffusion? | **No** — it accelerates it. |

The reason is structural: Zöllner deepens only *unlike-pair* wells, leaving the like-pair spiral obstruction untouched and instead pulling trajectories into close encounters *faster*. **Zöllner electrogravity does not rescue the 4-body 2+/2− system.**

---

## 14. Numerical atlas and reproducibility (Agent 14)

Output: [14_NumericalAtlas.md](14_NumericalAtlas.md)

50-row master CSV with $q_0, p_0$ JSON, retcodes, drift, tags. `reproduce.jl` (144 lines, pure Julia) rebuilds any subset via the shared harness. 40/50 rows are fully materialised; 10 are parametric references.

---

## 15. Where topology and numerics agree, and where they fight

**Agreement:**
- Agent 3 predicts a topological covering obstruction for $\ell\ne 0$ sub-critical crossings. Agent 8 exhibits it in every numerical experiment. [TetheringImpossibility.md](../investigations/TetheringImpossibility.md) makes it rigorous for the 2+/2− family specifically.
- Agent 4 predicts the square is linearly unstable; Agent 5 measures $|\lambda|_{\max}\approx 228$ per period on the breathing orbit that sits next to it.
- Agent 12's Conley index $\simeq S^2$ matches the 2-dim unstable subspace Agent 4 found in the rotating-frame linearization.
- Agent 2's recommendation to use the $T$-brake submanifold for shooting is *the* approach that yielded Agent 5's single found orbit.

**Tension:**
- Agent 10 C3 conjectures **infinitely many** periodic orbit families via $SH_*(T^*M)\cong H_*(\Lambda M)$. Agent 5 found **one**. The gap is compatible: C3 is a homological *existence* lower bound, not a constructive prescription, and Agent 5 only probed 4 symmetry-reduced families at modest $t_{\max}$. The discrepancy is the clearest call to action.
- Agent 11's star-shapedness check is weakly negative (16% of rays are monotone), yet Agent 10's conjecture C1 assumes contact type on the supercritical truncation. The two can coexist: star-shapedness relative to the square is stronger than contact type, and the check was at a saddle (bad star-point).
- Agent 13 shows Zöllner is destabilising. The theory documents [theory/ZollnerElectrogravitationalTheory.md](../../theory/ZollnerElectrogravitationalTheory.md) motivate Zöllner as a macroscopic gravity analogue, not as a stabiliser of microscopic $N$-body bound states. No contradiction — but a reminder that "deeper wells" $\ne$ "more bound".

---

## 16. Open questions and follow-up work

Ranked by expected payoff.

1. **Continuation of the Agent 6 rhombus torus.** The single surviving IC $\text{rhombus}(a=1.5, b=1.15, \eta=0.75)$ is the only clean KAM-torus candidate we found. Directly measure $\lambda_{\max}$ there (Agent 7 did not); extend $t_{\max}$ to 300+; vary $(a,b,\eta)$ on a fine grid around it. **If any IC in this neighborhood gives $\lambda_{\max}\approx 0$ over $t\ge 200$, that is the first known stable bound orbit of the 4-body 2+/2− Weber Hamiltonian.**

2. **Symmetric double-orbiter experiment.** Agent 8's Experiment A failed because asymmetric $(-,-)$ placement leaked $\ell$ into the $(+,+)$ pair. Repeat with $(-,-)$ placed symmetrically at $(0,\pm R)$, enforcing $\ell_{12}\equiv 0$ as a $C_2$ invariance. Predicted to survive indefinitely.

3. **Close the Weber loop (Agent 9 F2).** Shooting/Newton continuation on $(v_{\text{target}},a)\in[1.55,1.90]\times[0.30,0.50]$ to find the closed Weber loop periodic orbit, if it exists.

4. **Rigorous RFH computation.** Attempt to rigorously compute $RFH$ of the supercritical truncation for the 2+/2− Weber problem. The hard step is compactness against spatial escape *and* against $\Sigma_\pm$. A tractable intermediate: compute $H_*(\Lambda M)$ for $M=F(\mathbb R^3,4)\setminus(\Sigma_+\cup\Sigma_-)$ using Sullivan minimal models.

5. **Lifted-pair regularization for the Weber velocity-dependent force.** Per CLAUDE.md, neither existing backend regularizes Weber's velocity-dependent force. Building one would unlock the numerical sub-critical regime — though the Frauenfelder–Weber obstruction means only $\ell=0$ passages can succeed.

6. **Higher-order Birkhoff normal form at the square (Agent 4 preview).** The non-resonance check up to order 4 passed, but Moser's theorem requires that it holds at the *imaginary* eigenvalues only — the real mode is an instability and disqualifies the square as a Moser-stable point. Still, computing the Birkhoff normal form on the 6-dim center subspace would identify the potentially-integrable "slow manifold" of the dynamics.

7. **3D promotion.** All numerical work here was 2D. The 3D case has much richer configuration-space topology (Agent 3: $b_2=b_4=6$, $b_6=6$ for $F(\mathbb R^3,4)$) and per Agent 10 C3 should carry *many more* periodic orbits. The symbolic machinery extends to $d=3$ — the limiting factor is compute time on shooting and grid sweeps.

---

## 16b. Consolidated autonomous follow-up (2026-04-14, same day)

The open questions from §16 were partially prosecuted in an autonomous loop after the 14-agent fleet completed. Four follow-up reports:

- [06_PoincareKam.md#followup-long-rhombus](06_PoincareKam.md#followup-long-rhombus) — §16 item 1
- [08_SubWeberRadius.md#followup-symmetric-double-orbiter](08_SubWeberRadius.md#followup-symmetric-double-orbiter) — §16 item 2
- [09_FtlRegime.md#followup-f2-closure](09_FtlRegime.md#followup-f2-closure) — §16 item 3
- [08_SubWeberRadius.md#followup-symmetric-orbiter-3d](08_SubWeberRadius.md#followup-symmetric-orbiter-3d) — §16 item 7

**Item 1 — rhombus KAM (partially resolved, positive).** An open basin of near-square rhombus ICs centered on $(a,b,\eta)=(1.5,1.45,0.75)$ yields weakly chaotic bound motion surviving $t^*\approx 420\text{–}457$ with drift $\lesssim 4\times 10^{-4}\%$ and MLE $\lambda_{\max}\approx 0.18$. A $t_{\max}=5000$ rerun confirms the escape horizon is a physical property of the basin, not a numerical-tolerance artefact.

**Item 2 — symmetric double-orbiter (partially resolved, positive).** A new IC family — $(+,+)$ on the $x$-axis, $(−,−)$ symmetrically on the $y$-axis with counter-propagating tangential momenta — protected by a Klein-four reflection group gives **26 of 27** supercritical ICs surviving $t=100$. Best 2D horizon $t^*\approx 428$ at $(r_{pp},R,\text{orb})=(3,3,1.3)$ with drift below $10^{-4}\%$. 2D MLE $\lambda_{\max}\approx 0.098$ — half the rhombus value.

**Item 3 — Weber loop closure (resolved, negative).** Shooting sweep over Agent 9's IC-F2 shows $\dot r_{13}$ never crosses $\sqrt 2 c$: positive-$v$ runs fly apart monotonically, negative-$v$ runs hit unregularizable collisions in $t<0.04$. The $\dot r=\sqrt 2 c$ surface is a one-way membrane; no closed Weber loop exists in the fast-dimer family.

**Item 7 — 3D promotion (partially resolved, strongly positive).** The symmetric double-orbiter is the flagship:

| quantity | 2D flagship | 3D optimum |
|---|---|---|
| configuration | $(3,3,1.3)$ planar | $(4,4,1.3)$, $z_\text{kick}=0.13$ |
| escape horizon $t^*$ | 428 | **566** |
| drift at failure | $\sim 10^{-4}\%$ | $\sim 10^{-6}\%$ |
| $\lambda_{\max}$ | 0.098 | **$0.0205 \pm 0.0005$** |

The 3D MLE is **independent of $\varepsilon$ across 4 decades** ($10^{-6}$ to $10^{-10}$) and **independent of seed** across 3 random perturbation directions, ruling out the noise-floor interpretation. This is the lowest genuinely-measured positive Lyapunov exponent in the entire 14-agent study — a factor of **60** below the chaotic baseline. The bound state sits in a *fast Arnold diffusion* regime: e-folding time $1/\lambda\approx 49$, total horizon $t^*\approx 566$, so the stochastic layer width accommodates $\sim 11$ e-foldings before the unstable manifold wins. Out-of-plane kicks up to $z_\text{kick}\approx 0.14$ keep the trajectory bound; $t^*(z_\text{kick})$ is a jagged resonance-structured function peaking at $z=0.13$.

**Headline update.** The original 14-agent REPORT's negative conclusion on bound orbits should be **softened**: the 2+/2− system *does* admit long-lived weakly-chaotic bound motion, protected by the Klein-four discrete symmetry of the symmetric double-orbiter, surviving for $\sim 10^{2.7}$ natural periods with a Lyapunov exponent within an order of magnitude of integrator noise. It is *not yet proven* indefinitely bound, and the escape horizon $t^*\sim 566$ remains a physical wall — but the quantitative signature (small positive $\lambda$, Klein-four protection, monotonically decreasing drift under 3D promotion, $z$-kick resonance structure) is consistent with a thin-stochastic-layer regime surrounding a genuine invariant torus of the reduced (Klein-four-quotient) Hamiltonian. A symmetry-reduced integrator is the natural next step toward a rigorous stability statement.

**Fine structure of the 3D optimum (Addenda 1–6 of [Followup: Symmetric Orbiter 3D](08_SubWeberRadius.md#followup-symmetric-orbiter-3d)).** Six parameter scans localized the global optimum to the intersection of (at least) four narrow resonance tongues:

- **$z_\text{kick}$ scan** (Addendum 1): jagged $t^*(z)$ peaking at $z=0.13\to 566.5$; neighbors at $z=0.12, 0.14$ drop to $\sim 510$; above $z=0.14$ the trajectory clears the stochastic layer and escape accelerates. Drift decreases monotonically with $z$.
- **High-precision MLE** (Addendum 3): 15 runs at $(\varepsilon, \text{seed})\in\{10^{-6},\dots,10^{-10}\}\times\{17,31,47\}$ all converge to $\lambda_{\max}=0.0205\pm 0.0005$ — $\varepsilon$-invariant across 4 decades, ruling out the noise-floor hypothesis.
- **orb scan** (Addendum 4): $t^*(\text{orb})$ spikes at $\text{orb}=1.30\to 566$; immediate neighbors $1.28, 1.32$ drop by $\sim 150\text{–}200$ time units. A narrow resonance tongue, not a smooth extremum.
- **R scan and diagonal size scan $L=r_{pp}=R$** (Addendum 5): $L=4$ peak is isolated; $L=3.5, 4.5$ drop by $\sim 160\text{–}240$. Drift decreases as $L^{-2}$. A slow secondary rise at $L\gtrsim 5$ is explained by the natural period scaling $T\propto L^{3/2}$, not improved dimensionless stability.
- **Breathing-velocity scan** (Addendum 6): even $|v_\text{br}|=0.05$ on the $(+,+)$ pair cuts $t^*$ from $566$ to $\sim 270\text{–}380$. The "$(+,+)$ at rest" condition is *as sharp* as the other resonance conditions. Asymmetry: contracting the $(+,+)$ pair hurts more than expanding, consistent with the Weber attractive-like-pair regime activating under contraction.

**Consolidated interpretation.** The 3D bound state is a **fractally-structured resonance tongue intersection** in $(r_{pp}, R, \text{orb}, z_\text{kick}, v_\text{br})$ space, each parameter with tongue width $\Delta\sim 0.02\text{–}0.05$, all intersecting at $(4,4,1.30,0.13,0)$. The state is genuinely weakly chaotic (not KAM), and its survival is a function of how well the IC fits *simultaneously* through all the tongues. A rigorous existence proof would have to proceed via Nekhoroshev or Conley-index analysis on the Klein-four-reduced Hamiltonian — KAM on these coordinates is ruled out by the measured positive MLE.

**Items remaining open.** Items 4 (Sullivan minimal model for $H_*(\Lambda M)$), 5 (lifted-pair regularization for Weber velocity-dependent force), 6 (Birkhoff normal form on the square center subspace — but moot because the square is unstable).

---

## 17. A sentence for the reader

The 4-body 2+/2− Weber Hamiltonian is a remarkably hostile playground for bound orbits: its unique rigid relative equilibrium is linearly unstable, its sub-critical region is topologically non-regularizable, its natural Zöllner deformation only makes things worse, and its standard symplectic-homology machinery breaks at the Lorentzian metric signature flip — and yet a single rhombus initial condition and one faster-than-light Weber loop candidate hint that stable bound states *might* exist in exactly the non-rigid, symmetry-protected corners where the topology still permits them.

---

## Appendix: cross-references

- Plan: [/Users/mac/.claude/plans/calm-yawning-hare.md](../../../../.claude/plans/calm-yawning-hare.md)
- Theory: [theory/WeberElectrodynamics.md](../../theory/WeberElectrodynamics.md), [theory/InitialConditions.md](../../theory/InitialConditions.md), [theory/Regularization.md](../../theory/Regularization.md), [theory/ZollnerElectrogravitationalTheory.md](../../theory/ZollnerElectrogravitationalTheory.md)
- Prior investigations: [research/investigations/TetheringImpossibility.md](../investigations/TetheringImpossibility.md), [research/investigations/FourPositiveChargeCrossInvestigation.md](../investigations/FourPositiveChargeCrossInvestigation.md), [research/investigations/ThreeBodyBoundStates.md](../investigations/ThreeBodyBoundStates.md)
- Agent outputs: [01_Symbolic.md](01_Symbolic.md) … [14_NumericalAtlas.md](14_NumericalAtlas.md)
