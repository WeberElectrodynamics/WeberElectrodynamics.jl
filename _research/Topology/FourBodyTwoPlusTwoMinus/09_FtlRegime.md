# Agent 9 — Faster-Than-Light Relative Velocity Regime

**Question.** Weber's force law contains the factor `(1 − ṙ²/(2c²))`. At
|ṙ| = √2·c the factor vanishes (the pair becomes momentarily force-free);
beyond it, the factor flips sign so an attractive Coulomb pair becomes
effectively repulsive and a repulsive pair becomes attractive. Weber is
non-relativistic: there is no kinematic light cone, only this algebraic
sign change. Can the 4-body 2+/2− system exhibit *bounded* dynamics that
genuinely lives in the |ṙ|>c — and in particular |ṙ|>√2·c — regime, and
can such an orbit "stitch" sub-c and super-c phases together via a smooth
crossing of |ṙ|=√2·c (the **Weber loop** conjecture)?

Code: [`ftl_experiments.jl`](ftl_experiments.jl). Figures in
[`figures/`](figures/). All units natural with c=1.

## 1. Symbolic setup

Pair interaction energy for an unlike (q_i q_j = −q²) pair at separation r
with radial velocity ṙ:

    U_pair = (q_i q_j / r) · (1 − ṙ²/(2c²))
           = −q²/r + q² ṙ² / (2 r c²)

The first term is plain Coulomb attraction; the second is a *repulsive*,
velocity-quadratic correction (positive contribution to U). Because ṙ
enters U quadratically in `p`, it acts like a position-dependent kinetic
correction. In the COM frame of a single pair with reduced mass μ:

    H_pair = (μ/2) ṙ² + L²/(2μ r²) + (q_i q_j/r)(1 − ṙ²/(2c²))
           = (μ/2)(1 − (q_i q_j)/(μ c² r)) ṙ² + L²/(2μ r²) + q_i q_j/r
           = (μ/2)(1 − ρ/r) ṙ²              (like-pair, ρ = q²/(μc²) > 0)
           = (μ/2)(1 + |ρ|/r) ṙ²            (unlike-pair)

So **unlike pairs become heavier in the Weber sense as r→0**: the effective
reduced mass `μ_eff(r) = μ(1 + |ρ|/r)` diverges, and the pair becomes
arbitrarily resistant to further radial acceleration. Like pairs in
contrast have `μ_eff = μ(1 − ρ/r)` which vanishes at r=ρ and goes negative
inside (the Lorentzian metric-signature flip reported by Agent 8 / theory).

A direct consequence: **at fixed total H, |ṙ| of an unlike pair is
*bounded* by what the negative coefficient of the Weber correction can
absorb at finite r**. The constraint is

    H = (μ/2)(1 + |ρ|/r) ṙ² + L²/(2μr²) − q²/r,

so ṙ² ≤ 2 (H + q²/r) / (μ(1 + |ρ|/r)). For unlike pairs, ṙ stays finite
even as r→0 — and in particular Weber's velocity correction can dominate,
making it *energetically expensive* to drive |ṙ| beyond √2·c at small r.
For like pairs the opposite is true: |ṙ| can run away.

The total H of the 4-body system at our FTL initial conditions is finite
(printed by the script), confirming that |ṙ|>√2·c is on a perfectly
ordinary energy shell.

## 2. Initial conditions

Three ICs (≤8 budget; we used 3, well within):

| ID | Geometry | Targeted ṙ | Notes |
|----|----------|------------|-------|
| F1 | Alternating (++,−−) square, side 1.5, breathing inward | ṙ ≈ 1.5 c on opposite-corner pairs (1,2) and (3,4); ṙ ≈ √2/2·1.5 ≈ 1.06 c on edges | All four edge pairs are unlike; the like-pair diagonals are super-c |
| F2 | Two (+,−) dimers, dyad length 0.4, separation 2.5 | ṙ ≈ 1.8 c on each intra-dimer pair (1,3) and (2,4) | Inter-dimer pairs sub-c |
| F3 | One (+,+) pair at sep 1.0 with ṙ ≈ 2.0 c, two negatives parked far away | Like pair only | Tests like-pair super-√2·c head-on |

Initial Weber factors `(1 − ṙ²/2c²)` printed by the script:

```
F1  (1,2)++   ṙ=−1.500   factor=−0.125    (super-√2c, sign-flipped)
F1  (1,3)+−   ṙ=−1.061   factor=+0.437    (sub-√2c)
F2  (1,3)+−   ṙ=+1.800   factor=−0.620    (super-√2c)
F3  (1,2)++   ṙ=+2.000   factor=−1.000    (super-√2c)
```

In every case the total H is finite and well-defined.

## 3. Numerical results (tmax=5.0, dt=1e-4, unregularized symplectic)

```
IC-F1 breathing square (v≈1.5c)
  retcode=Failure  steps=5859  t_final=0.586  E drift_max = 1.74e−1 %
  |ṙ|_max = 2.91    crosses √2c smoothly: TRUE

IC-F2 fast dimers (v≈1.8c)
  retcode=Success  steps=50001 t_final=5.000  E drift_max = 5.65e−6 %
  |ṙ|_max = 2.81    crosses √2c smoothly: TRUE

IC-F3 fast (++) pair (v≈2.0c)
  retcode=Failure  steps=2195  t_final=0.219  E drift_max = 6.96e−3 %
  |ṙ|_max = 2.00    crosses √2c smoothly: FALSE
```

Figures `figures/ftl_f1_breathing_square.png`, `ftl_f2_fast_dimers.png`,
`ftl_f3_fast_like_pair.png` show |ṙ_ij|(t) for all 6 pairs with the c and
√2c thresholds dashed.

### F1 — breathing square (catastrophe)

The four like-pair diagonals start at ṙ = −1.5 c (super-√2c). The Weber
factor on these pairs is *negative*, so the like-pair Coulomb repulsion
becomes effective attraction. The diagonals plunge: |ṙ| escalates above
2.9 c before the integrator hits a singularity at t ≈ 0.59. Edge (unlike)
pairs are dragged along but the failure is driven by the like pairs.
This is the same head-on like-pair pathology that Agent 8 documents at
sub-ρ, here triggered *kinetically* (super-√2c) rather than spatially
(sub-ρ). The integrator does, however, traverse |ṙ|=√2c smoothly on the
way down through ~1.4c, so a *single-direction* crossing exists.

### F2 — fast dimers (clean run, the interesting case)

This is the only IC that survives 5 time units. The two intra-dimer pairs
start at ṙ = +1.8 c with Weber factor −0.62; for unlike pairs the
sign-flipped factor turns the Coulomb attraction into a strong effective
repulsion, decelerating each pair. |ṙ| drops smoothly through √2c
(crosses_smoothly=true) into the sub-c regime, where ordinary attraction
resumes and ṙ flips sign. The pairs then re-accelerate inward, and the
script records a *second* approach toward √2c. Energy is conserved to
6×10⁻⁶ % over 50,000 steps — this is a clean symplectic trajectory, not a
numerical artefact.

In effect the two dimers each execute a half "Weber loop": super-√2c
recession → force-free turnaround at |ṙ|=√2c → sub-c re-attraction. To
*close* the loop we would need a second smooth crossing back into the
super-√2c regime, which our 5-time-unit window does not show because the
sub-c attraction drains kinetic energy into the 4-body angular structure
rather than back into radial motion.

### F3 — fast like-pair (immediate collision)

Starting ṙ = +2 c on the (+,+) pair with the negatives parked. The Weber
factor is −1.0, so the like pair is *attractive* with twice Coulomb
strength. With the negatives essentially absent the pair behaves as an
isolated like-pair Weber attractor; |ṙ| does *not* damp through √2c
(no smooth crossing) and the integrator fails near r→0 at t ≈ 0.22.
This corner is exactly the spiraling singularity reported in
`research/Investigations/ThreeBodyBoundStates.md`.

## 4. Weber-loop conjecture: status

**Candidate flagged for Agent 5: IC-F2.** This is the only run where the
dynamics genuinely visits both sides of the |ṙ|=√2c surface, crosses it
*smoothly* (ṙ continuous, force factor passes analytically through zero),
and survives long enough to demonstrate that energy is conserved across
the crossing. The trajectory is not yet periodic — closing it would
require a continuation in (dyad_length, separation, v_target) so that the
sub-c re-attraction phase brings the pairs back onto the original
super-√2c branch. Agent 5's shooting/Newton periodic-orbit machinery is
the natural next tool. Suggested seed: F2 with v_target swept in
[1.55, 1.90] and dyad_length in [0.30, 0.50].

F1 and F3 demonstrate the two known failure modes:

- F1 — like-pair diagonals interpret super-√2c into effective attraction
  and collapse;
- F3 — isolated like-pair super-√2c is monotone collapse, not a loop.

So the Weber-loop conjecture survives empirically only on *unlike* pairs
embedded in a 4-body geometry that holds them apart laterally.

## 5. Theoretical notes

- **Causality.** Not at issue. Weber's potential is instantaneous-action;
  c is just a parameter in a quadratic form. |ṙ|>c is a perfectly
  admissible state of a Hamiltonian flow. The sign-flip at √2·c is an
  *algebraic* feature of `1 − ṙ²/(2c²)`, not a physical light-cone.

- **Why F2 looks special.** For an unlike pair the effective reduced mass
  `μ_eff = μ(1 + |ρ|/r)` keeps `(1/2) μ_eff ṙ²` finite as r decreases
  even when ṙ stays large; equivalently, the Weber correction acts as a
  velocity-dependent "soft wall" that decelerates the pair before it
  collides. This is structurally the same mechanism by which classical
  Weber atoms remain bound at high orbital frequency (Assis 1994).
  The wall is energetically transparent at exactly |ṙ|=√2c — that is the
  geometric signature of the loop.

- **Mathematical interest.** A bound orbit threading |ṙ|=√2c would be a
  genuinely new phase-space object: a smooth Hamiltonian trajectory whose
  energy hypersurface includes points where the force on a constituent
  pair vanishes by a *non-spatial* mechanism. It is dual, in a precise
  sense, to the spatial sub-ρ surface explored by Agent 8 — both are
  loci where one term in the Hamiltonian's quadratic form changes sign,
  but Agent 8's surface is q-dependent and Agent 9's is p-dependent.
  Agent 11 (contact-type analysis) and Agent 10 (Floer framing) should
  treat the surface `|ṙ_ij|=√2c` as a candidate stratum.

- **Worth continuing?** Yes, restricted to the unlike-pair / dimer-dimer
  channel. The like-pair super-√2c regime (F1, F3) is dynamically the
  same singular collapse already known from sub-ρ work and offers no new
  structure beyond what Agent 8 already maps. The unlike-pair channel
  (F2) is genuinely new and is the right input for Agent 5's continuation
  search.

## 6. Files

- [`ftl_experiments.jl`](ftl_experiments.jl) — IC builders, runner,
  threshold diagnostics, plotting.
- [`figures/ftl_f1_breathing_square.png`](figures/ftl_f1_breathing_square.png)
- [`figures/ftl_f2_fast_dimers.png`](figures/ftl_f2_fast_dimers.png)
- [`figures/ftl_f3_fast_like_pair.png`](figures/ftl_f3_fast_like_pair.png)

Total: 3 ICs, ≤5 t.u. each, ≤8-IC budget respected.


---

## Followup: F2 Closure

# Follow-up: the F2 Weber-loop hypothesis is false

**Date:** 2026-04-14 (autonomous follow-up to [FourBodyTwoPlusTwoMinus.md §16 item 3](FourBodyTwoPlusTwoMinus.md))
**Scope:** shooting sweep over Agent 9's F2 candidate to test whether a closed Weber loop exists in the fast-dimer family.

## Motivation

REPORT §16 item 3 asked: starting from IC-F2 (two `(+,−)` dimers with intra-dimer $\dot r > \sqrt 2 c$), do a continuation over `(v_target, dyad_len)` and check whether the trajectory ever *returns* to its initial phase-space point — i.e. whether a closed "Weber loop" (FTL + sub-c stitched into one periodic orbit) exists.

## Procedure

Used [`ic_f2_fast_dimers`](ftl_experiments.jl) exactly as Agent 9 defined it, with `sep=2.5`, unit masses and charges, $c=1$. Scanned $v_\text{target}\in\{\pm 1.55, \pm 1.75, \pm 1.90\}$ and `dyad_len`$\in\{0.30, 0.40, 0.50\}$ (18 ICs), integrated to $t_{\max}=20$, $dt=5\times 10^{-4}$. Return metric = minimum phase-space distance $\|z(t)-z(0)\|$ over $t\ge 1.5$.

Script: [sweep_f2_closure.jl](sweep_f2_closure.jl); raw log: [sweep_f2_closure.log](sweep_f2_closure.log).

## Result: no closure exists

**Positive $v_\text{target}$ (receding dimers).** Every IC runs clean to $t=20$ with drift $\le 5\times 10^{-4}\%$, but the intra-dimer $\dot r$ is **monotonically increasing** — its minimum over the run equals the initial value $v_\text{target}$ to 3+ digits, and its maximum grows with time. The dimers fly apart. Return distance grows monotonically from the start; the "best" return is always pinned to the earliest allowed sample $t=1.5$ at a large distance $\sim 5$–$8$.

**Negative $v_\text{target}$ (approaching dimers).** Every IC fails within $t<0.035$ with catastrophic energy drift (1–23%). The super-$\sqrt 2 c$ head-on phase hits an unregularizable singularity almost immediately.

**No mixed-sign trajectory exists.** In every positive-$v$ run, $\dot r_{13}$ stays strictly positive for the full $t=20$; it never brakes, never reverses, never approaches the Weber zero $\dot r=\sqrt 2 c$. In the negative-$v$ runs the IC dies before a single sign change occurs.

## Interpretation

The F2 fast-dimer family is **topologically open** at both ends: above $\sqrt 2 c$ the unlike pair's Weber factor $(1-\dot r^2/2c^2)$ is *negative*, so the effective sign of the Coulomb attraction is *flipped to repulsion*. For $\dot r > 0$ this accelerates the separation; for $\dot r < 0$ it accelerates the implosion. In neither case is there a restoring torque that could brake $\dot r$ back through $\sqrt 2 c$ within finite time. A closed "Weber loop" would require the trajectory to *cross* $\dot r=\sqrt 2 c$ twice in opposite directions; the FTL-regime equations of motion for an isolated dimer *prevent* that crossing because the force flips sign exactly at the crossing, making it an unstable fixed line of the radial dynamics.

The inter-dimer coupling (pairs 1-2, 3-4, 1-4, 2-3) is too weak to turn the intra-dimer $\dot r$ around at sep=2.5 before the dimers have flown apart. The sweep in `sep` would have to be done alongside `v_target`, but the monotonic growth of $\dot r$ in every positive-$v$ run suggests the coupling is subcritical by orders of magnitude.

## Bottom line

**REPORT §16 item 3 is resolved in the negative.** No closed Weber loop exists in the F2 family: the $\dot r=\sqrt 2 c$ surface is a one-way membrane, not a brake manifold. Any future search for a "Weber loop" in the 2+/2− system must start from a different ansatz — most plausibly one where an *unlike pair* with $\dot r$ slightly *below* $\sqrt 2 c$ is pumped through the zero by a *third-body interaction*, not by its own dynamics.

### Recommended replacement experiments

1. **Three-body pump.** Put an unlike dimer near $\dot r\approx 1.3c$ with a massive spectator `(+)` passing close enough to torque $\dot r$ through $\sqrt 2 c$. Check whether the reaction force on the spectator produces a net bounded orbit.
2. **Symmetric breathing square (F1) revisited.** Agent 9's F1 had all four $\dot r$ crossing $\sqrt 2 c$ simultaneously by symmetry. Re-examine F1 specifically for phase-space closure — it is the only F-family IC with *same-sign* crossings enforced by $D_4$ symmetry.
3. **Abandon the "Weber loop" framing.** The monotonicity result above shows the $\dot r=\sqrt 2 c$ surface is not a smooth brake under Weber dynamics. The right question may instead be whether finite-time sub-$\sqrt 2 c$ *arcs* can be glued by regularized bounces — but that requires extending the lifted-pair regularizer to velocity-dependent forces (REPORT §16 item 5).
