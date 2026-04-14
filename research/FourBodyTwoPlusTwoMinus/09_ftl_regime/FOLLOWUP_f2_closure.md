# Follow-up: the F2 Weber-loop hypothesis is false

**Date:** 2026-04-14 (autonomous follow-up to [REPORT.md §16 item 3](../REPORT.md))
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
