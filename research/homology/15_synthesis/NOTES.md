# Agent 15: Synthesis of 14-Agent Bound Orbit Study

## Overview

This document synthesizes results from Agents 01-14 of a 16-agent research study
on bound orbits in the Weber Hamiltonian. The study combined symplectic topology
(Rabinowitz-Floer homology, contact geometry, Conley-Zehnder indices), classical
perturbation theory (KAM, Nekhoroshev, Arnold diffusion), and large-scale numerical
computation (2062 two-body, 179 three-body, 66 four-body integrations, plus
variational orbit-finding and Poincare/frequency analysis).

Companion files in this directory:
- `master_orbit_catalog.csv` -- every orbit found, with metadata
- `existence_accounting.md` -- topology predictions vs numerical reality
- `open_conjectures.md` -- status of all conjectures (C1-C3, C11, C15a-e)

---

## 1. The Strongest Positive Results

### 1a. RFH theorem for 2-body unlike charges (Agent 02)

The single most important theoretical result of the study. For 2-body unlike charges
at any energy E < 0 and any c > 0, Rabinowitz-Floer homology is non-zero in
infinitely many degrees:

    RFH_k(Sigma_E) = H_{k+d}(Lambda S^d) != 0

This guarantees infinitely many geometrically distinct periodic orbits on every
negative-energy hypersurface. The proof chain (contact-type -> Moser regularization ->
c-continuation -> deformation invariance) is complete and rigorous for the 2-body case.

**Significance.** This upgrades the bound-orbit question for unlike charges from
"do they exist?" to "how many can we find?" The existence is now a theorem, not a
conjecture.

### 1b. First non-circular 2-body Weber periodic orbits (Agent 11)

Eight non-circular periodic orbits found by Newton shooting, with CZ indices 1-3 and
all elliptic (linearly stable). These are the first genuinely new periodic orbits of
the Weber Hamiltonian beyond the trivial circular family. They confirm the low-degree
RFH generators predicted by the theory.

### 1c. Three-body helium-like bound states (Agent 09)

The first systematic 3-body Weber survey found 27 bound orbits, with the 10 best
all in an asymmetric helium-like configuration: heavy nucleus (m=10, q=+2) with two
light satellites (m=1, q=-1). The best orbit (run 136) has energy drift effectively
zero and survives 10% random perturbations. This is the first evidence that multi-body
Weber systems can form robust bound states beyond the 2-body problem.

### 1d. Four-body triangular trap (Agent 10)

A new bound orbit family: three positive charges orbiting a central negative charge
(3+/1-) at high energy fraction. The R=4 configuration achieves energy drift ~1e-7%
(machine precision). This is only the second distinct bound geometry found for 4-body
Weber, after the 2+/2- breathing square.

### 1e. Contact-type verification (Agent 03)

Numerical confirmation that the 2-body unlike-charge Hill region is star-shaped at
all 50 negative energies tested (100% single-crossing). Combined with the analytical
proof (Agent 02, Prop 2.4), this establishes contact-type beyond doubt for 2-body.

### 1f. Complete FW2024 exegesis (Agent 01)

Definitive mapping between Frauenfelder-Weber's 2024 paper and the codebase.
Key clarifications: the "non-regularizability theorem" is NOT in FW2024 (it was
developed in the codebase); ell = 0 regularization is left explicitly open by FW;
no n-body theory exists in the literature.

---

## 2. The Strongest Negative Results

### 2a. No KAM tori in 4-body at c = 1 (Agent 12)

The most definitive negative result. All 138 tested ICs for 4-body 2+/2- showed
either escape or integrator failure. Frequency drift was 10-33% for the best
candidates (KAM requires < 0.01%). The "KAM basin" around the rhombus is not a
basin at all -- it is a region where dimer escape happens slowly.

### 2b. Double-orbiter is a transient, not bound (Agents 10, 12)

The 3D double-orbiter, previously the most promising 4-body candidate, is
definitively a long-lived transient. Agent 10 showed d_ratio grows to 32-76
(systematic expansion). Agent 12 revealed the ICs actually have E > 0 (unbound).
The "low Lyapunov exponent" was an artifact of near-free-flight.

### 2c. No pseudo-Riemannian Floer theory (Agent 08)

The subcritical (Lorentzian) region of the Weber metric lies genuinely outside all
current symplectic-topological machinery. The infinite Morse index of the Lorentzian
action prevents construction of a well-graded Floer chain complex. This is a
fundamental obstruction, not a technical gap.

### 2d. McGehee blow-up fails at critical radius (Agent 08)

The metric singularity at r = rho cannot be regularized by McGehee blow-up or any
standard technique. The singularity is in the kinetic energy (metric degeneracy),
not the potential, so Sundman/Levi-Civita absorption does not work.

### 2e. Nekhoroshev bounds are trivially short (Agent 07)

At c = 1, the Nekhoroshev stability time is T_N ~ exp(1) ~ 2.7 for the 4-body
system (6 effective DOF, exponent a = 1/12). Observed lifetimes (t* ~ 400-566)
exceed this by 200x, but this is in the fast Arnold diffusion regime, not the
Nekhoroshev regime. The perturbative framework is irrelevant at epsilon = 1.

### 2f. Elliptical orbits destabilize continuously (Agent 13)

No bifurcation or sharp transition: elliptical orbit instability grows continuously
as c decreases, with |lambda|_max scaling roughly as e^2/c^2. At c = 1, e = 0.7
orbits have |lambda| ~ 192. The Weber correction progressively destroys
integrability without creating new orbit families (for non-circular orbits).

---

## 3. Architecture of the Bound-Orbit Landscape

### 2-body unlike charges: fully understood

The bound-orbit structure is a continuous deformation of the Coulomb/Kepler problem.
Circular orbits are exact solutions (Weber vanishes), persisting at all c with
perfect stability. Non-circular orbits precess (quasi-periodic) and become
progressively unstable at low c, but remain energetically bound at all negative
energies. RFH guarantees infinitely many periodic orbits; 12 have been found
explicitly.

### 2-body like charges: theoretically classified, numerically inaccessible

FW2024's Theorem 2.1 completely classifies the radial motion. Subcritical bound
states exist for ell = 0 (bouncing) but are non-regularizable for ell != 0
(infinite winding). Supercritical bound states are unknown. The integrator cannot
access the subcritical regime (0/420 runs). This is a computational barrier, not
a physical one.

### 3-body: helium-like dominates

The only robust 3-body bound states are asymmetric helium-like configurations
with a heavy nucleus. This is physically intuitive: the heavy nucleus provides a
stable potential well, reducing the problem to an effective 2-body-plus-perturbation.
Equal-mass symmetric configurations are universally unstable. Same-sign triples
show zero bound states.

### 4-body: one periodic orbit, many transients

The breathing alternating square is the sole verified 4-body periodic orbit, but
it is violently unstable (CZ ~ 14, |lambda| ~ 229). The 3+/1- triangular trap is
a newly discovered bound family (3 orbits, good energy conservation) but has not
been verified as periodic. The 2+/2- double-orbiter and rhombus are long-lived
transients that eventually escape as dimers. No KAM tori exist at c = 1.

---

## 4. The Topology-Numerics Interface

### Where topology correctly predicted the numerics

1. **Infinitely many 2-body unlike orbits** (RFH, Agent 02) -- confirmed by 804
   bound orbits in the census (Agent 05) and 12 periodic orbits (Agent 11).

2. **Contact-type for 2-body unlike** (analytical, Agent 02; numerical, Agent 03) --
   100% single-crossing at all energies, consistent with the analytical proof.

3. **CZ index structure** (Agent 04) -- the found orbits have CZ 1-3, matching the
   predicted low-degree RFH generators. The breathing square at CZ 14 is a
   high-degree generator, consistent with its long period and high instability.

4. **c-continuation invariance** (Agent 02, Prop 3.2; Agent 13) -- circular orbits
   persist at all c, and RFH is constant along the continuation. The numerical
   continuation confirms no bifurcation.

### Where topology could not reach

1. **4-body contact-type** -- 90% single-crossing is suggestive but not sufficient.
   The topology-numerics gap is largest here.

2. **Like-charge dynamics** -- RFH is inapplicable (non-contact-type; Lorentzian
   metric). No topological tool exists for the subcritical regime.

3. **3-body existence** -- no RFH computation attempted. The helium-like bound states
   are a purely numerical discovery without topological backing.

4. **Low-CZ 4-body orbits** -- predicted by RFH (if C2 holds) but not found. This is
   the most important open computational challenge.

---

## 5. Literature Connections

Agent 06's survey of 28 Frauenfelder papers identified three priority research directions:

### Tier 1 (immediately actionable)

- **CZ-index numerical methods** (Frauenfelder-Koh-Moreno 2023): cell-mapping and
  CZ-tracking for Weber periodic orbit search. Agent 04 applied this successfully.
- **Doubly symmetric orbit constraints** (Frauenfelder-Moreno 2023): Weber's
  time-reversal symmetry constrains stability type. Applied implicitly by Agent 13.
- **Perturbation from Coulomb** (Moser 1970; FW2019): c-continuation. Implemented
  by Agents 11 and 13.

### Tier 2 (medium-term)

- **Convex embedding persistence** (Frauenfelder-van Koert-Zhao 2016): could prove
  4-body contact-type if the Weber perturbation preserves convexity.
- **Magnetic flow reduction** (Cieliebak-Frauenfelder-Paternain 2010): if Weber can
  be recast as a twisted cotangent bundle flow, the entire magnetic-flow toolkit
  applies. This is the most promising route to 3-body existence proofs.
- **Frozen planet variational approach** (Cieliebak-Frauenfelder-Volkov 2022): the
  nonlocal Levi-Civita regularization with multi-time-scale is directly relevant to
  the 3-body helium-like bound states found by Agent 09.

### Tier 3 (long-term foundational)

- **Pseudo-Riemannian Floer theory**: does not exist (Agent 08). Development would
  open the subcritical regime to topological methods.
- **Lagrangian RFH for noncompact surfaces** (Cieliebak-Frauenfelder-Miranda-
  Wisniewska 2024): handles the noncompactness near the critical radius.

---

## 6. The Biggest Obstruction

The **position-dependent kinetic metric** is the single thread connecting nearly every
open problem:

- It prevents Moser regularization near the critical radius (Agent 08).
- It prevents contact-type in the subcritical region (Agent 08).
- It prevents KAM/Nekhoroshev from giving useful bounds at c = 1 (Agent 07).
- It prevents standard Floer theory in the Lorentzian region (Agent 08).
- It makes the 4-body contact-type conjecture hard (Agent 03).
- It makes the integrator fail for like-charge subcritical orbits (Agent 05).

For unlike charges, the metric is positive definite everywhere (rho < 0), which is
why the 2-body theory is complete. For like charges, the metric degeneracy at r = rho
is the fundamental divide between what is accessible and what is not.

---

## 7. Summary Statistics

| Quantity | Value |
|----------|-------|
| Total integrations | 2307 (2062 + 179 + 66) |
| Total bound orbits found | 834 (804 2-body + 27 3-body + 3 4-body) |
| Verified periodic orbits | 13 (4 circular + 8 non-circular + 1 breathing square) |
| Like-charge bound orbits | 0 |
| KAM tori found | 0 |
| Papers surveyed | 28 (Agent 06) |
| Conjectures resolved | 1 (C15e confirmed: t* != c^2) |
| Conjectures still open | 4 (C1, C2/C11, C3) |
| New conjectures posed | 5 (C15a-e) |
| Strongest theorem | RFH != 0 for 2-body unlike at all E < 0, all c > 0 |
| Biggest gap | Low-CZ 4-body periodic orbits (predicted but not found) |
| Biggest obstruction | Position-dependent kinetic metric (Lorentzian signature) |

---

## References

Agents 01-14 NOTES.md files in `/research/homology/01_*/ through /14_*/`.
See individual files for detailed references to published literature.
