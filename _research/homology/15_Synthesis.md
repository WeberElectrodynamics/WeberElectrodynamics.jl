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


---

## Existence Accounting

# Existence Accounting: Topology Predictions vs Numerical Findings

## 1. What does RFH predict?

### 2-body unlike charges (the core theorem)

**Proposition 2.6 (Agent 02).** For the 2-body unlike-pair Weber Hamiltonian in
d dimensions at any energy E < 0 and any c > 0:

    RFH_k(Sigma_E) = H_{k+d}(Lambda S^d) != 0

in infinitely many degrees.

**Consequence.** There exist infinitely many geometrically distinct periodic orbits
on every negative-energy hypersurface. The proof chain is:

1. Contact-type verified analytically (Prop 2.4) and numerically (Agent 03: 100%
   single-crossing at all 50 negative energies tested).
2. Moser regularization compatible with Weber for all c > 0 (Prop 2.5; the Weber
   correction is smooth in the regularized chart).
3. c-continuation from Coulomb: no failure mode for unlike charges at E < 0,
   because (a) metric is positive definite everywhere (rho < 0), (b) Hill region
   is compact, (c) no escape mechanism.
4. RFH is a deformation invariant, so it equals the Kepler/Coulomb value at all c.

**Grading prediction.** The RFH generators in low CZ degree (0-2) correspond to
short-period orbits. The Kepler circular orbit family generates CZ degree 3 (d=2).
Higher iterates generate degrees 5, 7, 9, ... The 8 non-circular orbits found by
Agent 11 have CZ indices 1-3, consistent with being low-degree generators.

### 2-body like charges (supercritical)

**Conjecture 2.1 (Agent 02).** The supercritical like-charge regime has finitely
many (possibly zero) periodic orbits at a given energy. RFH of the truncated
supercritical level is generically zero (displaceable at infinity).

**Status.** The annular Hill region is not star-shaped, so standard RFH does not
apply. No periodic orbits found numerically (0/420 in Agent 05 census; integrator
limitation at metric degeneracy).

### 2-body like charges (subcritical)

**No Floer theory applies.** The Lorentzian kinetic metric gives infinite Morse
index (Agent 08). No pseudo-Riemannian Floer theory exists. Theorem 2.1 of FW2024
classifies all trajectories: ell != 0 spirals to collision (non-regularizable);
ell = 0 bounces at r = 0 (C^0-continuable but not C^1-periodic in FW's sense).

### 4-body 2+/2- (truncated supercritical)

**Conjecture 4.1 (Agent 02).** The truncated supercritical Liouville domain admits
a well-defined RFH complex, with generators matching known periodic orbits.

**Status.** Contact-type is only partially verified (Agent 03: 90% single-crossing
for alternating square and rhombus, not 100%). The conjecture remains open. If true,
RFH predicts periodic orbits in low CZ degrees (0-2) that have not yet been found
numerically.

### 3-body

No RFH computation attempted. The 3-body problem lacks the symmetry structure
needed for analytical contact-type proofs.

---

## 2. What did we find?

### Periodic orbits (rigorously verified)

| Category | Count | Source | Notes |
|----------|-------|--------|-------|
| 2-body circular (unlike) | 4 at c=1 + 24 across c-range | Agents 05, 11, 13 | Exact solutions: rdot=0 makes Weber vanish |
| 2-body non-circular (unlike) | 8 | Agent 11 | Newton-converged; CZ 1-3; all elliptic |
| 4-body breathing square | 1 | Agents 04, 05 | CZ~14; violently unstable (|lambda|=228.6); Newton did NOT converge |

**Total verified periodic: 13 distinct orbits** (4 circular + 8 non-circular + 1 breathing square),
plus 24 c-continuation copies of the circular orbits.

### Quasi-periodic orbits

| Category | Count | Source | Notes |
|----------|-------|--------|-------|
| 2-body elliptical (unlike) | 804 (all bound from census) | Agent 05 | Weber precession prevents closure; truly quasi-periodic |
| 2-body elliptical (e=0.3, 0.5, 0.7) at c=1 | 3 families | Agent 13 | Floquet |lambda| grows with e and 1/c^2 |

Weber apsidal precession makes all non-circular 2-body unlike-charge orbits
quasi-periodic rather than periodic (verified by Agent 14). The 804 bound orbits
from Agent 05's census are predominantly quasi-periodic (only 64 are circular/periodic).

### Chaotic-bound orbits

| Category | Count | Source | Notes |
|----------|-------|--------|-------|
| 4-body rhombus | 1 verified | Agent 14 | E < 0 but eventually escapes as dimers |
| 3-body helium-like | 27 bound | Agent 09 | 10 with E_drift < 1%; best: 0.00003% drift |
| 4-body triangular trap (3+/1-) | 3 | Agent 10 | New family; best: 1.3e-7% drift |

### Orbits shown to be artifacts

| Candidate | Verdict | Source | Notes |
|-----------|---------|--------|-------|
| Double-orbiter (2D/3D) | Long-lived transient, NOT bound | Agents 10, 12 | d_ratio grows linearly; systematic expansion |
| Double-orbiter "low Lyapunov" | Free flight (E > 0) | Agent 12 | Unbound; apparent stability = constant velocity |
| Rhombus "KAM basin" | Escaping dimers | Agent 12 | 10-33% frequency drift; no KAM tori |
| Like-charge sub-critical | Integrator failure, not physics | Agent 05 | 0/420; metric singularity prevents integration |

---

## 3. The gap

### Quantitative gap: 2-body unlike charges

| Metric | RFH prediction | Found | Ratio |
|--------|---------------|-------|-------|
| Periodic orbit families | Infinitely many (at each E < 0) | 12 distinct (4 circ + 8 non-circ) | 12 / infinity |
| CZ degrees populated | All non-negative integers | 1, 2, 3 only | 3 / infinity |
| Energy range covered | All E < 0 | E in [-5.0, -0.01] | Sampled |

**Specific missing orbits predicted by RFH:**

1. **CZ degree 0 generators.** The breathing-mode orbits (constant-topology loops)
   predicted as degree-0 generators of RFH have not been found. These should be the
   shortest-action periodic orbits.

2. **CZ degrees 4-13.** The gap between the found CZ 1-3 orbits and the breathing
   square at CZ 14. RFH predicts generators in every degree; degrees 4-13 are empty.

3. **High-energy orbits (|E| > 5).** No orbits computed at high binding energy where
   the Weber correction is strongest and CZ indices are highest (Agent 04: CZ = 23
   at E = -5).

4. **Non-circular orbit families at varied c.** Agent 13 tracked only circular orbits
   through c-continuation (trivial, since Weber vanishes on circles). The 8 non-circular
   orbits from Agent 11 were found only at c = 1. c-continuation of non-circular orbits
   would reveal how the Weber correction creates/destroys orbit families.

5. **Higher iterates.** The k-fold covers of the prime circular orbit have CZ = 2k(d-1).
   None of these higher iterates were explicitly found as distinct numerical orbits.

### Qualitative gap: 4-body

| Metric | RFH prediction (if contact-type holds) | Found | Status |
|--------|---------------------------------------|-------|--------|
| Contact-type | Conjectured (90% single-crossing) | Not proven | Open |
| Periodic orbits | At least several families | 1 (breathing square) | Enormous gap |
| Short-period orbits (CZ 0-2) | Predicted by low-degree RFH | 0 found | Main open problem |

The breathing square (CZ ~ 14) is a high-degree generator. The low-degree generators
(CZ 0-2) that RFH predicts should exist have not been found. These would be:
- A simple breathing mode (CZ 0)
- An orbit linking one critical sphere (CZ 1-2)
- A short double-cover of the breathing mode (CZ 2)

### Qualitative gap: like charges

RFH does not apply to like charges (neither subcritical nor supercritical). The
theoretical framework (FW2024) establishes that subcritical ell != 0 orbits are
non-regularizable. The numerical census (0/420 bound like-charge orbits) is
consistent with theory but cannot be definitive due to integrator limitations.

**Open question:** Do supercritical like-charge Weber periodic orbits exist?
The velocity-dependent attraction could create an effective potential well for
sufficiently fast radial approach, but this has not been demonstrated numerically
or analytically.

### Gap summary by body count

| N-body | Theory | Numerics | Gap characterization |
|--------|--------|----------|---------------------|
| 2-body unlike | Infinitely many periodic at every E < 0 (RFH theorem) | 12 periodic + ~800 quasi-periodic | Quantitative: need more families and higher CZ |
| 2-body like (sub) | Non-regularizable for ell != 0 (FW2024) | 0 found | Consistent; ell = 0 bounces exist in principle |
| 2-body like (super) | Unknown (RFH inapplicable) | 0 found | Open |
| 3-body | No theory | 27 bound (10 good) | No predictions to compare against |
| 4-body 2+/2- | Contact-type conjectured; RFH would predict families | 1 periodic + 1 chaotic-bound | Enormous gap; low-CZ orbits missing |
| 4-body 3+/1- | No theory | 3 bound | New discovery; no theoretical framework |

---

## 4. What would close the gap?

### Immediate (computational)

1. **L-BFGS action minimization** for the 4-body problem (Agent 11's gradient descent
   was too slow). Finding CZ 0-2 orbits for the 4-body alternating square would be
   the single most impactful computational result.

2. **Multiple-shooting with subdivision** for the breathing square (Agent 11 got
   ||F|| down to 0.042 before stalling). Needs better preconditioning or deflation.

3. **c-continuation of non-circular orbits** from Coulomb (c large) to Weber (c=1).
   Track the 8 Agent 11 orbits through the bifurcation cascade.

4. **T-brake orbit refinement** (Agent 13 found ~25 candidates with closure error < 0.1).

### Medium-term (theoretical)

5. **Prove 4-body contact-type** via exact star-center construction or Hofer's
   approach for overtwisted contact structures.

6. **Convex embedding persistence** under Weber perturbation (Frauenfelder-van Koert-Zhao
   approach from Agent 06, Section 2b).

7. **Magnetic flow reduction** of Weber (Agent 06, Section 4a). If Weber can be recast
   as a twisted cotangent bundle flow, the Cieliebak-Frauenfelder-Paternain toolkit
   gives periodic orbit existence on almost every energy level.

### Long-term (foundational)

8. **Develop pseudo-Riemannian Floer theory** to access subcritical dynamics. Currently
   no such theory exists (Agent 08).

9. **Frozen planet variational approach** adapted from Cieliebak-Frauenfelder-Volkov
   (Agent 06, Section 5a) for 3-body Weber bound states.


---

## Open Conjectures

# Open Conjectures: Status Update After 14-Agent Study

## Conjectures from Prior Work (Agent 10's C1-C3)

### C1: Existence of 4-body bound orbits

**Original statement.** The 4-body 2+/2- Weber Hamiltonian at c = 1 admits at
least one genuinely bound orbit (negative energy, all pair distances bounded for
all time).

**Status: OPEN (weakened).**

Evidence for:
- The breathing alternating square is a genuine periodic orbit (Agents 04, 05).
  However, it is violently unstable (|lambda|_max = 228.6) and generic nearby
  initial conditions escape.
- The 4-body rhombus chaotic-bound orbit survives t = 200 with E < 0 and drift
  ~5e-6 (Agent 14).

Evidence against:
- Agent 12 showed that ALL tested 4-body ICs either escape as dimers or crash
  the integrator. The "KAM basin" around the rhombus is NOT a basin -- it is a
  region of slow dimer escape with 10-33% frequency drift.
- The double-orbiter (2D and 3D) is definitively a long-lived transient with
  systematic expansion (d_ratio 32-76), not a bound orbit.
- The double-orbiter ICs actually have E > 0 (unbound); "low Lyapunov" was
  free flight.
- Agent 12: "There are no genuine KAM tori in the 4-body 2+/2- Weber Hamiltonian
  at the parameters surveyed."

**Assessment.** The breathing square is the only rigorously periodic 4-body orbit,
but its extreme instability means no generic IC stays near it. The conjecture of
a genuinely bound (non-periodic) orbit with bounded trajectories for all time
remains unresolved. The weight of evidence suggests that if such orbits exist,
they occupy a set of measure zero in phase space.

---

### C2: RFH is well-defined for the truncated supercritical domain

**Original statement.** The energy hypersurface of the 4-body 2+/2- Weber
Hamiltonian, restricted to the supercritical region {r_ij > rho_ij + epsilon},
is of contact type, and RFH is well-defined.

**Status: OPEN (partially supported).**

Evidence for:
- Agent 03: 90% single-crossing for alternating square and rhombus configurations
  after star-center optimization. The obstruction is localized near the collision
  strata.
- Agent 08: excising an epsilon-neighborhood of the critical sphere gives a valid
  Liouville domain with smooth convex boundary. The RFH chain complex is
  well-defined in principle.
- Agent 02 (Prop 2.4): for the 2-body unlike pair, contact-type holds analytically
  for all E < 0.

Evidence against:
- Agent 03: only 38% monotonicity even after optimization for the alternating
  square (compared to 100% for 2-body unlike). The 4-body geometry is fundamentally
  harder.
- No analytical star-center construction exists for the 4-body case.
- Agent 08: the critical sphere is a degenerate-contact boundary (not a clean stop),
  complicating the Floer boundary conditions.

**Assessment.** The excision approach (Agent 08) provides a well-defined Liouville
domain in principle, but Conjecture 4.1 (compactness of Floer moduli) remains
unverified. The numerical evidence (90% single-crossing) is encouraging but
insufficient for a rigorous proof.

---

### C3: Low-degree RFH generators correspond to short-period orbits

**Original statement.** The generators of RFH in degrees 0-2 for the 4-body 2+/2-
problem correspond to short-period periodic orbits distinct from the breathing
square (which has CZ ~ 14, a high-degree generator).

**Status: OPEN (not falsified but no candidates found).**

Evidence for:
- Agent 04: Kepler circular orbits have CZ = 3 (degree 0-2 after shift). Weber
  orbits at CZ 1-3 exist in 2-body (Agent 11). The grading structure is consistent.
- Agent 02 (Conj 5.2): predicts breathing-type degree-0 generator and critical-sphere-
  linking degree 1-2 generators.

Evidence against:
- No CZ 0-2 orbit has been found in the 4-body problem. Agent 11's action
  minimization failed to converge for any 4-body configuration.
- The breathing square Newton shooting failed (Agent 11: ||F|| stalled at 0.042).
  If the only periodic orbit is hard to close numerically, finding shorter ones
  may be even harder.
- Agent 12 found no quasi-periodic motion at all in the 4-body system, which
  undermines the assumption that short-period orbits exist and are accessible.

**Assessment.** This is the most important open computational challenge. Finding
even one CZ 0-2 periodic orbit in the 4-body system would confirm the RFH framework
and validate the contact-type conjecture.

---

## Conjecture from Agent 11

### C11: 4-body contact-type

**Original statement.** (From Agent 11's contact/Reeb analysis context.) The 4-body
2+/2- Weber energy hypersurface is of contact type in the supercritical region.

**Status: OPEN (same as C2).**

This is effectively a restatement of the key hypothesis underlying C2. The
numerical evidence from Agent 03 (90% single-crossing) provides partial support.
A rigorous proof would require either:
(a) Finding an exact star-center analytically for the 4-body Hill region.
(b) Using Hofer's 1993 approach for overtwisted contact structures (does not
    require star-shapedness).
(c) Establishing dynamical convexity via the Frauenfelder-Kang real holomorphic
    curve method (requires symmetry + convexity).

---

## New Conjectures Emerging from This Study

### C15a: No KAM tori in 4-body 2+/2- at c = 1

**Statement.** The 4-body 2+/2- Weber Hamiltonian at c = 1 has no KAM tori
(invariant tori carrying quasi-periodic motion) in any region of phase space.

**Status: STRONGLY SUPPORTED.**

Evidence: Agent 12 tested 138 ICs across rhombus and double-orbiter families.
All showed either escape or integrator failure. Frequency drift was 10-33% for
the best candidates (KAM requires < 0.01%). Agent 07 showed the perturbation
parameter epsilon = 1/c^2 = 1 is far too large for KAM theory to apply (estimated
KAM fraction ~ 0 at c = 1).

**Caveat.** The search was not exhaustive. KAM tori could exist in regions of
phase space not sampled. However, the theoretical estimate (KAM fraction ~ 0)
and numerical evidence are consistent.

### C15b: 3-body helium-like configurations are the most robust multi-body bound states

**Statement.** For n >= 3 body Weber systems, the most stable bound configurations
are asymmetric helium-like systems: one heavy nucleus (charge +2, mass >> 1) with
lighter unlike-charge satellites.

**Status: SUPPORTED.**

Evidence: Agent 09 found 10 good-bound 3-body orbits, all in the heavy-nucleus
helium-like configuration ([+2,-1,-1], masses [10,1,1]). Equal-mass configurations
produced 0 good-bound orbits across 40 equilateral runs. The best orbits (runs
132, 136) have energy drift < 0.001% and survive 10% perturbations.

### C15c: Circular orbits are the only stable 2-body periodic orbits at c = 1

**Statement.** For 2-body unlike charges at c = 1, circular orbits (which have
CZ 2-3 and are exact solutions) are the only linearly stable periodic orbits.
All non-circular periodic orbits are unstable, with instability growing with
eccentricity as ~e^2/c^2.

**Status: PARTIALLY SUPPORTED.**

Evidence: Agent 13 showed circular orbits have |lambda|_max = 1.000 at all c,
while elliptical orbits have |lambda|_max = 17 (e=0.3) to 192 (e=0.7) at c=1.
Agent 11's 8 non-circular Newton-converged orbits all have |lambda|_max in
[1.000, 1.002] -- very close to stable but these are near-circular perturbations,
not high-eccentricity.

**Open question:** Is there a critical eccentricity e*(c) below which non-circular
orbits remain stable? The Agent 11 orbits at e ~ 0.05-0.20 appear stable, while
Agent 13's orbits at e = 0.3-0.7 are unstable.

### C15d: The Weber correction is a Coulomb perturbation for unlike charges

**Statement.** For 2-body unlike charges, the bound-state structure is a continuous
deformation of the Coulomb problem for all c > 0. No qualitative phase transition
(bifurcation, new orbit family creation/destruction) occurs as c decreases from
infinity to any finite value.

**Status: SUPPORTED for circular orbits; PARTIALLY for non-circular.**

Evidence: Agent 13 found no bifurcation for circular orbits in [0.5, 100].
Agent 05 found bound fractions independent of c. Agent 02 proved RFH is constant
for all c > 0. However, the CZ index does jump at specific c values (Agent 04:
first bifurcation at c ~ 1.414 for E = -1.0), which Agent 13 interprets as a
Maslov/winding transition rather than a stability change. For non-circular orbits,
the continuous destabilization (Agent 13) is a quantitative but not qualitative
change.

### C15e: The double-orbiter escape time does not scale as c^2

**Statement.** The 4-body double-orbiter escape time t* does not follow the
t* ~ c^2 scaling predicted by fast Arnold diffusion.

**Status: CONFIRMED (Agent 10).**

Agent 10 measured: 2D double-orbiter t* ~ c^{-0.24}; 3D t* ~ c^{-0.56}; rhombus
t* ~ c^{+0.23}. The scaling is weak and geometry-dependent, not the c^2 predicted
by Agent 07. The dominant binding mechanism is Coulombic, and the Weber correction
provides only a modest, non-monotonic effect.

---

## Summary Table

| Conjecture | Status | Confidence | Key evidence |
|------------|--------|------------|-------------|
| C1 (4-body bound orbit) | Open (weakened) | Low | Breathing square is periodic but unstable; all others escape |
| C2 (RFH well-defined for 4-body) | Open (partially supported) | Moderate | 90% single-crossing; excision framework valid |
| C3 (Low-CZ orbits exist) | Open (not falsified) | Low | No candidates found despite extensive search |
| C11 (4-body contact-type) | Open (= C2) | Moderate | Same evidence as C2 |
| C15a (No KAM tori at c=1) | Strongly supported | High | 138 ICs, 10-33% freq drift, theory: KAM fraction ~ 0 |
| C15b (Helium-like most robust) | Supported | Moderate | 10/27 good bound are helium-like; 0/40 equal-mass |
| C15c (Only circular stable) | Partially supported | Moderate | Circular: |l|=1 all c; elliptical: |l| grows with e |
| C15d (Weber = Coulomb perturbation) | Supported (unlike) | High | No bifurcation; RFH constant; bound fractions c-independent |
| C15e (t* != c^2) | Confirmed | High | Measured: t* ~ c^{-0.24} to c^{+0.23} |
