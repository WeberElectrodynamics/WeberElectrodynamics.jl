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
