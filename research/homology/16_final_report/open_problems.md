# Open Problems in Weber Hamiltonian Bound Orbit Theory

Prioritized list emerging from the 16-agent RFH study (2026-04-16).

---

## Tier 1: Immediate computational targets (high payoff, feasible now)

### OP-1. Find low-CZ 4-body periodic orbits

**Priority: HIGHEST.** RFH predicts generators in CZ degrees 0-2 for the 4-body 2+/2- system (contingent on contact-type, Conjecture C2). The breathing alternating square at CZ ~14 is a high-degree generator; the low-degree generators -- which should be the shortest-action orbits -- have not been found despite extensive search.

**Recommended approach:** L-BFGS action minimization on a discretized loop space (Agent 11's gradient descent was too slow). Alternatively, deflated Newton methods to avoid reconverging to the breathing square. Multiple-shooting with subdivision may also succeed where single-shooting failed (Agent 11 stalled at ||F|| = 0.042).

**Why it matters:** Finding even one CZ 0-2 orbit would simultaneously (a) confirm the RFH framework for 4-body, (b) validate the contact-type conjecture, and (c) produce the first short-period 4-body Weber periodic orbit.

### OP-2. c-continuation of non-circular 2-body orbits

**Priority: HIGH.** Agent 11 found 8 non-circular periodic orbits at c = 1. Tracking them as c varies from 100 (near-Coulomb) to 0.5 (strong Weber) would reveal:
- Whether Weber creates new orbit families via bifurcation
- The critical eccentricity e*(c) separating stable from unstable non-circular orbits
- How CZ indices evolve under the continuation

Agent 13 tracked only circular orbits (trivial, since Weber vanishes on circles). The non-circular continuation is the natural next step.

### OP-3. T-brake orbit refinement

**Priority: MODERATE.** Agent 13 identified ~25 T-brake candidates (p = 0 initial conditions) with closure error < 0.1 but did not refine them with Newton iteration or compute Floquet multipliers. These are low-hanging fruit: each converged orbit adds to the periodic orbit census and populates new CZ degrees.

### OP-4. Verify 3-body helium-like periodicity

**Priority: MODERATE.** The 27 helium-like 3-body bound orbits (Agent 09) have not been tested for periodicity via Newton shooting. The best candidates (runs 132, 136) have near-zero energy drift and survive perturbations. If any are periodic (or near-periodic with small closure error), they would be the first 3-body Weber periodic orbits.

---

## Tier 2: Medium-term theoretical problems

### OP-5. Prove 4-body contact-type (Conjectures C2/C11)

**Status:** 90% single-crossing numerically (Agent 03); no analytical proof. The obstruction is near the collision strata, where the star-shapedness test fails.

**Approaches:**
- (a) Exact star-center construction for the 4-body Hill region
- (b) Hofer's 1993 approach for overtwisted contact structures (does not require star-shapedness)
- (c) Dynamical convexity via Frauenfelder-Kang real holomorphic curve method (requires symmetry + convexity)

### OP-6. Magnetic flow reduction of Weber

If the Weber Hamiltonian can be recast as a twisted cotangent bundle flow (magnetic Hamiltonian), then the Cieliebak-Frauenfelder-Paternain toolkit from "Symplectic topology of Mane's critical values" gives periodic orbit existence on almost every energy level. The Weber velocity-dependent correction has the structure of a magnetic term, making this a natural (but non-trivial) identification.

### OP-7. Frozen planet variational approach for 3-body

Adapt the Cieliebak-Frauenfelder-Volkov (2022) nonlocal Levi-Civita regularization with multi-timescale to the helium-like 3-body Weber bound states. The mass asymmetry (m_nucleus = 10, m_electron = 1) provides a natural separation of scales that the frozen-planet method exploits.

### OP-8. Convex embedding persistence

Apply Frauenfelder-van Koert-Zhao (2016) convex embedding theory to determine whether the Coulomb/Kepler convex embedding persists under the Weber perturbation. If so, this would prove contact-type for the 2-body problem by a second independent method, and potentially extend to restricted 3-body problems.

---

## Tier 3: Long-term foundational problems

### OP-9. Pseudo-Riemannian Floer theory

**The deepest open problem.** The subcritical region ($r < \rho$ for like charges) has a Lorentzian kinetic metric, giving infinite Morse index for the action functional. No well-graded Floer chain complex can be constructed with current tools (Agent 08). Developing a pseudo-Riemannian Floer theory would:
- Open the subcritical regime to topological methods
- Potentially classify like-charge bound states
- Have impact far beyond Weber electrodynamics (Lorentzian geometry, general relativity)

### OP-10. Like-charge regularization backend

The current integrator cannot handle the metric signature change at $r = \rho$, producing 0/420 bound orbits for like charges. A new regularization backend is needed. Key constraint: only $\ell = 0$ (head-on) passages can succeed; $\ell \neq 0$ spirals are non-regularizable (Frauenfelder-Weber 2024, Theorem 2.1). Possible approaches:
- Adaptive coordinate switching at $r = \rho + \epsilon$
- Implicit integration through the degenerate metric
- Symplectic splitting with separate treatment of the radial DOF

### OP-11. n-body RFH extension

Extend the 2-body RFH theorem (Prop 2.6, Agent 02) to 3+ bodies. Obstacles:
- No simple symmetry reduction for 3+ bodies
- Contact-type proofs require understanding the full Hill region geometry
- Collision strata become higher-codimension
- The free loop space homology $H_*(\Lambda M)$ is more complex

The 3-body helium-like configuration (Agent 09) is the most tractable target, due to the mass-ratio hierarchy.

### OP-12. Supercritical like-charge periodic orbits

Do supercritical ($r > \rho$) like-charge Weber periodic orbits exist? The annular Hill region is non-compact and not star-shaped, so RFH does not directly apply. The velocity-dependent attraction could create an effective potential well for sufficiently fast radial approach, but this has not been demonstrated. Agent 02 conjectures finitely many (possibly zero) periodic orbits at a given energy.

---

## Conjecture Status Summary

| ID | Statement | Status | Confidence |
|----|-----------|--------|------------|
| C1 | 4-body 2+/2- admits a genuinely bound orbit | Open (weakened) | Low |
| C2 | RFH well-defined for truncated supercritical 4-body | Open (partially supported) | Moderate |
| C3 | Low-CZ generators = short-period 4-body orbits | Open (not falsified) | Low |
| C11 | 4-body contact-type in supercritical region | Open (= C2) | Moderate |
| C15a | No KAM tori in 4-body at c = 1 | Strongly supported | High |
| C15b | Helium-like = most robust multi-body bound states | Supported | Moderate |
| C15c | Only circular orbits stable at c = 1 | Partially supported | Moderate |
| C15d | Weber = Coulomb perturbation for unlike charges | Supported (unlike) | High |
| C15e | $t^* \neq c^2$ scaling | **Confirmed** | High |
