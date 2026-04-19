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
