# 04 — Conley-Zehnder Indices for Weber Orbits

Agent 04 deliverable. Script: `cz_index.jl`. Data: `cz_table.csv`.

## 1. Background and conventions

The Conley-Zehnder index `mu_CZ(gamma)` is an integer assigned to each
non-degenerate closed Reeb orbit gamma that grades the Floer / RFH chain
complex. It counts (roughly) the number of half-rotations of the
linearized flow around gamma in the contact distribution.

**Convention.** We use the lower-semicontinuous CZ index of
Salamon-Zehnder 1992 / Robbin-Salamon 1993. For a symplectic path
Psi: [0,1] -> Sp(2n) with Psi(0) = I and Psi(1) = M (the monodromy),
mu_CZ is computed from the symplectic normal form of M:

- **Elliptic block** (eigenvalues `e^{+/- i theta}` on the unit circle,
  theta in (0, pi)): contribution = `2 * floor(theta / pi) + 1`.
  For the total rotation angle Theta = omega * T accumulated over the
  orbit: contribution = `2 * floor(Theta / pi) + 1`.
- **Positive hyperbolic** (real eigenvalues lambda, 1/lambda with
  lambda > 1): contribution = 0.
- **Negative hyperbolic** (real eigenvalues -lambda, -1/lambda with
  lambda > 1): contribution = 1.
- **Loxodromic** (complex eigenvalues off the unit circle): contribution
  = 0 per conjugate pair.
- **Trivial** (lambda = 1): degenerate; removed when computing the
  reduced CZ index (one pair per conserved quantity).

## 2. CZ index from Floquet multipliers

The monodromy matrix M = d(phi_T)/dz evaluated at a periodic orbit of
period T has Floquet multipliers as eigenvalues. For a Hamiltonian system
with n DOF, M is 2n x 2n symplectic. Conserved quantities (energy,
linear momentum, angular momentum) each contribute a trivial pair
(lambda = 1, 1). The reduced monodromy on the symplectically reduced
phase space gives the non-trivial CZ contributions.

**Implementation.** `cz_index.jl` provides:
- `conley_zehnder_index(M)` — from a full monodromy matrix
- `conley_zehnder_from_eigenvalues(evals)` — from Floquet multipliers
- `classify_symplectic_eigenvalue(lambda)` — block type classifier
- `cz_block_contribution(type, theta)` — per-block CZ formula

Sanity checks verified:
- Identity matrix -> mu_CZ = 0
- Rotation by theta < pi -> mu_CZ = 1
- Positive hyperbolic -> mu_CZ = 0
- Negative hyperbolic -> mu_CZ = 1
- Block-diagonal elliptic + hyperbolic -> correct sum

## 3. CZ index for the breathing alternating square

**Data from Agent 5.** The only periodic orbit found numerically:
- Period T = 11.78, Energy E = -0.646, max|lambda| = 228.6
- 4 particles, 2D, D_4 x T brake symmetry, L = 0
- Full phase space: 16-dim (8 DOF)
- Conserved: energy (1) + linear momentum (2) + angular momentum (1) = 4
- Reduced phase space: 8-dim -> 4 symplectic pairs

**Reconstructed reduced spectrum** (from Agent 4's equilibrium
eigenstructure, continued to the breathing orbit):

| Block | Type | Source | theta or log|lambda| | CZ contribution |
|-------|------|--------|----------------------|-----------------|
| 1 | pos_hyperbolic | dominant instability | log(228.6) = 5.43 | 0 |
| 2 | elliptic | omega = 0.6762 mode | theta = 7.97, theta/pi = 2.54 | 5 |
| 3 | elliptic | omega = 1.2290 mode | theta = 14.48, theta/pi = 4.61 | 9 |
| 4 | uncertain (Krein) | Krein quadruple | omega ~ 0.84, theta/pi ~ 3.15 | 0 or 7 |

**Result:**
- **mu_CZ = 14** (conservative: Krein pair loxodromic, contributing 0)
- **mu_CZ = 21** (upper bound: Krein pair elliptic, contributing 7)

The orbit contributes to Floer grading 14 (or 21). In either case this
is a high-index generator of the RFH complex, far from the lowest-action
generators in degrees 0-2.

## 4. CZ index for Kepler circular orbits (c -> infinity)

In the Kepler limit (c = infinity), all 2-body unlike-charge circular
orbits have mu_CZ = 3 (for the simple orbit). The n-th iterate has
mu_CZ = 2n + 1.

| orbit | n_iterate | mu_CZ |
|-------|-----------|-------|
| simple circular | 1 | 3 |
| 2nd iterate | 2 | 5 |
| 3rd iterate | 3 | 7 |

These are the "seed" indices for c-continuation to finite Weber coupling.

## 5. c-continuation of CZ index

As c decreases from infinity, the Weber correction introduces apsidal
precession. The rotation number changes as:

    nu(c) = 1 + k / (2 mu c^2 a)

where k = |q1 q2| and a is the orbital radius. The CZ index is locally
constant and jumps at bifurcation values of c where nu crosses an
integer:

    c_bif(n) = sqrt(k / (2 mu a (n - 1)))

**Key results at c = 1 (standard Weber):**

| Energy | Radius a | nu(c=1) | mu_CZ | Change from Kepler |
|--------|----------|---------|-------|--------------------|
| -0.25 | 2.0 | 1.5 | 3 | none |
| -0.50 | 1.0 | 2.0 | 5 | +2 (one bifurcation) |
| -1.00 | 0.5 | 3.0 | 7 | +4 (two bifurcations) |
| -2.00 | 0.25 | 5.0 | 11 | +8 (four bifurcations) |
| -5.00 | 0.1 | 11.0 | 23 | +20 (ten bifurcations) |

**Physical interpretation.** Tighter orbits (smaller a, more negative E)
experience stronger Weber correction and undergo more CZ index jumps.
Each jump corresponds to an additional half-rotation of the linearized
flow per period — physically, the Weber apsidal precession causes the
eccentricity vector to wind more rapidly.

The first bifurcation (CZ: 3 -> 5) occurs at:
- c = 0.707 for E = -0.25 (a = 2)
- c = 1.0 for E = -0.5 (a = 1)
- c = 1.414 for E = -1.0 (a = 0.5)

So for c = 1, orbits with E < -0.5 have already undergone at least one
CZ bifurcation from the Kepler value.

## 6. Predicted RFH generators

From Agent 10's conjectures C2-C3, the RFH complex should contain:

1. **Degree-0 generators**: short breathing orbits in each connected
   component of the loop space. These continue from Kepler circular
   orbits with mu_CZ = 3.

2. **Degree-2 generators**: linking-class-detecting orbits corresponding
   to the extra cohomology from excising the Weber critical spheres.
   Predicted mu_CZ = 5-7.

3. **Higher generators**: infinitely many, from the infinite total Betti
   number of the loop space Lambda(M) (Sullivan-Vigue-Poirrier).

The breathing square at mu_CZ = 14 (or 21) is a high-index generator,
consistent with being a long-period, highly unstable orbit that is not
among the shortest-action RFH generators.

## 7. Connection to other agents

- **Agent 4** (equilibria): eigenstructure used to reconstruct the
  breathing square's reduced monodromy spectrum.
- **Agent 5** (periodic orbits): Floquet data (max|lambda| = 228.6,
  T = 11.78) is the numerical input for the CZ computation.
- **Agent 10** (Floer/symplectic): RFH conjectures C1-C3 predict which
  CZ degrees should be populated.
- **Agent 11** (contact/Reeb): Weinstein contact-type criterion restricts
  CZ computations to the supercritical region.
- **Agent 13** (c-continuation): the c-bifurcation values predicted here
  are the loci where period-doubling or new orbit families should appear.

## 8. Caveats

1. The breathing square CZ is **estimated** from the equilibrium
   eigenstructure, not from a full numerical monodromy. The Krein pair
   fate (loxodromic vs elliptic) is unresolved without the full 16x16
   monodromy matrix.
2. The Weber apsidal precession formula is leading-order in 1/c^2 and
   breaks down when the orbit approaches the critical radius
   rho = q1 q2 / (mu c^2).
3. CZ indices at exact integer rotation numbers (nu = 2, 3, ...) are
   degenerate and require perturbation to resolve. The values in the
   table assume lower-semicontinuous convention (floor).
4. All results are in the supercritical region; CZ is undefined across
   the Weber critical spheres (Agent 11 section 3).

## Files

- `cz_index.jl` — CZ index calculator with sanity checks and table generator
- `cz_table.csv` — Index table for all orbits (14 rows)
