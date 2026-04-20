# Exegesis of Frauenfelder-Weber 2024: "A mathematical description of the Weber nucleus"

Agent 01 deliverable. Theorem-by-theorem analysis of FW2024 with
codebase mapping. Statements tagged **[thm]** (cited theorem/proposition),
**[obs]** (observation from this report), **[conj]** (conjecture from this report).

Reference: Frauenfelder, U., Weber, J. "A mathematical description of the
Weber nucleus as a classical and quantum mechanical system."
*Anal. Math. Phys.* **14**:31 (2024).
DOI: [10.1007/s13324-024-00891-5](https://doi.org/10.1007/s13324-024-00891-5).

---

## 0. Paper structure

| Section | Content |
|---------|---------|
| 1 | Introduction: Weber electrodynamics history, Lagrangian, metric rewriting, Theorem A statement |
| 2 | Classical motion: energy equation, complete trajectory classification (Theorem 2.1) |
| 3 | Derivation of the Weber-Schrodinger equation: Laplace-Beltrami on Weber plane, separation of variables |
| 4 | Inside critical radius: Sturm-Liouville analysis, Propositions 4.2 and 4.4 proving Theorem A |

The paper is 26 pages. It has exactly **one main theorem** (Theorem A), **one
classification theorem** (Theorem 2.1), **two propositions** (4.2, 4.4), **two
remarks** (4.1, 4.3, 4.5), and **one lemma** (3.1). There are no corollaries.

---

## 1. The Weber Lagrangian and metric rewriting

### 1.1 Lagrangian (Eq. 1.1)

FW2024 uses **normalized units**: two equal positive charges, reduced mass
mu = 1, charge product k = ee'/mu = 1. The Weber Lagrangian in polar
coordinates (r, phi) for relative motion is:

```
L_W(r, phi, v_r, v_phi) = (1/2)(v_r^2 + r^2 v_phi^2) - (1/r)(1 + v_r^2/(2c^2))
```

**[obs]** The sign convention: the potential term is `-1/r (1 + v_r^2/(2c^2))`.
Note the `+` inside the parentheses. This is for **like charges** (ee' > 0),
where the static Coulomb potential is `+1/r` (repulsive). The velocity-dependent
correction `v_r^2/(2c^2)` modifies the effective potential. The sign is opposite
to what one might expect from the standard Weber potential
`U = (ee'/r)(1 - rdot^2/(2c^2))` because FW absorb the static part into the
Lagrangian differently.

### 1.2 Metric rewriting (Eq. 1.1, second form)

By collecting `v_r^2` terms:

```
L_W = (1/2)(g_rr v_r^2 + g_phiphi v_phi^2) - 1/r
```

where:

```
g_rr = (r - r_c)/r,    g_phiphi = r^2,    g_rphi = 0
```

and the **critical (Weber) radius** is:

```
r_c := 1/c^2
```

**[obs]** In the codebase's physical units (CriticalRadiusAndLikeChargeAttraction.md),
rho = ee'/(mu c^2). With FW's normalization k = ee'/mu = 1, this gives
rho = 1/c^2 = r_c. The notations are consistent.

**[obs]** The metric tensor `(g_ij)` defines the **Weber plane**
`(R^2_x, g)` where `R^2_x = R^2 \ {0, x^2 + y^2 = r_c^2}`.

- For `r > r_c`: `g_rr > 0`, metric is **Riemannian** (positive definite)
- At `r = r_c`: `g_rr = 0`, metric is **degenerate**
- For `r < r_c`: `g_rr < 0`, metric is **Lorentzian** (indefinite)

**[obs]** This is exactly the metric signature flip documented in
`CriticalRadiusAndLikeChargeAttraction.md` Section 7, where it is noted
that `g_rr = 1 - rho/r` equals `mu_eff/mu`. The codebase derivation and
FW2024 are completely aligned.

### 1.3 Conjugate momenta (Section 2)

```
p_r = (r - r_c)/r * rdot,    p_phi = r^2 * phidot
```

Angular momentum: `ell := p_phi = r^2 phidot` is conserved.

### 1.4 Weber Hamiltonian (Eq. 2.2)

```
H(r, phi, p_r, p_phi) = (1/2)(r p_r^2/(r - r_c) + p_phi^2/r^2) + 1/r
```

**[obs]** The Hamiltonian has a **pole at `r = r_c`** in the `p_r^2` coefficient.
The Legendre transform `L -> H` is well-defined away from `r = r_c`, but
at `r = r_c` the coefficient `g_rr = 0` means the Legendre transform
degenerates: one cannot solve for `rdot` in terms of `p_r` at the critical
radius. This is a key subtlety not explicitly stated as a theorem but
noted implicitly. The Hamiltonian formulation is valid on `(0, r_c)` and
`(r_c, infinity)` separately.

**[obs]** Codebase connection: The symbolic Hamiltonian in
`src/weber_system.jl` (line 125) implements:
```julia
U_ij = k / r * (1 - r_dot^2 / (2 * c_squared))
```
This is the **Lagrangian-level** velocity-dependent potential. The codebase
does NOT use the FW metric rewriting; instead it works directly with the
velocity-dependent potential in the Hamiltonian via the standard Legendre
transform. The equivalence holds for `r != r_c`.

---

## 2. Theorem 2.1 -- Complete trajectory classification

### Statement

**[thm] Theorem 2.1** (Frauenfelder-Weber 2024). *The relative motion of two
equal charges in the plane under their mutual Weber force is as follows,
depending on their energy h compared to the critical energy h_c.*

The energy equation (Eq. 2.3):

```
rdot^2 = (ell^2 + 2r - 2hr^2) / (r(r_c - r))
```

The critical energy (Eq. 2.6):

```
h_c = V_eff(r_c) = ell^2/(2 r_c^2) + 1/r_c
```

where `V_eff(r) = ell^2/(2r^2) + 1/r`.

The classification:

| Case | Energy | Inside r_c | Outside r_c |
|------|--------|------------|-------------|
| **1** | h <= 0 | oscillation between 0 and r_c (finite time) | no solutions |
| **2b** | 0 < h < h_c | oscillation in (0, r_+), r_+ < r_c | escape to +/- infinity |
| **2a** | h > h_c | oscillation between 0 and r_c (both time dirs) | smooth reflection at r_+ > r_c, escape |
| **2c** | h = h_c | solutions reach 0 in one dir, infinity in other | solutions pass through r_c |

### Collision regularity within Theorem 2.1

FW2024 are very precise about the regularity at both endpoints:

**At r = r_c (critical radius):**
- Approach from `r < r_c`: speed `rdot ~ k/sqrt(r_c - r)` diverges (Eq. 2.4 neighborhood).
  Solution `r(t) = r_c - ((r_c - r_0)^{3/2} -/+ (3/2)kt)^{2/3}` reaches r_c in finite time.
  **Continuation**: solution can be continued `C^0` (continuously) but **NOT** `C^1` beyond r_c.
  FW call this "bouncing" back at r_c.
- Approach from `r > r_c`: the `r_+` turning point (Case 2a, 2b) gives `rdot ~ k sqrt(r_+ - r)`
  approaching r_+ at speed zero. Solution reflects **smoothly** (`C^infinity`).
- **Case 2c only** (`h = h_c`): solutions pass through r_c. This is the **unique energy**
  at which the critical radius can be crossed.

**At r = 0 (collision):**
- `ell = 0`: collision at speed `rdot -> +/- sqrt(2) c` (the Weber constant). `C^0`-continuable
  but not `C^1`. FW note this is a "bouncing" collision.
- `ell != 0`: collision at **infinite speed** (`rdot ~ k/sqrt(r)`, Eq. 2.4). This is the
  **spiraling collision** -- the angular velocity `phidot = ell/r^2` diverges and the
  total winding angle diverges. FW explicitly state (Section 1, "Interpretation"):
  "there are no periodic orbits" inside the Weber nucleus.

**[obs]** FW2024 do NOT state a non-regularizability theorem per se. Their
Theorem 2.1 is a **classification** of all radial motions. The statement
"not regularizable" appears only in Fig. 1 and the Interpretation section,
where they say the `ell != 0` singularity is "not regular" (in the ODE sense
of a regular singular point). The formal non-regularizability result was
developed in the codebase's `AngularMomentumRegularization.md` and
`TetheringImpossibility.md`, building on FW2024's classification.

**[obs]** FW2024's key observation about `ell = 0`: "it has yet not been studied
if there is a geometric regularization at the origin r = 0" (page 7).
This was subsequently resolved by the codebase: the collision bounce
approach in `solve.jl` (lines 1088-1141) implements precisely the
`C^0`-continuation that FW describe, reflecting `q_rel -> -q_rel` at a
finite bounce radius.

### The "no periodic orbits" statement

**[obs]** FW2024 state (Introduction, page 5): "inside the Weber nucleus there
are no periodic orbits." This is a consequence of Theorem 2.1: all interior
trajectories (`r < r_c`) oscillate between 0 and some `r_+ <= r_c`, but
for `ell != 0` they spiral into the origin in finite time and cannot be
continued. For `ell = 0`, the oscillation is purely radial (1D) and the
`C^0`-bouncing gives a periodic trajectory in the sense of a bouncing ball,
but FW do not count this as a periodic orbit because it lacks `C^1`
regularity.

**[obs]** This "no periodic orbits" statement is the classical counterpart of
the Gutzwiller trace formula issue they raise in the Interpretation: if
there are no classical periodic orbits, then Gutzwiller's semiclassical
approach to the quantum spectrum breaks down, making the quantum theory
genuinely different from the semiclassical approximation.

---

## 3. The Weber-Schrodinger equation

### 3.1 Laplace-Beltrami operator on Weber plane

**[thm] Lemma 3.1.** *The Laplace-Beltrami operator in polar coordinates
on the Weber plane `(R^2_x, g)` is:*

```
Delta f = (3/2 * 1/(r-r_c) - 1/2 * r/(r-r_c)^2) d_r f
          + r/(r-r_c) d_rr f + 1/r^2 d_phiphi f
```

(Eq. 3.8)

**[obs]** This is derived from the standard formula
`Delta f = (1/sqrt|g|) d_i(sqrt|g| g^{ij} d_j f)` with `|g| = r |r - r_c|`
and the cometric `g^{rr} = r/(r-r_c)`, `g^{phiphi} = 1/r^2`.

### 3.2 The Schrodinger equation (Eq. 3.10)

```
-(1/2) Delta_g psi + (1/r) psi = E psi
```

for `psi: R^2_x -> C`.

### 3.3 Separation of variables (Section 3.2)

Ansatz `psi(r, phi) = R(r) Y(phi)`.

Angular equation: `Y''(phi) = -2 ell Y(phi)` with solutions
`Y(phi) = c e^{ik phi}` where `k = sqrt(2 ell)`, so `ell = k^2/2`
for `k in N_0`.

**[obs]** The angular quantum number `k` is related to the angular momentum
quantum number by `ell = k^2/2`. This differs from the standard hydrogen atom
where `ell = l(l+1)`. The factor of 2 arises from the 2D (not 3D) geometry.

### 3.4 Radial equation

After separation, the radial equation inside the critical radius
`r in (0, r_c)` takes the Sturm-Liouville normal form (Eq. 4.14):

```
(p Rdot)' + q R = w E R
```

with coefficients depending on whether `ell = 0` or `ell != 0`.

---

## 4. Theorem A -- The main result

### Statement

**[thm] Theorem A** (Frauenfelder-Weber 2024). *The radial part of
Schrodinger's equation of the Weber nucleus is limit circle at both ends
of `(0, r_c)`, namely at the origin `r = 0` and at the critical radius
`r = r_c`.*

### Weyl's dichotomy

The theorem uses Weyl's classification of singular Sturm-Liouville problems:

- **Limit point**: only one solution is `L^2` near the singularity. The
  Schrodinger operator is essentially self-adjoint (no boundary condition
  needed).
- **Limit circle**: all solutions are `L^2` near the singularity. The
  operator is NOT essentially self-adjoint; a boundary condition must be
  chosen.

Since the Weber-Schrodinger equation is limit circle at BOTH endpoints,
**two boundary conditions** must be imposed -- one at `r = 0` and one at
`r = r_c` -- to select a self-adjoint extension and hence a discrete
spectrum. The spectrum depends on the choice of boundary conditions.

### Proof structure

Theorem A is proved via two propositions:

**[thm] Proposition 4.2** (Zero angular momentum, ell = 0). *The singular
Sturm-Liouville problem given by the 1-dimensional Weber Schrodinger
equation (4.17) on the interval (0, r_c) is:*
*(a) limit circle at the left origin boundary singularity 0;*
*(b) limit circle at the right critical radius boundary singularity r_c.*

**[thm] Proposition 4.4** (Non-zero angular momentum, ell != 0). *The singular
Sturm-Liouville problem given by the 1-dimensional Weber Schrodinger
equation (4.27) on the interval (0, r_c) is:*
*(a) limit circle at the left origin boundary singularity 0;*
*(b) limit circle at the right critical radius boundary singularity r_c.*

### Key difference: regular vs irregular singularity

**[thm] Remark 4.1.** FW note that the two cases differ fundamentally in
the ODE-theoretic classification of the singularity at `r = 0`:

- `ell = 0`: the singularity at `r = 0` is a **regular** singular point
  (Fuchsian). The indicial equation gives power-law solutions `r^{k_1}`,
  `r^{k_2}` with exponents:
  ```
  k_1 = -1/4 - sqrt(1/16 - 2r_c),    k_2 = -1/4 + sqrt(1/16 - 2r_c)
  ```
  Both are `> -1`, hence `L^2` near 0. The classical solutions (collisions)
  are regularizable, and natural boundary conditions exist in the quantum theory.

- `ell != 0`: the singularity at `r = 0` is an **irregular** singular point.
  Solutions oscillate wildly near 0 (asymptotic behavior involves Bessel
  functions `J_nu`, `Y_nu` with the argument going to infinity). Despite
  the wild oscillation, FW prove all solutions are still `L^2` (the key
  estimate uses the `1/rho^{1/4}` decay of Bessel functions at infinity).
  But no natural boundary condition exists. The classical non-regularizability
  is mirrored by this lack of natural BCs.

### The parallel (Fig. 1 of FW2024)

| Angular momentum | Classical | Quantum |
|---|---|---|
| ell = 0 | collisions (regularizable) | non-oscillating (natural BCs exist) |
| ell != 0 | spiraling (not regularizable) | oscillating (no natural BCs) |

**[obs]** This parallel is the deepest conceptual result of FW2024. It
establishes that the classical regularization obstruction has a precise
quantum counterpart: the failure of essential self-adjointness at the
irregular singular point.

### Proof of limit circle at r = 0 (ell != 0) -- the hard case

The proof (Section 4.2.1) is the technical heart of the paper. Key steps:

1. **Change of variable**: `rho = 1/r` maps `(0, r_c) -> (c^2, infinity)`.
   The ODE becomes a Bessel-type equation (Eq. 4.28-4.29):
   ```
   w'' + (lambda^2/(4rho) + (1 - nu^2)/(4rho^2)) w = 0
   ```
   where `lambda^2 = 8 ell r_c = (2k/c)^2` and `nu^2 = (1 - 32 r_c)/4`.

2. **Homogeneous solutions** (Eq. 4.30):
   ```
   u(rho) = rho^{1/4} J_nu(lambda rho^{1/2})
   v(rho) = rho^{1/4} Y_nu(lambda rho^{1/2})
   ```
   These are both `L^2` on `(c^2, infinity)` because Bessel functions decay
   as `rho^{-1/4}` for large argument, giving `|u|, |v| ~ rho^{-1/4} * rho^{1/4} = O(1)` --
   but after the change of variable the `L^2` weight compensates.

3. **Variation of parameters**: the inhomogeneous solution is controlled via
   Gronwall's lemma, using the Bessel decay bounds (Remark 4.5).

4. **Remark 4.5** (Bessel function bounds):
   - `|J_nu| <= 1` on `[0, infinity)`, `|Y_nu| <= 1` on `(1/r_c, infinity)`.
   - The critical decay estimate: `J_nu(rho) ~ sqrt(2/(pi rho)) cos(...)` as
     `rho -> infinity`, giving `1/rho^{1/4}` decay. This exponent is
     **smaller** than `1/2`, which is what makes the `L^2` integral converge.
   - `c_Y := cot(nu pi) + 1/(c^{2nu} Gamma(1-nu))` is a bounded constant.

**[obs]** The Bessel function analysis is where the paper's analytical difficulty
lies. The finiteness of the `L^2` norm at `rho -> infinity` (i.e. `r -> 0`)
hinges on the `1/rho^{1/4}` decay of Bessel functions, which just barely
makes the integral converge. A slightly different exponent would break
the argument.

### Proof of limit circle at r = r_c -- both cases

The proof at the critical radius (Sections 4.1.2 and 4.2.2) is comparatively
straightforward. The key ODE near `r = r_c` (for `ell = 0`, Eq. 4.21):

```
Rddot + (1/(2(r_c - r))) Rdot + (3/(2r)) Rdot = -(r_c - r)/r^2 * 2R
```

Homogeneous solutions: `u(r) = 1` and `v(r) = integral of s^{-3/2} sqrt(r_c - s) ds`.
The function `v` is bounded on `[r_0, r_c]` (the integrand is integrable),
so both solutions are `L^2`. Gronwall's lemma gives uniform boundedness
of all solutions.

For `ell != 0` at `r = r_c` (Section 4.2.2), the same approach works with
essentially the same homogeneous solutions and Gronwall argument, since the
`ell`-dependent terms in the ODE are smooth at `r = r_c`.

---

## 5. What FW2024 says about specific topics

### 5.1 Non-regularizability: precise scope

**[obs]** FW2024 does NOT contain a formal "non-regularizability theorem."
The paper classifies trajectories (Theorem 2.1) and notes that:

- `ell = 0` collisions reach `r = 0` at finite speed `sqrt(2) c` and are
  `C^0`-continuable (bouncing). FW leave open whether a "geometric
  regularization" exists (page 7).
- `ell != 0` collisions reach `r = 0` at infinite speed with infinite
  winding. FW call this singularity "not regular" (in ODE sense: irregular
  singular point).

The **formal non-regularizability theorem** (infinite winding number is a
topological invariant preserved by smooth coordinate changes) was developed
in the codebase (`AngularMomentumRegularization.md`), not in FW2024 itself.

**Boundary cases**: The ell = 0 case IS regularizable. FW confirm this
implicitly by noting the finite collision speed and `C^0` continuability.

### 5.2 The ell = 0 question

FW2024 explicitly states (page 7): "it has yet not been studied if there
is a geometric regularization at the origin r = 0."

**[obs]** The codebase has answered this question affirmatively for practical
purposes: the collision bounce (`CollisionBounce(r)` callback passed to
`solve`) implements a `C^0` continuation. The `AngularMomentumRegularization.md`
confirms that standard Sundman/Levi-Civita regularization works for `ell = 0`
(Experiment 3 in that document).

### 5.3 Quantum-mechanical treatment and bound states

**[obs]** FW2024's quantum analysis reveals that the Weber nucleus supports
a **discrete spectrum** inside `(0, r_c)`, but the spectrum is
**non-unique**: it depends on the choice of boundary conditions at both
endpoints. This is fundamentally different from the hydrogen atom, where
the spectrum is uniquely determined.

FW2024 does NOT compute the spectrum explicitly. They establish:
1. The problem is limit circle at both ends (Theorem A)
2. Therefore boundary conditions must be chosen (not uniquely determined)
3. The choice of BCs determines the energy levels

The Gutzwiller trace formula connection (Introduction, "Interpretation"):
since there are no classical periodic orbits inside the Weber nucleus,
the standard semiclassical approach via Gutzwiller fails. FW suggest
finding "a semi-classical interpretation" as an open problem.

**[obs]** The earlier paper Frauenfelder-Weber 2019 ("The fine structure of
Weber's hydrogen atom: Bohr-Sommerfeld approach", Z. Angew. Math. Phys.
70(4), 105-116) treated the **unlike-charge** (hydrogen-like) case using
Bohr-Sommerfeld quantization. That paper is referenced as [6] in FW2024
and is the origin of the Lorentzian metric interpretation.

### 5.4 n-body extensions

**[obs]** FW2024 treats **only the 2-body problem** (two equal charges).
There is no discussion of n-body extensions whatsoever. The paper's
results are:
- Specific to the 2-body relative motion in the plane
- Dependent on separation of variables (angular/radial), which is
  unavailable for n > 2
- Dependent on conservation of angular momentum ell, which generalizes
  but complicates for n > 2

The codebase's 4-body investigations (FourBodyTwoPlusTwoMinus study)
go far beyond FW2024's scope by considering multi-pair interactions
where the critical radius constraint applies to each like-charge pair
independently.

### 5.5 Legendre transform subtleties

**[obs]** FW2024 implicitly handles the Legendre transform issue by working
on `(0, r_c)` and `(r_c, infinity)` separately. The transform:

```
p_r = (r - r_c)/r * rdot    =>    rdot = r/(r - r_c) * p_r
```

is well-defined for `r != r_c` but singular at `r = r_c`. The Hamiltonian
(Eq. 2.2) has a pole `r p_r^2/(r - r_c)` at the critical radius.

**Key point**: The Legendre transform maps the Lagrangian system to the
Hamiltonian system on each connected component of the Weber half-line
`R_x = (0, infinity) \ {r_c}`. The two components `(0, r_c)` and
`(r_c, infinity)` are **dynamically disconnected** except at the unique
energy `h = h_c` (Theorem 2.1, Case 2c).

The codebase implements the Hamiltonian globally (including both sides
of the critical radius), using the velocity-dependent potential form
rather than the metric form, which avoids the Legendre transform issue
at `r = r_c`.

### 5.6 The "Weber nucleus" concept

FW2024 defines the Weber nucleus as: two equal (like) charges bound inside
the critical radius `r_c = 1/c^2`. The binding mechanism is the Lorentzian
metric: inside `r_c`, the effective mass is negative, inverting the dynamical
response to the Coulomb repulsion so that it acts as attraction.

**At what energies can a like-charge pair be bound?**

From Theorem 2.1:
- `h <= 0`: always bound inside `r_c` (oscillation between 0 and r_c)
- `0 < h < h_c`: bound inside `r_c` with turning point `r_+ < r_c`
- `h = h_c`: transitional (can cross r_c)
- `h > h_c`: interior trajectories reach 0 in both time directions

So the Weber nucleus exists for all energies `h < h_c`, where
`h_c = ell^2/(2 r_c^2) + 1/r_c`. For `ell = 0`, `h_c = 1/r_c = c^2`.

**At what angular momenta?**

All values of `ell` permit binding inside `r_c`. But:
- `ell = 0`: trajectories bounce at `r = 0` (regularizable)
- `ell != 0`: trajectories spiral into `r = 0` (non-regularizable, finite-time singularity)

So the Weber nucleus is classically well-defined only for `ell = 0`.
For `ell != 0`, the classical motion terminates in a finite-time
singularity. Whether the quantum Weber nucleus with `ell != 0` is
well-defined depends on the choice of boundary conditions (Theorem A).

---

## 6. What FW2024 does NOT address

### 6.1 Gaps relevant to bound orbit searches

1. **No stability analysis.** Theorem 2.1 classifies trajectories but says
   nothing about orbital stability, Lyapunov exponents, or KAM-type
   persistence under perturbation.

2. **No n-body theory.** The entire paper is 2-body. The 4-body (2+/2-)
   problem investigated in the codebase has no theoretical support from
   FW2024 beyond the per-pair critical radius constraint.

3. **No regularization theorem.** FW classify the `ell = 0` collision as
   `C^0`-continuable and explicitly leave geometric regularization as an
   open question.

4. **No explicit quantum spectrum.** Theorem A establishes limit-circle
   at both endpoints but does not compute eigenvalues for any choice of
   boundary conditions. The spectrum remains unknown.

5. **No Bohr-Sommerfeld quantization inside r_c.** The 2019 paper [6]
   treats Bohr-Sommerfeld for the hydrogen-like (unlike-charge) Weber atom.
   No analogous quantization is done for the like-charge Weber nucleus.

6. **No discussion of the Zollner extension.** The kappa parameter
   (mismatch factor) from Zollner's electrogravitational theory is not
   mentioned. FW treat pure Weber electrodynamics only.

7. **No collision-bounce energy conservation analysis.** FW note the `C^0`
   bouncing at `r = 0` and `r = r_c` but do not analyze whether energy
   (or any other quantity) is conserved across the bounce. The codebase's
   numerical work (`CollisionBounceRegularization.md`) shows that the
   symplectic integrator with bounce reflection achieves ~0.01% energy
   conservation.

8. **No phase-space structure analysis.** McGehee blow-up, Conley index,
   and Morse theory of the effective potential are not discussed.
   (These were developed by Agent 12 in `12_homology_morse/NOTES.md`.)

9. **No Floer/symplectic homology.** The Rabinowitz-Floer homology
   approach to proving existence of periodic orbits (Agent 10,
   `10_floer_symplectic/NOTES.md`) is not mentioned.

### 6.2 Open questions explicitly posed by FW2024

1. "It has yet not been studied if there is a geometric regularization
   at the origin r = 0" (page 7). [Partially answered by codebase:
   collision bounce works practically; Sundman/LC work for ell = 0.]

2. "It would be interesting to find a semi-classical interpretation"
   of the spectrum (page 5). [Open. The absence of classical periodic
   orbits inside the nucleus blocks standard Gutzwiller.]

---

## 7. Connections to the 14-agent study

### 7.1 Confirmed results

- **Critical radius as metric signature flip**: FW2024 Section 1 = CriticalRadiusAndLikeChargeAttraction.md Section 7. Exact match.
- **Energy equation (2.3)**: Used extensively in TetheringImpossibility.md as equation (FW).
- **ell = 0 regularizability**: FW2024 confirms finite-speed collision. Codebase implements collision bounce.
- **ell != 0 non-regularizability**: FW2024's classification (Theorem 2.1) is the foundation for AngularMomentumRegularization.md and TetheringImpossibility.md.
- **No periodic orbits inside r_c**: Confirmed numerically by the 4-body study. Only the breathing square (super-critical) was found as a periodic orbit.

### 7.2 Extended results

- **Tethering impossibility** (TetheringImpossibility.md): Goes beyond FW2024 by proving external charges cannot stabilize sub-critical ell != 0 pairs. Uses FW's energy equation as starting point.
- **Seven failed regularizations** (AngularMomentumRegularization.md): Systematically tests approaches that FW's classification predicts should fail.
- **Floer homology framing** (Agent 10): Uses FW's critical-radius obstruction as the fundamental barrier preventing Floer-theoretic orbit counting across Sigma_+/-.
- **Morse/Conley analysis** (Agent 12): Operates in the super-critical regime where FW's theory does not apply.

### 7.3 Contradicted results

None. All codebase results are consistent with FW2024.

---

## 8. Technical errata and clarifications

**[obs]** FW2024 uses the convention `r_c = 1/c^2` throughout, with
normalized units `k = ee'/mu = 1`. To convert to the codebase convention:
`rho = ee'/(mu c^2) = k/c^2 = 1/c^2 = r_c`. The two conventions are
identical under the normalization.

**[obs]** FW2024's energy `h` is the total Hamiltonian value, not a
reduced energy. The critical energy `h_c = ell^2/(2 r_c^2) + 1/r_c`
is the value of the effective potential at the critical radius. In
physical units: `h_c = ell^2 mu c^4 / (2 (ee')^2) + mu c^2 / (ee')`.

**[obs]** The "Weber constant" `sqrt(2) c` appearing as the `ell = 0`
collision speed is `sqrt(2)` times the speed of light. In Weber's
original notation using `c_W = sqrt(2) c`, this is simply `c_W`.
The codebase (`CriticalRadiusAndLikeChargeAttraction.md` Section 5)
correctly identifies this as `rdot^2 -> 2c^2` at collision.

---

## References

- Frauenfelder, U., Weber, J. *Anal. Math. Phys.* **14**:31 (2024).
- Frauenfelder, U., Weber, J. *Z. Angew. Math. Phys.* **70**(4), 105-116 (2019).
- `theory/Regularization.md`
- `research/investigations/CriticalRadiusAndLikeChargeAttraction.md`
- `research/investigations/AngularMomentumRegularization.md`
- `research/investigations/TetheringImpossibility.md`
- `research/FourBodyTwoPlusTwoMinus/10_floer_symplectic/NOTES.md`
- `research/FourBodyTwoPlusTwoMinus/12_homology_morse/NOTES.md`
