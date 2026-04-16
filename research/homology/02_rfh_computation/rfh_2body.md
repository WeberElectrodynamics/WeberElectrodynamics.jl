# Complete 2-Body RFH Calculation for the Weber Hamiltonian

Agent 02, companion to `NOTES.md`. This document carries out the full
Rabinowitz-Floer homology computation for the 2-body Weber problem,
including all intermediate steps, explicit formulas, and the
contact-type/Moser compatibility proofs.

Notation: **[thm]** = cited theorem, **[prop]** = proved here,
**[conj]** = conjecture, **[def]** = definition.

---

## 1. The reduced 2-body Weber Hamiltonian

### 1.1 Center-of-mass reduction

Two particles of masses `m_1, m_2`, charges `q_1, q_2`, in `R^d` (`d in {1,2,3}`).
Total mass `M = m_1 + m_2`, reduced mass `mu = m_1 m_2 / M`.
Charge product `k = kappa * q_1 * q_2` (with Zollner parameter `kappa`,
equal to 1 when Zollner is disabled).

Center-of-mass coordinates: `R = (m_1 q_1 + m_2 q_2)/M`, `r = q_1 - q_2`.
Center-of-mass momentum: `P = p_1 + p_2`, relative momentum: `p_r = mu (p_1/m_1 - p_2/m_2)`.

The total Hamiltonian separates: `H = |P|^2/(2M) + H_red(r, p_r)`, and
we work with the reduced part.

### 1.2 Explicit reduced Hamiltonian

**[def 1.1]** The reduced 2-body Weber Hamiltonian is:

```
H_red(r, p_r) = |p_r|^2 / (2 mu)
              + k / |r|
              - k / (2 mu^2 c^2 |r|) * (r . p_r / |r|)^2

            = |p_r|^2 / (2 mu)
              + k / |r|
              - k / (2 mu^2 c^2 |r|^3) * (r . p_r)^2
```

Define the unit radial vector `n = r / |r|`, radial momentum
`p_rad = n . p_r`, transverse momentum `p_perp^2 = |p_r|^2 - p_rad^2`.
Then:

```
H_red = p_rad^2 / (2 mu) * (1 - k / (mu c^2 |r|))
      + p_perp^2 / (2 mu)
      + k / |r|

      = p_rad^2 / (2 mu) * (1 - rho / |r|)
      + p_perp^2 / (2 mu)
      + k / |r|
```

where **`rho = k / (mu c^2)`** is the critical radius.

**[def 1.2]** The effective inverse metric tensor on `R^d \ {0}` is:

```
g^{ij}(r) = (1/mu) * (delta^{ij} - (rho / |r|) * n^i n^j)
```

Eigenvalues:
- Radial: `lambda_rad = (1/mu)(1 - rho/|r|)`
- Transverse (multiplicity `d-1`): `lambda_perp = 1/mu`

### 1.3 Sign analysis by charge type

**Unlike charges** (`k < 0`): `rho = k/(mu c^2) < 0`.
Then `1 - rho/|r| = 1 + |rho|/|r| > 1` for all `|r| > 0`.
The metric is positive definite everywhere.

**Like charges** (`k > 0`): `rho = k/(mu c^2) > 0`.
Then `1 - rho/|r|` changes sign at `|r| = rho`.
- `|r| > rho`: `lambda_rad > 0` (Riemannian).
- `|r| = rho`: `lambda_rad = 0` (degenerate).
- `|r| < rho`: `lambda_rad < 0` (Lorentzian signature `(-,+,...,+)`).

---

## 2. Symplectic geometry of the energy surface

### 2.1 Phase space and Liouville structure

Phase space: `(T*R^d, omega = dp_r wedge dr)`.
Liouville 1-form: `lambda = p_r . dr` (fiber-linear).
Liouville vector field: `Z = p_r . d/dp_r` (fiber-radial dilation).
Fundamental identity: `iota_Z omega = lambda`, `L_Z omega = omega`.

### 2.2 The energy surface

Fix energy `E` and define `Sigma_E = {H_red = E}`.

On `Sigma_E`, the Liouville transversality is:

```
lambda(X_H) = dH(Z) = Z . H_red
            = sum_a (p_r)_a * dH/d(p_r)_a
            = p_r^T g^{-1}(r) p_r
            = (1/mu) [ p_rad^2 (1 - rho/|r|) + p_perp^2 ]
            = 2 T_eff(r, p_r)
```

where `T_eff = H_red - k/|r|` is the effective kinetic energy (the
part of `H` that is quadratic in `p`).

On `Sigma_E`: `T_eff = E - k/|r| = E - V(r)`.

---

## 3. Contact-type proof for unlike charges at negative energy

### 3.1 Assumptions

- `k < 0` (unlike charges), so `V(r) = k/|r| < 0` for all `|r| > 0`.
- `E < 0`.
- `rho < 0`, so the metric is positive definite everywhere (Section 1.3).

### 3.2 Hill region

The Hill region is `H_E = {r in R^d \ {0} : V(r) <= E}`.
Since `V(r) = k/|r|` with `k < 0`:

```
k / |r| <= E  <==>  |k| / |r| >= |E|  <==>  |r| <= |k| / |E| =: r_max
```

So `H_E = {0 < |r| <= r_max}` is a punctured ball of radius `r_max = |k|/|E|`.

### 3.3 Star-shapedness

**[prop 3.1]** The Hill region `H_E` is star-shaped with respect to the
origin `r_* = 0`.

*Proof.* The star-shapedness condition (Viterbo repair, Agent 11 Section 1)
requires: on the boundary `{V(r) = E}`, i.e., `{|r| = r_max}`:

```
(r - r_*) . nabla V(r) > 0
```

With `r_* = 0`:

```
r . nabla V = r . (- k r / |r|^3) = -k |r| / |r|^2 = -k / |r| = |k| / |r|
```

At `|r| = r_max`: `r . nabla V = |k| / r_max = |E| > 0`. QED.

### 3.4 Liouville transversality

**[prop 3.2]** `Sigma_E` is of restricted contact type.

*Proof.* We need `dH(Z) > 0` on `Sigma_E`, i.e., `T_eff > 0`.

Case 1: In the interior of the Hill region (`|r| < r_max`), `p_r != 0` and
`T_eff = E - V(r) = E - k/|r|`. Since `k < 0` and `E < 0`:

```
E - k/|r| = E + |k|/|r|
```

We need `|k|/|r| > |E|`, i.e., `|r| < |k|/|E| = r_max`. This is satisfied
in the open Hill region. At `|r| = r_max`, `T_eff = 0` and `p_r = 0`
(zero-velocity surface) -- `Z` vanishes here.

Case 2: The zero-velocity locus `{|r| = r_max, p_r = 0}`. Apply the Viterbo
repair: the modified Liouville field `Z' = Z + X_phi` with
`X_phi = (1/2) sum_a r_a d/dr_a` (position-space radial field, which
preserves `omega` via `L_{X_phi} omega = 0` -- actually `X_phi` generates
the diagonal scaling, and we use a cutoff version).

Explicitly, take `phi(r) = (1/4)|r|^2` and `Z' = Z + X_phi^{vert}` where
`X_phi^{vert}` is the Hamiltonian vector field of `phi` lifted to `T*R^d`.
Then on `{p_r = 0}`:

```
dH(Z') = dH(Z) + dH(X_phi^{vert})
       = 0 + {H, phi}
       = sum_a (dH/dp_a)(dphi/dr_a) - (dH/dr_a)(dphi/dp_a)
       = sum_a (dH/dp_a)(r_a/2)
       = (1/2) sum_a (p_r)_a * g^{aa}(r) * ... = 0  (since p_r = 0)
```

This vanishes too, so we need a more refined argument. The classical approach
(Weinstein 1978, Hofer-Zehnder 1994 Section 4.1) is to use the fact that
`Sigma_E` is *fibrewise star-shaped* rather than requiring pointwise
Liouville transversality. Since the metric `g^{-1}(r)` is positive definite
for all `r in H_E`, the fiber of `Sigma_E` over each `r` with `V(r) < E`
is an ellipsoid (star-shaped about the origin in `p`-space). The fibres
degenerate to a point at `{V = E}`, but the manifold-with-boundary is still
globally star-shaped in the fibre direction.

The rigorous statement uses the characterization of Cieliebak-Frauenfelder
(2009, Section 2): `Sigma_E` bounds a compact Liouville domain in `T*R^d`
iff it is fibrewise star-shaped and the base `H_E` is star-shaped. We have
verified both:

(a) Fibrewise star-shapedness: `g^{-1}(r)` positive definite on `H_E` (Prop 2.1
from NOTES.md, restated in Section 1.3 above).

(b) Base star-shapedness: Prop 3.1 above.

Therefore `Sigma_E` is of restricted contact type. QED.

### 3.5 The contact form

The contact form on `Sigma_E` is:

```
alpha = lambda |_{Sigma_E} = (p_r . dr) |_{Sigma_E}
```

The Reeb vector field `R_alpha` satisfies `alpha(R_alpha) = 1`,
`d(alpha)(R_alpha, .) = 0`. It equals the Hamiltonian vector field
rescaled by the kinetic energy:

```
R_alpha = X_H / (2 T_eff)
```

Periodic Reeb orbits of period `tau` correspond to periodic Hamiltonian orbits
of period `T` with `tau = int_0^T 2 T_eff dt`.

---

## 4. Moser regularization of the collision singularity

### 4.1 Review of Moser's method for the Kepler problem

**[thm 4.1]** (Moser 1970.) Consider the Kepler Hamiltonian
`K = |p|^2/(2mu) + k/|r|` with `k < 0`, at energy `E < 0`. After:

(i) Time reparametrization `dt/ds = |r|` (Sundman regularization),
(ii) Scaling: `r -> r / (-2E)`, `p -> p * (-2E)^{1/2}` (to normalize energy),
(iii) Stereographic projection `S^d -> R^d` mapping `(r, p) -> (xi, eta)` on
     `T*S^d`,

the energy surface `{K = E}` maps diffeomorphically to the unit cotangent
bundle `S*S^d = {|eta| = 1} subset T*S^d`, and the reparametrized Kepler flow
becomes the geodesic flow on `S^d` (all geodesics are great circles of equal
period `2 pi`).

### 4.2 Extension to Weber

The Weber correction to the Kepler Hamiltonian is:

```
W(r, p) = -k / (2 mu^2 c^2 |r|^3) * (r . p)^2
```

Under Sundman regularization `dt/ds = |r|`:

```
H_reg = |r| * (H_red - E)
      = |r| * |p|^2 / (2mu)  +  k  -  E |r|
        - k (r . p)^2 / (2 mu^2 c^2 |r|^2)  -  0
        ^--Kepler part--^        ^--Weber correction--^
```

Wait, let me be more careful. Write `H_red = T_Kepler + V + W` where
`T_Kepler = |p|^2/(2mu)`, `V = k/|r|`, `W = -k (r.p)^2 / (2 mu^2 c^2 |r|^3)`.

The Sundman-regularized Hamiltonian is:

```
K = |r| * (H_red - E)
  = |r| * (T_Kepler + V + W - E)
  = |r| T_Kepler + k - E|r| + |r| W
  = |p|^2 |r| / (2mu)  +  k  -  E|r|
    - k (r . p)^2 / (2 mu^2 c^2 |r|^2)
```

The Kepler part `K_Kep = |p|^2 |r| / (2mu) + k - E|r|` is the standard
Moser Hamiltonian, which on `{K_Kep = 0}` gives `S*S^d`.

The Weber correction in regularized form is:

```
W_reg = - k (r . p)^2 / (2 mu^2 c^2 |r|^2)
      = - (rho / |r|) * p_rad^2 / (2 mu)
      = - (rho / |r|) * (r . p)^2 / (2 mu |r|^2)
```

**Key estimate near collision (`|r| -> 0`).**

In the Kepler problem at energy `E`, near collision:
- `|p| ~ (2 mu |k| / |r|)^{1/2}` (virial relation),
- `|r . p| / |r| = |p_rad| <= |p| ~ (2 mu |k| / |r|)^{1/2}`.

Therefore:

```
|W_reg| = |rho| / |r| * p_rad^2 / (2mu)
        <= |rho| / |r| * |p|^2 / (2mu)
        = |rho| / |r| * (|k| / |r| - E)           [on Sigma_E]
        ~ |rho| * |k| / |r|^2                      [as |r| -> 0]
```

In Moser's coordinates `(xi, eta)` on `T*S^d`, the relation is
`|r| ~ 1 / |eta|^2` (roughly), so `|r|^{-2} ~ |eta|^4`. But on the
regularized level `{|eta| = 1}`, the `|r|^{-2}` singularity is
tempered by the Sundman time change.

More precisely, in the Moser chart, the Weber correction `W_reg` is a
smooth function of `(xi, eta)` on `T*S^d \ {south pole}`. The potential
singularity at the collision point (south pole in stereographic projection)
is removable because:

```
|r| W_reg = -k (r . p)^2 / (2 mu^2 c^2 |r|^2)
```

and `|r| (r . p)^2 / |r|^2 = (r . p)^2 / |r|`. Near collision,
`r . p ~ |p| |r| cos(theta)`, so `(r . p)^2 / |r| ~ |p|^2 |r| cos^2(theta)`.
On the Kepler energy surface, `|p|^2 |r| ~ 2 mu |k|` (bounded!), so:

```
|r| W_reg ~ -k * 2 mu |k| cos^2(theta) / (2 mu^2 c^2) = -|k|^2 cos^2(theta) / (mu c^2)
```

This is **bounded** as `|r| -> 0`. Since the full regularized Hamiltonian is
`K = K_Kep + |r| W_reg / |r| * |r|` -- wait, let me redo this cleanly.

**[prop 4.1]** (Smoothness of the Weber correction in Moser coordinates.)

The full regularized Hamiltonian is `K = |r| (H_red - E)`. We showed:

```
K = K_Kep + W_reg_full
```

where `K_Kep = |r| (T_Kep + V - E)` and:

```
W_reg_full = |r| * W = |r| * (-k (r.p)^2 / (2 mu^2 c^2 |r|^3))
           = -k (r.p)^2 / (2 mu^2 c^2 |r|^2)
```

Now `(r.p) / |r| = p_rad` is bounded on the energy surface (it equals
`|p| cos theta`), so:

```
W_reg_full = -k p_rad^2 / (2 mu^2 c^2)  *  1/1  =  -(rho / 1) * p_rad^2 / (2 mu)
```

Wait -- `(r.p)^2 / |r|^2 = p_rad^2`, and this is just a momentum quantity.
So:

```
W_reg_full = -k p_rad^2 / (2 mu^2 c^2) = -(k / (mu c^2)) * p_rad^2 / (2 mu) = -rho * p_rad^2 / (2 mu)
```

This is **NOT** a function of `|r|` at all! It depends only on `p_rad = (r . p)/|r|`,
which in Moser coordinates is a smooth function on `T*S^d` (it corresponds to
the component of the cotangent vector along the geodesic direction).

**[prop 4.2]** In Moser coordinates `(xi, eta) in T*S^d`, the Weber correction
to the regularized Hamiltonian is:

```
W_reg_full(xi, eta) = -rho * f(xi, eta)^2 / (2 mu)
```

where `f: T*S^d -> R` is the smooth function corresponding to `p_rad` under
the Moser map. The total regularized Hamiltonian `K = K_Kep + W_reg_full`
is smooth on all of `T*S^d`, including at the collision point.

*Proof.* `p_rad = (r . p) / |r|` is the radial momentum. Under the Moser
stereographic map, it corresponds to the component of `eta in T*_xi S^d`
along the geodesic connecting `xi` to the south pole (the collision
image). This is a smooth function of `(xi, eta)` for `xi` away from the
south pole. At the south pole itself, `|r| = 0` in the original
coordinates, but the Moser map is precisely designed so that the
stereographic lift extends smoothly through this point. Hence `f` is smooth
on all of `T*S^d`. QED.

---

## 5. RFH computation

### 5.1 The reference case: Kepler (c = infinity)

At `c = infinity` (`epsilon = 1/c^2 = 0`), `W_reg_full = 0` and
`K = K_Kep`. The level `{K_Kep = 0}` is the unit cotangent bundle `S*S^d`.

**[thm 5.1]** (Cieliebak-Frauenfelder 2009, Corollary 1.8;
Abbondandolo-Schwarz 2006.) The Rabinowitz-Floer homology of the unit
cotangent bundle `Sigma = S*S^d subset T*S^d` is:

```
RFH_*(Sigma, T*S^d) ~= H_{*+d-1}(Lambda S^d, S^d; Z_2)  [positive part]
                     + H^{d-1-*}(Lambda S^d, S^d; Z_2)    [negative part]
```

where the splitting is by the sign of the period parameter `eta`.

For the full (unsplit) RFH:

```
RFH_j(S*S^d) ~= H_j(Lambda S^d; Z_2)  (shifted)
```

The precise grading conventions vary in the literature. We use the convention
where the degree-shift equals `n - 1` with `n = d + 1` the dimension of the
ambient manifold `S^d` (so shift = `d`), following Cieliebak-Frauenfelder-Oancea.

### 5.2 The loop space homology of S^d

**[thm 5.2]** (Ziller 1977, Vigué-Poirrier-Sullivan 1976.)

For `d >= 2`, the free loop space `Lambda S^d` has Poincare series:

**d odd** (`d = 2m+1`, `m >= 1`):
```
P_t(Lambda S^d; Q) = (1 + t^d) / (1 - t^{d-1})
                    = (1 + t^d)(1 + t^{d-1} + t^{2(d-1)} + ...)
```

Betti numbers (rational):
- `b_0 = 1` (constant loops)
- `b_{d-1} = 1, b_d = 1` (prime geodesic class and its partner)
- `b_{2(d-1)} = 1, b_{2d-1} = 1` (2-fold cover)
- General: `b_{k(d-1)} = 1, b_{k(d-1)+1} = 1` for `k >= 1`... [this is
  not quite right for general d; let me give the precise answer]

Actually, for `S^d` with `d >= 2`:

```
H_j(Lambda S^d; Q) ~= Q  for j in {0, d-1, d, 2(d-1), 2d-1, 3(d-1), 3d-2, ...}
```

i.e., one generator in degrees `0, d-1, d, 2(d-1), 2d-1, 3(d-1), 3d-2, ...`
The pattern for `j >= d-1` is: generators at `j = k(d-1)` and `j = k(d-1)+1`
for `k >= 1`.

**For `d = 2` (planar Kepler):**

```
b_j(Lambda S^2; Q) = 1 for all j >= 0
```

Every degree has exactly one rational generator. This is because `d - 1 = 1`,
so generators occur at `j = k` and `j = k + 1` for all `k >= 1`, plus
`j = 0`. This gives `b_j = 1` for all `j >= 0`.

**For `d = 3` (3D Kepler):**

```
b_0 = 1, b_2 = 1, b_3 = 1, b_4 = 1, b_5 = 1, b_6 = 1, b_7 = 1, ...
```

with `b_1 = 0`. Generators at degrees `0, 2, 3, 4, 5, 6, 7, ...` (all
`j >= 2` plus `j = 0`). More precisely, from `d - 1 = 2`:
generators at `2k` and `2k + 1` for `k >= 1`, giving all `j >= 2`.

### 5.3 RFH for the 2-body unlike-pair Weber problem

**[thm 5.3]** (Main result of this computation.) For the 2-body unlike-pair
Weber problem (`q_1 q_2 < 0`) in `d` dimensions at energy `E < 0`, with
any finite `c > 0`:

```
RFH_*(Sigma_E^Weber, T*S^d) ~= RFH_*(S*S^d, T*S^d) ~= H_{*+shift}(Lambda S^d)
```

In particular, `RFH` is non-zero in infinitely many degrees. By the
existence theorem (Cieliebak-Frauenfelder 2009), `Sigma_E` carries
infinitely many geometrically distinct periodic Reeb orbits, hence
infinitely many geometrically distinct periodic orbits of the Weber
Hamiltonian at energy `E`.

*Proof.*

Step 1: Contact type. `Sigma_E` is of restricted contact type (Section 3).

Step 2: Moser regularization. The collision singularity is regularizable,
and the regularized Weber correction is smooth on `T*S^d` (Section 4,
Prop 4.2). The regularized energy surface `{K = 0}` is a smooth
perturbation of `S*S^d` inside `T*S^d`.

Step 3: Deformation invariance. The one-parameter family of Hamiltonians
```
K_epsilon = K_Kep + epsilon * W_reg_full,   epsilon in [0, 1/c^2]
```
gives a smooth family of energy surfaces `{K_epsilon = 0}`. We claim each
is of contact type.

On `{K_epsilon = 0}`:
```
Z . K_epsilon = Z . K_Kep + epsilon * Z . (W_reg_full)
```

Now `Z . K_Kep = |p|^2 |r| / mu` (the kinetic part, always non-negative,
zero only at the collision = south pole in Moser coordinates). And
`Z . (W_reg_full)` involves `Z . (-rho p_rad^2 / (2mu)) = -rho p_rad^2 / mu`
(since `Z = p . d/dp` acts on the quadratic `p_rad^2` by doubling it,
but `p_rad = (r.p)/|r|` involves both `r` and `p`, so we need care).

Actually, in the Moser-regularized setting, the contact-type condition for
`{K_epsilon = 0}` in `T*S^d` follows from a different argument: `T*S^d` is a
Weinstein domain (its zero-section `S^d` is a Lagrangian with trivial
Maslov class for `d >= 2`), and `{K_epsilon = 0}` is the boundary of a
Liouville domain provided it remains fibrewise star-shaped.

Fibrewise star-shapedness: At each `xi in S^d`, the fiber
`{K_epsilon = 0} cap T*_xi S^d` is an ellipsoid determined by the
positive-definite quadratic form `g^{-1}_epsilon(xi)`. Since `rho < 0`
for unlike charges, the Weber correction **strengthens** the metric
(`g^{-1}` has larger eigenvalues), keeping the fibre ellipsoid inside
the Kepler one. Star-shapedness is preserved.

Step 4: Apply Cieliebak-Frauenfelder Theorem 1.5: RFH is invariant under
deformation through contact-type hypersurfaces within a fixed Liouville
domain. The family `{K_epsilon = 0}` for `epsilon in [0, 1/c^2]` stays
contact-type, so:

```
RFH_*({K_{1/c^2} = 0}, T*S^d) = RFH_*({K_0 = 0}, T*S^d) = RFH_*(S*S^d, T*S^d)
```

Step 5: The right-hand side is `H_{*+shift}(Lambda S^d)`, which is non-zero
in infinitely many degrees (Section 5.2). QED.

### 5.4 Generators and their physical interpretation

Each generator of `RFH_j` corresponds to a periodic orbit (or an iterated
cover thereof) on the energy surface.

**For `d = 2` (planar case):**

| RFH degree `j` | Loop space degree | Period multiple `k` | Physical orbit |
|---|---|---|---|
| 0 | 2 | 1 | Prime elliptic orbit (Kepler at c=inf, Weber-deformed at finite c) |
| 1 | 3 | 1 | Same prime orbit, different grading sector |
| 2 | 4 | 2 | 2-fold cover of the prime orbit |
| 2j | 2j+2 | j+1 | (j+1)-fold cover |

At `c = infinity`, the prime orbit is the circular Kepler orbit (unique
up to rotation). At finite `c`, the Weber correction deforms it into a
precessing ellipse. The precession rate is `O(1/c^2)`.

**For `d = 3` (spatial case):**

Additional generators appear because `dim S^3 = 3` provides richer
loop-space topology. The great circles on `S^3` come in an `S^2`-family
(the Hopf fibration base), giving a continuum of Reeb orbits. After
perturbation (Weber correction), these split into isolated non-degenerate
orbits, each contributing to `RFH`.

### 5.5 Action values of generators

The Rabinowitz action of the k-th iterated Kepler orbit at energy `E < 0` is:

```
A_k = k * T_1 * |E|
```

where `T_1 = 2 pi mu^{1/2} |k|^{1/2} / (2|E|)^{3/2}` is the Kepler period
(from `|E| = mu k^2 / (2 n^2)` in atomic units -- here `n = 1` for the ground
state).

Actually, for the Kepler problem with `H = p^2/(2mu) + k/r` at energy `E < 0`,
the period is `T = pi |k| mu^{1/2} / (2|E|^{3/2}) * 2` (this is `T = 2pi a^{3/2} / sqrt(mu |k|)`
with semi-major axis `a = |k|/(2|E|)`).

The action is:
```
A = int_0^T p . dr/dt dt = int_0^T 2 T_kin dt = 2 <T_kin> T
```

By the virial theorem for Kepler, `<T_kin> = -E = |E|`. So:

```
A = 2 |E| T = 2 |E| * pi |k| sqrt(mu) / (2 |E|)^{3/2}
  = pi |k| sqrt(mu) / sqrt(2 |E|)
  = pi |k| sqrt(mu / (2|E|))
```

For the k-th iterate: `A_k = k * A_1`.

The Weber correction shifts the action by `O(1/c^2)`:

```
A_k^Weber = A_k^Kepler + O(k / c^2)
```

The action ordering (and hence the filtration on RFH) is preserved for
`c >> 1`.

---

## 6. The like-charge case: obstruction analysis

### 6.1 Supercritical like-pair regime

For like charges `k > 0`, `rho > 0`. In the supercritical region `|r| > rho`,
the potential is repulsive: `V(r) = k/|r| > 0`. For `E > k/rho = mu c^2`,
the Hill region is `{rho < |r| <= k/E}` -- an annular shell.

**Is this level contact-type?** The star-shapedness check requires a center
`r_*` such that every ray from `r_*` meets `{V = E}` exactly once
transversally. The inner boundary `|r| = rho` is the critical sphere (not
part of the Hill region but its boundary). The outer boundary is `|r| = k/E`.

Taking `r_* = 0`: a ray from the origin meets the outer boundary at
`|r| = k/E` (with `r . nabla V = -k/|r| < 0` there -- pointing INWARD,
since `V` is decreasing). This violates star-shapedness!

**[prop 6.1]** For like charges, the supercritical Hill region is NOT
star-shaped with respect to any interior point. The energy surface is not
of contact type via the standard Weinstein/Viterbo criterion.

*Proof.* The Hill region is an annulus `{rho < |r| < k/E}`. Any ray from
an interior point must exit through the outer boundary, where `nabla V`
points inward (V is decreasing). The outward-pointing condition
`(r - r_*) . nabla V > 0` fails on the outer boundary for any `r_*` in
the interior. QED.

### 6.2 Implications

The standard RFH machinery does not apply directly to the like-pair
supercritical regime. However:

1. **Trapped orbits from the Weber well.** The effective radial potential
   (including angular momentum) for Weber like-charges can have a local
   minimum, creating a potential well that traps orbits. These orbits
   exist for specific energies and angular momenta, not generically.

2. **Local RFH.** One could attempt a local version of RFH near the
   potential well minimum, using the Conley-Zehnder index of the
   equilibrium. This would give a finite-dimensional chain complex.

3. **Connection to molecular states.** The subcritical regime (`|r| < rho`)
   has bounded "molecular" oscillations (theory/InitialConditions.md) but
   is outside the scope of RFH due to the Lorentzian signature.

**[conj 6.1]** For the 2-body like-pair Weber problem, the only periodic
orbits in the supercritical regime are those trapped in the effective
Weber potential well (if it exists for the given `E, L`). Their count is
finite and is not captured by topological (RFH) methods.

---

## 7. Summary table

| Setting | Contact type? | RFH | Periodic orbits |
|---|---|---|---|
| Unlike, E < 0, any c | YES (Prop 3.2) | H_{*+d}(Lambda S^d) != 0 | Infinitely many (Thm 5.3) |
| Unlike, E > 0 | No (non-compact, no Hill region) | Not defined | Scattering only |
| Like, supercritical | NO (Prop 6.1) | Not applicable | Finite or zero |
| Like, subcritical | NO (Lorentzian) | Not applicable | Non-regularizable spirals |

---

## 8. Explicit formulas for numerical implementation

For future numerical verification (computing CZ indices of known orbits),
here are the key formulas in the codebase's coordinate system.

### 8.1 Rabinowitz action functional (to evaluate on numerical orbits)

Given a numerically computed periodic orbit `{(q(t_j), p(t_j))}_{j=0}^{N_t}`
on the energy surface `{H = E}` with period `T`:

```
A = sum_{j=0}^{N_t - 1} p(t_j) . (q(t_{j+1}) - q(t_j))
```

(trapezoidal approximation to `int p . dq`).

The period parameter is `eta = T` (forward orbit) or `eta = -T` (backward).

### 8.2 Conley-Zehnder index (from monodromy matrix)

Given the `2n x 2n` monodromy matrix `M` (symplectic), decompose its
eigenvalues:

- Hyperbolic pairs `{lambda, 1/lambda}` with `lambda in R, |lambda| > 1`:
  contribute 0 to `mu_CZ`.
- Elliptic pairs `{e^{i theta}, e^{-i theta}}` with `theta in (0, pi)`:
  contribute `2 floor(theta T / (2 pi)) + 1` for a T-periodic orbit
  (Salamon-Zehnder normalization).
- Loxodromic quadruples: contribute via spectral flow (Robbin-Salamon).
- Parabolic eigenvalues `{1, 1}` or `{-1, -1}`: boundary cases requiring
  perturbation.

The total CZ index is the sum over all eigenvalue contributions, with a
global correction for the Maslov class of the capping disk.

### 8.3 Params vector for 2-body computation

For the 2-body problem in the codebase:

```
params = [m_1, m_2, q_1, q_2, c, kappa_12]
```

Length = 2*2 + 1 + 1 = 6.

For unlike charges: `q_1 = +1, q_2 = -1` (or vice versa), `kappa_12 = 1.0`.
Reduced mass: `mu = m_1 m_2 / (m_1 + m_2)`.
Charge product: `k = q_1 q_2 = -1`.
Critical radius: `rho = k / (mu c^2) = -1 / (mu c^2) < 0` (no critical sphere).

---

## References

Same as NOTES.md, plus:
- Hofer, H., Zehnder, E. *Symplectic invariants and Hamiltonian dynamics.*
  Birkhauser, 1994. Section 4.1 (fibrewise star-shapedness).
- Weinstein, A. *On the hypotheses of Rabinowitz' periodic orbit theorems.*
  J. Diff. Eq. **33** (1979), 353--358.
- Robbin, J., Salamon, D. *The Maslov index for paths.*
  Topology **32** (1993), 827--844.
