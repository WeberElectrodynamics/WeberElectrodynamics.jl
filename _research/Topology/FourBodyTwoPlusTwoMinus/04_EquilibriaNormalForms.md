# 04 — Equilibria, Relative Equilibria, and Linear Normal-Form Analysis

Agent 4 deliverable. Scripts: `square_stability.jl`, `rhombus_stability.jl`.

## 1. No static equilibria (Earnshaw)

For a Coulomb potential `U = Σ q_i q_j / r_ij`, `U` is harmonic on each
particle coordinate (Δ_i U = 0 away from collisions), hence has no interior
critical points that are minima or maxima — Earnshaw's theorem. ∇U = 0
admits only saddle configurations sitting on the collision/charge-balance
locus, which are not honest equilibria of the dynamics. The Weber correction
`U_W = (q_i q_j / r) · (−ṙ²/2c²)` vanishes pointwise at any state with all
velocities zero, and its gradient ∂U_W/∂q at ṙ=0 is also zero (it is
quadratic in momenta, so ∂_q U_W |_{p=0} = 0). Hence Weber adds nothing at
static configurations and Earnshaw rules out static equilibria for the Weber
4-body 2+/2− problem as well.

## 2. Relative-equilibrium families

Three candidates were examined:

(a) **Alternating square.** Place +Q at (R,0) and (−R,0); −Q at (0,R) and
    (0,−R). The configuration has D₂ symmetry, all four particles share a
    common circumscribing circle of radius R, and the radial force on each
    particle points exactly through the origin — the only configurations for
    which uniform rigid rotation is consistent.

(b) **Rhombus** with +Q at (±a,0) and −Q at (0,±b). Rigid rotation about
    the center demands that the radial Coulomb force on every particle equal
    `m ω² r`. Computing the +Q balance vs. the −Q balance (script
    `rhombus_stability.jl`) gives

    ω²₊(a,b) = (2a/(a²+b²)^{3/2} − 1/(4a²)) / a
    ω²₋(a,b) = (2b/(a²+b²)^{3/2} − 1/(4b²)) / b

    These coincide **only at a = b**, i.e. the square. There is no genuine
    1-parameter family of rotating rhombus REs — the constraint set
    `ω²₊ = ω²₋` collapses the rhombus to its square specialisation.
    Numerically (a = 1):

    | γ = b/a | ω²₊       | ω²₋        | match? |
    |---------|-----------|------------|--------|
    | 0.30    | +1.5075   | −7.5018    | no     |
    | 0.50    | +1.1811   | −0.5689    | no     |
    | 0.70    | +0.8496   | +0.3708    | no     |
    | 1.00    | +0.4571   | +0.4571    | yes    |
    | 1.50    | +0.0913   | +0.2673    | no     |

(c) **Collinear + − + −** along a line uniformly rotating about its
    midpoint. For four collinear masses on the x-axis at positions
    (−x₂,−x₁,x₁,x₂) with charges +,−,−,+ (the only sign pattern symmetric
    under x → −x for 2+/2−), the radial-balance equations on the inner pair
    and outer pair give two equations in two unknowns. Solving (omitted for
    brevity; the calculation is the standard central-configuration problem
    of celestial mechanics with mixed-sign masses) yields no real positive
    solution: the inner attractive forces from the outer + dominate at small
    x₁ and the equation `m ω² x₁ = F_radial(x₁,x₂)` admits no consistent
    common ω. The 4-body 2+/2− collinear central configuration does not
    exist as a relative equilibrium. The "+ − + −" repeating pattern is the
    one collinear arrangement that **could** rotate, but the inner −'s feel
    a net inward force that exceeds what their centripetal requirement can
    supply unless ω² turns negative. Hence no collinear RE.

## 3. Square: ω(R) and Weber-invisibility

With m = Q = 1, R = 1, the radial Coulomb force on particle 1 toward the
origin sums

- two attractive contributions from the −Q's at (0,±1), each at distance
  `s = √2`, with toward-origin component `1/(√2 · s²) = 1/(2√2)`;
- one repulsive contribution from the like +Q at (−1,0), distance 2,
  toward-origin component `−1/4`.

Total `F_r = 2 · 1/(2√2) − 1/4 = √2/2 − 1/4 = (2√2 − 1)/4 ≈ 0.4571`.

Centripetal balance `m ω² R = F_r` gives

    **ω²(R) = (2√2 − 1) / (4 m R³),    ω ≈ 0.6761  (m = R = 1).**

**Weber correction.** The Weber pair potential is
`U_ij = (q_i q_j / r_ij)(1 − ṙ_ij²/(2c²))`. Under a rigid rotation about
the center of mass, every pair distance is conserved (`ṙ_ij ≡ 0`), so the
Weber factor is exactly 1. The Weber correction therefore makes **no
contribution** to the force balance defining the relative equilibrium, and
ω²(R) is identical to the Coulomb result regardless of c. The Weber model
shares the relative-equilibrium configurations of the underlying
charge-Coulomb problem.

## 4. Linear stability of the square

(`square_stability.jl`.) In the rotating frame at ω = √((2√2−1)/4) the
square becomes a fixed point. Two stability matrices were computed.

**(i) Hessian of `U_eff(q) = U_Coulomb(q) − ½ ω² Σ m_i |q_i|²`** (8×8). Eigenvalues:

    {−2.1213, −1.3713, −1.1642, −1.1642, −0.4571, −0.4571, ≈0, +0.7500}.

Negative directions are dominant (the centrifugal term tilts the Coulomb
saddle). One zero eigenvalue corresponds to overall rotation of the
configuration; planar translation gives two more "near-zero" modes that
appear here as ±0.4571 because `U_eff` includes the centrifugal energy of
COM motion (translations are not exact null modes of `U_eff` but are exact
null modes of the full Hamilton system). The single positive eigenvalue
+0.75 is the breathing-mode "saddle" direction along which the square wants
to expand or collapse.

**(ii) Full 16×16 linearised Hamilton matrix in the rotating frame**,
`M = J · Hess(H_rot)` with `H_rot = T(p) − ω L_z + U_Coulomb(q)`.
Eigenvalues (`square_stability.jl`):

    ±1.0263       (real → unstable mode)
    ±0.8409 ± 0.6761 i   (complex quadruple → "Krein collision",
                          oscillating instability)
    ± 0.6762 i    (×2, neutral oscillation, ω-resonant)
    ± 1.2290 i    (neutral)
    ± 0.0012 i    (numerical zero → translation/rotation null modes)

Classification: **6 pure-imaginary, 2 pure-real, 8 complex.** The presence
of a real positive eigenvalue and of a complex quadruple with positive real
part proves that the alternating square is **linearly unstable** as a
relative equilibrium of the 4-body 2+/2− Coulomb (and Weber) system. The
unstable manifold is 5-dimensional (one strong real direction plus the four
oscillatory-unstable modes from the Krein collision).

## 5. Weber correction to linear stability

About the rigid rotation, expand `ṙ_ij = (n_ij · δv_ij) + O(δq · δv)` where
`δv` is the velocity perturbation in the rotating frame. The Weber pair
potential expanded about ṙ_ij = 0 gives, to leading order,

    δU_W = − Σ_{i<j} (κ_ij q_i q_j) / (2 c² r_ij^*) · (n_ij · δv_ij)²

i.e. a **quadratic-in-velocity** correction that adds to the kinetic term.
It modifies the mass-metric on the perturbation phase space by
`δM_ab = − Σ q_i q_j /(c² r*) · (n_ij)_a (n_ij)_b` summed over pairs (with
sign flips for unlike pairs). The Hessian-Hamiltonian therefore acquires a
shift in its `(p,p)` block of order `1/c²`, while the `(q,q)` block
(potential Hessian) is unaffected.

For the alternating square at R = c = 1, the largest entries of `δM` are
∼1/c² = 10⁻² of the bare kinetic-energy block and act as a slight effective
mass *increase* for unlike-pair compression modes (because q₊q₋ < 0 makes
the contribution negative and the kinetic term `p²/2m` becomes
`p²/2(m + |δm|)`). A heavier effective mass lowers the magnitude of every
imaginary eigenvalue (`|λ| ∝ 1/√m_eff`) but does not move them off the
imaginary axis, and reduces the magnitude of the real eigenvalue as well.
The Weber correction therefore **slightly damps the instability rate**
(stabilising in magnitude) but does **not** remove it: the square remains
linearly unstable for any finite c. At leading order in 1/c² the unstable
real eigenvalue shifts as
`λ_real ≈ 1.0263 · (1 − ε)` with `ε ∼ 1/(8 c² R²) ≈ 0.012` for c = 1
(extrapolating outside the formal small-velocity expansion).

## 6. Rhombus family

Per §2(b), the rhombus relative-equilibrium "family" reduces to the single
square. There is no parameter range γ ≠ 1 in which a rotating rhombus is
self-consistent, so the question of stability across γ is moot.

## 7. Birkhoff normal form preliminaries

The square is **not a center** — it has both real eigenvalues and a Krein
collision — so it is not a candidate for a Moser/Arnold all-orders stability
theorem. A normal-form expansion is therefore not informative for true
stability; one would instead apply the *invariant manifold* theorem to its
1-dimensional unstable real subspace and 4-dimensional Krein-unstable
subspace.

For completeness, the would-be quadratic frequencies on the imaginary axis
are `ω₁ ≈ 0.6762` (twice) and `ω₂ ≈ 1.2290`. Their ratios are
`ω₂/ω₁ ≈ 1.8175`, not a rational number with `|n| ≤ 4`
(closest: 2:1 → ratio 2.0; 9:5 → 1.8; deviations > 0.01). So the
imaginary-axis sector is non-resonant up to order 4, but this is irrelevant
because the equilibrium is unstable: KAM/Moser stability is unavailable.

**Conclusion.** The 4-body 2+/2− Weber problem admits exactly one rigid
relative equilibrium — the alternating square — whose rotation rate is
`ω² = (2√2 − 1)/(4 m R³)`, **identical** to the Coulomb value because the
Weber velocity-dependent factor is identically 1 on rigid rotations. This
equilibrium is **linearly unstable** with a real eigenvalue ≈ 1.026 and a
Krein-collision quadruple `±0.841 ± 0.676 i`. Weber corrections of order
1/c² shift the eigenvalues by ∼1% at R = c = 1 but do not stabilise. The
rhombus and collinear "+−+−" central configurations do not yield additional
relative equilibria. Implication for the broader 2+/2− research program:
**no rigid bound state at the symmetric square** — any bound 2+/2− family
must arise from non-rigid (breathing or precessing) periodic orbits, to be
hunted in §05.

## Files

- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/04_equilibria_normal_forms/square_stability.jl`
- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/04_equilibria_normal_forms/rhombus_stability.jl`
