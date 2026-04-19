# 11 — Contact Geometry and Reeb Dynamics for the 4-Body 2+/2− Weber Hamiltonian

Agent 11 deliverable. Script: `star_shaped_check.jl`. Refs: Weinstein 1979,
Viterbo 1987, Hofer 1993. Cross-refs: Agent 4 (square eigenstructure),
Agent 3 (config topology), Agent 10 (RFH — pending). Non-theorems: **[conj]**.

## 1. Weinstein's contact-type criterion for natural Hamiltonians

Let `H : T*ℝⁿ → ℝ`, `H = ½ Σ pᵢ²/mᵢ + V(q)`, `E` a regular value, `ω = dλ`
with `λ = Σ pᵢ dqᵢ`. The fibre Liouville vector `Y = ½ Σ pᵢ ∂/∂pᵢ` satisfies
`L_Y ω = ω`. `Σ_E = {H = E}` is **of contact type** iff some such `Y` is
transverse (`Y·H ≠ 0` on `Σ_E`); the contact form is `α = ι_Y ω|_{Σ_E}`.

```
Y · H = ½ Σ pᵢ ∂H/∂pᵢ = ½ Σ pᵢ²/mᵢ = T = E − V(q).
```

On the interior `{V < E}`, `Y·H = E − V > 0` automatically. The subtlety is
the boundary `{V = E, p = 0}` where `Y` vanishes. The Viterbo (1987) repair
replaces `Y` by `Y' = Y + X(q)·∂_q` with `L_{Y'}ω = ω` (`X = ½∇φ`), giving
the cleaner equivalent criterion:

> **Star-shaped Hill region.** `Σ_E` is of contact type iff the Hill
> region `H_E = {V ≤ E}` is star-shaped with respect to some center
> `q* ∈ H_E°`: every ray from `q*` meets the level set `{V = E}` exactly
> once and transversally (i.e. `(q − q*) · ∇V > 0` on `∂H_E`).

Equivalently `dV·(q − q*) − 2(V − E) > 0` on `H_E`. When this holds, the
Hamiltonian flow on `Σ_E` is, after positive time-reparameterisation, the
**Reeb flow** of `α`.

## 2. Caveat: Weber is not a natural Hamiltonian

The Weber pair Hamiltonian for particles `i,j` reads

```
H_ij = pᵢ²/(2mᵢ) + pⱼ²/(2mⱼ) + (qᵢqⱼ/rᵢⱼ) − (qᵢqⱼ/rᵢⱼ)·(ṙᵢⱼ²)/(2c²),
```

with `ṙᵢⱼ = nᵢⱼ·(vᵢ − vⱼ)` linear in `p`, `nᵢⱼ = (qᵢ−qⱼ)/rᵢⱼ`. The
velocity term is therefore *quadratic in `p`* with position-dependent
coefficient. Writing `ṙᵢⱼ = (nᵢⱼ·pᵢ)/mᵢ − (nᵢⱼ·pⱼ)/mⱼ`:

```
−(qᵢqⱼ)/(2c² rᵢⱼ) · ((nᵢⱼ·pᵢ)/mᵢ − (nᵢⱼ·pⱼ)/mⱼ)².
```

a rank-one quadratic form in `(pᵢ,pⱼ)` with projector
`Π_{ij}^{ab} = (nᵢⱼ)_a (nᵢⱼ)_b`. The whole Weber Hamiltonian rewrites as
`H = ½ p^T G(q) p + V_C(q)` with effective inverse metric `G(q)` block
matrix on `(ℝᵈ)⁴`.

For **N=2** with reduced mass `μ` and relative momentum `p_r`:

```
H_{2-body} = p_r²/(2μ) − (q₁q₂)/(c² r) · (p_r·n)²/(2μ²) + q₁q₂/r,
          = (p_r²/2μ) · (1 − (q₁q₂)/(μc²r) · cos²θ_p) + q₁q₂/r,
```

so `g⁻¹(q) = (1/μ)(𝟙 − (q₁q₂/(μc²r)) n⊗n)`.

For **N=4**, sum the rank-one corrections over all six pairs:
`G_{ij}(q) = (δ_{ij}/mᵢ)𝟙 − Σ_k κ_{ik}(qᵢqₖ)/(c²r_{ik}mᵢmₖ) Π_{ik}·(δᵢⱼ−δₖⱼ)` —
a small position-dependent perturbation of `diag(1/mᵢ 𝟙)`.

## 3. Lorentzian signature in the sub-critical regime — a fundamental obstruction

For a single pair, the eigenvalue of `g^{−1}` along `n` is
`(1/μ)(1 − q₁q₂/(μc²r))`. Define the **critical radius**
`ρ = q₁q₂/(μc²)` (well-defined and positive for like pairs). Then:

- `r > ρ` (super-critical): all eigenvalues positive → **Riemannian**.
- `r < ρ` (sub-critical, like-pair only): the radial eigenvalue **flips
  sign** while the transverse eigenvalues stay positive →
  **Lorentzian** signature `(−,+,…,+)`.

Standard symplectic/Floer machinery for cotangent bundles assumes a
**Riemannian/Finsler** kinetic structure: convexity makes `Σ_E` fibrewise
star-shaped and supplies the energy estimates for Floer trajectories. In
the Lorentzian regime the fibres of `Σ_E ⊂ T*ℝⁿ` are hyperboloids, not
ellipsoids, and `Y·H = ½ p^T G p` can be zero or negative.

> **Obstruction.** The sub-critical regime is not just topologically
> pathological (Frauenfelder–Weber 2024 spirals) but **geometrically
> outside** standard symplectic homology / RFH. Contact-type / Reeb
> statements must be restricted to the supercritical region. Cf.
> theory/Regularization.md §"Metric signature".

## 4. Supercritical region: the standard machinery applies

Define the supercritical configuration set

```
M_+ = { q ∈ ℝ⁸ \ Δ_collision  :  r_{ij}(q) > ρ_{ij}  ∀ i<j }.
```

On `M_+` the inverse metric `G(q)` is positive definite — a small
perturbation of `diag(1/mᵢ)` controllably below `1/m` whenever `r > ρ`.
`H` on `T*M_+` is a genuine Riemannian kinetic + potential system; the
Weinstein criterion applies with `V(q) = Σ qᵢqⱼ/rᵢⱼ` (the Weber
correction vanishes at `p = 0`).

## 5. Star-shapedness check at the alternating-square energy

Take Agent 4's alternating square `q*`, side `s = 1`, charges `±1`,
masses `1`, `c = 1`. The static potential energy at `q*` is

```
V(q*) = 2 · (+1·−1)/1 + 2 · (+1·−1)/1 + (+1)(+1)/√2 + (−1)(−1)/√2
     = −4 + √2  ≈ −2.5858,
```

matching the script output `E_sq = −2.585786…`. Sample 50 unit-norm rays
`v ∈ S⁷ ⊂ ℝ⁸` and walk outward `t ∈ (0, 50]`, evaluating
`V(q* + t v)` via `sys.hamiltonian_compiled(·, 0, params)`. Two metrics:

| metric                                       | result   |
|----------------------------------------------|----------|
| outward-monotone rays (V non-decreasing)     | 8 / 50 = **16%** |
| rays crossing `{V = E_sq}` at most once      | 40 / 50 = **80%** |

**Interpretation.** The square is a **saddle** of `V` (Agent 4 Hessian:
7 negative/near-zero, 1 positive), so most rays are not monotone. But
80% cross the level set at most once — the Hill region is locally
connected and the obstruction is mild. The 20% multi-crossings trace
through cul-de-sacs near the diagonal-collision stratum (Agent 3).

The square is a **bad** star-point. **[conj]** A better center — e.g. an
inflated configuration `s ≫ 1` strictly in `M_+` — should give a
star-shaped Hill region for `E ∈ (E_coll, E_escape)`. Follow-up should
pick `q*` as a strict local max of the restricted potential, not a saddle.

## 6. Reeb orbits ↔ periodic orbits

On a contact-type `Σ_E` with `α = ι_Y ω`, the Reeb field `R_α` (defined by
`α(R_α)=1`, `dα(R_α,·)=0`) satisfies `R_α = X_H / (Y·H)`: the Hamiltonian
flow reparameterised by `Y·H = T = E − V > 0`. **Periodic orbits of `H`
on `Σ_E` are in bijection with closed Reeb orbits of `α`.** Any orbit
Agent 5 finds — breathing square, dimer-dimer, figure-8 analogue — is, on
its supercritical level, a closed Reeb orbit.

## 7. Conley–Zehnder index preview

For a non-degenerate closed Reeb orbit `γ`, the linearised flow on the
contact distribution `ξ = ker α` is a path of symplectic matrices and
yields the **Conley–Zehnder index** `μ_CZ(γ) ∈ ℤ`, a Maslov-type integer
counting half-rotations. In Floer/RFH, `μ_CZ` grades the chain complex.

Agent 4 reports for the alternating square in the rotating frame:

- one real pair `±1.0263` (pure unstable),
- a Krein-collision quadruple `±0.8409 ± 0.6761 i`,
- imaginary pairs `±0.6762 i, ±1.2290 i`.

**[conj]** Treating the rigidly rotating square as a formal Reeb orbit
(its linearised data are intrinsic even though the orbit itself is
unstable), each imaginary pair contributes `+1` to `μ_CZ` per
quarter-rotation, while the real unstable pair contributes via
Robbin–Salamon spectral flow as a transverse hyperbolic direction. The
orbit is then "bad" in the SFT sense (orientation cancellation, zero
contribution to the differential) rather than a stable bound generator —
matching the dynamical fact that the square is not long-lived.

## 8. Conjecture C11

> **[conj]** **C11.** Let `H` be the 4-body 2+/2− Weber Hamiltonian on
> `T*ℝ⁸`, and let `M_+` be the supercritical configuration region of
> §4. Then there exists an open energy interval
> `I = (E_coll, E_escape) ⊂ ℝ`, with `E_coll < 0 < E_escape`, such
> that for every `E ∈ I` the level set `Σ_E ∩ T*M_+` is of contact
> type, and the induced Reeb flow admits **at least one closed Reeb
> orbit** (equivalently, `H` admits a periodic orbit at energy `E`).

Existence follows from proven cases of the Weinstein conjecture: Viterbo
1987 (star-shaped in `ℝ²ⁿ`) and Hofer 1993 (tight contact 3-manifolds).
The Weber correction is `O(1/c²)` on `M_+`, so for `c` large the level
is a small perturbation of the (star-shaped) Coulomb level and the
contact-type property is open. This is **weaker** than Agent 10's RFH
prediction (which counts orbits) but **more robust**: it needs only a
transversality argument, no Floer compactness in the Lorentzian region.

## References

- Weinstein, *On the hypotheses of Rabinowitz' periodic orbit theorems*,
  J. Diff. Eq. 33 (1979), 353–358.
- Viterbo, *A proof of Weinstein's conjecture in ℝ²ⁿ*, Ann. IHP Anal.
  Non Linéaire 4 (1987), 337–356.
- Hofer, *Pseudoholomorphic curves in symplectisations…*, Invent. Math.
  114 (1993), 515–563.
- Cieliebak–Frauenfelder, *Rabinowitz Floer homology*, JEMS 2009.
- Frauenfelder–Weber 2024, sub-critical spiral non-regularisability.

## Files

- `/Users/mac/dev/Weber/WeberElectrodynamics/research/FourBodyTwoPlusTwoMinus/11_contact_reeb/NOTES.md`
- `/Users/mac/dev/Weber/WeberElectrodynamics/research/FourBodyTwoPlusTwoMinus/11_contact_reeb/star_shaped_check.jl`
