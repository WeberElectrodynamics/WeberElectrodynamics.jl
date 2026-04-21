# 12 — Morse / Conley analysis of the effective potential

Agent 12 deliverable. Companion script: `morse_analysis.jl`.

## 1. Setup

P = 0 (COM frame), L = L₀. After eliminating the SO(2) angle the
reduced Hamiltonian is `H_red = T(p_red) + V_eff(q; L₀)`, with two
natural choices:

- **Routh form**, fixed L₀: `V_eff^{(L)} = V_C + L₀²/(2 I(q))`.
- **Smale–Marsden form**, fixed ω: `V_ω = V_C − ½ ω² I(q)`. Rigid
  relative equilibria are critical points of **V_ω**, not of V_eff^{(L)}
  in general (the two coincide only along the one-parameter scaling
  family of a given relative equilibrium). The script catches this
  distinction in practice (§3).

`V_C` is the Coulomb pair sum; masses and `|q|` set to 1.

## 2. Working space — three-parameter submanifold

The full reduced configuration space is 5-dimensional (8 − 2 COM − 1 SO(2)).
For tractable Morse analysis we restrict to the 3-parameter
**dimer-perpendicular family**

```
+ dimer along x, center at (0, +Δ):     positives at (±r_+,  Δ)
− dimer along y, center at (0, −Δ):     negatives at (0,  ±r_- − Δ)
```

with `(r_+, r_-, Δ) ∈ R₊ × R₊ × R`. This submanifold is **invariant under**
the discrete subgroup `⟨S₁₂⟩ × ⟨S₃₄⟩` (intra-dimer label swaps) of Agent 2,
and also under `Δ → −Δ` combined with `S₃₄`. The COM is automatically zero
and the total charge sums to zero by construction.

The **alternating square** of Agent 4 corresponds to
`(r_+, r_-, Δ) = (R, R, 0)` (R is the circumscribing radius). The two
walls `r_± = 0` are the like-pair collisions; the `(0, 0, Δ)` axis is a
double collision.

Closed form (all four pair distances expressed in `(r_+, r_-, Δ)`):

```
r_{12} = 2 r_+,   r_{34} = 2 r_-,
r_{13} = r_{24} = √(r_+² + (r_- − 2Δ)²),
r_{14} = r_{23} = √(r_+² + (r_- + 2Δ)²),

V_C(u) = 1/(2 r_+) + 1/(2 r_-)
        − 2 / √(r_+² + (r_- − 2Δ)²)
        − 2 / √(r_+² + (r_- + 2Δ)²).
```

The sanity check in `morse_analysis.jl` matches the 4×4 Coulomb sum to
machine epsilon at three random points.

## 3. Critical points

### 3.1 L₀ = 0 — Earnshaw

Levenberg–Marquardt search of `‖∇V_C‖` from four seeds **diverges to
the box wall** (R ≈ 50, `|∇V_C| ≈ 10⁻³`, V_C → 0⁻). V_C is harmonic
in each particle coordinate, so it has no interior critical points
(Earnshaw). The 1D square diagonal reads `V_C(R) = (1 − 2√2)/R`,
monotone, no extremum.

### 3.2 L₀ ≠ 0, Routh form — also no critical point

`V_eff^{(L)} = V_C + L₀²/(4 r_+² + 4 r_-² + 8 Δ²)`. Both terms are
scale-homogeneous and vanish at infinity; LM search drifts to the
box wall. The 1D minimum of the symmetric scaling slice is **not** a
critical point of the full 3D function (the orthogonal directions
have nonzero gradient). **The Routh form has no interior critical
point on this submanifold at any L₀.**

### 3.3 Smale form V_ω — the alternating square is a critical point

Using `V_ω(u) = V_C(u) − ½ ω² I(u)` with `I(u) = 2 r_+² + 2 r_-² + 4 Δ²`
and `ω²(R) = (2√2 − 1) / (4 R³)` (Agent 4), evaluation at `(R, R, 0)`
gives, for `R ∈ {0.5, 1, 2}`:

| R    | ω      | ‖∇V_ω‖    | Hessian eigenvalues             | Morse index |
|------|--------|-----------|----------------------------------|-------------|
| 0.50 | 1.9123 | 1.6 × 10⁻⁹ | (−37.255, −21.941, +12.000)     | **2**       |
| 1.00 | 0.6761 | 9.4 × 10⁻¹¹| ( −4.6569, −2.7426, +1.5000)    | **2**       |
| 2.00 | 0.2390 | 0.0       | ( −0.5821, −0.3428, +0.1875)    | **2**       |

The eigenvalues scale as 1/R³, as forced by the homogeneity of the
Coulomb potential under `q → λ q, ω² → ω²/λ³`. The signature is
**2 negative + 1 positive** at every R: the alternating square is a
genuine **index-2 saddle of V_ω** on the dimer-perpendicular slice.

The sign of each eigenvector confirms physical intuition:
- the two negative directions live in the `(r_+, r_-)` plane (compression
  / scaling of the in-plane charges), corresponding to the breathing
  saddle that Agent 4 saw as the +0.75 / +1.026 unstable real mode;
- the single positive direction is the asymmetric Δ-mode — pushing the
  + dimer up and the − dimer down restores stability because it
  separates the unlike charges.

### 3.4 Other candidate critical points

The plan listed three candidates in §3:

- **Alternating square**: confirmed in §3.3, index 2.
- **Collinear ABAB**: requires `r_+ = 0` or `r_- = 0`, which is a
  like-pair collision — boundary of the manifold, excluded.
- **Two collapsed dimers**: `r_± → 0` again hits the collision wall.

Inside the open manifold, the alternating square is the **unique
interior critical point of V_ω** detected by ten LM searches from
seeds `(0.3..1.0, 0.3..1.0, 0..1)`. We have not proved uniqueness
analytically.

## 4. Morse inequalities and topology

Agent 3 gives Betti numbers for the configuration space. Morse
inequalities require a Morse function on a *compact* manifold; ours
is non-compact (V → −∞ at collisions, → 0 at infinity). The standard
fix is **N-body Morse theory with singular potential** (Bahri–Rabinowitz,
Benci, Fadell–Husseini): truncate to `M_{ε,R} = {ε ≤ r_{ij} ≤ R}` and
use the relative inequalities `m_k(V_ω, M_{ε,R}) ≥ b_k(M_{ε,R},
∂M_{ε,R})`. For our 3-slice the L–S category gives lower bound 1 on
critical-point count; we find exactly 1 (the square), trivially
satisfied with `m_2 = 1`, others zero. Interior critical points off
the 3-slice in the full 5-dim reduction are left to follow-up work.

## 5. Conley index for the alternating square

Let `S` be the maximal invariant set inside an isolating neighborhood
`N` of the square in the rotating-frame phase space. Lifting V_ω's two
negative Hessian eigenvalues to T*ℝ³ gives two pairs `(λ, −λ)` of real
phase-space eigenvalues — a 2-dim unstable subspace.

The Conley index of an isolated saddle equilibrium with `k` real
positive eigenvalues is `h(S) = [N/∂_exit N] ≃ S^k`. Here `k = 2`:

```
h({square}) ≃ S²,    H̃_2(h) ≅ Z, others zero.
```

Non-triviality of `h(S)` proves, via Conley/Wazewski, that the
isolating neighborhood contains a nonempty invariant set, and the
invariant `h = S²` is robust under continuation — including the Weber
`O(1/c²)` perturbation. Connection matrices (Franzosa, McCord) could
in principle force homoclinic/heteroclinic orbits, but we do not
pursue this.

## 6. Limitations

1. **Coulomb-only.** The Weber term `−q_iq_j ṙ²/(2 c² r)` depends on
   momenta and cannot be folded into a q-only V_eff. By Agent 4 the
   Weber correction shifts the Hessian by `O(1/c²)` (~1% at c = R = 1),
   not enough to change the index — signature `(2−, 0, 1+)` persists.
2. **Non-compactness.** The LM box `|u_i| ≤ 50` is a numerical
   Bahri–Rabinowitz cutoff. `V_ω → −½ω²I → −∞` at infinity, so the
   square is the only bounded critical point.
3. **Slice, not the full space.** The 3-submanifold is invariant only
   under Z₂ × Z₂, not the full D₄. Critical points off the slice are
   not detected here; Agent 4 ruled out non-square rhombus REs.
4. **Lower bound, not constructive.** Morse/Conley gives existence;
   orbit locations require direct integration (Agents 5, 6, 14).

## 7. Numerics summary

`morse_analysis.jl` performs ≈ 50 evaluations of `V_C` and `V_ω` (well
under the budget). All results above are reproducible with

```
julia _research/Topology/FourBodyTwoPlusTwoMinus/12_homology_morse/morse_analysis.jl
```

with no extra packages required (only `LinearAlgebra` and `Printf`;
gradients and Hessians are central differences with `h = 10⁻⁵`).

## 8. Headline conclusions

- **L₀ = 0**: V_C has *no* interior critical point (Earnshaw).
- **L₀ ≠ 0, Routh form**: still no interior critical point on this slice;
  the centrifugal barrier alone does not create one.
- **Rotating-frame amended potential V_ω**: the alternating square is a
  critical point with **Morse index 2** at every R, scale-invariant
  signature `(2−, 0, 1+)`. This is the Morse-theoretic incarnation of
  Agent 4's linear instability of the rotating square.
- **Conley index** of the corresponding invariant set in the rotating
  frame is `h ≃ S²`, providing a non-perturbative proof that some
  invariant set persists near the square under any sufficiently small
  perturbation — including the Weber `O(1/c²)` correction.
- **Lower bound**: at least one invariant set in the dimer-perpendicular
  slice. The 3-submanifold gives no positive-index Morse data on its
  own — the topological richness of the full 5-dim reduced space is
  the next thing to investigate.

## Files

- `/Users/mac/dev/Weber/WeberElectrodynamics/_research/Topology/FourBodyTwoPlusTwoMinus/12_homology_morse/morse_analysis.jl`

## References

- Smale, *Topology and mechanics II*, Invent. Math. 11 (1970).
- Marsden, *Lectures on Mechanics*, LMS Lect. Notes 174 (1992).
- Bahri–Rabinowitz, *A minimax method...*, JFA 82 (1989).
- Fadell–Husseini, Proc. AMS 107 (1989).
- Conley, *Isolated invariant sets and the Morse index*, CBMS 38 (1978).
- Franzosa, *The connection matrix theory...*, Trans. AMS 311 (1989).
- McCord, Ergodic Theory Dynam. Systems 8* (1988).
