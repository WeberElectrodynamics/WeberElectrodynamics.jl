# Reduced Hamiltonian — 4-body 2+/2− Weber in Jacobi coordinates

Equal masses `m₁=m₂=m₃=m₄=m`, charges `q₁=q₂=+q`, `q₃=q₄=−q`, all `κ_ij=1`.

## Coordinates

```
r₊  = x₂ − x₁                       (intra-positive,  μ₊ = m/2)
r₋  = x₄ − x₃                       (intra-negative,  μ₋ = m/2)
R   = (x₃+x₄)/2 − (x₁+x₂)/2         (inter-dimer,    μ_R = m)
R_c = (x₁+x₂+x₃+x₄)/4               (COM, removed)
```

Conjugate momenta:
```
P₊  = (p₂ − p₁)/2
P₋  = (p₄ − p₃)/2
P   = (p₃+p₄)/2 − (p₁+p₂)/2
P_c = p₁+p₂+p₃+p₄  ≡ 0  (COM frame)
```

Inverse (in the COM frame, R_c=0, P_c=0):
```
x₁ = −R/2 − r₊/2,    p₁ = −P/2 − P₊
x₂ = −R/2 + r₊/2,    p₂ = −P/2 + P₊
x₃ = +R/2 − r₋/2,    p₃ = +P/2 − P₋
x₄ = +R/2 + r₋/2,    p₄ = +P/2 + P₋
```

(Check Σpᵢ = 0 ✓ ; one verifies `{r_α, P_β} = δ_{αβ}` for α,β ∈ {+,−,R}.)

## Pair separations

```
r₁₂ = r₊
r₃₄ = r₋
r₁₃ = R + (r₊ − r₋)/2
r₁₄ = R + (r₊ + r₋)/2
r₂₃ = R − (r₊ + r₋)/2
r₂₄ = R − (r₊ − r₋)/2
```
with magnitudes `r_ij = |r_ij|` (Euclidean norm in `d` dimensions).

## Pairwise relative velocities

In the equal-mass case, `ẋ_i = p_i/m`, so `ẋ_j − ẋ_i = (p_j − p_i)/m`. Substituting:

```
ṙ₁₂ = 2 P₊ / m
ṙ₃₄ = 2 P₋ / m
ṙ₁₃ = (P + P₊ − P₋) / m
ṙ₁₄ = (P + P₊ + P₋) / m
ṙ₂₃ = (P − P₊ − P₋) / m
ṙ₂₄ = (P − P₊ + P₋) / m
```
(These are vector relative velocities. The scalar Weber `ṙ_ij = (r_ij · ṙ_ij)/|r_ij|`.)

## Kinetic term

```
T = P²/(2m) + P₊²/m + P₋²/m
```
(Using μ_R = m, μ₊ = μ₋ = m/2.)

## Potential — Coulomb limit (c → ∞)

```
U_Coul = +q²/r₁₂ + q²/r₃₄                  (like-pairs, repulsive)
       − q²/r₁₃ − q²/r₁₄ − q²/r₂₃ − q²/r₂₄  (unlike pairs, attractive)
```
Substitute the Jacobi pair expressions:

```
U_Coul(r₊,r₋,R) =
   q² [ 1/|r₊| + 1/|r₋|
        − 1/|R + (r₊−r₋)/2|
        − 1/|R + (r₊+r₋)/2|
        − 1/|R − (r₊+r₋)/2|
        − 1/|R − (r₊−r₋)/2| ]
```

The full Coulomb-limit reduced Hamiltonian:
```
H_Coul = P²/(2m) + P₊²/m + P₋²/m + U_Coul(r₊, r₋, R).
```
This Hamiltonian is invariant under (i) `r₊ → −r₊`, (ii) `r₋ → −r₋`, (iii) `r₊ ↔ r₋, R → −R`, (iv) parity, (v) time reversal — exactly the discrete group enumerated in `NOTES.md`.

## Weber correction (full reduced H)

Each Coulomb pair gets multiplied by `(1 − ṙ_ij²/(2c²))`:
```
H_Weber = T + Σ_{i<j} (q_i q_j / r_ij) · (1 − ṙ_ij² / (2c²))
```
where `q_iq_j` is `+q²` for (12),(34) and `−q²` for the four cross pairs, the `r_ij` are the Jacobi expressions above, and `ṙ_ij` is the scalar radial speed:
```
ṙ_ij = (r_ij · ṙ_ij^{vec}) / r_ij
```
with the six `ṙ^{vec}` already listed above.

Schematically:
```
H_Weber = P²/(2m) + P₊²/m + P₋²/m
        + q² [ (1/r₁₂)(1 − (2P₊/m)·r̂₁₂)²/(2c²)
              + (1/r₃₄)(1 − (2P₋/m)·r̂₃₄)²/(2c²)  ]
        − q² Σ_{cross} (1/r_ij)(1 − (ṙ_ij)²/(2c²)).
```
Expanding fully would produce a long expression with `O(20)` rational+square-root terms — the symbolic version is available directly from `WeberSystem(4,2).hamiltonian_symbolic` and a follow-up Jacobi substitution; the algebraic simplification was found to be expensive (Symbolics `simplify` did not converge in 15 min on this input). For analytic work it is most useful to keep `H_Weber` in the schematic form above and substitute Jacobi expressions only in the specific submanifold being studied.

## Symmetry-restricted simplifications

### Square-symmetric subspace (`C`-fixed: r₊ = r₋ ≡ s, P₊ = P₋ ≡ Π, R = 0, P = 0)

In this 4D submanifold the four cross pairs collapse to two pair lengths and the kinetic part collapses too:
```
r₁₂ = r₃₄ = s
r₁₃ = r₂₄ = (r₊ − r₋)/2 = 0    (when r₊ = r₋)
r₁₄ = r₂₃ = r₊ = s
```
Wait — when R=0 and r₊=r₋=s, the cross pairs degenerate (r₁₃=0). This is the **collinear coincidence** singularity, not a regular submanifold. The proper "square" relative equilibrium uses **r₊ ⟂ r₋** with equal magnitudes, not r₊ = r₋ as vectors. Let `r₊ = s ê_x`, `r₋ = s ê_y`, `R = 0`. Then:
```
r₁₂ = s, r₃₄ = s,
r₁₃ = (s ê_x − s ê_y)/2  →  |r₁₃| = s/√2
r₁₄ = (s ê_x + s ê_y)/2  →  |r₁₄| = s/√2
r₂₃ = (−s ê_x − s ê_y)/2 →  |r₂₃| = s/√2
r₂₄ = (−s ê_x + s ê_y)/2 →  |r₂₄| = s/√2
```
This is the alternating square. Coulomb potential:
```
U_□(s) = 2·(q²/s) − 4·(q²/(s/√2)) = q²/s · (2 − 4√2).
```
Negative ⇒ bound; the symmetric square is a Coulomb relative equilibrium with rotation rate ω determined by `2 m ω² s² = −s · dU/ds` — Agent 4 will compute the Weber correction.

### Brake submanifold (`T`-fixed: P = P₊ = P₋ = 0)

```
H|_T = U_Coul(r₊, r₋, R)        (and the Weber correction vanishes because every ṙ_ij = 0)
```
Hence on the brake submanifold the Weber problem coincides with the Coulomb problem. This is a **rigorous** simplification, and means brake-orbit boundary conditions for periodic orbits can be searched using only the Coulomb potential — the Weber correction enters only off the brake.

### Two-dimer asymptotic (|R| ≫ |r₊|, |r₋|)

Multipole-expand the four cross potentials in (r₊, r₋)/R. To leading order the (+−) dipole moments d₊ = q r₊, d₋ = q r₋ (or rather the centers-of-charge difference) couple as a dipole-dipole term ~ (d₊·d₋ − 3(d₊·R̂)(d₋·R̂)) / R³. This is the regime where the "dimer-dimer orbital" candidate orbits (Agent 5) live. The internal dimer dynamics are 2-body Coulomb (Kepler) plus Weber, giving rise to the standard precessing ellipses studied in `examples/two_body_reference.ipynb`.

## Cross-references

- `NOTES.md` — symmetry tables, reduction counts, invariant submanifolds.
- `verify_symmetries.jl` — numerical proof of every symmetry listed above.
- `theory/WeberElectrodynamics.md` — canonical Hamiltonian.
- Agent 4 — uses these reduced coordinates to find relative equilibria.
- Agent 5 — uses the `T`-fixed and `C`-fixed submanifolds for periodic-orbit shooting.
- Agent 11 — uses the schematic `H_Weber` to evaluate the contact-type condition on energy hypersurfaces.
