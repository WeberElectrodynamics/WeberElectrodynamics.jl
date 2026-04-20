# Agent 1 — Symbolic Hamiltonian Anatomy (4-body, 2+/2−, 2D)

**Scope.** Build the symbolic Weber Hamiltonian for N=4 in 2D, substitute the
canonical 2+/2− parameters (m=1, q=(+1,+1,−1,−1), c=1, κ≡1), decompose it, and
verify continuous/discrete conservation laws.

Artefacts in this directory:
- `run_analysis.jl` — reproducible driver.
- `hamiltonian_expanded.txt`, `hamiltonian.tex` — full expanded H.
- `pair_table.md` — per-pair Coulomb/Weber/ρ summary.
- `virial.txt` — q- and p-virial identities (numerically verified).
- `conserved_quantities.md` — Poisson-bracket table and discrete symmetries.

All non-trivial symbolic identities below were **verified numerically** at 6 random
phase-space points (seed 42). Using `Symbolics.simplify(...; expand=true)` directly on
the full H chokes inside `MultivariatePolynomials.gcd` on the square-root rational
terms, so the pipeline is: symbolic `substitute` → symbolic derivatives → numerical
evaluation at random points to certify identities to machine precision (<1e-10).

## 1. Hamiltonian form

```
H = Σ_i p_i² /(2 m_i)  +  Σ_{i<j} (q_i q_j / r_ij) [1 − ṙ_ij² /(2 c²)]
```

with `r_ij = ‖r_i − r_j‖` and `ṙ_ij = ((r_i − r_j)·(v_i − v_j))/r_ij`,
`v_i = p_i/m_i`. For the canonical substitution m=c=κ=1, q=(+1,+1,−1,−1):

- `T_kin = ½ Σ_i (px_i² + py_i²)`
- `U_c   = +1/r₁₂ + 1/r₃₄  − 1/r₁₃ − 1/r₁₄ − 1/r₂₃ − 1/r₂₄`
- `U_w   = − Σ_{i<j} (q_i q_j) ṙ_ij² /(2 r_ij)`

The numerical identity `H_sub − (T + U_c + U_w) = 0` is verified (max residual
3.6e-15 over 5 random samples). **DECOMP_OK = true.**

## 2. Pair anatomy (see `pair_table.md`)

| pair | q_i q_j | Coulomb | Weber sign | ρ_ij = q_i q_j /(μ c²), μ=½ |
|---|---|---|---|---|
| (1,2) | +1 | repulsive  | − lowers H (weakens repulsion) | **+2** |
| (1,3) | −1 | attractive | + raises H  (weakens attraction) | −2 |
| (1,4) | −1 | attractive | + raises H  (weakens attraction) | −2 |
| (2,3) | −1 | attractive | + raises H  (weakens attraction) | −2 |
| (2,4) | −1 | attractive | + raises H  (weakens attraction) | −2 |
| (3,4) | +1 | repulsive  | − lowers H (weakens repulsion) | **+2** |

The Weber correction `−q_i q_j ṙ²/(2c² r)` always carries the *opposite* sign of
the static Coulomb interaction. For like-charge pairs (1,2) and (3,4) this is
the "attractive sub-Weber-radius" channel: at ṙ² → 2c² the effective Coulomb
coefficient `q_i q_j (1 − ṙ²/(2c²))` flips sign. The critical radius ρ_ij is
positive (+2) for like-charge pairs and negative for the four cross pairs — a
negative ρ means the unlike-sign pair has no radial flip; its Weber term only
provides a velocity-dependent "weakening" of the attraction (and becomes
repulsive for `ṙ² > 2c²`, which would place the relative motion above the
Weber radial light-cone).

## 3. Poisson-bracket conservation laws

| F | definition | max \|{F,H}\| | conserved? |
|---|---|---|---|
| P_x | Σ p_{x,i} | 4.4e-16 | **yes** |
| P_y | Σ p_{y,i} | 1.9e-16 | **yes** |
| L   | Σ (x_i p_{y,i} − y_i p_{x,i}) | 1.6e-15 | **yes** |
| H   | H itself (sanity) | 0 | **yes** |
| D   | Σ q_k p_k (dilation) | 4.58e+1 | **no** |
| d_x | Σ charge_i · x_i  (charge dipole, x) | 3.75e+1 | **no** |
| d_y | Σ charge_i · y_i  (charge dipole, y) | 5.83e+0 | **no** |

Translational (Px, Py), rotational (L) and time-translation (H) invariances all
survive the Weber velocity coupling. This is expected: H is manifestly written
in terms of pair-relative scalars, so the Galilean action of translations and 2D
rotations acts trivially on each summand.

Closed forms for the non-vanishing residuals (derived in §5 below, verified
numerically):

- `{D,H} = 2T + U_c + 3 U_w`.
- `{d_x,H} = Σ charge_i · ∂H/∂p_{x,i}`, i.e. the x-component of the
  charge-weighted *velocity* (electric current) — non-zero whenever the total
  polarisation current is non-zero, which is generic for 2+/2−.

## 4. Discrete symmetries

| symmetry | action on (q,p) | residual max \|H∘g − H\| | invariant |
|---|---|---|---|
| C (charge swap) | (1↔2),(3↔4) | 1.1e-16 | **yes** |
| P (parity)      | q→−q, p→−p  | 1.1e-16 | **yes** |
| T (time reverse)| p→−p        | 2.2e-16 | **yes** |

All three are exact. Reasoning:

- **T.** H depends on p only through `p²` (kinetic term) and `ṙ²` (Weber term),
  both quadratic and even in p; the flip `p→−p` leaves them invariant.
- **P.** Under `q→−q, p→−p`, every pair separation `r_ij` is even (depends on
  differences squared), and `ṙ_ij² = ((Δr·Δv)/r)²` is even in both the r and v
  simultaneous flip. `p²` is obviously even.
- **C.** The two + particles (1,2) share identical mass (1) and charge (+1),
  and similarly (3,4). Swapping (1↔2) permutes the pair list
  `{(1,3),(2,3)}→{(2,3),(1,3)}`, etc., leaving each pair summand unchanged.

Additional symmetries not tested here but implied: full charge conjugation
`(1,2)↔(3,4)` combined with parity, and the D2 spatial dihedral subgroup that
acts on the square configuration — deferred to Agent 2.

## 5. Virial theorem (Weber-modified)

Both `U_coulomb` and `U_weber` are homogeneous of degree **−1** in q (the latter
because `ṙ² = ((Δr·Δv)/r)²` supplies a factor `r⁰ · v²` and `U_w ∝ v²/r`). So
Euler's theorem gives

    Σ_k q_k ∂U/∂q_k = −U   (residual 3.9e-15).

In contrast, `U_weber` is homogeneous of degree **+2** in p, while `T` is also
degree +2 in p, so

    Σ_k p_k ∂H/∂p_k = 2 T + 2 U_weber   (residual 7.1e-15).

Combining via `d/dt(Σ q·p) = Σ p·∂H/∂p − Σ q·∂H/∂q`:

    d/dt (Σ q·p) = (2T + 2 U_w) − (−U_c − U_w) = **2T + U_c + 3 U_w**.

Time-averaging along any **bounded** trajectory kills the LHS, giving the
Weber-modified virial identity:

> **2⟨T⟩ + ⟨U_c⟩ + 3⟨U_w⟩ = 0    ⇔    2⟨T⟩ = −⟨U_c⟩ − 3⟨U_w⟩.**

Comments:
- Reduces to the classical `2⟨T⟩ = −⟨U⟩` in the `c→∞` limit where `U_w → 0`.
- The factor-of-3 on `U_w` comes from the fact that `U_w` is simultaneously a
  potential (degree −1 in q) **and** a mass-correction kinetic term (degree +2
  in p) — it is counted once on each side and adds up.
- For bound 2+/2− states `U_c < 0` (unlike-sign pairs dominate) and `U_w` has
  mixed sign; the modified virial constrains ⟨T⟩ non-trivially and tightens
  the classical bound whenever `⟨U_w⟩ > 0`, which is the generic unlike-pair-
  dominated regime.
- `D = Σ q·p` therefore acts like a "virial charge": its non-conservation
  `{D,H} = 2T + U_c + 3 U_w` is precisely the left-hand side of the modified
  virial identity, so on a time average it vanishes — **D is a bounded but not
  conserved quantity on every bounded trajectory.**

## 6. Key takeaways for the wave-1 cohort

1. **Good news for Agent 2 (symmetry reduction).** Continuous (P,L,H) and
   discrete (C,P,T) symmetries are all exact; no anomalies from the Weber
   velocity coupling. The full symmetry group is at least
   `(R² ⋊ SO(2)) × Z₂^C × Z₂^P × Z₂^T`.
2. **Good news for Agent 3 (config-space topology).** `ρ_{ij}` has a clean
   sign dichotomy: +2 for the two like-charge pairs, −2 for the four unlike
   pairs. The critical-radius strata are therefore only two real spheres
   (`r₁₂ = √2`, `r₃₄ = √2` in canonical units) — the unlike pairs contribute
   *no* critical-radius surface.
3. **Modified virial.** Any periodic-orbit or KAM work (Agents 5,6) should use
   `2⟨T⟩ = −⟨U_c⟩ − 3⟨U_w⟩`, not the Coulomb version.
4. **Dilation charge D.** Bounded but not conserved; useful as a detector of
   secular drift / unboundedness in numerical experiments (Agent 14).

## 7. Reproducibility

```bash
julia --project=. research/FourBodyTwoPlusTwoMinus/01_symbolic/run_analysis.jl
```

The script rebuilds `HamiltonianSystem(4,2)`, performs the substitutions, dumps the
expanded H, computes all Poisson brackets symbolically, and verifies each
identity numerically (6 samples). Total wall clock: ~90 s on a warm Julia.
