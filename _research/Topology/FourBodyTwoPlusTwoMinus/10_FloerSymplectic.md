# Floer / Symplectic Homology Framing of the 2+/2− Weber Problem

Agent 10 deliverable. Purely expository. No numerics.
Statements are tagged **[thm]** (cited theorem) or **[conj]** (this report).

## 1. Cotangent setup

Let `M` denote the reduced configuration space of the 2+/2− 4-body Weber
problem in dimension `d ∈ {2,3}`, after removing collisions and the two
Weber critical hypersurfaces:
```
M = F(R^d, 4) \ (Σ_+ ∪ Σ_−),
```
with Agent 3's notation: `Σ_± = {r_{12,34} = ρ_{12,34}}`. Translation
reduction kills one `R^d` factor; rotation/momentum reduction is deferred to
Agent 2's coordinates. For the framing in this note we treat `M` as a smooth
open manifold; symmetry quotients descend formally.

The phase space is the cotangent bundle
```
T*M,    ω = dp ∧ dq = -dλ,    λ = p dq.
```
`T*M` is an **exact** symplectic (Liouville) manifold; the Liouville vector
field `Z = p ∂_p` satisfies `ι_Z ω = λ` and generates the radial-fiber
dilation. The 4-body Weber Hamiltonian
```
H(q, p) = T(p) + V(q) + W(q, p)
```
splits into kinetic `T = Σ |p_i|^2/(2 m_i)`, Coulomb `V`, and the velocity-
dependent Weber correction `W` (per-pair `−q_i q_j ṙ_{ij}^2 /(2 c^2 r_{ij})`,
made fiber-quadratic via Legendre transform). `H` is smooth on `T*M`.

Agent 3 established: `H^*(F(R^3,4))` already has rank 24, and removing
`Σ_+ ∪ Σ_−` adds two linking classes in degree `d−1`. So `M` has *richer*
cohomology than the Coulomb-only configuration space — a fact that will be
crucial when we invoke loop-space arguments below.

## 2. Energy hypersurfaces and contact type

For an exact symplectic `(W, dλ)` a regular level `Σ = H^{-1}(E) ⊂ W` is of
**restricted contact type** if `λ|_Σ` is a contact form, equivalently if the
Liouville vector field `Z` is transverse to `Σ` and points outward along it.
For a mechanical Hamiltonian `H = |p|^2/(2m) + U(q)` the radial fiber field
`Z = p ∂_p` satisfies
```
dH(Z) = |p|^2/m = 2 (E − U(q))   on   Σ.
```
This is positive throughout `Σ` exactly when `E > U(q)` strictly on the
projection `π(Σ) = {U ≤ E}`, i.e. when the Hill region has no
configuration-space boundary touching the level (no zero-velocity locus
inside the singular set). Then `Σ` is contact and the closed Hamiltonian
orbits on `Σ` coincide (up to reparametrization) with the closed Reeb orbits
of `(Σ, λ|_Σ)`.

For Weber:
- `V` has Coulomb singularities `−∞` at unlike collisions and `+∞` at like
  collisions (above the critical radius), so the Hill region is open and
  *unbounded*: particles can escape to spatial infinity.
- The like-pair critical spheres `Σ_±` produce *finite* divergences in
  effective masses (metric signature flip), not in `V` itself; nevertheless
  the dynamics is incomplete across them in the spiraling case
  (Frauenfelder–Weber 2024; Agent 3 §6).
- The Weber correction `W` is fiber-quadratic; it deforms the kinetic metric
  but does not destroy the Liouville-radial argument provided the deformed
  metric is positive definite — which fails precisely on `Σ_±`. **This is
  why we excised them from `M`.**

**Discussion of admissible energies.**
- For `E` very negative (deeply bound), the Hill region in `M` is a union of
  small neighborhoods of unlike-pair clusters; collisions of unlike pairs
  bring `V → −∞` so `H ≤ E` is satisfied on a non-compact slice of `T*M`
  reaching arbitrarily close to the unlike-pair diagonal. After Moser-style
  collision regularization (§4) these slices compactify to a smooth contact
  level.
- For `E` near zero, the level extends to spatial infinity in escape
  directions; standard symplectic-homology truncation handles this.
- For `E > 0` the system is unbounded both at infinity and at the unlike
  collisions; the level is non-compact and contact in the bulk but not
  Liouville-fillable in any obvious way.

## 3. Rabinowitz–Floer homology

**[thm]** (Cieliebak–Frauenfelder 2009.) For a compact restricted-contact-
type hypersurface `Σ ⊂ (W, dλ)` in an exact symplectic manifold, the
Rabinowitz action functional
```
A^H(γ, η) = ∫_0^1 γ*λ − η ∫_0^1 H(γ(t)) dt
```
on `LW × R` has critical points exactly the closed Reeb orbits on `Σ` (with
period `η`). Its Morse–Floer chain complex is the **Rabinowitz–Floer
homology** `RFH(Σ, W)`. It is invariant under Hamiltonian deformations of
`Σ` that stay within the contact-type category, depends only on the filling
data, and vanishes if `Σ` is *displaceable* in `W` by a compactly supported
Hamiltonian isotopy.

**Corollary (existence).** If `RFH(Σ, W) ≠ 0` then `Σ` carries at least one
closed Reeb orbit — i.e. at least one closed Hamiltonian orbit of `H` at
energy `E`. This is the homological existence theorem we want to apply to
the Weber level sets.

Two further properties we will need:
- **[thm]** (Abbondandolo–Schwarz 2006.) For `W = T*N` with `N` closed, the
  symplectic homology `SH_*(T*N)` is isomorphic to `H_*(ΛN)`, the homology
  of the free loop space.
- **[thm]** (Viterbo 1999; Salamon–Zehnder for the index correspondence.)
  `SH_*(T*N)` is computed by closed Hamilton orbits on starshaped
  hypersurfaces, with index = Morse index in the loop space. In particular,
  `dim SH_k = b_k(ΛN)` and these are nontrivial for infinitely many `k`
  whenever `N` has nontrivial topology (Sullivan; Vigué-Poirrier–Sullivan).

## 4. Obstructions specific to Weber

**(a) Collision regularization.** For the pure Coulomb 2-body problem at
negative energy, Moser (1970) and Ligon–Schaaf compactify the energy level
to the unit cotangent bundle of a sphere: `H^{-1}(E)/Reeb ≅ S*S^d`. RFH
(equivalently `SH^+`) of `T*S^d` is non-zero, recovering the existence of
closed Kepler orbits homologically. This is the bedrock case our framing
generalizes.

**(b) The Frauenfelder–Weber obstruction.** Agent 3's Conjecture 6.1
(building on Frauenfelder–Weber 2024 Thm 2.1) says: ℓ ≠ 0 sub-critical
inspirals reach `r = 0` with infinite winding in finite time. **Floer
interpretation**: any compactification of `T*M` that includes the
sub-critical stratum `W_± = {r < ρ}` would have to add cells whose link is a
covering of `S^1` of infinite degree — there is no such finite CW
compactification. Therefore RFH-type theories *cannot cross* the critical
stratum: each connected component of the supercritical region must be
treated as a separate Liouville domain. Numerical Floer counts are valid
only on those super-critical components.

**(c) Non-compactness in escape directions.** Standard fix: replace `SH_*`
with `SH^+_*`, the cone of the natural map from singular cohomology
(constants) into `SH_*`. `SH^+` is generated by *non-constant* closed orbits
and is the right object to count actual periodic Reeb orbits on
non-compact contact-type levels of cotangent bundles (Cieliebak–Oancea).

## 5. Conjectures

**[conj C1] (Contact type, intermediate `E`).** There exists a non-empty
open interval `(E_*, 0) ⊂ R` such that for every `E` in it, the truncated
energy level
```
Σ_E := H^{-1}(E) ∩ {r_{12} > ρ_{12} + δ, r_{34} > ρ_{34} + δ},
```
once compactified by Moser regularization across the unlike-pair collision
diagonal and bounded above by a Liouville-finite cutoff in the spatial
escape direction, is of restricted contact type in (the corresponding
Liouville completion of) `T*M`. The argument: on this region `V` is smooth
modulo Coulomb singularities of *opposite-sign* type (Moser-regularizable),
and the kinetic + Weber metric remains positive definite (we are bounded
away from `Σ_±`).

**[conj C2] (Non-vanishing RFH ⇒ periodic orbits).** With `Σ_E` as in C1,
`RFH_*(Σ_E) ≠ 0`. Hence there exists at least one closed Reeb orbit on
`Σ_E`, i.e. at least one periodic orbit of the 4-body 2+/2− Weber system at
every such energy in the supercritical region. Justification sketch: the
truncated `Σ_E` is Hamiltonian-isotopic, through contact-type levels, to a
slight perturbation of the corresponding Coulomb level (turning Weber off,
`c → ∞`); and the Coulomb level at the same energy is a deformation of
several copies of `S*S^{d−1}` glued along Moser-regularized diagonals, with
`RFH ≠ 0` as in (a). Invariance of RFH under contact-type deformation
transports non-vanishing.

**[conj C3] (Topology lower-bounds the orbit count).** The number of
geometrically distinct periodic-orbit families on `Σ_E` is bounded below by
```
∑_k dim SH^+_k(T*M)  =  ∑_k b_k(ΛM) − b_k(M),
```
where `M` is the supercritical reduced configuration space of §1. By
Abbondandolo–Schwarz / Viterbo this is infinite as soon as `M` has
nontrivial rational homotopy in degree ≥ 2 (Sullivan dichotomy), which is
the case for `d = 3` already at the level of `F(R^3,4)`: Agent 3 gives
`b_0=1, b_2=6, b_4=11, b_6=6` for `F`, plus two more `b_2` classes for
`F\Σ`. Hence the loop space `ΛM` has infinite total Betti number and
`Σ_E` carries **infinitely many distinct periodic-orbit families** in the
supercritical region.

## 6. Caveats

These are conjectures assembled from standard machinery, **not theorems**.
Several non-trivial items are unresolved for the Weber potential:

1. **Transversality / Fredholm theory.** The Weber Hamiltonian is
   fiber-quadratic but with a non-positive-definite metric near `Σ_±`. The
   standard Floer cylinder operator and its compactness require uniformly
   convex generating Hamiltonians; we have these only after excising a
   neighborhood of `Σ_±` and only for `E < 0`.
2. **Compactness against escape.** Maximum principle / Hofer-energy
   estimates need the Liouville completion to be of finite type. The
   spatial-infinity end of `T*M` is fine (cone-like), but the unlike-pair
   collision end requires Moser regularization, which is classical for
   Coulomb but not yet checked to commute with the Weber `(1 − ṙ^2/2c^2)`
   correction. Plausibly it does because `W` vanishes at the regularized
   collision (`ṙ` is bounded and `r → 0` in the kinetic term dominates),
   but this needs verification.
3. **Compactness against critical-sphere spiraling.** This is the hard
   obstruction (point (b) above). It blocks any RFH argument that would
   try to count orbits crossing `Σ_±`. We treat it as an axiom: **the
   theory only sees supercritical periodic orbits.**
4. **Symmetry quotient.** All of the above descends to the quotient by
   `S_{12} × S_{34}` (Agent 2/3) provided the action lifts to a strict
   exact-symplectic action on `T*M`, which it does (it acts by linear
   symplectomorphisms on each fiber).

We offer C1–C3 as **research directions** for a future analyst willing to
carry out the Floer-theoretic work for this specific potential.

## 7. Connection to Agent 5 (numerical periodic-orbit search)

Agent 5 hunts periodic orbits numerically. If the RFH framing of §3–§5 is
correct in spirit, then for `d = 3` Agent 5 should find at least *some*
short-period orbits in the supercritical region at moderate negative
energies — the topology guarantees their existence (in fact, infinitely
many). Concrete predictions from the framing:

- **Energy bands**: orbits should populate the entire interval `(E_*, 0)`,
  not just one isolated value. Failure to find any across a coarse energy
  scan is evidence either that `E_*` is closer to `0` than expected, or
  that the candidate orbits are very long-period and the integrator-time
  budget is too tight.
- **Action filtration**: `SH^+` is filtered by Hamiltonian action ≈ period
  × `E`. Short-period orbits correspond to the lowest-action generators
  of `H_*(ΛM)`, which by §5 live in degree `0` (the constants — meaning
  short orbits in each connected component of `ΛM`) and degree `2`
  (linking-class-detecting orbits). Geometrically, these are the
  *breathing* and *rotating* short orbits of the dimer-dimer family.
- **Negative result is also informative.** If Agent 5 finds nothing,
  the most likely culprits are (i) all Hill regions in our energy window
  are sub-critical (so the integrator dies on the spiral obstruction),
  or (ii) periods are >> tmax. Both are theoretically diagnosable from
  the contact-type interval `(E_*, 0)` of C1.

## 8. References

- Cieliebak, K., Frauenfelder, U. *A Floer homology for exact contact
  embeddings.* Pacific J. Math. **239** (2009), 251–316.
- Cieliebak, K., Oancea, A. *Symplectic homology and the Eilenberg–Steenrod
  axioms.* Algebr. Geom. Topol. **18** (2018).
- Frauenfelder, U., Schlenk, F. *Hamiltonian dynamics on convex symplectic
  manifolds.* Israel J. Math. **159** (2007), 1–56.
- Abbondandolo, A., Schwarz, M. *On the Floer homology of cotangent bundles.*
  Comm. Pure Appl. Math. **59** (2006), 254–316.
- Viterbo, C. *Functors and computations in Floer homology with applications
  I.* Geom. Funct. Anal. **9** (1999), 985–1033.
- Salamon, D., Zehnder, E. *Morse theory for periodic solutions of
  Hamiltonian systems and the Maslov index.* Comm. Pure Appl. Math. **45**
  (1992).
- Moser, J. *Regularization of Kepler's problem and the averaging method on
  a manifold.* Comm. Pure Appl. Math. **23** (1970), 609–636.
- Ligon, T., Schaaf, M. *On the global symmetry of the classical Kepler
  problem.* Rep. Math. Phys. **9** (1976).
- Vigué-Poirrier, M., Sullivan, D. *The homology theory of the closed
  geodesic problem.* J. Differential Geom. **11** (1976), 633–644.
- Frauenfelder, U., Weber, J. *A mathematical description of the Weber
  nucleus as a classical and quantum mechanical system.* Anal. Math. Phys.
  **14**:31 (2024).
- `theory/Regularization.md`, `_research/Topology/FourBodyTwoPlusTwoMinus/03_config_space_topology/NOTES.md`.
