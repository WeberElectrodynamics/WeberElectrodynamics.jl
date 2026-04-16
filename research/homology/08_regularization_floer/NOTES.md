# 08 -- Regularization of the Weber Critical-Radius Singularity in the Context of Floer Theory

Agent 08 deliverable. Cross-refs: Agent 11 (contact/Reeb, section 3 on Lorentzian
obstruction), Agent 10 (RFH), AngularMomentumRegularization.md (7 failed attempts),
TetheringImpossibility.md (external-charge stabilization impossible).

Files: `mcgehee_weber.md`, `l0_regularization.jl`.

## 1. McGehee blow-up at the critical sphere

### Setup

The 2-body Weber Hamiltonian in relative polar coordinates:

```
H = (r / (2 mu (r - rho))) p_r^2  +  ell^2 / (2 mu r^2)  +  e^2 / r
```

with critical radius `rho = e^2 / (mu c^2)`. The radial eigenvalue of the
inverse metric is `(1/mu)(1 - rho/r)`, which vanishes at `r = rho`.

### Blow-up substitution

Set `r = rho + s^2` (super-critical) or `r = rho - s^2` (sub-critical), `s >= 0`.
With canonical momentum `P_s = 2 s p_r`:

```
H = P_s^2 (rho + s^2) / (8 mu s^2)  +  ell^2/(2 mu (rho+s^2)^2)  +  e^2/(rho+s^2)
```

**Result: the blow-up does NOT regularize.** The kinetic term diverges as `s^{-2}`
at `s = 0`. A Sundman factor `dt = s^2 d(tau)` absorbs this, but the coordinate
change `s = sqrt(r - rho)` reintroduces a divergence of the same order. The net
effect is no improvement.

**Why this differs from collision regularization.** At `r = 0`, the singularity
is in the potential (`1/r`); Sundman + Levi-Civita absorbs it by multiplying `r`
into both sides. At `r = rho`, the singularity is in the *metric* (the kinetic
coefficient diverges); no analogous absorption is possible because the divergence
multiplies `p^2`, not the potential.

### Boundary manifold

The blown-up boundary at `s = 0` is `S^1 x R` (parameterized by angle `phi`
and momentum `P_s`). The Hamiltonian does not extend continuously to this
boundary -- it is a *degenerate* boundary where the metric tensor has a zero
eigenvalue. See `mcgehee_weber.md` for full calculation.

## 2. ell = 0 regularization through the critical sphere

### Analytical result

For radial (ell = 0) trajectories, the Frauenfelder-Weber energy equation reduces to:

```
r_dot^2 = (2(e^2/mu) - 2 h r) / (rho - r)
```

At `r = rho`, the denominator vanishes. The numerator `2(e^2/mu) - 2 h rho`
also vanishes if and only if `h = h* = e^2/(mu rho) = c^2`. By L'Hopital:

```
lim_{r -> rho} r_dot^2 = 2 h* = 2 c^2     (at the critical energy)
```

So `|r_dot| -> sqrt(2) c` at `r = rho` -- **finite and nonzero**.

### Numerical verification (l0_regularization.jl)

The script confirms:

1. **Super-critical infall** from `r = 1.5` at `h = h*`: reaches `rho` at
   `t = 0.354` with `|r_dot| = 1.414 = sqrt(2)` (exact match to theory).

2. **Through-critical integration**: passes smoothly through `rho`, continues
   to `r = 0`, arriving at `t = 1.061` total. Velocity is continuous at the
   crossing.

3. **Sub-critical ell=0 collision**: from `r = 0.99 rho`, reaches `r = 0` at
   finite time with `|r_dot| -> sqrt(2) c`, confirming the known result.

4. **Energy scan**: only `h = h*` allows crossing. For `h < h*`, `r_dot^2 -> +inf`
   (momentum blowup). For `h > h*`, turning point before `rho`.

### Key finding

The critical sphere is traversable for ell = 0, but only at a single energy
value `h* = c^2`. This is a **codimension-1 condition** (one constraint on a
1-parameter family of radial orbits), so the set of through-critical orbits
has measure zero in the radial phase space.

## 3. Floer-theoretic interpretation of the critical sphere

### Is the critical sphere a "stop"?

In the stopped Liouville sector framework (Sylvan 2019, Ganatra-Pardon-Shende 2020),
a **stop** is a codimension-2 isotropic submanifold of the contact boundary where
the Liouville structure degenerates. The critical sphere `Sigma_rho = {r = rho}`
does NOT fit this description:

- `Sigma_rho` is **codimension 1** in configuration space (and codimension 1 in
  the energy hypersurface), not codimension 2.
- It is a **metric degeneracy** (zero eigenvalue of the kinetic metric), not a
  singularity of the Liouville form.
- The Liouville vector field `Y = sum p_i d/dp_i` is well-defined at `Sigma_rho`;
  it is `Y . H = T` (kinetic energy) that can vanish there.

**Classification.** The critical sphere is a **signature-change hypersurface**
of the kinetic metric. In the symplectic/contact language:

- On the super-critical side `{r > rho}`, the energy hypersurface `Sigma_E` is
  of contact type (fibrewise star-shaped, Agent 11 section 4).
- On the sub-critical side `{r < rho}`, `Sigma_E` is NOT of contact type (the
  fibres are hyperboloids, not ellipsoids; `Y . H` can be negative).
- At `r = rho`, the boundary is **degenerate contact**: the contact form
  `alpha = iota_Y omega` degenerates along the radial direction.

This is closest to a **concave boundary** (Liouville vector field tangent or
inward-pointing along the degenerate direction). Excising a neighborhood of
`Sigma_rho` from the super-critical side yields a Liouville domain with
smooth convex boundary `{r = rho + epsilon}` for any `epsilon > 0`.

### Consequences for the RFH chain complex

**Excision approach.** Define the super-critical Liouville domain:

```
W_+ = T*{q : r_{ij} > rho_{ij} + epsilon  for all i<j}  intersect  {H <= E}
```

For `epsilon > 0`, this is a compact Liouville domain with smooth contact-type
boundary (Agent 11 section 4). The RFH complex `RFH_*(W_+)` is well-defined
and computes Rabinowitz Floer homology of the super-critical region.

**What is lost.** The excision discards:
- All sub-critical orbits (the Lorentzian region)
- Orbits that approach `Sigma_rho` asymptotically
- The measure-zero ell=0 orbits that cross `Sigma_rho`

The RFH generators are periodic orbits (Reeb orbits) on `{H = E} cap T*M_+`.
These are super-critical orbits that never approach `rho`. The critical sphere
acts as a **hard wall** for the Floer equation: pseudo-holomorphic curves
cannot cross `Sigma_rho` because the energy estimate fails in the Lorentzian
region (the `L^2` norm of `du/ds` can be negative).

**Algebraic consequence.** The RFH complex for the excised domain is the
"best available" Floer theory. It counts super-critical periodic orbits with
their Conley-Zehnder indices. It does NOT detect:
- Sub-critical bound states (ell=0 collision-bounce oscillations)
- Through-critical orbits (the measure-zero ell=0 family at h=h*)
- Any ell != 0 sub-critical dynamics (non-regularizable)

## 4. Pseudo-Riemannian Floer theory

### Literature survey

The sub-critical region has Lorentzian kinetic metric signature `(-,+,...,+)`.
Standard Floer theory assumes Riemannian kinetic structure for:
- Energy estimates (bounding the `L^2` gradient of the action functional)
- Compactness (Gromov-type bounds on pseudo-holomorphic curves)
- The Weinstein conjecture (existence of Reeb orbits on contact-type levels)

**Known results on indefinite/Lorentzian Floer theory:**

1. **Frauenfelder (2004), "The Arnold-Givental conjecture and moment Floer
   homology."** Handles group-invariant Lagrangians with mixed-signature
   quadratic forms, but only in the compact setting (loop spaces of compact
   Lie groups). Does not apply to the Weber non-compact case.

2. **Abbondandolo-Schwarz (2006), "On the Floer homology of cotangent bundles."**
   The foundational result for `T*M` Floer homology requires a **Riemannian**
   metric on `M`. No extension to pseudo-Riemannian base metrics is known.

3. **Abbondandolo-Majer (2006), "Lectures on the Morse complex for infinite-
   dimensional manifolds."** The Morse-theoretic approach to Floer theory
   requires the Hessian of the action to be Fredholm with finite Morse index.
   For Lorentzian kinetic terms, the Hessian is **indefinite** (infinite Morse
   index), so the standard Morse complex is not defined.

4. **No Weinstein conjecture in Lorentzian geometry.** The Weinstein conjecture
   ("every compact contact-type hypersurface carries a closed Reeb orbit")
   has no known analog for hypersurfaces in pseudo-Riemannian cotangent bundles.
   The fundamental issue is that the fibrewise star-shapedness (which makes
   `Y . H > 0`) fails when the metric has indefinite signature.

5. **Suhr (2013-2018), work on Lorentzian geodesic theory.** Stefan Suhr has
   developed variational methods for Lorentzian geodesics using the Fermat
   principle and causal structure. These give existence results for causal
   geodesics but do not provide a Floer-type chain complex.

6. **Giambio-de Luce-Masiello (2004), "Bernstein-Nagumo type conditions for
   geodesics on Lorentzian manifolds."** Variational existence results for
   Lorentzian geodesics, but the methods are finite-dimensional (minimax on
   finite-dimensional approximations) and do not yield a homological invariant.

**Conclusion.** There is NO established Floer theory for pseudo-Riemannian
Hamiltonians. The sub-critical Weber region lies genuinely outside current
symplectic-topological machinery. This is not merely a technical gap but a
fundamental obstruction: the infinite Morse index of the Lorentzian action
prevents the construction of a well-graded chain complex.

## 5. Practical regularization proposal for the codebase

### What works now

- **Collision bounce** (r ~ 0, ell = 0): elastic reflection `q_rel -> -q_rel`
  at a small radius. Implemented in the codebase. Works well with the
  unregularized symplectic integrator (bounded energy error).

- **Levi-Civita/KS regularization** (r ~ 0, Coulomb/Kepler only): handles
  close encounters in the super-critical regime. Does NOT regularize Weber's
  velocity-dependent force.

### Proposed additions

**1. Weber-aware radial regularization for ell = 0 through rho.**

For purely radial trajectories, the energy equation `r_dot^2 = f(r)` can be
integrated directly without the Hamiltonian formulation. The singularity at
`r = rho` is removable (0/0 by L'Hopital at `h = h*`).

Implementation sketch:
```julia
function weber_radial_step!(r, v, dt, rho, h)
    if abs(r - rho) < epsilon
        # Use L'Hopital limit: r_dot^2 = 2h
        v_new = sign(v) * sqrt(2h)
        r_new = r + v_new * dt
    else
        # Normal integration
        v2 = (2*e2_mu - 2*h*r) / (rho - r)
        v_new = sign(v) * sqrt(abs(v2))
        r_new = r + v_new * dt
    end
    return r_new, v_new
end
```

**Caveat:** This only works for ell = 0. Any angular momentum, however small,
makes the critical sphere non-traversable (the ell^2 term in the numerator
prevents the 0/0 cancellation).

**2. Metric perturbation (approximate regularization).**

Replace the Weber inverse metric eigenvalue `(1/mu)(1 - rho/r)` with a
regularized version that stays positive:

```
g_reg(r) = (1/mu)(1 - rho/r + delta^2)
```

for small `delta > 0`. This keeps the metric Riemannian everywhere, at the
cost of introducing an `O(delta^2)` error in the dynamics. The error is
localized near `r = rho` and vanishes as `delta -> 0`.

**Pros:** All existing machinery (symplectic integrator, LC/KS regularization,
Floer theory) works unchanged. Contact-type holds globally.

**Cons:** The regularized system is NOT Weber electrodynamics. Sub-critical
bound states (which depend on the metric sign flip) do not exist in the
regularized system. The perturbation changes the qualitative dynamics.

**Recommendation:** The metric perturbation is useful for studying
super-critical dynamics near `rho` (avoiding numerical stiffness) but should
NOT be used to study sub-critical physics. For sub-critical ell = 0, the
collision bounce remains the correct approach. For sub-critical ell != 0,
no regularization exists (Frauenfelder-Weber 2024, Theorem 2.1).

**3. Detection and early warning.**

The most practical addition to the codebase: detect when a like-charge pair
approaches `rho` and estimate the time to reach it:

```julia
function time_to_critical(r, rdot, rho)
    if r > rho && rdot < 0
        # Linear estimate
        return (r - rho) / abs(rdot)
    end
    return Inf
end
```

This allows the integrator to:
- Reduce the timestep near `rho` (stiffness handling)
- Warn the user that the metric is about to degenerate
- Switch to the collision-bounce regime if `r` crosses `rho` with ell ~ 0

## 6. Summary of findings

| Question | Answer |
|----------|--------|
| McGehee blow-up at rho? | Does NOT regularize (metric singularity, not potential) |
| ell = 0 through rho? | YES, but only at critical energy h* = c^2 |
| ell != 0 through rho? | NO (infinite winding, topological obstruction) |
| Critical sphere as "stop"? | No -- codimension 1, not codimension 2 |
| Critical sphere classification | Signature-change hypersurface (degenerate contact) |
| Excised RFH well-defined? | YES, for super-critical domain with epsilon-buffer |
| Pseudo-Riemannian Floer theory? | Does not exist (infinite Morse index) |
| Practical regularization? | Collision bounce (ell=0 sub-critical), metric perturbation (approximate), early warning (detection) |

## References

- Frauenfelder, U., Weber, J. "A mathematical description of the Weber nucleus."
  Anal. Math. Phys. 14:31 (2024).
- McGehee, R. "Triple collision in the collinear three-body problem."
  Inventiones math. 27, 191-227 (1974).
- Sylvan, Z. "On partially wrapped Fukaya categories." J. Topol. 12 (2019), 372-441.
- Ganatra, S., Pardon, J., Shende, V. "Covariantly functorial wrapped Floer
  theory on Liouville sectors." Publ. Math. IHES 131 (2020), 73-200.
- Abbondandolo, A., Schwarz, M. "On the Floer homology of cotangent bundles."
  Comm. Pure Appl. Math. 59 (2006), 254-316.
- Suhr, S. "Theory of geodesics on Lorentzian manifolds." Various 2013-2018.
- Cieliebak, K., Frauenfelder, U. "A Floer homology for exact contact embeddings."
  Pacific J. Math. 239 (2009), 251-316.
- Viterbo, C. "A proof of Weinstein's conjecture in R^{2n}." Ann. IHP 4 (1987).
- Hofer, H. "Pseudoholomorphic curves in symplectisations." Invent. Math. 114 (1993).

## Files

- `NOTES.md` -- this file
- `mcgehee_weber.md` -- detailed McGehee blow-up calculation
- `l0_regularization.jl` -- Julia script testing ell=0 regularization
