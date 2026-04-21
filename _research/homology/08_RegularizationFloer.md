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


---

## McGehee-Weber

# McGehee Blow-Up at the Weber Critical Sphere

Agent 08 deliverable. Cross-refs: Agent 11 (contact/Reeb, section 3),
AngularMomentumRegularization.md (Experiment 4), TetheringImpossibility.md.

## 1. Setup: 2-body Weber Hamiltonian in relative coordinates

Two like charges with reduced mass mu, charges e1 = e2 = e > 0,
critical radius rho = e^2 / (mu c^2). In polar relative coordinates (r, phi):

```
H = (r / (2 mu (r - rho))) p_r^2  +  ell^2 / (2 mu r^2)  +  e^2 / r
```

The radial kinetic coefficient is `1/(2 mu) * r/(r - rho)`, which has the
eigenvalue `(1/mu)(1 - rho/r)` along the radial direction. This changes sign
at r = rho.

The Frauenfelder-Weber energy equation:

```
r_dot^2 = (ell^2 + 2(e^2/mu) r - 2 h r^2) / (r (rho - r))     for r < rho
```

## 2. McGehee blow-up of the collision singularity (r = 0)

The classical McGehee (1974) blow-up addresses the r = 0 singularity. Define:

```
s = sqrt(r),     r = s^2,     s >= 0
w = s * s_dot = s * (r_dot / (2s)) = r_dot / 2
```

Then `r_dot = 2w` and the energy equation becomes:

```
4 w^2 = (ell^2 + 2(e^2/mu) s^2 - 2 h s^4) / (s^2 (rho - s^2))
```

For ell = 0 as s -> 0: `4 w^2 -> 2(e^2/mu) / rho = 2 c^2`, so w -> +/- c/sqrt(2).
The collision manifold {s = 0} is invariant, and the flow extends smoothly.

For ell != 0 as s -> 0: `4 w^2 ~ ell^2 / (rho s^2) -> infinity`. Both w and
the angular velocity omega = ell/s^2 diverge. The collision manifold is NOT
invariant -- this is the topological obstruction (infinite winding).

## 3. McGehee blow-up of the critical sphere (r = rho)

This is the NEW calculation. We blow up the metric degeneracy at r = rho,
not the collision at r = 0. Define:

```
s = sqrt(r - rho),     r = rho + s^2,     s >= 0  (super-critical side)
s = sqrt(rho - r),     r = rho - s^2,     s >= 0  (sub-critical side)
```

### Super-critical side (r = rho + s^2, s >= 0)

The radial eigenvalue:

```
(1/mu)(1 - rho/r) = (1/mu)(s^2 / (rho + s^2))
```

The effective radial kinetic energy:

```
T_r = p_r^2 / (2 mu) * (r / (r - rho)) = p_r^2 / (2 mu) * (rho + s^2) / s^2
```

This diverges as s -> 0. To remove the divergence, define the McGehee-type
momentum:

```
P_s = 2 s p_r     (from the canonical relation p_r dr = P_s ds, since dr = 2s ds)
```

so `p_r = P_s / (2s)`. Then:

```
T_r = P_s^2 / (8 mu s^2) * (rho + s^2) / s^2 = P_s^2 (rho + s^2) / (8 mu s^4)
```

This is WORSE -- it diverges as s^{-4} rather than s^{-2}. The blow-up
does not remove the kinetic singularity.

### Sundman time rescaling

Apply dt = s^2 d(tau) to absorb the s^{-2} from the original metric. The
extended Hamiltonian at energy h:

```
K = s^2 (H - h) = s^2 [ p_r^2 (rho + s^2) / (2 mu s^2) + ell^2/(2 mu (rho+s^2)^2) + e^2/(rho+s^2) - h ]
  = p_r^2 (rho + s^2) / (2 mu) + s^2 [ ell^2/(2 mu (rho+s^2)^2) + e^2/(rho+s^2) - h ]
```

In terms of P_s = 2 s p_r:

```
K = P_s^2 (rho + s^2) / (8 mu s^2) + s^2 [ ... ]
```

The first term still diverges as s -> 0. The Sundman factor s^2 cancels one
power of s^2 from the metric, but the coordinate change s = sqrt(r - rho)
introduces another. Net result: the singularity persists.

### Why the critical sphere is worse than a collision

At r = 0 (collision), the kinetic energy scales as p_r^2 / (2 mu) with
a Coulomb 1/r potential. The Sundman factor g = r (equivalently u^2 in
LC coordinates) absorbs the 1/r singularity, giving a regular Hamiltonian.

At r = rho (critical sphere), the kinetic coefficient itself diverges:
`r/(r - rho) -> infinity`. This is a singularity IN the metric, not in
the potential. The metric divergence is of order (r - rho)^{-1}, which
is the same order as a simple pole, but it multiplies p_r^2 rather than
appearing additively. No Sundman factor can simultaneously:
1. Absorb the metric pole (requires g ~ (r - rho))
2. Keep the angular dynamics regular (requires g ~ r^2 for ell != 0)
3. Reach the critical sphere in finite fictitious time

### Boundary manifold structure

Despite the blow-up failing to regularize, it reveals useful structure.
The blown-up phase space at s = 0 (the critical sphere) has:

```
Boundary manifold: S^1 x R  (parameterized by phi and P_s)
```

On the super-critical side, as s -> 0+:
- The angular dynamics (phi, ell) remain regular (ell is conserved, phi evolves smoothly)
- The radial momentum p_r has a definite limit (can be zero)
- The radial velocity r_dot = s^2 p_r / (mu * rho) -> 0 (the particle takes infinite time to reach rho)

On the sub-critical side, as s -> 0+ (r -> rho-):
- Same boundary manifold S^1 x R
- But the SIGN of the radial kinetic energy is reversed
- Trajectories approaching rho from below have p_r -> infinity (since the effective mass -> 0)

The critical sphere is a **degenerate** boundary: the metric tensor has
a zero eigenvalue there. It is neither a regular boundary nor a collision
singularity -- it is a **signature change hypersurface**.

## 4. Explicit formulas for the 2-body case

### Hamiltonian in blown-up coordinates (super-critical side)

With r = rho + s^2, phi, and conjugate momenta P_s, ell:

```
H(s, phi, P_s, ell) = P_s^2 (rho + s^2) / (8 mu s^2)
                     + ell^2 / (2 mu (rho + s^2)^2)
                     + e^2 / (rho + s^2)
```

Hamilton's equations:

```
ds/dt     = dH/dP_s = P_s (rho + s^2) / (4 mu s^2)
dP_s/dt   = -dH/ds  = P_s^2 rho / (4 mu s^3) - P_s^2/(4 mu s)
                     + ell^2 s / (mu (rho + s^2)^3)
                     + e^2 s / ((rho + s^2)^2) * 2
dphi/dt   = dH/dell = ell / (mu (rho + s^2)^2)
dell/dt   = -dH/dphi = 0
```

Note: ds/dt diverges as s -> 0 unless P_s -> 0 fast enough. Specifically,
ds/dt ~ P_s rho / (4 mu s^2), so for ds/dt to remain bounded we need
P_s = O(s^2), i.e., p_r = P_s/(2s) = O(s) -> 0.

This means: **only trajectories with vanishing radial momentum can reach
the critical sphere in finite time from the super-critical side.**

### Sub-critical Hamiltonian

With r = rho - s^2 (s >= 0, sub-critical):

```
H(s, phi, P_s, ell) = -P_s^2 (rho - s^2) / (8 mu s^2)
                     + ell^2 / (2 mu (rho - s^2)^2)
                     + e^2 / (rho - s^2)
```

The kinetic term is NEGATIVE (Lorentzian signature). The boundary
at s = 0 again has the metric eigenvalue vanishing.

## 5. Comparison with the AngularMomentumRegularization.md McGehee experiment

The earlier Experiment 4 (McGehee blow-up at r = 0) found:
- ell = 0: w -> -sqrt(2)/2 (bounded), collision manifold invariant
- ell != 0: both w and omega diverge, manifold not invariant

The McGehee blow-up at r = rho (this document) finds an analogous but
distinct structure:
- The critical sphere is not a collision but a metric degeneracy
- Even for ell = 0, the blow-up does not regularize the Hamiltonian
  (the kinetic singularity persists)
- The critical sphere acts as an accumulation boundary: trajectories
  approach it asymptotically but (generically) do not cross it

## 6. Conclusion

The McGehee blow-up at the critical sphere r = rho does NOT yield a
regularized Hamiltonian. The fundamental reason is that the singularity
is in the kinetic metric (a zero eigenvalue), not in the potential
(a pole). Standard regularization techniques (Sundman, Levi-Civita, KS)
are designed for potential singularities and fail here.

The blown-up phase space has a boundary manifold S^1 x R at s = 0,
but the Hamiltonian does not extend continuously to this boundary
(the kinetic term diverges). This confirms the picture from Agent 11
section 3: the critical sphere is a true obstruction, not removable
by coordinate changes.

The one exception is the ell = 0 radial case, where the problematic
kinetic term is the only degree of freedom and can be treated by
passing to a different canonical structure (see l0_regularization.jl).
