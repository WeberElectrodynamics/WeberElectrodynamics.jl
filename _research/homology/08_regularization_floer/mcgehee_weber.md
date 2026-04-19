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
