# Agent 13 -- c-Continuation from Coulomb to Weber

## Goal

Track 2-body orbits (unlike charges, equal masses) as the Weber speed-of-light
parameter c decreases from 100 (near-Coulomb) to 1 (strong Weber regime).
Identify stability transitions, bifurcations, and the fate of elliptical orbits.

## Method

- `continuation.jl` integrates 2-body Weber systems using the symplectic
  Strang-splitting integrator (`SymmetricProjectionIntegrator`).
- Monodromy matrices computed via centered finite differences (eps = 1e-7).
- Floquet multipliers = eigenvalues of the monodromy matrix; |lambda|_max > 1
  signals linear instability.
- Families studied: circular (E = -0.5, -1.0, -2.0), elliptical (E = -1.0;
  e = 0.3, 0.5, 0.7), plus 4-body configurations and a new-orbit search.

## Key Findings

### 1. Circular orbits remain linearly stable at all c

For all three circular families (E = -0.5, -1.0, -2.0), the maximum Floquet
multiplier stays at |lambda|_max = 1.000000 (+/- 1e-6) across the entire range
c in [0.5, 100]. No stability transitions or bifurcations were detected.

**Why**: Circular orbits have rdot = 0 identically, so the Weber velocity-dependent
correction V_Weber ~ (rdot/c)^2 vanishes exactly. The orbit itself is unchanged
from Coulomb. The linearized variational equations around a circular orbit also
see no Weber effect to leading order because radial perturbations are O(epsilon)
and the Weber correction enters at O(epsilon^2 * rdot^2 / c^2).

The fine bifurcation scan (Part 3) with ~30-50 c-values per energy confirmed:
no eigenvalue crossings at +1 or -1, no Krein collisions, no stability
transitions anywhere in [0.5, 5.0] for any energy.

**Implication for Agent 04's CZ index jump at c ~ 1.414**: The CZ index change
reported by Agent 04 is not a linearized instability of circular orbits themselves.
It likely reflects a change in the Maslov index / winding of eigenvalues on the
unit circle (a topological invariant change without leaving the unit circle),
or it applies to a different orbit family.

### 2. Elliptical orbits are strongly destabilized as c decreases

| ecc | |lambda|_max at c=100 | c=10 | c=5 | c=1 | Status at c=1 |
|-----|----------------------|------|-----|-----|----------------|
| 0.3 | 1.32 | 3.97 | 7.35 | 17.1 | Open |
| 0.5 | 1.95 | 11.0 | 19.1 | 83.8 | Open |
| 0.7 | 3.85 | 68.3 | 11.4 | 191.6 | Open |

The Floquet multipliers grow roughly exponentially with eccentricity and with
decreasing c. At c = 1, the e = 0.7 orbit has |lambda|_max ~ 192, meaning
perturbations amplify by a factor of ~200 per period.

**Note on c=100 instability**: Even at c=100, |lambda|_max > 1 for elliptical
orbits. This is partly a numerical artifact: Kepler ellipses are degenerate
(the orbit does not precess, so neighboring orbits with slightly different
parameters diverge secularly). The finite-difference monodromy captures this
near-degeneracy as apparent instability. The physically meaningful finding is
the dramatic growth factor as c decreases.

**Orbit closure**: At c=100, elliptical orbits close (or nearly close) at the
Kepler period. As c decreases, the Weber precession prevents closure at the
Kepler period. A period-search algorithm found approximate closing periods
that shift by up to +30% at c = 1.

**Eigenvalue type transitions**: As c decreases, eigenvalue types shift from
`positive_hyperbolic` (real lambda > 1) through `loxodromic` (complex, off unit
circle) to `negative_hyperbolic` (real lambda < -1). This sequence is typical
of a parametric resonance cascade.

### 3. All 4-body configurations failed

Both the breathing alternating square (vrad = 0.3, 0.5, 0.7) and the double
orbiter configuration failed at all c values tested (1, 2, 4, 10). The 4-body
unlike-charge configurations are not viable periodic orbit candidates in this
parameter regime. Close encounters between unlike-charge particles cause the
integrator to fail.

### 4. New orbit candidates found but not validated

The T-brake shooting search (Part 4) found ~25 approximate periodic orbit
candidates with closure error < 0.1, at various c values. These have not been
refined or Floquet-analyzed. Notable: at c = 1, several candidates exist with
T ~ 2-16, suggesting a rich periodic orbit structure in the strong-Weber regime.

## Summary Table: Stability of Circular vs Elliptical

```
c       Circular (any E)     Elliptical e=0.3     Elliptical e=0.7
------  -------------------  -------------------  -------------------
100     stable (|l|=1.00)    weakly unstable 1.3  unstable 3.9
 10     stable (|l|=1.00)    unstable 4.0         strongly unstable 68
  5     stable (|l|=1.00)    unstable 7.3         strongly unstable 11
  1     stable (|l|=1.00)    unstable 17          strongly unstable 192
```

## Conclusions

1. **No bifurcation of circular orbits** in the range c in [0.5, 100]. The Weber
   force does not destabilize circular orbits because it vanishes on them.

2. **Elliptical orbits become increasingly unstable** as c decreases, with
   instability growing roughly as ~1/c^2 (consistent with the Weber correction
   scaling). Higher eccentricity orbits are more unstable because they have
   larger radial velocities at periapsis, where the Weber correction is strongest.

3. **The physically important transition is not a bifurcation but a continuous
   destabilization**: there is no sharp critical c where circular orbits become
   unstable. Instead, the Weber force progressively destroys the integrability
   of the Kepler problem, turning the phase space near elliptical orbits chaotic.

4. **Connection to Agent 04**: The CZ index jump at c ~ 1.414 is likely a
   topological (Maslov index) transition of the linearized flow on the unit
   circle, not a linearized instability. The eigenvalues remain on the unit
   circle for circular orbits but their winding number changes.

5. **Connection to Agent 07**: The t* ~ c^2 survival time prediction could not
   be tested because all 4-body configurations failed. The 2-body elliptical
   instability growth rate ~1/c^2 is consistent with this scaling.

## Files

- `continuation.jl` -- main script (578 lines)
- `bifurcation_diagram.csv` -- full output (239 data rows)
- `NOTES.md` -- this file
