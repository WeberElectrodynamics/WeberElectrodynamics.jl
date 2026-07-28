# The Weber Hamiltonian

This page states the canonical formulation the package implements, and the one
rule that follows from it: **canonical momentum is not mass times velocity.**

## Canonical momentum

The Weber Lagrangian in absolute (Gauss–Weber) units is

```math
L(\vec r, \vec v) = \sum_i \tfrac12 m_i \lVert\vec v_i\rVert^2
- \sum_{i<j} \frac{q_i q_j}{r_{ij}}\left(1 + \frac{\dot r_{ij}^2}{2c^2}\right)
```

where ``m_i`` and ``q_i`` are the mass and charge of particle ``i``, ``r_{ij}``
is the pair separation, ``\dot r_{ij}`` its time derivative, and ``c`` the speed
of light. It depends on velocity through both terms, so

```math
\vec p_i = \frac{\partial L}{\partial \vec v_i}
= m_i \vec v_i - \sum_{j\ne i} \frac{q_i q_j}{c^2}
  \frac{\dot r_{ij}}{r_{ij}^2}(\vec r_i - \vec r_j).
```

The correction is radial and depends on the pair's physical radial velocity. It
vanishes only when every ``\dot r_{ij} = 0``, or in the Coulomb limit
``c \to \infty``. Everywhere else ``\vec p_i \ne m_i \vec v_i``.

## Recovering physical velocities

Writing ``k_{ij} = q_i q_j/(c^2 r_{ij})`` and

```math
s_{ij} = \hat r_{ij}\cdot\left(\frac{\vec p_i}{m_i} - \frac{\vec p_j}{m_j}\right)
```

for the *naive* radial rate (what ``\dot r_{ij}`` would be if momentum were
kinetic), the physical radial velocities solve an ``n(n-1)/2`` dimensional
linear system coupling every pair that shares a particle:

```math
\dot r_{ab} = s_{ab} + \sum_{(c,d)} G_{ab,cd}\, k_{cd}\, \dot r_{cd},
\qquad
G_{ab,cd} = (\hat r_{ab}\cdot\hat r_{cd})
\left[\frac{\delta_{ac}-\delta_{ad}}{m_a} - \frac{\delta_{bc}-\delta_{bd}}{m_b}\right].
```

The velocities then follow explicitly:

```math
\vec v_i = \frac{\vec p_i}{m_i}
+ \frac{1}{m_i}\sum_{j\ne i} k_{ij}\,\dot r_{ij}\,\hat r_{ij}.
```

For two particles the system is one-dimensional and collapses to a division. With
``\mu = m_1m_2/(m_1+m_2)`` and ``p_r = \mu\, \hat r \cdot (\vec p_1/m_1 - \vec p_2/m_2)``:

```math
p_r = \left(\mu - \frac{q_1q_2}{r c^2}\right)\dot r.
```

Use [`physical_velocities`](@ref) — never `p ./ masses`:

```julia
sol = solve(prob, SymmetricProjectionIntegrator())
v   = physical_velocities(sol.prob, sol.q[end], sol.p[end])
```

## The canonical Hamiltonian

The Legendre transform gives

```math
H(\vec r, \vec p) = \tfrac12 \sum_i \vec p_i \cdot \vec v_i + \sum_{i<j}\frac{q_iq_j}{r_{ij}}
= \sum_i \frac{\lVert\vec p_i\rVert^2}{2m_i}
+ \tfrac12\sum_{i<j} k_{ij}\,\dot r_{ij}\,s_{ij}
+ \sum_{i<j}\frac{q_iq_j}{r_{ij}}
```

which equals the velocity-space energy
``\sum_i \tfrac12 m_i\lVert\vec v_i\rVert^2 + \sum_{i<j} U_{ij}``. The canonical
equations are

```math
\dot{\vec r}_i = \vec v_i(\vec r, \vec p),
\qquad
\dot{\vec p}_i = \sum_{j\ne i}\frac{q_iq_j}{r_{ij}^2}
\left[\hat r_{ij}\left(1 + \frac{3\dot r_{ij}^2}{2c^2}\right)
- \frac{\dot r_{ij}}{c^2}(\vec v_i - \vec v_j)\right].
```

Canonical momentum rate is not the mechanical Weber force. Since ``\vec p_i``
carries the velocity-dependent correction ``\vec\alpha_i``, the acceleration
obeys ``m_i \dot{\vec v}_i = \dot{\vec p}_i + \mathrm{d}\vec\alpha_i/\mathrm{d}t``,
which reproduces Weber's force law.

Full derivations: [`theory/WeberElectrodynamics.md`](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/WeberElectrodynamics.md).
Independent symbolic verification of every equation above:
[`papers/Computational-Weber-Electrodynamics/verify_formulas.py`](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/papers/Computational-Weber-Electrodynamics/verify_formulas.py).

## Weber's critical radius

For a like-charge pair the effective radial inertia ``\mu - q_1q_2/(rc^2)``
vanishes at

```math
\rho = \frac{q_1 q_2}{\mu c^2}
```

where the canonical mass matrix is singular and a canonical momentum no longer
determines a finite physical velocity. Evaluating there raises
[`WeberCriticalRadiusError`](@ref) rather than returning silent garbage.

Below ``\rho`` the effective inertia is negative but finite: the pair responds
to the still-repulsive Coulomb force as if it had negative inertia and
accelerates *towards* the other charge. This sub-critical regime is well defined
and integrable; only the crossing itself is singular.

## Practical consequences

| Quantity | Correct | Wrong |
|---|---|---|
| Physical velocity | `physical_velocities(prob, q, p)` | `p ./ masses(prob)` |
| Kinetic energy | ``\sum_i \tfrac12 m_i\lVert\vec v_i\rVert^2`` | ``\sum_i \lVert\vec p_i\rVert^2/(2m_i)`` |
| Pair radial velocity | ``\dot r_{ij}`` from the solve | ``s_{ij}`` |
| Radial initial condition | ``p_r = (\mu - q_1q_2/(rc^2))\dot r`` | ``p_r = \mu\dot r`` |

`EnergyData.kinetic_energy`, `PairEnergyData.radial_velocity`, and the
`PairForceData` velocity/acceleration decompositions all report the physical
quantities. Pass `c` to [`two_body_initial_conditions`](@ref) whenever
`radial_velocity` is nonzero.

## API

[`physical_velocities`](@ref), [`WeberCriticalRadiusError`](@ref), and
[`has_symbolic_hamiltonian`](@ref) are documented on the
[System](@ref) API page.
