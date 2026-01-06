# Computational Examples Specification

Three numerical examples for the Weber electrodynamics paper.

---

## Example 1: Two-Body Elliptical Orbit with Precession

Weber's velocity-dependent corrections cause apsidal precession in bound orbits.

### Parameters

| Parameter | Value |
|-----------|-------|
| $m_1$ | 1.0 |
| $m_2$ | 0.1 |
| $q_1$ | +1.0 |
| $q_2$ | −1.0 |
| $k$ (coupling) | 0.1 |
| $r_0$ | 2.0 |
| $c$ | 4.0 |

### Initial Conditions (center-of-mass frame)
With $M = m_1 + m_2 = 1.1$ and $v_{\text{circ}} = \sqrt{k|q_1 q_2| M / (m_1 m_2 r_0)} \approx 0.74$:

- Particle 1 at $(-m_2 r_0/M, 0) \approx (-0.18, 0)$, velocity $(0, -m_2 v_{\text{scale}} v_{\text{circ}}/M)$
- Particle 2 at $(+m_1 r_0/M, 0) \approx (+1.82, 0)$, velocity $(0, +m_1 v_{\text{scale}} v_{\text{circ}}/M)$
- $v_{\text{scale}} = 0.7$ (elliptical orbit)

### Integration
- Time span: 50 orbital periods (~500 time units)
- Step size: $\Delta t = 0.01$

### Output
1. Trajectory plot showing rosette pattern
2. Comparison with Coulomb ($c \to \infty$) showing closed ellipse

---

## Example 2: Two-Body Near-Relativistic Dynamics

Explore Weber dynamics when $v/c \to 1$.

### Parameters

| Parameter | Value |
|-----------|-------|
| $m_1, m_2$ | 1.0 |
| $q_1$ | +1.0 |
| $q_2$ | −1.0 |
| $k$ (coupling) | 1.0 |
| $r_0$ | 1.0 |
| $c$ | 1.0 |

### Initial Conditions (center-of-mass frame, three configurations)
For equal masses, each particle is at $(\pm r_0/2, 0) = (\pm 0.5, 0)$. Circular orbit: $v_{\text{rel}} = \sqrt{2k|q_1 q_2|/(m \cdot r_0)} \approx 1.41$, so each particle has $v \approx 0.71$.

1. **Circular attempt**: Tangential velocity $v = 0.71$ each (relative $v_{\text{rel}} = 1.41$)
2. **Radial approach**: Inward radial velocity $v = 0.25$ each (relative $0.5c$)
3. **High velocity**: Tangential velocity $v = 0.45$ each (relative $0.9c$)

### Integration
- Time span: 10–50 time units (adjust per configuration)
- Step size: $\Delta t = 0.001$
- Max Newton iterations: 200

### Output
1. Trajectory plots for each configuration
2. Phase space portraits $(r, \dot{r})$
3. Energy time series

---

## Example 3: Standing Wave in 2D Particle Grid

Longitudinal oscillations in a lattice of charged particles.

### Physical Setup

```
Column:  1    2    3    4    5    6    7    8    9   10   11   12
Row 1:   ●────○────○────○────○────○────○────○────○────○────○────●
Row 2:   ●────○────○────○────○────○────○────○────○────○────○────●
Row 3:   ●────○────○────○────○────○────○────○────○────○────○────●

● = Fixed particle (columns 1 and 12)
○ = Mobile particle (columns 2–11)
```

- Total: 36 particles (6 fixed, 30 mobile)

### Parameters

| Parameter | Value |
|-----------|-------|
| Lattice spacing $a$ | 1.0 |
| $q_i$ | +1.0 (all positive, repulsive) |
| $m_i$ | 1.0 |
| $k$ (coupling) | 1.0 |
| $c$ | $\infty$ (Coulomb limit) |
| Perturbation $\delta$ | 0.05 |

### Initial Conditions

**Equilibrium positions**:
```julia
for row in 1:3, col in 1:12
    x[row, col] = (col - 1) * a
    y[row, col] = (row - 1) * a
end
```

**Perturbation** (fundamental mode):
```julia
δ = 0.05
for row in 1:3
    x[row, 2] += δ    # Column 2: displaced right
    x[row, 11] -= δ   # Column 11: displaced left
end
```

**Initial momenta**: Zero.

**Fixed particles**: Include columns 1 and 12 in force calculation but do not update their positions.

### Integration
- Time span: 100 time units (~15 oscillation periods)
- Step size: $\Delta t = 0.01$

### Output
1. Snapshots at $t = 0, T/4, T/2, 3T/4, T$
2. Displacement $x_i(t) - x_i^{\text{eq}}$ vs $t$ for select particles
3. Energy vs time

---

## Summary

| Example | Bodies | Key Parameter | Phenomenon |
|---------|--------|---------------|------------|
| 1. Precession | 2 | $v/c \approx 0.2$ | Apsidal precession |
| 2. Relativistic | 2 | $v/c \to 1$ | High-velocity dynamics |
| 3. Standing wave | 36 | Lattice grid | Collective oscillation |
