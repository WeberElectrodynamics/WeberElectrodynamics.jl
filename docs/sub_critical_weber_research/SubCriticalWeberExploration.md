# Sub-Critical Weber Exploration Report

_Generated: 2026-03-14 21:59_

This report systematically analyzes the angular momentum regularization barrier
in Weber electrodynamics below the critical radius ρ. Two like charges inside ρ
exhibit negative effective inertial mass μ_eff = μ(1 − ρ/r), creating bound
oscillatory states. For angular momentum ℓ ≠ 0, particles spiral into collision
at infinite speed with infinite winding — a topological obstruction that no smooth
coordinate transform can remove.

**Goal**: Identify the most promising research directions to overcome or circumvent
this barrier, using symbolic analysis, numerical experiments, and literature survey.


## The Weber Hamiltonian


### Cartesian Form (Primary — Relational Mechanics)

Following Assis's relational mechanics, we retain individual particle masses m₁, m₂
and charges q₁, q₂:

```
H = p₁²/(2m₁) + p₂²/(2m₂) + q₁q₂/r · (1 − ṙ²/(2c²))
```

where ṙ = (r⃗ · v⃗)/r with r⃗ = r⃗₁ − r⃗₂ and v⃗ = p⃗₁/m₁ − p⃗₂/m₂.

**Hamilton's equations (2D):**

Velocities (∂H/∂p):
```
  d(x1)/dt = (m1*px2*q1*q2*(x1^2) - (2//1)*m1*px2*q1*q2*x1*x2 + m1*px2*q1*q2*(x2^2) + m1*py2*q1*q2*x1*y1 - m1*py2*q1*q2*x1*y2 - m1*py2*q1*q2*x2*y1 + m1*py2*q1*q2*x2*y2 - m2*px1*q1*q2*(x1^2) + (2//1)*m2*px1*q1*q2*x1*x2 - m2*px1*q1*q2*(x2^2) - m2*py1*q1*q2*x1*y1 + m2*py1*q1*q2*x1*y2 + m2*py1*q1*q2*x2*y1 - m2*py1*q1*q2*x2*y2 + (c^2)*m1*m2*px1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3)) / ((c^2)*(m1^2)*m2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3))
  d(y1)/dt = (m1*px2*q1*q2*x1*y1 - m1*px2*q1*q2*x1*y2 - m1*px2*q1*q2*x2*y1 + m1*px2*q1*q2*x2*y2 + m1*py2*q1*q2*(y1^2) - (2//1)*m1*py2*q1*q2*y1*y2 + m1*py2*q1*q2*(y2^2) - m2*px1*q1*q2*x1*y1 + m2*px1*q1*q2*x1*y2 + m2*px1*q1*q2*x2*y1 - m2*px1*q1*q2*x2*y2 - m2*py1*q1*q2*(y1^2) + (2//1)*m2*py1*q1*q2*y1*y2 - m2*py1*q1*q2*(y2^2) + (c^2)*m1*m2*py1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3)) / ((c^2)*(m1^2)*m2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3))
  d(x2)/dt = (-m1*px2*q1*q2*(x1^2) + (2//1)*m1*px2*q1*q2*x1*x2 - m1*px2*q1*q2*(x2^2) - m1*py2*q1*q2*x1*y1 + m1*py2*q1*q2*x1*y2 + m1*py2*q1*q2*x2*y1 - m1*py2*q1*q2*x2*y2 + m2*px1*q1*q2*(x1^2) - (2//1)*m2*px1*q1*q2*x1*x2 + m2*px1*q1*q2*(x2^2) + m2*py1*q1*q2*x1*y1 - m2*py1*q1*q2*x1*y2 - m2*py1*q1*q2*x2*y1 + m2*py1*q1*q2*x2*y2 + (c^2)*m1*m2*px2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3)) / ((c^2)*m1*(m2^2)*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3))
  d(y2)/dt = (-m1*px2*q1*q2*x1*y1 + m1*px2*q1*q2*x1*y2 + m1*px2*q1*q2*x2*y1 - m1*px2*q1*q2*x2*y2 - m1*py2*q1*q2*(y1^2) + (2//1)*m1*py2*q1*q2*y1*y2 - m1*py2*q1*q2*(y2^2) + m2*px1*q1*q2*x1*y1 - m2*px1*q1*q2*x1*y2 - m2*px1*q1*q2*x2*y1 + m2*px1*q1*q2*x2*y2 + m2*py1*q1*q2*(y1^2) - (2//1)*m2*py1*q1*q2*y1*y2 + m2*py1*q1*q2*(y2^2) + (c^2)*m1*m2*py2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3)) / ((c^2)*m1*(m2^2)*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^3))
```

Forces (−∂H/∂q):
```
  d(px1)/dt = (-(2//1)*(m1^2)*(px2^2)*q1*q2*(x1^3) + (6//1)*(m1^2)*(px2^2)*q1*q2*(x1^2)*x2 - (6//1)*(m1^2)*(px2^2)*q1*q2*x1*(x2^2) + (2//1)*(m1^2)*(px2^2)*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m1^2)*(px2^2)*q1*q2*(x2^3) - (2//1)*(m1^2)*(px2^2)*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*(m1^2)*px2*py2*q1*q2*(x1^2)*y1 + (4//1)*(m1^2)*px2*py2*q1*q2*(x1^2)*y2 + (8//1)*(m1^2)*px2*py2*q1*q2*x1*x2*y1 - (8//1)*(m1^2)*px2*py2*q1*q2*x1*x2*y2 - (4//1)*(m1^2)*px2*py2*q1*q2*(x2^2)*y1 + (4//1)*(m1^2)*px2*py2*q1*q2*(x2^2)*y2 + (2//1)*(m1^2)*px2*py2*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m1^2)*px2*py2*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m1^2)*(py2^2)*q1*q2*x1*(y1^2) + (4//1)*(m1^2)*(py2^2)*q1*q2*x1*y1*y2 - (2//1)*(m1^2)*(py2^2)*q1*q2*x1*(y2^2) + (2//1)*(m1^2)*(py2^2)*q1*q2*x2*(y1^2) - (4//1)*(m1^2)*(py2^2)*q1*q2*x2*y1*y2 + (2//1)*(m1^2)*(py2^2)*q1*q2*x2*(y2^2) + (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*x1 - (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*x2 + (4//1)*m1*m2*px1*px2*q1*q2*(x1^3) - (12//1)*m1*m2*px1*px2*q1*q2*(x1^2)*x2 + (12//1)*m1*m2*px1*px2*q1*q2*x1*(x2^2) - (4//1)*m1*m2*px1*px2*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*px1*px2*q1*q2*(x2^3) + (4//1)*m1*m2*px1*px2*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*px1*py2*q1*q2*(x1^2)*y1 - (4//1)*m1*m2*px1*py2*q1*q2*(x1^2)*y2 - (8//1)*m1*m2*px1*py2*q1*q2*x1*x2*y1 + (8//1)*m1*m2*px1*py2*q1*q2*x1*x2*y2 + (4//1)*m1*m2*px1*py2*q1*q2*(x2^2)*y1 - (4//1)*m1*m2*px1*py2*q1*q2*(x2^2)*y2 - (2//1)*m1*m2*px1*py2*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*m1*m2*px1*py2*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*px2*py1*q1*q2*(x1^2)*y1 - (4//1)*m1*m2*px2*py1*q1*q2*(x1^2)*y2 - (8//1)*m1*m2*px2*py1*q1*q2*x1*x2*y1 + (8//1)*m1*m2*px2*py1*q1*q2*x1*x2*y2 + (4//1)*m1*m2*px2*py1*q1*q2*(x2^2)*y1 - (4//1)*m1*m2*px2*py1*q1*q2*(x2^2)*y2 - (2//1)*m1*m2*px2*py1*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*m1*m2*px2*py1*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*py1*py2*q1*q2*x1*(y1^2) - (8//1)*m1*m2*py1*py2*q1*q2*x1*y1*y2 + (4//1)*m1*m2*py1*py2*q1*q2*x1*(y2^2) - (4//1)*m1*m2*py1*py2*q1*q2*x2*(y1^2) + (8//1)*m1*m2*py1*py2*q1*q2*x2*y1*y2 - (4//1)*m1*m2*py1*py2*q1*q2*x2*(y2^2) - (2//1)*(m2^2)*(px1^2)*q1*q2*(x1^3) + (6//1)*(m2^2)*(px1^2)*q1*q2*(x1^2)*x2 - (6//1)*(m2^2)*(px1^2)*q1*q2*x1*(x2^2) + (2//1)*(m2^2)*(px1^2)*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m2^2)*(px1^2)*q1*q2*(x2^3) - (2//1)*(m2^2)*(px1^2)*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*(m2^2)*px1*py1*q1*q2*(x1^2)*y1 + (4//1)*(m2^2)*px1*py1*q1*q2*(x1^2)*y2 + (8//1)*(m2^2)*px1*py1*q1*q2*x1*x2*y1 - (8//1)*(m2^2)*px1*py1*q1*q2*x1*x2*y2 - (4//1)*(m2^2)*px1*py1*q1*q2*(x2^2)*y1 + (4//1)*(m2^2)*px1*py1*q1*q2*(x2^2)*y2 + (2//1)*(m2^2)*px1*py1*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m2^2)*px1*py1*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m2^2)*(py1^2)*q1*q2*x1*(y1^2) + (4//1)*(m2^2)*(py1^2)*q1*q2*x1*y1*y2 - (2//1)*(m2^2)*(py1^2)*q1*q2*x1*(y2^2) + (2//1)*(m2^2)*(py1^2)*q1*q2*x2*(y1^2) - (4//1)*(m2^2)*(py1^2)*q1*q2*x2*y1*y2 + (2//1)*(m2^2)*(py1^2)*q1*q2*x2*(y2^2) + (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2)) / (2(c^2)*(m1^2)*(m2^2)*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^5))
  d(py1)/dt = (-(2//1)*(m1^2)*(px2^2)*q1*q2*(x1^2)*y1 + (2//1)*(m1^2)*(px2^2)*q1*q2*(x1^2)*y2 + (4//1)*(m1^2)*(px2^2)*q1*q2*x1*x2*y1 - (4//1)*(m1^2)*(px2^2)*q1*q2*x1*x2*y2 - (2//1)*(m1^2)*(px2^2)*q1*q2*(x2^2)*y1 + (2//1)*(m1^2)*(px2^2)*q1*q2*(x2^2)*y2 - (4//1)*(m1^2)*px2*py2*q1*q2*x1*(y1^2) + (8//1)*(m1^2)*px2*py2*q1*q2*x1*y1*y2 - (4//1)*(m1^2)*px2*py2*q1*q2*x1*(y2^2) + (2//1)*(m1^2)*px2*py2*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*(m1^2)*px2*py2*q1*q2*x2*(y1^2) - (8//1)*(m1^2)*px2*py2*q1*q2*x2*y1*y2 + (4//1)*(m1^2)*px2*py2*q1*q2*x2*(y2^2) - (2//1)*(m1^2)*px2*py2*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m1^2)*(py2^2)*q1*q2*(y1^3) + (6//1)*(m1^2)*(py2^2)*q1*q2*(y1^2)*y2 - (6//1)*(m1^2)*(py2^2)*q1*q2*y1*(y2^2) + (2//1)*(m1^2)*(py2^2)*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m1^2)*(py2^2)*q1*q2*(y2^3) - (2//1)*(m1^2)*(py2^2)*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*y1 - (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*y2 + (4//1)*m1*m2*px1*px2*q1*q2*(x1^2)*y1 - (4//1)*m1*m2*px1*px2*q1*q2*(x1^2)*y2 - (8//1)*m1*m2*px1*px2*q1*q2*x1*x2*y1 + (8//1)*m1*m2*px1*px2*q1*q2*x1*x2*y2 + (4//1)*m1*m2*px1*px2*q1*q2*(x2^2)*y1 - (4//1)*m1*m2*px1*px2*q1*q2*(x2^2)*y2 + (4//1)*m1*m2*px1*py2*q1*q2*x1*(y1^2) - (8//1)*m1*m2*px1*py2*q1*q2*x1*y1*y2 + (4//1)*m1*m2*px1*py2*q1*q2*x1*(y2^2) - (2//1)*m1*m2*px1*py2*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*px1*py2*q1*q2*x2*(y1^2) + (8//1)*m1*m2*px1*py2*q1*q2*x2*y1*y2 - (4//1)*m1*m2*px1*py2*q1*q2*x2*(y2^2) + (2//1)*m1*m2*px1*py2*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*px2*py1*q1*q2*x1*(y1^2) - (8//1)*m1*m2*px2*py1*q1*q2*x1*y1*y2 + (4//1)*m1*m2*px2*py1*q1*q2*x1*(y2^2) - (2//1)*m1*m2*px2*py1*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*px2*py1*q1*q2*x2*(y1^2) + (8//1)*m1*m2*px2*py1*q1*q2*x2*y1*y2 - (4//1)*m1*m2*px2*py1*q1*q2*x2*(y2^2) + (2//1)*m1*m2*px2*py1*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*py1*py2*q1*q2*(y1^3) - (12//1)*m1*m2*py1*py2*q1*q2*(y1^2)*y2 + (12//1)*m1*m2*py1*py2*q1*q2*y1*(y2^2) - (4//1)*m1*m2*py1*py2*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*py1*py2*q1*q2*(y2^3) + (4//1)*m1*m2*py1*py2*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m2^2)*(px1^2)*q1*q2*(x1^2)*y1 + (2//1)*(m2^2)*(px1^2)*q1*q2*(x1^2)*y2 + (4//1)*(m2^2)*(px1^2)*q1*q2*x1*x2*y1 - (4//1)*(m2^2)*(px1^2)*q1*q2*x1*x2*y2 - (2//1)*(m2^2)*(px1^2)*q1*q2*(x2^2)*y1 + (2//1)*(m2^2)*(px1^2)*q1*q2*(x2^2)*y2 - (4//1)*(m2^2)*px1*py1*q1*q2*x1*(y1^2) + (8//1)*(m2^2)*px1*py1*q1*q2*x1*y1*y2 - (4//1)*(m2^2)*px1*py1*q1*q2*x1*(y2^2) + (2//1)*(m2^2)*px1*py1*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*(m2^2)*px1*py1*q1*q2*x2*(y1^2) - (8//1)*(m2^2)*px1*py1*q1*q2*x2*y1*y2 + (4//1)*(m2^2)*px1*py1*q1*q2*x2*(y2^2) - (2//1)*(m2^2)*px1*py1*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m2^2)*(py1^2)*q1*q2*(y1^3) + (6//1)*(m2^2)*(py1^2)*q1*q2*(y1^2)*y2 - (6//1)*(m2^2)*(py1^2)*q1*q2*y1*(y2^2) + (2//1)*(m2^2)*(py1^2)*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m2^2)*(py1^2)*q1*q2*(y2^3) - (2//1)*(m2^2)*(py1^2)*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2)) / (2(c^2)*(m1^2)*(m2^2)*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^5))
  d(px2)/dt = ((2//1)*(m1^2)*(px2^2)*q1*q2*(x1^3) - (6//1)*(m1^2)*(px2^2)*q1*q2*(x1^2)*x2 + (6//1)*(m1^2)*(px2^2)*q1*q2*x1*(x2^2) - (2//1)*(m1^2)*(px2^2)*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m1^2)*(px2^2)*q1*q2*(x2^3) + (2//1)*(m1^2)*(px2^2)*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*(m1^2)*px2*py2*q1*q2*(x1^2)*y1 - (4//1)*(m1^2)*px2*py2*q1*q2*(x1^2)*y2 - (8//1)*(m1^2)*px2*py2*q1*q2*x1*x2*y1 + (8//1)*(m1^2)*px2*py2*q1*q2*x1*x2*y2 + (4//1)*(m1^2)*px2*py2*q1*q2*(x2^2)*y1 - (4//1)*(m1^2)*px2*py2*q1*q2*(x2^2)*y2 - (2//1)*(m1^2)*px2*py2*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m1^2)*px2*py2*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m1^2)*(py2^2)*q1*q2*x1*(y1^2) - (4//1)*(m1^2)*(py2^2)*q1*q2*x1*y1*y2 + (2//1)*(m1^2)*(py2^2)*q1*q2*x1*(y2^2) - (2//1)*(m1^2)*(py2^2)*q1*q2*x2*(y1^2) + (4//1)*(m1^2)*(py2^2)*q1*q2*x2*y1*y2 - (2//1)*(m1^2)*(py2^2)*q1*q2*x2*(y2^2) - (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*x1 + (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*x2 - (4//1)*m1*m2*px1*px2*q1*q2*(x1^3) + (12//1)*m1*m2*px1*px2*q1*q2*(x1^2)*x2 - (12//1)*m1*m2*px1*px2*q1*q2*x1*(x2^2) + (4//1)*m1*m2*px1*px2*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*px1*px2*q1*q2*(x2^3) - (4//1)*m1*m2*px1*px2*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*px1*py2*q1*q2*(x1^2)*y1 + (4//1)*m1*m2*px1*py2*q1*q2*(x1^2)*y2 + (8//1)*m1*m2*px1*py2*q1*q2*x1*x2*y1 - (8//1)*m1*m2*px1*py2*q1*q2*x1*x2*y2 - (4//1)*m1*m2*px1*py2*q1*q2*(x2^2)*y1 + (4//1)*m1*m2*px1*py2*q1*q2*(x2^2)*y2 + (2//1)*m1*m2*px1*py2*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*m1*m2*px1*py2*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*px2*py1*q1*q2*(x1^2)*y1 + (4//1)*m1*m2*px2*py1*q1*q2*(x1^2)*y2 + (8//1)*m1*m2*px2*py1*q1*q2*x1*x2*y1 - (8//1)*m1*m2*px2*py1*q1*q2*x1*x2*y2 - (4//1)*m1*m2*px2*py1*q1*q2*(x2^2)*y1 + (4//1)*m1*m2*px2*py1*q1*q2*(x2^2)*y2 + (2//1)*m1*m2*px2*py1*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*m1*m2*px2*py1*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*py1*py2*q1*q2*x1*(y1^2) + (8//1)*m1*m2*py1*py2*q1*q2*x1*y1*y2 - (4//1)*m1*m2*py1*py2*q1*q2*x1*(y2^2) + (4//1)*m1*m2*py1*py2*q1*q2*x2*(y1^2) - (8//1)*m1*m2*py1*py2*q1*q2*x2*y1*y2 + (4//1)*m1*m2*py1*py2*q1*q2*x2*(y2^2) + (2//1)*(m2^2)*(px1^2)*q1*q2*(x1^3) - (6//1)*(m2^2)*(px1^2)*q1*q2*(x1^2)*x2 + (6//1)*(m2^2)*(px1^2)*q1*q2*x1*(x2^2) - (2//1)*(m2^2)*(px1^2)*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m2^2)*(px1^2)*q1*q2*(x2^3) + (2//1)*(m2^2)*(px1^2)*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*(m2^2)*px1*py1*q1*q2*(x1^2)*y1 - (4//1)*(m2^2)*px1*py1*q1*q2*(x1^2)*y2 - (8//1)*(m2^2)*px1*py1*q1*q2*x1*x2*y1 + (8//1)*(m2^2)*px1*py1*q1*q2*x1*x2*y2 + (4//1)*(m2^2)*px1*py1*q1*q2*(x2^2)*y1 - (4//1)*(m2^2)*px1*py1*q1*q2*(x2^2)*y2 - (2//1)*(m2^2)*px1*py1*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m2^2)*px1*py1*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m2^2)*(py1^2)*q1*q2*x1*(y1^2) - (4//1)*(m2^2)*(py1^2)*q1*q2*x1*y1*y2 + (2//1)*(m2^2)*(py1^2)*q1*q2*x1*(y2^2) - (2//1)*(m2^2)*(py1^2)*q1*q2*x2*(y1^2) + (4//1)*(m2^2)*(py1^2)*q1*q2*x2*y1*y2 - (2//1)*(m2^2)*(py1^2)*q1*q2*x2*(y2^2) - (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2)) / (2(c^2)*(m1^2)*(m2^2)*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^5))
  d(py2)/dt = ((2//1)*(m1^2)*(px2^2)*q1*q2*(x1^2)*y1 - (2//1)*(m1^2)*(px2^2)*q1*q2*(x1^2)*y2 - (4//1)*(m1^2)*(px2^2)*q1*q2*x1*x2*y1 + (4//1)*(m1^2)*(px2^2)*q1*q2*x1*x2*y2 + (2//1)*(m1^2)*(px2^2)*q1*q2*(x2^2)*y1 - (2//1)*(m1^2)*(px2^2)*q1*q2*(x2^2)*y2 + (4//1)*(m1^2)*px2*py2*q1*q2*x1*(y1^2) - (8//1)*(m1^2)*px2*py2*q1*q2*x1*y1*y2 + (4//1)*(m1^2)*px2*py2*q1*q2*x1*(y2^2) - (2//1)*(m1^2)*px2*py2*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*(m1^2)*px2*py2*q1*q2*x2*(y1^2) + (8//1)*(m1^2)*px2*py2*q1*q2*x2*y1*y2 - (4//1)*(m1^2)*px2*py2*q1*q2*x2*(y2^2) + (2//1)*(m1^2)*px2*py2*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m1^2)*(py2^2)*q1*q2*(y1^3) - (6//1)*(m1^2)*(py2^2)*q1*q2*(y1^2)*y2 + (6//1)*(m1^2)*(py2^2)*q1*q2*y1*(y2^2) - (2//1)*(m1^2)*(py2^2)*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m1^2)*(py2^2)*q1*q2*(y2^3) + (2//1)*(m1^2)*(py2^2)*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*y1 + (m1*px2*x1 - m1*px2*x2 + m1*py2*y1 - m1*py2*y2 - m2*px1*x1 + m2*px1*x2 - m2*py1*y1 + m2*py1*y2)*(-m1*px2*x1 + m1*px2*x2 - m1*py2*y1 + m1*py2*y2 + m2*px1*x1 - m2*px1*x2 + m2*py1*y1 - m2*py1*y2)*q1*q2*y2 - (4//1)*m1*m2*px1*px2*q1*q2*(x1^2)*y1 + (4//1)*m1*m2*px1*px2*q1*q2*(x1^2)*y2 + (8//1)*m1*m2*px1*px2*q1*q2*x1*x2*y1 - (8//1)*m1*m2*px1*px2*q1*q2*x1*x2*y2 - (4//1)*m1*m2*px1*px2*q1*q2*(x2^2)*y1 + (4//1)*m1*m2*px1*px2*q1*q2*(x2^2)*y2 - (4//1)*m1*m2*px1*py2*q1*q2*x1*(y1^2) + (8//1)*m1*m2*px1*py2*q1*q2*x1*y1*y2 - (4//1)*m1*m2*px1*py2*q1*q2*x1*(y2^2) + (2//1)*m1*m2*px1*py2*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*px1*py2*q1*q2*x2*(y1^2) - (8//1)*m1*m2*px1*py2*q1*q2*x2*y1*y2 + (4//1)*m1*m2*px1*py2*q1*q2*x2*(y2^2) - (2//1)*m1*m2*px1*py2*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*px2*py1*q1*q2*x1*(y1^2) + (8//1)*m1*m2*px2*py1*q1*q2*x1*y1*y2 - (4//1)*m1*m2*px2*py1*q1*q2*x1*(y2^2) + (2//1)*m1*m2*px2*py1*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*px2*py1*q1*q2*x2*(y1^2) - (8//1)*m1*m2*px2*py1*q1*q2*x2*y1*y2 + (4//1)*m1*m2*px2*py1*q1*q2*x2*(y2^2) - (2//1)*m1*m2*px2*py1*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*m1*m2*py1*py2*q1*q2*(y1^3) + (12//1)*m1*m2*py1*py2*q1*q2*(y1^2)*y2 - (12//1)*m1*m2*py1*py2*q1*q2*y1*(y2^2) + (4//1)*m1*m2*py1*py2*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (4//1)*m1*m2*py1*py2*q1*q2*(y2^3) - (4//1)*m1*m2*py1*py2*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m2^2)*(px1^2)*q1*q2*(x1^2)*y1 - (2//1)*(m2^2)*(px1^2)*q1*q2*(x1^2)*y2 - (4//1)*(m2^2)*(px1^2)*q1*q2*x1*x2*y1 + (4//1)*(m2^2)*(px1^2)*q1*q2*x1*x2*y2 + (2//1)*(m2^2)*(px1^2)*q1*q2*(x2^2)*y1 - (2//1)*(m2^2)*(px1^2)*q1*q2*(x2^2)*y2 + (4//1)*(m2^2)*px1*py1*q1*q2*x1*(y1^2) - (8//1)*(m2^2)*px1*py1*q1*q2*x1*y1*y2 + (4//1)*(m2^2)*px1*py1*q1*q2*x1*(y2^2) - (2//1)*(m2^2)*px1*py1*q1*q2*x1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (4//1)*(m2^2)*px1*py1*q1*q2*x2*(y1^2) + (8//1)*(m2^2)*px1*py1*q1*q2*x2*y1*y2 - (4//1)*(m2^2)*px1*py1*q1*q2*x2*(y2^2) + (2//1)*(m2^2)*px1*py1*q1*q2*x2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(m2^2)*(py1^2)*q1*q2*(y1^3) - (6//1)*(m2^2)*(py1^2)*q1*q2*(y1^2)*y2 + (6//1)*(m2^2)*(py1^2)*q1*q2*y1*(y2^2) - (2//1)*(m2^2)*(py1^2)*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(m2^2)*(py1^2)*q1*q2*(y2^3) + (2//1)*(m2^2)*(py1^2)*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) - (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*y1*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2) + (2//1)*(c^2)*(m1^2)*(m2^2)*q1*q2*y2*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^2)) / (2(c^2)*(m1^2)*(m2^2)*(sqrt((x1 - x2)^2 + (y1 - y2)^2)^5))
```


### Polar Reduced Form (Theoretical Analysis)

In dimensionless units (ρ = 1, k = 1, μ = 1):

**Standard form:** H = (1//2)*(p_r^2)*(1 + -1 / r) + (ℓ^2) / (2(r^2)) + 1 / r

**Metric form:** H = (ℓ^2) / (2(r^2)) + ((p_r^2)*r) / (2(-1 + r)) + 1 / r

**Energy equation:** ṙ² = (ℓ² + 2r − 2hr²) / (r(1 − r))

**Weber plane metric:** g_rr = (r − 1)/r, g_φφ = r²
- r > 1: Riemannian (normal dynamics)
- r = 1: degenerate (critical radius barrier)
- r < 1: Lorentzian (inverted inertia → like-charge binding)

**Kepler decomposition:** H = H_kepler + ρ · H_correction

- H_kepler = (ℓ^2) / (2(r^2)) + 1 / r + (1//2)*(p_r^2)
- H_correction = (-(p_r^2)) / (2r)


## The Collision Dichotomy

The central result: collision behavior splits fundamentally by angular momentum.
We verify all analytical predictions numerically.


### Trajectory Comparison: ℓ = 0 vs ℓ ≠ 0


![Collision trajectories for different angular momenta](sub_critical_weber/01_trajectories.png)


| ℓ | min r reached | winding (revolutions) | collision time |
| --- | --- | --- | --- |
| 0.0 | 1.00e-12 | 0.0 | 0.3789 |
| 0.01 | 1.60e-11 | 79413.9 | 0.3787 |
| 0.05 | 4.70e-11 | 46422.1 | 0.3751 |
| 0.1 | 7.48e-11 | 36788.3 | 0.3662 |
| 0.2 | 1.18e-10 | 29285.7 | 0.3400 |


### Power-Law Verification

**Analytical predictions:**
- ℓ = 0: r ~ (T_coll − t)^1, collision speed → √2·c
- ℓ ≠ 0: r ~ (T_coll − t)^(2/3), collision speed → ∞


![Power-law scaling near collision](sub_critical_weber/02_power_law.png)


| ℓ | measured exponent β | theoretical β | R² |
| --- | --- | --- | --- |
| 0.0 | 0.998 | 1.000 | 0.999993 |
| 0.1 | 0.701 | 0.667 | 0.998376 |


### Infinite Winding: The Topological Obstruction

The total angle accumulated from r₀ to r_min diverges as r_min → 0:

    φ = ∫ (ℓ/r²) / |ṙ| dr ~ ∫ r^(-3/2) dr → ∞

This infinite winding number is a **topological invariant** preserved by all smooth
coordinate changes. It is the fundamental reason why 7 regularization approaches fail.


![Winding number divergence as r_min → 0](sub_critical_weber/03_winding_growth.png)


**Winding growth exponent:** N_wind ~ r_min^-0.567 (theory: r_min^(-0.5))

| r_min | revolutions |
| --- | --- |
| 1e-02 | 1.0 |
| 1e-03 | 6.6 |
| 1e-04 | 27.7 |
| 1e-05 | 96.4 |
| 1e-06 | 314.2 |
| 1e-07 | 1003.2 |
| 1e-08 | 3182.7 |


### Collision Time (Finite for Both Cases)


![Collision time is finite for both ℓ=0 and ℓ≠0](sub_critical_weber/04_collision_times.png)

Both collision types reach r = 0 in **finite time**. The difference is that ℓ ≠ 0 accumulates infinite winding during that finite time.


### Cartesian Spiral Trajectory


![Cartesian spiraling collision: lab frame (left) and relative frame (right)](sub_critical_weber/05_cartesian_spiral.png)


**Energy conservation:** max relative error = 8.22e-01


![Energy conservation of the Cartesian ODE integrator](sub_critical_weber/06_energy_conservation.png)


## Effective Potential and Weber Metric


![Effective potential outside (left) and inside (right) the critical radius](sub_critical_weber/07_effective_potential.png)

**Key insight:** In Kepler, the centrifugal barrier ℓ²/(2r²) prevents collision for ℓ ≠ 0.
In Weber's sub-critical regime, the **negative effective mass** overcomes this barrier,
allowing the particle to reach r = 0 despite having angular momentum. This is the
origin of the non-regularizable spiral.


![Weber metric signature change at the critical radius](sub_critical_weber/08_metric_signature.png)


## Coordinate Transform Analysis

12 canonical transforms cataloged, each tracking singularity effects and smoothness.


### Singularity Summary

- **r = 0** (pole, order 2): Collision singularity. For ℓ=0: regularizable (finite speed √2c). For ℓ≠0: non-regularizable (infinite winding, topological obstruction).
- **r = 1** (metric_degeneracy, order 1): Metric signature change. Coordinate singularity removable by Kruskal-Szekeres. g_rr changes sign: Riemannian (r>ρ) → Lorentzian (r<ρ).
- **r = 0 (centrifugal)** (pole, order 2): Centrifugal term ℓ²/(2r²). In Kepler, creates a barrier preventing collision. In Weber (sub-critical), the negative effective mass overcomes this barrier, ALLOWING collision.


### Transform Catalog

| Transform | Type | Smooth | r=0 effect | r=ρ effect |
| --- | --- | --- | --- | --- |
| birkhoff_inversion | point | ✓ | moved_to_infinity | preserved |
| kruskal_szekeres | non_smooth | ✗ | preserved | removed |
| levi_civita | point | ✓ | removed | preserved |
| log_map | point | ✗ | moved_to_neg_infinity | preserved |
| mcgehee_blowup | point | ✓ | preserved | preserved |
| midpoint_shift | point | ✓ | shifted | shifted |
| momentum_space | generating_function | ✓ | transformed | transformed |
| radial_flattening | point | ✓ | preserved | removed |
| sigmoid | generating_function | ✓ | moved_to_neg_infinity | moved_to_pos_infinity |
| sundman_α=1 | time_scaling | ✓ | partially_regularized | preserved |
| sundman_α=2 | time_scaling | ✓ | partially_regularized | preserved |
| tortoise | point | ✓ | preserved | moved_to_neg_infinity |


### Winding Number Preservation

The winding number (total revolutions around origin) is a **topological invariant**
under smooth diffeomorphisms. Only non-smooth maps can potentially alter it.

| Transform | Smooth | Winding preserved? | Implication |
| --- | --- | --- | --- |
| birkhoff_inversion | ✓ | Yes | Cannot escape barrier |
| kruskal_szekeres | ✗ | **No** ★ | Candidate for barrier escape |
| levi_civita | ✓ | Yes | Cannot escape barrier |
| log_map | ✗ | **No** ★ | Candidate for barrier escape |
| mcgehee_blowup | ✓ | Yes | Cannot escape barrier |
| midpoint_shift | ✓ | Yes | Cannot escape barrier |
| momentum_space | ✓ | Yes | Cannot escape barrier |
| radial_flattening | ✓ | Yes | Cannot escape barrier |
| sigmoid | ✓ | Yes | Cannot escape barrier |
| sundman_α=1 | ✓ | Yes | Cannot escape barrier |
| sundman_α=2 | ✓ | Yes | Cannot escape barrier |
| tortoise | ✓ | Yes | Cannot escape barrier |

**★ Non-smooth transforms** (Kruskal-Szekeres, log map) are candidates for escaping the
topological barrier. However, log map was tested (Experiment 7) and failed due to velocity
divergence. Kruskal-Szekeres removes the r = ρ singularity but not r = 0.


### Interesting Pairwise Compositions

Compositions where the r = 0 singularity is affected AND smoothness might be broken:

| Transform 1 | Transform 2 | r=0 effect | Smooth? |
| --- | --- | --- | --- |
| birkhoff_inversion | log_map | moved_to_neg_infinity | ✗ (non-smooth) |
| kruskal_szekeres | levi_civita | removed | ✗ (non-smooth) |
| kruskal_szekeres | log_map | moved_to_neg_infinity | ✗ (non-smooth) |
| kruskal_szekeres | sigmoid | moved_to_neg_infinity | ✗ (non-smooth) |
| levi_civita | log_map | removed | ✗ (non-smooth) |
| log_map | sigmoid | moved_to_neg_infinity | ✗ (non-smooth) |


## Research Direction Ranking


![Research direction ranking (blue=beyond smooth, green=redefine, orange=Lorentzian, purple=algebraic, brown=computational)](sub_critical_weber/09_ranking.png)


| # | Direction | Score | Category | Status | Key Question |
| --- | --- | --- | --- | --- | --- |
| 1 | ○ Blow-up and desingularization (Bendixson-Dumortier) | 4.4 | A | untested | Does the iterated blow-up terminate? What is the desingularized phase portrait? |
| 2 | ◐ Angular averaging / homogenization | 4.2 | C | partially_tested | Does the angular-averaged dynamics preserve energy? Is it unique? |
| 3 | ○ Distributional / measure-valued continuation | 4.2 | A | untested | What is the correct measure-valued continuation that matches the quantum L² solution? |
| 4 | ○ Renormalization group approach | 3.6 | C | untested | What are the relevant and irrelevant directions at the RG fixed point? |
| 5 | ○ Painlevé analysis | 3.6 | D | untested | Does the Weber ODE pass the Painlevé test? If not, what type of movable singularity? |
| 6 | ○ Picard-Fuchs monodromy | 3.4 | D | untested | Is the monodromy group at r=0 finite or infinite? |
| 7 | ○ Piecewise-smooth regularization | 3.4 | A | untested | Can the stitching maps be chosen so that the dynamics is continuous across cuts? |
| 8 | ◐ Partial radial regularization with angular bookkeeping | 3.4 | E | partially_tested | Is 'modulo infinite winding' a physically meaningful definition? |
| 9 | ○ Stochastic regularization | 3.4 | C | untested | What is the correct noise strength that matches the quantum distribution? |
| 10 | ○ Penrose diagram / conformal compactification | 3.2 | B | untested | What does the Penrose diagram of the sub-critical Weber plane look like? |
| 11 | ◐ Covering space / Riemann surface unwinding | 3.2 | A | partially_tested | The log map was tested (Experiment 7) and failed. Can a COMPACTIFICATION of the log-image work? |
| 12 | ○ Characteristic (null) coordinates | 3.0 | B | untested | Do the null coordinates simplify the angular dynamics? |
| 13 | ○ Galois theory of variational equations | 2.8 | D | untested | Is the identity component of the Galois group abelian (→ integrable) or not? |
| 14 | ○ Newman-Penrose / spin-weighted formalism | 2.4 | B | untested | Is the Weber plane sufficiently rich to benefit from NP formalism? |
| 15 | ○ Adaptive angular mesh integration | 2.4 | E | untested | Can adaptive refinement track the spiral accurately enough to be useful? |


### Top 3 Detailed Scores


| Direction | Promise | Feasibility | Novelty | Connection | Scope | Composite |
| --- | --- | --- | --- | --- | --- | --- |
| Blow-up and desingularization (Bendixson-Dumortier) | 5 | 4 | 5 | 4 | 4 | 4.4 |
| Angular averaging / homogenization | 4 | 5 | 4 | 5 | 3 | 4.2 |
| Distributional / measure-valued continuation | 4 | 3 | 5 | 5 | 4 | 4.2 |


### Category Analysis

| Category | Directions | Avg Score | Best Direction |
| --- | --- | --- | --- |
| A beyond smooth | 4 | 3.80 | Blow-up and desingularization (Bendixson-Dumortier) |
| C redefine | 3 | 3.73 | Angular averaging / homogenization |
| D algebraic | 3 | 3.27 | Painlevé analysis |
| E computational | 2 | 2.90 | Partial radial regularization with angular bookkeeping |
| B lorentzian | 3 | 2.87 | Penrose diagram / conformal compactification |


## Top Research Directions — Detailed Analysis


### Blow-up and desingularization (Bendixson-Dumortier)

Apply iterated blow-ups at the singular point, as in algebraic geometry resolution of singularities. Each blow-up replaces a point by a projective line. For ODEs, the Bendixson-Dumortier procedure desingularizes non-elementary singular points. The Weber spiral singularity IS a non-elementary singular point.

**Key question:** Does the iterated blow-up terminate? What is the desingularized phase portrait?

**Notes:** TOP RANKED. This is a constructive algorithm specifically designed for this class of singularity.

**Related literature:** Bendixson-Dumortier, Dumortier-Llibre-Artes, resolution of singularities


### Angular averaging / homogenization

Near collision, φ cycles infinitely faster than r changes (scale separation). Average over the angle to obtain effective radial dynamics. The averaged radial dynamics is known to be regular (Experiment 5). Question: is this a physically meaningful continuation?

**Key question:** Does the angular-averaged dynamics preserve energy? Is it unique?

**Notes:** TOP RANKED. Experiment 5 showed radial regularization works; this formalizes it.

**Related literature:** homogenization theory, averaging for fast-slow systems, Arnol'd averaging


### Distributional / measure-valued continuation

Replace pointwise trajectory by probability measure that spreads over angles near collision. Connect to quantum limit-circle analysis: the Schrödinger solution exists in L² even when the classical trajectory fails.

**Key question:** What is the correct measure-valued continuation that matches the quantum L² solution?

**Related literature:** Frauenfelder-Weber 2024, limit circle Sturm-Liouville, measure-valued solutions ODE


### Renormalization group approach

Treat approach to r=0 as a scaling limit. The self-similar asymptotics r~t^(2/3) suggest an RG fixed point. Classify all continuations by their behavior near this fixed point.

**Key question:** What are the relevant and irrelevant directions at the RG fixed point?

**Related literature:** self-similar solutions ODE, renormalization group dynamical systems


### Painlevé analysis

Apply the Painlevé test (Laurent series near movable singularity). Classify the singularity type: algebraic, logarithmic, or essential. This is a standard algorithmic test that provides definitive information about the singularity.

**Key question:** Does the Weber ODE pass the Painlevé test? If not, what type of movable singularity?

**Notes:** TOP RANKED for feasibility. Painlevé analysis is algorithmic and can be done symbolically.

**Related literature:** ARS algorithm, Painlevé analysis Hamiltonian systems, movable singularities


## Literature Survey


### averaging (69 papers)

Queries: angular averaging homogenization fast rotation; fast-slow system averaging singular limit; homogenization rapidly oscillating ODE; adiabatic invariant fast angular variable

1. **Chern Theorem and Topological Matter in Fast-Rotating Atomic Nuclei.** (2026) — M. Guidry; Yang Sun [link](https://www.semanticscholar.org/paper/9b38508034c8f6d80b37303b4dde9257fa22ff01)
2. **exoALMA. VI. Rotating under Pressure: Rotation Curves, Azimuthal Velocity Substructures, and Gas Pressure Variations** (2025) — J. Stadler; M. Benisty; A. Winter; A. Izquierdo; C. Longarini; Maria Galloway-Sp... [link](https://www.semanticscholar.org/paper/a5241742086d739eb11b467f79b0c604c415bf46)
3. **The tidal interaction of an orbiting giant planet with a star near the Kraft break: the excitation of r-modes and the retention of orbital and spin angular momenta misalignment** (2025) — J. Papaloizou; G. Savonije [link](https://www.semanticscholar.org/paper/01b2c999319830a4512f2167e58ade5d5767484a)
4. **Homogenization of attractors to reaction - diffusion equations in domains with rapidly oscillating boundary: supercritical case** (2025) — Gaziz Faizullaevich Azhmoldayev; K. Bekmaganbetov; G. Chechkin; V. Chepyzhov [link](https://www.semanticscholar.org/paper/99a3bed1b256082ad038bd7dfb0dcbea9b369a10)
5. **Homogenization of Attractors to Reaction–Diffusion Equations in Domains with Rapidly Oscillating Boundary: Subcritical Case** (2025) — G. Azhmoldaev; K. Bekmaganbetov; G. Chechkin; V. Chepyzhov [link](https://www.semanticscholar.org/paper/287235e0661a0981f3d1d507083b2d8fb53e8d09)
6. **Post-adiabatic dynamics and waveform generation in self-force theory: an invariant pseudo-Hamiltonian framework** (2025) — Jack Lewis; Takafumi Kakehi; Adam Pound; Takahiro Tanaka [link](https://www.semanticscholar.org/paper/00608d267affa9290ba2e425d7391a3a004871eb)
7. **Angular resampling-sparse representation classification (AR-SRC): A new method for bearing fault diagnosis in non-stationary conditions** (2025) — Mohamed Sekini; Bilal El Yousfi; T. Benkedjouh; Kamal Medjaher; A. Lourari [link](https://www.semanticscholar.org/paper/3d3695c40183db36e123f9ab76daa9b726f5dbcf)
8. **Inversion for Inferring Solar Meridional Circulation: The Case with Constraints on Angular Momentum Transport inside the Sun** (2024) — Y. Hatta; Hideyuki Hotta; T. Sekii [link](https://www.semanticscholar.org/paper/49798bd1a74b9ba14ff41aded1c0122feab86bd9)
9. **Renormalization group and elliptic homogenization in high contrast** (2024) — Scott Armstrong; Tuomo Kuusi [link](http://arxiv.org/abs/2405.10732v3)
10. **Asymptotic analysis of a coupled ODE‐PDE system arising from heterogeneous diffusion‐reaction kinetics** (2024) — Michal Beneš; Michael Eden; Adrian Muntean [link](https://www.semanticscholar.org/paper/499b3926d586d2a89c9f6aedc06c130ea658cbef)


### limit_circle (53 papers)

Queries: limit circle Sturm-Liouville oscillatory endpoint; singular Sturm-Liouville self-adjoint extension; Weyl limit circle limit point classification; wildly oscillating solutions Sturm-Liouville

1. **On Sesquilinear Forms for Lower Semibounded (Singular) Sturm-Liouville Operators** (2025) — Jussi Behrndt; Fritz Gesztesy; Seppo Hassi; Roger Nichols; Henk de Snoo [link](http://arxiv.org/abs/2509.07205v1)
2. **Self-adjoint extensions of singular Sturm-Liouville operators on graphs and Weyl's law** (2025) — Elisha Falbel [link](http://arxiv.org/abs/2510.18922v1)
3. **ON THE SQUARE-INTEGRABLE SOLUTIONS OF SINGULAR IMPULSIVE DYNAMIC STURM-LIOUVILLE EQUATIONS** (2025) — B. Allahverdiev; H. Tuna; Hamlet A. Isayev [link](https://www.semanticscholar.org/paper/6d7f240b2e0e2838d9436253b16d43d1a551e644)
4. **Finer limit circle/limit point classification for Sturm-Liouville operators** (2024) — Mateusz Piorkowski; Jonathan Stanfill [link](http://arxiv.org/abs/2407.04847v1)
5. **On the non-existence of oscillation numbers in Sturm-Liouville theory** (2024) — Angelo B. Mingarelli [link](http://arxiv.org/abs/2404.19575v1)
6. **Wintner-type nonoscillation theorems for conformable linear Sturm-Liouville differential equations** (2024) — Kazuki Ishibashi [link](https://www.semanticscholar.org/paper/ab41e98ab3c50b8ef01c70eb69f3f184c1ab07ef)
7. **Limit point and limit circle trichotomy for Sturm-Liouville problems with complex potentials** (2023) — Florian Leben; Edison Leguizamón; Carsten Trunk; Monika Winklmeier [link](http://arxiv.org/abs/2310.16128v2)
8. **The limit-point/limit-circle classification for ordinary differential equations with distributional coefficients** (2023) — Varun Bhardwaj; Rudi Weikard [link](http://arxiv.org/abs/2309.05068v1)
9. **A Fixed-Point Approach to Non-Commutative Central Limit Theorems** (2023) — Jad Hamdan [link](http://arxiv.org/abs/2305.06960v6)
10. **Titchmarsh–Weyl Theory for Impulsive Dynamic Dirac System** (2023) — B. Allahverdiev; H. Tuna [link](https://www.semanticscholar.org/paper/390431248819bc5fb0bac1ed0279168da3644af9)


### lorentzian (58 papers)

Queries: Lorentzian signature change classical mechanics; pseudo-Riemannian particle dynamics geodesic; metric signature change horizon classical; Penrose diagram classical mechanics

1. **Signature Change in $f(R, T_\phi)$ Theory** (2026) — Serkan Doruk Hazinedar; Y. Heydarzade [link](https://www.semanticscholar.org/paper/515f458af8647021e5918151203681b534c6e22b)
2. **Spontaneous Emergence of Lorentzian Signature from Curvature-Minimizing Geometry** (2025) — Miguel Bermudez [link](http://arxiv.org/abs/2510.07891v2)
3. **Lorentzian-Euclidean black holes and Lorentzian to Riemannian metric transitions** (2025) — Rossella Bartolo; E. Caponio; A. Germinario; Miguel Sánchez [link](https://www.semanticscholar.org/paper/11fdb21fdbe78321fd0ea3efb9cbeb556faba869)
4. **A Blockchain-Assisted Fair Exchange Signature Protocol Using Quantum Key Distribution for Metaverse Environment** (2025) — Sunil Prajapat; Pankaj Kumar; Goutham Reddy Alavalapati [link](https://www.semanticscholar.org/paper/d8a449b010f69d063543a5fce156a0a9602ad352)
5. **Geodesics of charged particle in electromagnetic field** (2025) — Nitish Yadav; Seema Jangir [link](http://arxiv.org/abs/2501.14814v2)
6. **Horizon-scale variability of M87* from 2017--2021 EHT observations** (2025) — The Event Horizon Telescope Collaboration [link](http://arxiv.org/abs/2509.24593v1)
7. **A Conceptual Introduction To Signature Change Through a Natural Extension of Kaluza-Klein Theory** (2025) — Vincent Moncrief; N. E. Rieger [link](https://www.semanticscholar.org/paper/653682362a2940bd410eff61f063401ceb79311f)
8. **The Construction and Application of Penrose Diagrams, with a Focus on the Maximally Analytically Extended Schwarzschild Spacetime** (2025) — Christian Röken [link](http://arxiv.org/abs/2507.23514v2)
9. **A New Space-Time Theory Unravels the Origins of Classical Mechanics for the Dirac Equation** (2025) — Wei Wen [link](https://www.semanticscholar.org/paper/bda5b44afaa3f501ba775666b0194ea52b3d5d6f)
10. **Penrose's eight-conic theorem** (2024) — Russell Arnold; Albert Chern; Morten Eide; Charles Gunn; Thomas Neukirchner; Rog... [link](http://arxiv.org/abs/2409.17150v8)


### regularization (84 papers)

Queries: spiral singularity regularization ODE; infinite winding collision regularization; collision regularization celestial mechanics; triple collision regularization n-body; Sundman regularization singular ODE; Levi-Civita Kustaanheimo-Stiefel regularization

1. **A Freefall-based Switching Criterion for P3T Hybrid N-body Methods in Collisional Stellar Systems** (2026) — Long 龙 Wang 王; David M. Hernandez; Zepeng 泽鹏 Zheng 郑; Wanhao 万豪 Huang 黄 [link](https://www.semanticscholar.org/paper/9085b08f4ae056696d22250a3c07c9ca0b22d871)
2. **SPIRAL: Self-Play on Zero-Sum Games Incentivizes Reasoning via Multi-Agent Multi-Turn Reinforcement Learning** (2025) — Bo Liu; Leon Guertler; Simon Yu; Zichen Liu; Penghui Qi; Daniel Balcells; Mickel... [link](http://arxiv.org/abs/2506.24119v3)
3. **Regular celestial amplitudes** (2025) — Reiko Liu; Wen-Jie Ma [link](http://arxiv.org/abs/2512.05882v2)
4. **Global Regularity for Navier-Stokes on T3 via Bounded Vorticity-Response Functionals** (2025) — Jeffrey Camlin [link](https://www.semanticscholar.org/paper/531368e58e7b665b0b5447d510e696b4dcd46053)
5. **Zero noise limit for singular ODE regularized by fractional noise** (2024) — Łukasz Mądry; Paul Gassiat [link](http://arxiv.org/abs/2401.09970v1)
6. **One-Shot Method for Computing Generalized Winding Numbers** (2024) — Cedric Martens; Mikhail Bessmeltsev [link](http://arxiv.org/abs/2408.04466v2)
7. **An FFT based chemo-mechanical framework with fracture: application to mesoscopic electrode degradation** (2024) — Gabriel Zarzoso; Eduardo Roque; Francisco Montero-Chacón; Javier Segurado [link](http://arxiv.org/abs/2411.16583v1)
8. **Non-linear Triple Changes Estimator for Targeted Policies** (2024) — Sina Akbari; Negar Kiyavash [link](http://arxiv.org/abs/2402.12583v1)
9. **Direct N-body Simulations of Satellite Formation around Small Asteroids: Insights from DART’s Encounter with the Didymos System** (2024) — H. Agrusa; Yun Zhang; D. C. Richardson; P. Pravec; M. Ćuk; P. Michel; R. Ballouz... [link](https://www.semanticscholar.org/paper/80de01fe148401772b4d86064373f064ecddf37e)
10. **Regularizing fuel-optimal multi-impulse trajectories** (2024) — Kenta Oshima [link](https://www.semanticscholar.org/paper/17e4223bcef27c8a6eb92d177aaad1b13d268ea6)


### relational_mechanics (48 papers)

Queries: Assis relational mechanics inertia; Mach principle Weber electrodynamics; Weber force relational inertia gravitational; relational mechanics n-body Weber gravitational

1. **Deep Biomechanically-Guided Interpolation for Keypoint-Based Brain Shift Registration** (2025) — Tiago Assis; Ines P. Machado; Benjamin Zwick; Nuno C. Garcia; Reuben Dorent [link](http://arxiv.org/abs/2508.13762v1)
2. **THE LAGRANGE - JACOBI EQUATION AND ITS APPICATION TO THE N - BODY PROBLEM** (2025) — G.Т. Оmarova; Zh.Т. Оmarova [link](https://www.semanticscholar.org/paper/67b19cf026e1e3602fb0b7df961377eb38497f16)
3. **An FFT based chemo-mechanical framework with fracture: application to mesoscopic electrode degradation** (2024) — Gabriel Zarzoso; Eduardo Roque; Francisco Montero-Chacón; Javier Segurado [link](http://arxiv.org/abs/2411.16583v1)
4. **Kerr Geodesics in horizon-penetrating Kerr coordinates: description in terms of Weierstrass functions** (2024) — Zuzanna Bakun; Angelika Łukanty; Anastasiia Untilova; Adam Cieślik; Patryk Mach [link](http://arxiv.org/abs/2409.03722v2)
5. **Unfairly Splitting Separable Necklaces** (2024) — Patrick Schnider; Linus Stalder; Simon Weber [link](http://arxiv.org/abs/2408.17126v1)
6. **Equivariant optimisation for the gravitational n-body problem: A computational factory of symmetric orbits** (2024) — V. Barutello; Gian Marco Canneori; Roberto Ciccarelli; Susanna Terracini; M. Ber... [link](https://www.semanticscholar.org/paper/8aec21415006988228e166bdbd5529fd10991695)
7. **A Concrete Model for the Quantum Permutation Group on 4 Points** (2023) — Nicolas Faroß; Moritz Weber [link](http://arxiv.org/abs/2304.09124v1)
8. **Framework for the full N-body problem in SE(3) and its reduction to the circular restricted full three-body problem** (2023) — Morad Nazari; D. Canales; Brennan S. McCann; Eric A. Butcher; K. Howell [link](https://www.semanticscholar.org/paper/a037e27c113d019813929643f95506cb0ca2c873)
9. **Sizing of Energy Storage System for Virtual Inertia Emulation** (2022) — Mohamed Abuagreb; Ahmed Abuhussein; Saif alZahir [link](http://arxiv.org/abs/2201.06566v2)
10. **Nonradial stability of expanding Goldreich-Weber stars** (2022) — Mahir Hadžić; Juhi Jang; King Ming Lam [link](http://arxiv.org/abs/2212.11420v2)


### singularity_theory (76 papers)

Queries: blow-up desingularization planar ODE Bendixson; Dumortier desingularization singular point; resolution of singularities dynamical system; Painleve analysis movable singularity Hamiltonian; movable singularity classification ODE

1. **Group-circulant singularities and partial desingularization preserving normal crossings** (2026) — André Belotto da Silva; Edward Bierstone [link](http://arxiv.org/abs/2602.09114v1)
2. **Secure Multiuser Beamforming With Movable Antenna Arrays** (2026) — Zhenqiao Cheng; Chongjun Ouyang; Boqun Zhao; Xingqi Zhang [link](http://arxiv.org/abs/2601.05686v1)
3. **ECH capacities of concave singular toric domains** (2025) — Jonathan Trejos [link](https://www.semanticscholar.org/paper/9e2b0bff65890331a06b156cf38821804a5eeca2)
4. **On entry-exit formulas for degenerate turning point problems in planar slow-fast systems** (2025) — Renato Huzak; K. U. Kristiansen [link](https://www.semanticscholar.org/paper/e9d545a7f010d4b6d0e19e3e77748ba747d8b354)
5. **Cosmological study of finite-time singularities under dynamical systems survey in covariant modified teleparallel gravity** (2025) — M. Ganiou; M. Toure; C. Aïnamon; S. I. V. Hontinfinde; M. Houndjo [link](https://www.semanticscholar.org/paper/026f89298b282cbdc32b331c6c7447e6c55d2143)
6. **On the geometry of a 4-dimensional extension of a q-Painlevé I equation with symmetry type A1(1)** (2025) — Alexander Stokes; T. Takenawa; Adrian Stefan Carstea [link](https://www.semanticscholar.org/paper/3cd44252d1e60b63ddbdd66c1380cc1c9d756571)
7. **Center and Degenerate Hopf Bifurcation Cyclicity of high-order Singularity in a Class of Three-Dimensional Systems** (2025) — Jingping Lu; Jie Yao; Qinlong Wang [link](https://www.semanticscholar.org/paper/da4fad903e2bf2fbefce21d908218776e2eb3b0a)
8. **On 2-Movable Total Domination in the Join and Corona of Graphs** (2025) — Ariel C. Pedrano; Rolando N. Paluga [link](http://arxiv.org/abs/2508.10952v1)
9. **ROMA: ROtary and Movable Antenna** (2025) — Jiayi Zhang; Wenhui Yi; Bokai Xu; Zhe Wang; Huahua Xiao; Bo Ai [link](http://arxiv.org/abs/2501.13403v2)
10. **On the blow-up formula of weighted stability for polarized toric manifolds** (2024) — King Leung Lee; Naoto Yotsutani [link](http://arxiv.org/abs/2407.10082v2)


### weber (63 papers)

Queries: Weber electrodynamics velocity-dependent potential; Weber force law classical mechanics; Frauenfelder Weber nucleus atomic; Weber electrodynamics critical radius binding; Assis relational mechanics Weber

1. **Atomic-resolution imaging reveals nucleus-free crystallization in two-dimensional amorphous ice on graphite** (2025) — Zifeng Yuan; Ye Tian; Binze Tang; Tiancheng Liang; Chon-Hei Lo; Zixiang Yan; Don... [link](https://www.semanticscholar.org/paper/13c2675fa40799e6f800e9a49a356b13b6596510)
2. **Interactive 3D Visualization of Bohr's Atomic Model: Enhancing Educational Tools with WebGL and Force-Directed Algorithms** (2025) — Zaid Kraitem; Hamza Alhaj; Mohamad Taky [link](https://www.semanticscholar.org/paper/55f9aec973d9e2c82862f74fbf488973edb78998)
3. **Bound Deuteron-Antideuteron System (Deuteronium): Leading Radiative and Internal-Structure Corrections to Bound-State Energies** (2025) — G. S. Adkins; U. Jentschura [link](https://www.semanticscholar.org/paper/66ca8bdc2bd63f51713e59a2cf26117ad0d713e7)
4. **Application of Hard X-Ray and Gamma-Ray TES Microcalorimeter at Accelerator Facility** (2025) — T. Saito; Shinji Okada; Y. Toyama; T. Azuma; Baptista Goncalo; D. Becker; D. Ben... [link](https://www.semanticscholar.org/paper/98435cba6465e040fab67cf63c3f4165f7ac4b62)
5. **Deep Biomechanically-Guided Interpolation for Keypoint-Based Brain Shift Registration** (2025) — Tiago Assis; Ines P. Machado; Benjamin Zwick; Nuno C. Garcia; Reuben Dorent [link](http://arxiv.org/abs/2508.13762v1)
6. **Unfairly Splitting Separable Necklaces** (2024) — Patrick Schnider; Linus Stalder; Simon Weber [link](http://arxiv.org/abs/2408.17126v1)
7. **An FFT based chemo-mechanical framework with fracture: application to mesoscopic electrode degradation** (2024) — Gabriel Zarzoso; Eduardo Roque; Francisco Montero-Chacón; Javier Segurado [link](http://arxiv.org/abs/2411.16583v1)
8. **A mathematical description of the Weber nucleus as a classical and quantum mechanical system** (2024) — Urs Frauenfelder; Joa Weber [link](https://www.semanticscholar.org/paper/ce5f0301cf507e17d4cd86f7cf369a67e961a2b2)
9. **Isomorphic Fluorescent Nucleosides** (2024) — Y. Tor [link](https://www.semanticscholar.org/paper/c8439d655e29f5f42a3a15ee1aab136361aecc6b)
10. **Relational Quantum Dynamics (RQD): An Informational Ontology** (2024) — A. Zaghi [link](https://www.semanticscholar.org/paper/0111c7d7c03d8400a25b41d28f84a47375340e15)

**Total papers found:** 451


## Asymptotic Exponents Reference

| Quantity | ℓ = 0 (head-on) | ℓ ≠ 0 (spiraling) |
| --- | --- | --- |
| r(t) near collision | ~ (T−t)¹ | ~ (T−t)^(2/3) |
| ṙ at r → 0 | → √2·c (finite) | ~ r^(−1/2) → ∞ |
| φ̇ at r → 0 | = 0 | ~ ℓ/r² → ∞ |
| dφ/dr near r = 0 | = 0 | ~ r^(−3/2) → ∞ |
| Total winding | 0 | ∞ (divergent) |
| Collision time | Finite | Finite |
| Regularizable? | **Yes** (collision bounce) | **No** (topological obstruction) |
| Quantum analog | Non-oscillating, natural BCs | Wildly oscillating, no natural BCs |


## Conclusions and Next Steps

### Key Finding

The angular momentum regularization barrier is **topological**: the infinite winding number
of spiraling collisions is preserved by all smooth coordinate changes. Seven standard
regularization approaches have been tested and all fail. Progress requires fundamentally
different mathematical tools.

### Most Promising Directions

1. **Blow-up/desingularization (Bendixson-Dumortier)** [Score: 4.4] — A constructive
   algorithm from algebraic geometry specifically designed for non-elementary singular
   points of planar ODEs. This has NOT been tried. It is the top-ranked direction.

2. **Angular averaging / homogenization** [Score: 4.2] — Near collision, the angle
   cycles infinitely faster than the radius changes. Averaging over the fast angular
   variable gives regular effective radial dynamics. This is partially validated
   (Experiment 5 showed radial regularization works).

3. **Distributional continuation** [Score: 4.2] — Replace the pointwise trajectory
   by a probability measure near collision, connecting to the quantum L² solutions
   (Frauenfelder-Weber 2024 limit-circle analysis).

4. **Painlevé analysis** [Score: 3.6] — An algorithmic test that definitively classifies
   the movable singularity type. High feasibility (can be done symbolically).

### Recommended First Actions

1. Implement Painlevé analysis (most feasible, provides definitive classification)
2. Construct the Bendixson-Dumortier blow-up sequence (highest theoretical promise)
3. Formalize the angular averaging as a rigorous homogenization (partially validated)
4. Solve the Weber-Schrödinger equation numerically to understand the quantum continuation

