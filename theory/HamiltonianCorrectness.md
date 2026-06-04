# Hamiltonian Correctness Note: `H_naive` vs. the Weber Hamiltonian

> **Status: review finding, not yet acted on.** This note documents a
> physics-correctness issue identified in a v0.5.1 whole-package review and
> traced to the companion paper. No code has been changed in response to it. It
> is recorded here so the finding is not lost and so anyone can independently
> reproduce it. See the "Reproduction" section for runnable scripts.

## Summary

The default system built by `HamiltonianSystem(n, dims)` integrates

$$
H_{\text{naive}} = \sum_i \frac{\lVert p_i \rVert^2}{2 m_i}
+ \sum_{i<j} \frac{q_i q_j}{r_{ij}}\left(1 - \frac{\dot r_{ij}^2}{2c^2}\right),
\qquad
\dot r_{ij} \equiv \hat r_{ij} \cdot \left(\frac{p_i}{m_i} - \frac{p_j}{m_j}\right)
$$

— that is, it substitutes $v = p/m$ directly into the **velocity-dependent**
Weber potential (see [`src/hamiltonian/builders/weber.jl`](../src/hamiltonian/builders/weber.jl)).

Because the Weber potential depends on velocity, the canonical momentum is
**not** $m\,v$:

$$
p_i = \frac{\partial L}{\partial v_i} = m_i v_i - \frac{\partial S}{\partial v_i},
\qquad
S = \sum_{i<j} \frac{q_i q_j}{r}\left(1 + \frac{\dot r^2}{2c^2}\right)
\quad \text{(Lagrangian sign)}
$$

So the naive $v \to p/m$ substitution does **not** produce the Legendre transform
of the Weber Lagrangian. $H_{\text{naive}}$ is a *different dynamical system*. It
agrees with Weber electrodynamics only in the Coulomb limit $c \to \infty$; at
finite $c$ it deviates at $O(1/c^2)$ — which is exactly the order of the Weber
correction itself (the entire reason the theory is interesting).

This is consistent with the repo's own theory doc
[`theory/WeberElectrodynamics.md`](WeberElectrodynamics.md): the explicit
equations of motion there (the section "Equations of motion", e.g. the
$\dot x_1$ line) carry a $+$ velocity–momentum coupling — the correct Legendre
transform — whereas the compiled `dq_dt` carries the opposite $-$ sign. The
doc's prose $H = \sum_i T_i + \sum_{i<j} U_{ij}$ reads like $H_{\text{naive}}$,
but its explicit EOM correctly assume $p \neq m v$. The same cannot be said of
the **companion paper**: the paper's equations of motion match the code (the $-$
sign), not this doc — the error originates in the paper's derivation, which the
code then faithfully implements (see *Origin: the paper's derivation* below). So
this is not a code-only slip: among the repo's three statements of the dynamics,
only this theory markdown's velocity equation is correct.

## The correct Hamiltonian

$S = U_{\text{Coulomb}} + \sum_{i<j} (q_i q_j/r)\,\dot r^2/(2c^2)$ is **quadratic
in the velocities**, so the Legendre transform is closed-form and clean:

$$
H_{\text{Weber}}(q, p) = \tfrac{1}{2}\, p^{\top} M(q)^{-1} p + \sum_{i<j} \frac{q_i q_j}{r}
$$

- The **potential is pure Coulomb**.
- The **entire** velocity correction moves into a configuration-dependent
  effective-mass matrix $M(q) = \operatorname{diag}(m_1 I, \dots, m_N I) - K(q)$,
  where $K(q)$ is assembled from the per-pair blocks
  $\dfrac{q_i q_j}{r c^2}\,\hat r \otimes \hat r$ acting on the relative
  coordinate of each pair.
- This $H$ is **non-separable** — precisely the class the package's
  semi-explicit (Tao / Jayawardana–Ohsawa) integrator is designed for, so the
  existing integrator could integrate it unchanged.
- Expanding
  $M(q)^{-1} = \operatorname{diag}(1/m) + \operatorname{diag}(1/m)\,K\,\operatorname{diag}(1/m) + O(1/c^4)$
  reproduces $v_i = \dfrac{p_i}{m_i} + \dfrac{q_i q_j}{m_i r c^2}\,\dot r\,\hat r + \dots$
  — i.e. the $+$ coupling of the theory doc.

For $N = 2$ the inverse is closed-form (a $2\times 2$ system along $\hat r$ in
relative coordinates). For general $N$ it is an $Nd \times Nd$ solve evaluated
inside the EOM.

## Origin: the paper's derivation

The defect is not introduced in the code — the code faithfully implements the
companion paper
[`Computational-Weber-Electrodynamics.tex`](../papers/Computational-Weber-Electrodynamics/Computational-Weber-Electrodynamics.tex).
Two specific steps there produce $H_{\text{naive}}$ rather than the Weber
Hamiltonian:

1. **Naive substitution.** The paper's eq (H) writes
   $H = \tfrac{1}{2} m_1 \vec v_1^2 + \tfrac{1}{2} m_2 \vec v_2^2
   + \dfrac{q_1 q_2}{r}\left(1 - \dfrac{\dot r^2}{2c^2}\right)$, then instructs
   that *"each velocity component $\dot x_i$ is replaced by $p_{x_i}/m_i$."* For
   a velocity-dependent potential this is not the Legendre transform — the
   canonical momentum is not $m\,v$.
2. **Kinetic vs. canonical momentum.** The paper's appendix defines
   $p_{x_i} = m_i \dot x_i$ (kinetic momentum, hence $\dot x_i = p_{x_i}/m_i$),
   but Hamilton's equations $\dot q = \partial H/\partial p$,
   $\dot p = -\partial H/\partial q$ hold only for the **canonical** momentum
   $p_{x_i} = \partial L/\partial \dot x_i = m_i \dot x_i - \alpha_{x_i}$, with
   $\alpha_{x_i} = \dfrac{q_1 q_2}{c^2}\dfrac{\dot r\,(x_i - x_j)}{r^2}$. The two
   momenta differ at $O(1/c^2)$.

The two readings give opposite signs for the velocity equation:

$$
\dot x_1 = \frac{1}{m_1}\left(p_{x_1} - \alpha_x\right) \quad \text{(paper, matches code)}
\qquad \text{vs.} \qquad
\dot x_1 = \frac{1}{m_1}\left(p_{x_1} + \alpha_x\right) \quad \text{(correct Legendre transform)}
$$

The paper is therefore internally inconsistent in two ways: its appendix
$p = m v$ would give $\dot x_1 = p_{x_1}/m_1$ (no $\alpha$ at all), contradicting
its own $\dot x_1 = (p_{x_1} - \alpha_x)/m_1$; and integrating its equations of
motion does **not** reproduce its own stated Weber force law (the paper's eq 2),
diverging at $O(1/c^2)$ (see the RK4 comparison below). The paper labels the
naive substitution itself "the Legendre transformation", which is the conceptual
slip.

A one-point check ($m_1 = m_2 = 1$, $q_1 q_2 = -1$, $c = 1$, $q = (1,0,0,0)$,
$p = (0.3, 0.1, -0.2, 0)$, so $\dot r = 0.5$ and $\alpha_x = -0.5$):

| source | $\dot x_1$ |
|---|---|
| compiled code (`dq_dt[1]`) | $0.8$ |
| paper, $(p_{x_1} - \alpha_x)/m_1$ | $0.8$ (matches code) |
| `theory/WeberElectrodynamics.md`, $(p_{x_1} + \alpha_x)/m_1$ | $-0.2$ (correct) |
| paper appendix $p = mv \Rightarrow p_{x_1}/m_1$ | $0.3$ (contradicts the paper's own EOM) |

> **Certified vs. not.** The velocity equation / substitution error above is
> rigorously verified (analytically and numerically). The momentum equation
> $\dot p$ is written identically in the paper and in
> `theory/WeberElectrodynamics.md`; this note does **not** separately certify
> that $\dot p$ formula equals the exact $-\partial H_{\text{Weber}}/\partial q$.

## Why the current test suite does not catch it

- `test/test_physics.jl` "Energy conservation – Weber" conserves
  `weber_energy_2body_2d` (`test/test_utils.jl`), which **is** $H_{\text{naive}}$
  (it uses $v = p/m$). The integrator is built to conserve $H_{\text{naive}}$, so
  the test is self-referential: it verifies symplecticity, not Weber correctness.
- The "Weber force decomposition – consistency" test checks only that the
  *vector form* and *radial form* of the post-hoc force diagnostic agree with
  each other (an algebraic identity). It never compares $m\,a$ along the
  trajectory to the Weber force.
- `compute_pair_force_timeseries` (`src/statistics/forces.jl`) reconstructs the
  *true* Weber force from $v = p/m$ and finite-difference accelerations — which
  is **not** the force the integrator actually applied (`dp_dt` of
  $H_{\text{naive}}$). For finite $c$ with $\dot r \neq 0$ they differ.

A meaningful regression would compare `solve(...)` against an independent
integration of the Weber **force law** at finite $c$ (see below) and assert
convergence under `dt` refinement — not just $H$ conservation.

## Reproduction

### 1. Equation-of-motion sign check

```julia
# weber_verify.jl
using WeberElectrodynamics

sys = HamiltonianSystem(2, 2)
m1, m2 = 1.0, 1.0
Q1, Q2 = 1.0, -1.0
c = 1.0
q = [1.0, 0.0, 0.0, 0.0]        # r1=(1,0), r2=(0,0) => R=(1,0), r=1
p = [0.3, 0.1, -0.2, 0.0]
params = [m1, m2, Q1, Q2, c]
kappas = [1.0]

dq = zeros(4); dp = zeros(4)
sys.dq_dt_compiled(dq, q, p, 0.0, params, kappas)

x1,y1,x2,y2 = q
dx = x1-x2; r = abs(dx)
rdot_p = dx*((p[1]/m1)-(p[3]/m2))/r          # ṙ using v = p/m

pred_naive_minus = p[1]/m1 - (Q1*Q2/(m1*c^2))*rdot_p*dx/r^2   # H_naive (code)
pred_doc_plus    = (1/m1)*(p[1] + (Q1*Q2/c^2)*rdot_p*dx/r^2)  # theory doc EOM

@show dq[1] pred_naive_minus pred_doc_plus
# dq[1] = 0.8  (matches H_naive);  doc formula predicts -0.2
```

### 2. Trajectory vs. the true Weber force law ($1/c^2$ scaling)

Integrate the standard Weber force

$$
\vec F = \frac{Q_1 Q_2}{r^2}\,\hat r
\left(1 + \frac{\vec v \cdot \vec v + \vec r \cdot \vec a - \tfrac{3}{2}(\hat r \cdot \vec v)^2}{c^2}\right)
$$

with an independent RK4 (the force is central, so the implicit $\vec r \cdot \vec a$
term is solved algebraically), starting at an apsis ($\dot r = 0$, so the package
and the reference share identical initial states), and compare the final
positions.

```julia
# weber_ref.jl
using WeberElectrodynamics, Printf

function weber_accel(s, m1,m2,Q1,Q2,c)
    x1,y1,x2,y2, vx1,vy1,vx2,vy2 = s
    rx,ry = x1-x2, y1-y2
    r = hypot(rx,ry); rhx,rhy = rx/r, ry/r
    vx,vy = vx1-vx2, vy1-vy2
    vdotv = vx^2+vy^2; rhdotv = rhx*vx+rhy*vy
    mu = m1*m2/(m1+m2)
    A = (Q1*Q2/r^2)*(1 + (vdotv - 1.5*rhdotv^2)/c^2)
    B = Q1*Q2/(r*mu*c^2)
    G = A/(1-B)                       # solve implicit r·a (central force)
    F1x,F1y = G*rhx, G*rhy
    return (F1x/m1, F1y/m1, -F1x/m2, -F1y/m2)
end
deriv(s,p...) = (s[5],s[6],s[7],s[8], weber_accel(s,p...)...)
function rk4(s0, dt, n, p)
    s = collect(Float64, s0)
    for _ in 1:n
        k1=collect(deriv(s,p...)); k2=collect(deriv(s .+ dt/2 .*k1,p...))
        k3=collect(deriv(s .+ dt/2 .*k2,p...)); k4=collect(deriv(s .+ dt .*k3,p...))
        s = s .+ (dt/6).*(k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    s
end

m1,m2,Q1,Q2 = 1.0,1.0,1.0,-1.0; w = 0.3; T = 2.0
for c in (3.0, 30.0)
    p = (m1,m2,Q1,Q2,c)
    s0 = (1.0,0.0,-1.0,0.0, 0.0,w,0.0,-w)            # apsis, ṙ=0
    sref = rk4(s0, 1e-5, round(Int,T/1e-5), p)
    sys = HamiltonianSystem(2,2)
    prob = HamiltonianProblem(sys,(0.0,T),[1.0,0.0,-1.0,0.0],[0.0,m1*w,0.0,-m2*w];
                              masses=[m1,m2],charges=[Q1,Q2],c=c,dt=1e-4)
    qend = solve(prob).q[end]
    d = sqrt((qend[1]-sref[1])^2+(qend[2]-sref[2])^2+(qend[3]-sref[3])^2+(qend[4]-sref[4])^2)
    @printf("c=%5.1f  |Δr_final| pkg-vs-WeberForceLaw = %.3e\n", c, d)
end
```

Observed:

| $c$ | $\lvert \Delta r_{\text{final}} \rvert$ (package vs. true Weber) |
|---|---|
| 3.0 | 1.46e-01  (≈15% of the orbit scale) |
| 30.0 | 1.53e-03 |

The discrepancy scales as $1/c^2$ (10× larger $c$ ⇒ ~95× smaller error) at
fixed `dt`. That is an $O(1/c^2)$ *physics* discrepancy — the Weber-correction
order — not integrator truncation (truncation is governed by `dt`, identical in
both runs). Both runs converge to the Coulomb/Kepler orbit as $c \to \infty$,
confirming the reference is correct and the divergence is the Weber term itself.

## Scope of impact

- Affects the package's central claim and any finite-$c$ quantitative result
  (orbital precession, eccentric orbits), including the showcased "precessing
  ellipse" at $c=4$.
- The defect **originates in the companion paper's derivation** (the eq (H)
  naive substitution and the kinetic-vs-canonical momentum conflation), not
  merely the code — the code faithfully implements the paper. Among the repo's
  three statements of the dynamics, only `theory/WeberElectrodynamics.md`'s
  velocity equation is correct; the paper and the code share the error.
- The **Zöllner** extension (`src/hamiltonian/builders/zollner.jl`) is built on
  the same naive potential and inherits the same issue.
- Regularization lifts the same $H_{\text{naive}}$ pair Hamiltonian.
- **Not affected:** the Coulomb/Kepler limit ($c \to \infty$) is correct; the
  semi-explicit integrator is a valid symplectic method *for the Hamiltonian it
  is given*; the software architecture, regularization machinery, statistics,
  and tests are sound as software.

## If/when a fix is desired

Add a correct `weber_hamiltonian` builder producing
$H = \tfrac{1}{2} p^{\top} M(q)^{-1} p + U_{\text{Coulomb}}$, route the default
`HamiltonianSystem(n, dims)` through it (keeping $H_{\text{naive}}$ available only
if explicitly requested), re-derive the Zöllner correction on the corrected
potential, and add the force-law comparison above as a regression test asserting
`dt`-convergence (not just $H$ conservation).
