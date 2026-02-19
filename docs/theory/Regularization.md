# Collision Regularization in 1D, 2D, and 3D $N$-Body Hamiltonian Systems

## Goal

Regularization removes collision singularities in pairwise $1/r$-type interactions while keeping a Hamiltonian structure. The core idea is:

1. Replace singular relative coordinates by collision-resolving coordinates.
2. Rescale physical time so trajectories slow down near collisions in physical time but remain smooth in a fictitious time.

This document is self-contained and independent of any specific implementation.

## Singular Structure of the $N$-Body Problem

Consider a Hamiltonian of the form

$$
H(q,p,t)=\sum_{i=1}^N \frac{\|p_i\|^2}{2m_i}+\sum_{1\le i\lt j\le N}\frac{\kappa_{ij}}{r_{ij}}+R(q,p,t),
$$

where

$$
r_{ij}=\|q_i-q_j\|,
$$

$\kappa_{ij}$ may have either sign, and $R$ is the remaining non-collision-singular part. Binary collision singularities occur at $r_{ij}=0$.

For each close pair $(i,j)$, use center-of-mass and relative variables:

$$
R_{ij}=\frac{m_i q_i+m_j q_j}{m_i+m_j},\quad r_{ij}=q_i-q_j,
$$

$$
P_{ij}=p_i+p_j,\quad p_{ij}=\mu_{ij}\left(\frac{p_i}{m_i}-\frac{p_j}{m_j}\right),\quad
\mu_{ij}=\frac{m_i m_j}{m_i+m_j}.
$$

The singular part is always in the relative coordinate $r_{ij}$, so regularization is applied there.

## Time Regularization in Extended Phase Space

Introduce conjugate pair $(t,p_t)$ and a positive monitor function $g(q,p,t)>0$. Define

$$
\mathcal{K}(q,p,t,p_t)=g(q,p,t)\,\big(H(q,p,t)+p_t\big).
$$

Hamilton's equations in fictitious time $\tau$ are

$$
q'=\frac{\partial \mathcal{K}}{\partial p},\quad
p'=-\frac{\partial \mathcal{K}}{\partial q},\quad
t'=\frac{\partial \mathcal{K}}{\partial p_t}=g,\quad
p_t'=-\frac{\partial \mathcal{K}}{\partial t},
$$

where $'\equiv d/d\tau$. On the physical manifold $H+p_t=0$:

$$
q'=g\,\frac{\partial H}{\partial p},\quad
p'=-g\,\frac{\partial H}{\partial q},\quad
\frac{dt}{d\tau}=g.
$$

Hence

$$
dt=g\,d\tau.
$$

Choosing $g\to 0$ near collision slows physical time and removes blow-ups in the transformed vector field.

For discrete midpoint-type schemes, a common practical choice is to evaluate $g$ at substep start and keep it fixed over that substep. This preserves a consistent physical-time increment while retaining close-encounter slowing from the $g$-dependence between substeps.

## 1D Binary Regularization

For a 1D relative coordinate $x$, define $r=|x|$ and local chart sign $s\in\{+1,-1\}$:

$$
x=s\,u^2,\quad r=u^2.
$$

Canonical one-form preservation $p_x\,dx=p_u\,du$ gives

$$
p_u=2su\,p_x,\quad
p_x=\frac{p_u}{2su}.
$$

For pair Hamiltonian

$$
H_{\text{pair}}=\frac{p_x^2}{2\mu}+\frac{\kappa}{|x|}+R(x,p_x,t),
$$

the kinetic term transforms to

$$
\frac{p_x^2}{2\mu}=\frac{p_u^2}{8\mu u^2}.
$$

With Sundman scaling $dt=r\,d\tau=u^2 d\tau$ and fixed energy $E$:

$$
K_{1D}=u^2\big(H_{\text{pair}}-E\big)
=\frac{p_u^2}{8\mu}+\kappa+u^2\big(R-E\big),
$$

which is finite at $u=0$. Crossing $x=0$ is handled by switching chart sign $s$.

## 2D Levi-Civita Regularization

Let $\mathbf{r}=(x,y)$, $u=(u_1,u_2)$, and define

$$
x=u_1^2-u_2^2,\quad y=2u_1u_2.
$$

Equivalent complex form: $z=x+iy=(u_1+i u_2)^2$. Then

$$
r=\sqrt{x^2+y^2}=u_1^2+u_2^2=\|u\|^2.
$$

Jacobian:

$$
J=\frac{\partial (x,y)}{\partial (u_1,u_2)}
=2\begin{bmatrix}u_1&-u_2\\u_2&u_1\end{bmatrix},
\quad
JJ^T=4r\,I_2.
$$

Canonical momenta from $p\cdot d\mathbf{r}=U\cdot du$:

$$
U=J^T p,\quad
U_1=2(u_1p_x+u_2p_y),\quad
U_2=2(-u_2p_x+u_1p_y).
$$

Therefore

$$
\|p\|^2=\frac{\|U\|^2}{4r},\quad
\frac{\|p\|^2}{2\mu}=\frac{\|U\|^2}{8\mu r}.
$$

With $dt=r\,d\tau$:

$$
K_{2D}=r\big(H_{\text{pair}}-E\big)
=\frac{\|U\|^2}{8\mu}+\kappa+r\big(R-E\big),
$$

regular at $r=0$. The map is two-to-one: $(u_1,u_2)$ and $(-u_1,-u_2)$ represent the same physical point.

## 3D Kustaanheimo-Stiefel (KS) Regularization

For 3D relative coordinate $\mathbf{r}=(x_1,x_2,x_3)$, introduce $u=(u_1,u_2,u_3,u_4)\in\mathbb{R}^4$:

$$
x_1=u_1^2-u_2^2-u_3^2+u_4^2,
$$

$$
x_2=2(u_1u_2-u_3u_4),
$$

$$
x_3=2(u_1u_3+u_2u_4).
$$

Then

$$
r=\sqrt{x_1^2+x_2^2+x_3^2}=u_1^2+u_2^2+u_3^2+u_4^2=\|u\|^2.
$$

Let $J=\partial x/\partial u\in\mathbb{R}^{3\times 4}$. It satisfies

$$
JJ^T=4r\,I_3.
$$

Define 4D momentum $U$ by

$$
U=J^T p+\lambda\,n(u),\quad n(u)=(u_4,-u_3,u_2,-u_1)^T.
$$

The KS bilinear constraint removes gauge freedom:

$$
\Psi(u,U)=u_4U_1-u_3U_2+u_2U_3-u_1U_4=0.
$$

On $\Psi=0$, $U=J^T p$ and

$$
\|p\|^2=\frac{\|U\|^2}{4r},\quad
\frac{\|p\|^2}{2\mu}=\frac{\|U\|^2}{8\mu r}.
$$

With $dt=r\,d\tau$:

$$
K_{3D}=r\big(H_{\text{pair}}-E\big)
=\frac{\|U\|^2}{8\mu}+\kappa+r\big(R-E\big),
$$

with additional constraint $\Psi=0$. Binary collisions map to $u=0$ without singular derivatives in $\tau$.

## Lifting to $N$-Body Systems

For $N>2$, regularization is applied to close binary subsystems in relative coordinates, then coupled back to the full system.

Global time scaling uses a positive $g$ depending on all pair distances. Common families:

$$
dt=g\,d\tau,\quad
g=r_{ab}\ \text{(single active pair)},
$$

$$
g=\left(\sum_{i\lt j}\frac{w_{ij}}{r_{ij}}\right)^{-1}\ \text{(global close-encounter monitor)},
$$

$$
g=\frac{1}{\alpha T+\beta\Omega+\gamma},\quad
\Omega=\sum_{i\lt j}\frac{w_{ij}}{r_{ij}},\quad
\alpha,\beta,\gamma\ge 0.
$$

In transformed canonical variables $(Q,P)$:

$$
Q'=g\,\frac{\partial H}{\partial P},\quad
P'=-g\,\frac{\partial H}{\partial Q},\quad
t'=g.
$$

When multiple KS pairs are active, constraints $\Psi_a(Q,P)=0$ are enforced for each pair.

## Practical Activation and Chain Fallback

In numerical $N$-body work, regularization is typically activated only for close encounters, not globally:

1. Detect candidate close pairs with an activation radius $r_{\text{on}}$.
2. Keep regularization active until all relevant separations exceed a larger radius $r_{\text{off}}>r_{\text{on}}$ (hysteresis).
3. If the active encounter involves one close binary, apply pair regularization directly.
4. If several close pairs share particles (overlapping encounter graph), switch to a local chain-style subsystem and evolve it with a global monitor $g$ (for example $g=1/\Omega$).

This practical policy preserves the main benefit of regularization near singular events while keeping the far-field integration in ordinary coordinates.

## Multiple Close Encounters and Chain Coordinates

When several nearby bodies interact simultaneously, direct pairwise relative coordinates can become ill-conditioned. A standard remedy is chain coordinates:

1. Choose an ordered chain of bodies $i_1,\dots,i_M$ through the closest separations.
2. Use link vectors $s_k=q_{i_{k+1}}-q_{i_k}$ as primary relative coordinates.
3. Regularize each short link (1D/2D/3D map as appropriate) and evolve with a global $g$.

This reduces subtraction of large nearby numbers and improves robustness in few-body subsystems embedded in large $N$-body dynamics.

## Key Mathematical Outcomes

- Binary collision singularities are removed from the transformed vector field.
- Physical time is recovered by integrating $t'=g$.
- Canonical structure is preserved by coordinate maps plus extended-phase-space Hamiltonian flow.
- 1D, 2D, and 3D binary regularizations share the same pattern:
  singular radius $r$ becomes quadratic in regularized coordinates and is canceled by $dt=r\,d\tau$.
- In $N$-body systems, regularization is exact for isolated binary collisions and remains effective for close-encounter subsystems with suitable global $g$ and coordinate organization.

## References

- Levi-Civita regularization and Sundman/Poincare time transform:
  [A Family of Canonical Transformations for Regularizing Binary Collisions in the Spatial Kepler Problem](https://www.mdpi.com/2227-7390/10/11/1852)
- KS transformation formulas, bilinear relation, and Sundman map:
  [Kustaanheimo-Stiefel Transformation with an Arbitrary Defining Vector](https://www.mdpi.com/2227-7390/10/20/3787)
- Few-body regularization context (KS, chain, and logarithmic Hamiltonian approaches):
  [A New Class of Symplectic Integrators for the Gravitational N-body Problem](https://arxiv.org/abs/0803.4441)
- Time-transformed Hamiltonian form for many-singularity few-body integration:
  [Algorithmic Regularization with Velocity-Dependent Forces](https://arxiv.org/abs/1306.0197)
