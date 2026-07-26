# Semi-Explicit Symplectic Integrator for Non-Separable Hamiltonians

A symplectic numerical integrator for non-separable Hamiltonians $H(q,p)$ combining extended phase space with symmetric projection.

## Hamilton's Equations

$$
\dot{q} = \frac{\partial H}{\partial p}, \quad \dot{p} = -\frac{\partial H}{\partial q}
$$

Note: The momentum derivative function $\dot{p}$ should already include the minus sign.

## Extended Phase Space

- Original system: $(q,p) \in \mathbb{R}^{2d}$ where $q, p \in \mathbb{R}^d$
- Extended system: $(q,x,p,y) \in \mathbb{R}^{4d}$ where $q, x, p, y \in \mathbb{R}^d$
- Invariant subspace $\mathcal{N}$: points where $x = q$ and $y = p$

## Extended System Dynamics

The extended Hamiltonian $\hat{H}(q,x,p,y) = H(q,y) + H(x,p)$ generates the system:

$$
\dot{q} = D_2 H(x,p), \quad \dot{p} = -D_1 H(q,y)
$$

$$
\dot{x} = D_2 H(q,y), \quad \dot{y} = -D_1 H(x,p)
$$

where $D_i$ denotes the partial derivative with respect to the $i$-th argument.

Starting from initial conditions $(q_0, q_0, p_0, p_0) \in \mathcal{N}$, the solution satisfies $(q(t), p(t)) = (x(t), y(t))$ for all $t$, and this pair solves the original Hamilton's equations.

## Flow Maps

### Flow A

Hamiltonian:

$$
\hat{H}_A(q,x,p,y) = H(q,y)
$$

Flow map:

$$
\hat{\Phi}^A_t(q_0, x_0, p_0, y_0) = \left(q_0,\, x_0 + t\,\dot{q}(q_0, y_0),\, p_0 + t\,\dot{p}(q_0, y_0),\, y_0\right)
$$

Variables $q$ and $y$ are frozen; $x$ and $p$ evolve.

### Flow B

Hamiltonian:

$$
\hat{H}_B(q,x,p,y) = H(x,p)
$$

Flow map:

$$
\hat{\Phi}^B_t(q_0, x_0, p_0, y_0) = \left(q_0 + t\,\dot{q}(x_0, p_0),\, x_0,\, p_0,\, y_0 + t\,\dot{p}(x_0, p_0)\right)
$$

Variables $x$ and $p$ are frozen; $q$ and $y$ evolve.

## Second-Order Integrator

Strang splitting:

$$
\hat{\Phi}_{\Delta t} = \hat{\Phi}^A_{\Delta t/2} \circ \hat{\Phi}^B_{\Delta t} \circ \hat{\Phi}^A_{\Delta t/2}
$$

## Higher-Order Integrators

Higher-order methods are constructed by composing the 2nd-order integrator with carefully chosen time steps.

### Triple Jump Composition

Given a 2nd-order method $\hat{\Phi}^{(2)}_{\Delta t}$, construct order $n$ (even) recursively:

$$
\hat{\Phi}^{(n)}_{\Delta t} := \hat{\Phi}^{(n-2)}_{\gamma_3 \Delta t} \circ \hat{\Phi}^{(n-2)}_{\gamma_2 \Delta t} \circ \hat{\Phi}^{(n-2)}_{\gamma_1 \Delta t}
$$

where:

$$
\gamma_1 = \gamma_3 = \frac{1}{2 - 2^{1/(n-1)}}, \quad \gamma_2 = -\frac{2^{1/(n-1)}}{2 - 2^{1/(n-1)}}
$$

This yields a $3^{n/2}$-stage method:
- 4th order: 9 stages
- 6th order: 27 stages

### Suzuki Composition (5-stage)

Alternative composition for order $n$:

$$
\hat{\Phi}^{(n)}_{\Delta t} := \hat{\Phi}^{(n-2)}_{\gamma_5 \Delta t} \circ \hat{\Phi}^{(n-2)}_{\gamma_4 \Delta t} \circ \hat{\Phi}^{(n-2)}_{\gamma_3 \Delta t} \circ \hat{\Phi}^{(n-2)}_{\gamma_2 \Delta t} \circ \hat{\Phi}^{(n-2)}_{\gamma_1 \Delta t}
$$

where:

$$
\gamma_1 = \gamma_2 = \gamma_4 = \gamma_5 = \frac{1}{4 - 4^{1/(n-1)}}, \quad \gamma_3 = -\frac{4^{1/(n-1)}}{4 - 4^{1/(n-1)}}
$$

### Yoshida 6th-Order (7-stage)

A 7-stage 6th-order symmetric composition (Yoshida 1990).

## Projection Matrix

The projection matrix $A \in \mathbb{R}^{2d \times 4d}$ and its transpose $A^T \in \mathbb{R}^{4d \times 2d}$ are:

$$
A = \begin{bmatrix} I_d & -I_d & 0 & 0 \\ 0 & 0 & I_d & -I_d \end{bmatrix}, \quad A^T = \begin{bmatrix} I_d & 0 \\ -I_d & 0 \\ 0 & I_d \\ 0 & -I_d \end{bmatrix}
$$

where $I_d$ is the $d \times d$ identity matrix.

## Algorithm Steps

**Input:** $z_n = (q_n, p_n) \in \mathbb{R}^{2d}$, time step $\Delta t$

1. **Embed to extended space:**

   $$
   Z_n = (q_n, q_n, p_n, p_n) \in \mathbb{R}^{4d}
   $$

2. **Solve for shift vector** $\mu \in \mathbb{R}^{2d}$ such that:

   $$
   \hat{\Phi}_{\Delta t}(Z_n + A^T \mu) + A^T \mu \in \mathcal{N}
   $$

   Equivalently, solve:

   $$
   A\left(\hat{\Phi}_{\Delta t}(Z_n + A^T \mu) + A^T \mu\right) = 0
   $$

3. **Shift starting point:**

   $$
   \hat{Z}_n = Z_n + A^T \mu
   $$

4. **Explicit evolution step:**

   $$
   \hat{Z}_{n+1} = \hat{\Phi}_{\Delta t}(\hat{Z}_n)
   $$

5. **Shift endpoint back to subspace:**

   $$
   Z_{n+1} = \hat{Z}_{n+1} + A^T \mu
   $$

6. **Extract result:**

   $$
   z_{n+1} = (q_{n+1}, p_{n+1})
   $$

   where $Z_{n+1} = (q_{n+1}, q_{n+1}, p_{n+1}, p_{n+1})$

**Output:** $z_{n+1} = (q_{n+1}, p_{n+1})$

## Projection Step Iteration

Define the nonlinear function $f: \mathbb{R}^{2d} \to \mathbb{R}^{2d}$:

$$
f(\mu) = A\left(\hat{\Phi}_{\Delta t}(Z_n + A^T \mu) + A^T \mu\right)
$$

**Initial guess:**

$$
\mu^{(0)} = 0
$$

> **Implementation note:** Using the converged $\mu$ from the previous time step as a warm start is permitted and typically reduces iterations from 2-3 to 1-2.

**Simplified Newton iteration:**

$$
\mu^{(k+1)} = \mu^{(k)} - \frac{1}{4}f(\mu^{(k)})
$$

This is equivalent to a relaxed fixed-point iteration with relaxation factor $\omega = 0.25$, since $Df(0) = 2AA^T = 4I_{2d}$.

### Derivation of the Factor 1/4

The Jacobian of $f$ at $\mu$ is:

$$
Df_{\Delta t}(\mu) = A \cdot D\hat{\Phi}_{\Delta t}(Z_n + A^T \mu) \cdot A^T + AA^T
$$

At $\Delta t = 0$, we have $D\hat{\Phi}_0 = I_{4d}$, so:

$$
Df_0(\mu) = A \cdot I_{4d} \cdot A^T + AA^T = 2AA^T
$$

Computing $AA^T$:

$$
AA^T = \begin{bmatrix} I_d & -I_d & 0 & 0 \\ 0 & 0 & I_d & -I_d \end{bmatrix} \begin{bmatrix} I_d & 0 \\ -I_d & 0 \\ 0 & I_d \\ 0 & -I_d \end{bmatrix} = \begin{bmatrix} 2I_d & 0 \\ 0 & 2I_d \end{bmatrix} = 2I_{2d}
$$

Therefore $Df_0 = 2 \cdot 2I_{2d} = 4I_{2d}$, giving the Newton step $(Df_0)^{-1} = \frac{1}{4}I_{2d}$.

**Stopping criterion:**

$$
\left\|\mu^{(N+1)} - \mu^{(N)}\right\| < \varepsilon
$$

and set $\mu = \mu^{(N)}$.

## Broyden's Method (Alternative)

More robust quasi-Newton approach using rank-1 updates.

### Sherman-Morrison Formula

For an invertible matrix $M \in \mathbb{R}^{n \times n}$ and vectors $u, v \in \mathbb{R}^n$, the inverse of the rank-1 update $M + uv^T$ is:

$$
(M + uv^T)^{-1} = M^{-1} - \frac{M^{-1}uv^T M^{-1}}{1 + v^T M^{-1}u}
$$

This allows efficient computation of $J_k^{-1}$ from $J_{k-1}^{-1}$ in $O(n^2)$ time instead of $O(n^3)$ for full matrix inversion.

### Broyden's Iteration

$$
\mu^{(k+1)} = \mu^{(k)} - J_k^{-1} f(\mu^{(k)})
$$

where $J_0 = 4I_{2d}$ and $J_k$ is maintained via the rank-1 update:

$$
J_k = J_{k-1} + \frac{(\Delta f^{(k)} - J_{k-1}\Delta\mu^{(k)})(\Delta\mu^{(k)})^T}{\|\Delta\mu^{(k)}\|^2}
$$

with:
- $\Delta f^{(k)} = f(\mu^{(k)}) - f(\mu^{(k-1)}) \in \mathbb{R}^{2d}$ (column vector)
- $\Delta\mu^{(k)} = \mu^{(k)} - \mu^{(k-1)} \in \mathbb{R}^{2d}$ (column vector)
- $(\Delta\mu^{(k)})^T \in \mathbb{R}^{1 \times 2d}$ (row vector)
- The numerator is an **outer product** producing a $2d \times 2d$ matrix

Applying the Sherman-Morrison formula to this rank-1 update yields the direct inverse update (maintained iteratively, never inverting $J_k$ explicitly):

$$
J_k^{-1} = J_{k-1}^{-1} + \frac{\Delta\mu^{(k)} - J_{k-1}^{-1}\Delta f^{(k)}}{(\Delta\mu^{(k)})^T J_{k-1}^{-1}\Delta f^{(k)}} \left(\Delta\mu^{(k)}\right)^T J_{k-1}^{-1}
$$

## Key Properties

- Symplectic in original phase space (preserves $dq \wedge dp$)
- Same order of accuracy as underlying extended integrator $\hat{\Phi}$
- Symmetric if $\hat{\Phi}$ is symmetric
- Phase space defect $\|(q,p) - (x,y)\| \sim \varepsilon$ (tolerance)
- Requires few iterations (typically 1-3) for small $\Delta t$
- Projection step overhead per Newton iteration is $O(d)$ (one $\hat{\Phi}$ evaluation plus sparse matrix-vector products with $A$, $A^T$)

## Reference

Jayawardana, B. & Ohsawa, T. (2021). *Semiexplicit Symplectic Integrators for Non-separable Hamiltonian Systems*. [arXiv:2111.10915](https://arxiv.org/abs/2111.10915)

See also: [Regularization.md](Regularization.md), [RegularizedIntegrationDesign.md](RegularizedIntegrationDesign.md), [WeberElectrodynamics.md](WeberElectrodynamics.md).
