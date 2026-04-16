# Nekhoroshev stability time estimates -- explicit formulas and numbers

## 1. General Nekhoroshev bound

For a near-integrable Hamiltonian $H = H_0(I) + \varepsilon H_1(I, \phi)$ with
$H_0$ satisfying a steepness condition (convexity suffices), Nekhoroshev's theorem
guarantees:

$$
|I(t) - I(0)| \leq R_0 \varepsilon^b, \quad |t| \leq T_0 \exp\!\left(\left(\frac{\varepsilon_0}{\varepsilon}\right)^a\right)
$$

### Exponents

| Source | $a$ | $b$ | Notes |
|--------|-----|-----|-------|
| Nekhoroshev (1977) | $1/(2n)$ | $1/(2n)$ | Original, $n$ = DOF |
| Lochak (1992) | $1/(2n)$ | $1/(2n)$ | Simultaneous approx. method |
| Poschel (1993) | $1/(2(n-1))$ | $1/(2n)$ | Improved for convex $H_0$ |
| Niederman (2004) | $1/(2n)$ | $1/(2n)$ | Optimal for steep systems |
| Bounemoura-Niederman (2012) | $1/(2\lfloor n/2\rfloor)$ | -- | Quasi-convex, best known |

For the Weber system at $c = 1$, the perturbation parameter is $\varepsilon = 1/c^2 = 1$.

### Pre-factor constants

$T_0$ and $\varepsilon_0$ depend on:
- The analyticity width $\sigma$ of $H_1$ (distance to nearest singularity in complexified angle space)
- The steepness constants of $H_0$
- The size of the action domain

For order-of-magnitude estimates, we take $T_0 \sim \varepsilon_0 \sim 1$ in natural units,
recognizing this underestimates the true pre-factors.

---

## 2. Application to the 2+/2- quasi-bound candidates

### 2.1. Effective degrees of freedom

| Configuration | Raw DOF (3D, after COM) | Discrete symmetry reduction | $n_{\text{eff}}$ |
|--------------|-------------------------|---------------------------|-----------------|
| Double-orbiter (3D) | 9 | Klein-four ($\mathbb{Z}_2 \times \mathbb{Z}_2$): 2 generators | 6 |
| Rhombus (2D) | 6 | $D_2$ in-plane | 4 |
| Alternating square (2D) | 6 | $D_4$ broken to $\mathbb{Z}_4$ | 4 |

Note: discrete symmetries constrain the motion to a lower-dimensional invariant
subspace but do not reduce the number of action variables in the same way as
continuous symmetries. The values above are estimates; the true effective DOF
for Nekhoroshev purposes depends on the rank of the frequency map restricted
to the invariant subspace.

### 2.2. Nekhoroshev times for the double-orbiter ($n_\text{eff} = 6$)

Using $a = 1/(2n) = 1/12$ (Nekhoroshev/Lochak):

$$
T_N(c) = T_0 \exp\!\left(\left(\frac{\varepsilon_0}{\varepsilon}\right)^{1/12}\right)
       = T_0 \exp\!\left(\left(\varepsilon_0 \cdot c^2\right)^{1/12}\right)
$$

| $c$ | $\varepsilon$ | $(1/\varepsilon)^{1/12}$ | $T_N / T_0$ |
|-----|--------------|------------------------|-------------|
| 1 | 1.000 | 1.000 | 2.72 |
| 2 | 0.250 | 1.122 | 3.07 |
| 4 | 0.0625 | 1.259 | 3.52 |
| 10 | 0.010 | 1.468 | 4.34 |
| 100 | $10^{-4}$ | 2.154 | 8.62 |
| $10^3$ | $10^{-6}$ | 3.162 | 23.6 |
| $10^6$ | $10^{-12}$ | 10.00 | $2.2 \times 10^4$ |

The bound becomes non-trivial (say $T_N > 100$) only for $c \gtrsim 10^4$.

Using the improved Poschel exponent $a = 1/(2(n-1)) = 1/10$:

| $c$ | $(1/\varepsilon)^{1/10}$ | $T_N / T_0$ |
|-----|------------------------|-------------|
| 1 | 1.000 | 2.72 |
| 4 | 1.320 | 3.74 |
| 100 | 2.512 | 12.3 |
| $10^6$ | 15.85 | $7.7 \times 10^6$ |

Still very weak. The $1/(2n)$ and $1/(2(n-1))$ exponents differ by less than
an order of magnitude for practical $c$ values.

### 2.3. Nekhoroshev times for the rhombus ($n_\text{eff} = 4$)

Using $a = 1/8$:

| $c$ | $(1/\varepsilon)^{1/8}$ | $T_N / T_0$ |
|-----|------------------------|-------------|
| 1 | 1.000 | 2.72 |
| 4 | 1.414 | 4.11 |
| 100 | 3.162 | 23.6 |

Slightly better than the 6-DOF case due to the larger exponent, but still
trivially short at $c = 1$.

---

## 3. Comparison with observations

### 3.1. Table

| Candidate | $n_\text{eff}$ | $a$ | $T_N/T_0$ at $c=1$ | $t^*$ (observed) | Ratio $t^*/T_N$ |
|-----------|----------------|------|---------------------|-------------------|----------------|
| Double-orbiter (3D) | 6 | 1/12 | 2.72 | 566 | 208 |
| Rhombus | 4 | 1/8 | 2.72 | 440 | 162 |

The observed lifetimes exceed the Nekhoroshev lower bounds by two orders of
magnitude. This is **normal** -- Nekhoroshev bounds are worst-case estimates
over all initial conditions and all perturbations of the given size. The actual
stability time for a specific IC in a specific system can be (and usually is)
much longer.

### 3.2. What Nekhoroshev bounds DO tell us

1. **At $c = 1$, the Nekhoroshev machinery provides no useful constraint.** The
   bound $T_N \sim 3$ is smaller than a single orbital period.

2. **At $c = 100$, useful bounds begin to emerge.** $T_N \sim 10$ orbital periods
   for 6 DOF. This is the threshold where the perturbative framework starts to
   bite.

3. **The observed $t^* \sim 500$ is not explained by Nekhoroshev theory.** It
   must arise from specific geometric properties of the phase space (basin
   structure, resonance web topology) rather than from generic perturbative
   estimates.

---

## 4. Arnold diffusion rate estimates

### 4.1. Diffusion coefficient

In the Arnold-diffusion regime, the characteristic diffusion coefficient in action
space is

$$
D_I \sim \varepsilon^2 \cdot \omega^2 \cdot \exp\!\left(-\frac{C}{\sqrt{\varepsilon}}\right)
$$

for analytic perturbations (exponentially small splitting). At $\varepsilon = 1$
this gives $D_I \sim e^{-C}$ -- no exponential suppression.

From the measured Lyapunov exponent, a simpler estimate:

$$
D_I \sim \lambda_\max \cdot (\Delta I_\text{layer})^2
$$

where $\Delta I_\text{layer}$ is the action-space width of the stochastic layer.

### 4.2. Escape time from diffusion

$$
t_\text{escape} \sim \frac{(\Delta I_\text{basin})^2}{D_I}
$$

For the double-orbiter:
- $\Delta I_\text{basin} / I_0 \sim 0.03$ (from parameter scans: orb varies by $\pm 0.02$ around 1.30)
- $\lambda_\max = 0.0205$
- Assuming $\Delta I_\text{layer} \sim \Delta I_\text{basin}$ (fast diffusion regime):
  $D_I \sim 0.0205 \times (0.03)^2 \sim 1.8 \times 10^{-5}$ (normalized)
- $t_\text{escape} \sim (0.03)^2 / (1.8 \times 10^{-5}) \sim 50$

This underestimates the observed $t^* \approx 566$ by an order of magnitude,
suggesting the effective basin in action space is wider than the parameter-scan
width, or the diffusion has a sub-linear (sticky) component near island remnants.

For the rhombus:
- $\Delta I_\text{basin} / I_0 \sim 0.10$ (broader basin from the $b$-scan)
- $\lambda_\max = 0.18$
- $D_I \sim 0.18 \times (0.10)^2 \sim 1.8 \times 10^{-3}$
- $t_\text{escape} \sim (0.10)^2 / (1.8 \times 10^{-3}) \sim 5.6$

Again underestimates (observed $\sim 440$), by nearly two orders of magnitude.
The discrepancy suggests the basin structure has fractal/sticky boundaries
that slow transport significantly beyond the simple random-walk estimate.

---

## 5. Predictions for $c$-continuation experiments

| $c$ | $\varepsilon$ | Predicted $\lambda_\max$ (power law $\propto 1/c^2$) | Predicted $T_N$ (Nekh, $a=1/12$) | Predicted $t^*$ (fast diff, $\propto c^2$) |
|-----|--------------|------------------------------------------------------|-----------------------------------|--------------------------------------------|
| 1 | 1.0 | 0.020 (measured) | 2.7 | 566 (measured) |
| 2 | 0.25 | 0.005 | 3.1 | 2264 |
| 4 | 0.0625 | 0.00125 | 3.5 | 9056 |
| 10 | 0.01 | $2\times 10^{-4}$ | 4.3 | 56600 |
| 100 | $10^{-4}$ | $2\times 10^{-6}$ | 8.6 | $5.66\times 10^6$ |

If the fast-diffusion scaling $t^* \propto c^2$ holds, the double-orbiter at $c = 4$
should survive to $t^* \sim 9000$. If instead Nekhoroshev scaling takes over at some
$c_\text{crit}$, $t^*$ would grow much faster beyond that point.

**Recommended test**: run the 3D double-orbiter at $c = 2$ and $c = 4$ with
$t_\max = 10000$ and measure $t^*$. If $t^*(c=2) \approx 2300$ and
$t^*(c=4) \approx 9000$, the fast-diffusion scaling is confirmed. If $t^*$
grows faster than $c^2$, Nekhoroshev exponential stability has been reached.
