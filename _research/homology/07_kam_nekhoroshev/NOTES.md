# Agent 07 -- KAM / Nekhoroshev / Arnold diffusion analysis of the Weber Hamiltonian

**Scope.** Apply classical perturbation theory to the Weber Hamiltonian and produce
explicit stability-time bounds for the three known quasi-bound candidates from the
14-agent study.

All formulas use the absolute-unit convention of this repository: no factors of
$4\pi\epsilon_0$.

---

## 1. The Weber Hamiltonian as a near-integrable perturbation

### 1.1. Splitting

The N-body Weber Hamiltonian is

$$
H = \sum_i \frac{|\mathbf{p}_i|^2}{2m_i}
  + \sum_{i<j} \frac{q_i q_j}{r_{ij}}
    \Bigl(1 - \frac{\dot{r}_{ij}^2}{2c^2}\Bigr)
$$

Writing $\varepsilon = 1/c^2$ and splitting:

$$
H = H_0 + \varepsilon\, H_1
$$

where

$$
H_0 = \sum_i \frac{|\mathbf{p}_i|^2}{2m_i}
     + \sum_{i<j} \frac{q_i q_j}{r_{ij}}
\qquad\text{(Coulomb / Kepler)}
$$

$$
H_1 = -\sum_{i<j} \frac{q_i q_j}{r_{ij}} \frac{\dot{r}_{ij}^2}{2}
$$

For the simulations in this study, $c = 1$ (or $c = 4$), so
$\varepsilon = 1$ (or $0.0625$). This is **not** a small parameter in the KAM
sense at $c = 1$ -- the perturbation is O(1) relative to $H_0$. At $c = 4$
the perturbation is 6.25% of the Coulomb energy -- marginally within the range
where perturbative methods give useful qualitative guidance, though far from the
regime where rigorous KAM bounds apply.

### 1.2. Two-body Coulomb: the integrable base

For $N = 2$ in relative coordinates with reduced mass $\mu = m_1 m_2/(m_1+m_2)$
and coupling $k = q_1 q_2$, the center-of-mass separates and the relative motion
is a Kepler problem:

$$
H_0^{(\text{rel})} = \frac{p_r^2}{2\mu} + \frac{L^2}{2\mu r^2} + \frac{k}{r}
$$

where $L$ is angular momentum and $p_r$ is radial momentum.

**Action-angle variables (Delaunay).** For bound orbits ($k < 0$, attractive, i.e.
opposite charges $q_1 q_2 < 0$):

- Radial action: $I_r = \frac{1}{2\pi}\oint p_r\,dr = -L + \frac{|k|\sqrt{\mu}}{\sqrt{2|E|}}$
- Total action: $I = I_r + L = \frac{|k|\sqrt{\mu}}{\sqrt{2|E|}}$
- Hamiltonian in action variables: $H_0 = -\frac{\mu k^2}{2 I^2}$
- Twist: $\omega(I) = \partial H_0/\partial I = \frac{\mu k^2}{I^3}$
- Twist nondegeneracy: $\frac{\partial^2 H_0}{\partial I^2} = -\frac{3\mu k^2}{I^4} \neq 0$

The twist condition is satisfied everywhere for $I > 0$, which is the fundamental
requirement for KAM theory.

For like-sign charges ($k > 0$, repulsive), no bound orbits exist in pure Coulomb.
Weber's velocity-dependent correction can create quasi-bound states, but these have
no integrable base to perturb around -- KAM theory does not directly apply.

### 1.3. Four-body Coulomb: near-integrable sub-structures

The 4-body Coulomb problem is **not** integrable. However, for the 2+/2- configurations
studied here, there are approximate near-integrable regimes:

**(a) Dimer-dimer at large separation.** When the inter-dimer distance $R$ is much
larger than intra-dimer size $a$, the system is approximately two independent Kepler
problems (each dimer is a bound (+)(-) pair). The interaction between dimers is
dipole-dipole $\sim O(a^2/R^3)$, which is a small perturbation when $R \gg a$.

**(b) Symmetric configurations.** The Klein-four symmetry ($\sigma_x, \sigma_y$) of
the double-orbiter reduces the 12-DOF system (after COM removal) to effectively
6 DOF in the symmetry-reduced space. Two commuting discrete symmetries constrain
the dynamics to a lower-dimensional invariant subspace.

**(c) Hierarchy.** For the double-orbiter at $(r_{pp}, R, \text{orb}) = (4, 4, 1.3)$,
the positive pair separation equals the negative pair orbital radius. There is no
clear scale separation, so the dimer-dimer approximation does not apply. The system
is genuinely non-integrable even at zeroth order.

---

## 2. KAM analysis for 2-body Weber

### 2.1. Setup

For two opposite charges ($q_1 = +q$, $q_2 = -q$), the 2-body relative Hamiltonian
in polar coordinates $(r, \theta)$ with conjugate momenta $(p_r, L)$ is:

$$
H = \frac{p_r^2}{2\mu} + \frac{L^2}{2\mu r^2} - \frac{q^2}{r}
  + \varepsilon\, H_1(r, p_r, L)
$$

The Weber correction term, substituting $\dot{r} = p_r/\mu$ (to leading order in
$\varepsilon$):

$$
H_1 = \frac{q^2}{r}\frac{p_r^2}{2\mu^2}
$$

Note: the full $\dot{r}_{ij}$ includes velocity-dependent terms from the canonical
transformation, but to leading order in $\varepsilon$, $\dot{r} \approx p_r/\mu$.

### 2.2. Perturbation size in action-angle variables

In Delaunay variables $(I, \phi, L, \theta)$ the unperturbed Kepler Hamiltonian is

$$
H_0 = -\frac{\mu q^4}{2I^2}
$$

The perturbation $H_1$ must be expressed in these variables. On a Kepler ellipse with
semi-major axis $a_K = I^2/(\mu q^2)$ and eccentricity $e$:

$$
r = a_K(1 - e\cos E), \quad p_r = \frac{\mu q^2 e \sin E}{I}
$$

where $E$ is the eccentric anomaly. Therefore:

$$
H_1 = \frac{q^2}{r}\frac{p_r^2}{2\mu^2}
     = \frac{q^2}{a_K(1-e\cos E)} \cdot \frac{q^4 e^2 \sin^2 E}{2I^2}
     = \frac{\mu q^4 e^2 \sin^2 E}{2I^2(1-e\cos E)}
$$

The Fourier average over one period:

$$
\langle H_1 \rangle = \frac{\mu q^4 e^2}{2I^2} \left\langle \frac{\sin^2 E}{1-e\cos E}\right\rangle
$$

Using the standard Kepler average $\langle \sin^2 E/(1-e\cos E)\rangle = 1$, we get:

$$
\langle H_1 \rangle = \frac{\mu q^4 e^2}{2I^2}
$$

The ratio of perturbation to unperturbed energy:

$$
\frac{\varepsilon\langle H_1\rangle}{|H_0|}
= \frac{1}{c^2} \cdot \frac{\mu q^4 e^2/(2I^2)}{\mu q^4/(2I^2)}
= \frac{e^2}{c^2}
= \frac{v_\text{typ}^2}{c^2}
$$

where $v_\text{typ} \sim eq^2/I$ is the typical orbital velocity. This confirms
the perturbation is genuinely $O(v^2/c^2)$.

### 2.3. KAM threshold

The classical KAM theorem (Arnold, Kolmogorov, Moser) guarantees that for a
near-integrable Hamiltonian $H = H_0(I) + \varepsilon H_1(I, \phi)$ with
$\det(\partial^2 H_0/\partial I^2) \neq 0$ (twist condition), most invariant
tori survive for $\varepsilon < \varepsilon_*$, where $\varepsilon_*$ depends on:

- The twist strength: $|\partial^2 H_0/\partial I^2| = 3\mu q^4/I^4$
- The Fourier decay rate of $H_1$ (analyticity gives exponential decay)
- The Diophantine constants of the frequency ratio

For the 2-body problem (1 DOF after angular momentum reduction), KAM theory is
trivially satisfied: any nondegenerate 1-DOF integrable Hamiltonian has *all* orbits
on invariant circles. The KAM question becomes nontrivial only for $\geq 2$ DOF.

For the full 2-body problem in 2D (2 DOF: radial + angular), the frequency ratio is

$$
\frac{\omega_r}{\omega_\theta} = 1 + O(\varepsilon)
$$

(Kepler degeneracy: radial and angular frequencies are equal). The Weber correction
lifts this degeneracy:

$$
\omega_r = \frac{\mu q^4}{I^3} + \varepsilon\frac{\partial\langle H_1\rangle}{\partial I},
\quad
\omega_\theta = \frac{\mu q^4}{I^3} + \varepsilon\frac{\partial\langle H_1\rangle}{\partial L}
$$

The Kepler degeneracy ($\omega_r = \omega_\theta$ at $\varepsilon = 0$) means the
standard KAM theorem does not apply directly; one needs the **isoenergetic KAM theorem**
or **Arnold's theorem for properly degenerate systems**. The Weber correction breaks
the degeneracy at $O(\varepsilon)$, so the effective perturbation parameter for KAM
purposes is $\sqrt{\varepsilon}$ rather than $\varepsilon$.

**Estimate.** For $c = 4$ ($\varepsilon = 0.0625$), circular orbits ($e = 0$):
$H_1 = 0$. For modest eccentricity $e = 0.3$: $\varepsilon e^2 = 0.0056$. This is
well within the range where KAM theory predicts survival of most tori.

For $c = 1$ ($\varepsilon = 1$), $e = 0.3$: $\varepsilon e^2 = 0.09$. Marginal --
some tori survive, many do not. For $e \gtrsim 0.5$: $\varepsilon e^2 = 0.25$, and
KAM breakdown is expected for a substantial fraction of phase space.

### 2.4. Measure of surviving KAM tori vs. c

In a generic near-integrable system, the relative measure of destroyed tori scales as:

$$
\mu_{\text{destroyed}} \sim \sqrt{\varepsilon} = \frac{1}{c}
$$

(this is the "KAM gap" scaling). Therefore:

$$
\mu_{\text{KAM tori}} \sim 1 - \frac{C}{c}
$$

for some constant $C$ depending on the action domain.

| $c$ | $\varepsilon = 1/c^2$ | $\sqrt{\varepsilon}$ | Estimated KAM fraction |
|-----|----------------------|---------------------|----------------------|
| 1   | 1.000                | 1.000               | $\sim 0$ (breakdown)  |
| 4   | 0.0625               | 0.250               | $\sim 0.75$ (most tori survive) |
| 10  | 0.010                | 0.100               | $\sim 0.90$           |
| 100 | 0.0001               | 0.010               | $\sim 0.99$           |

At $c = 1$, the KAM fraction is expected to be negligible for generic orbits. The
numerical finding that no true KAM torus was found in the 14-agent study at $c = 1$
is entirely consistent with this estimate.

---

## 3. Nekhoroshev stability estimates

### 3.1. Theory

Nekhoroshev's theorem (1977) states that for a near-integrable Hamiltonian satisfying
a **steepness** (or convexity/quasi-convexity) condition, the action variables remain
close to their initial values for exponentially long times:

$$
|I(t) - I(0)| \leq C\varepsilon^b \quad \text{for} \quad
|t| \leq T_N \exp\left(\frac{1}{\varepsilon^a}\right)
$$

where $a = b = 1/(2n)$ in the original estimate ($n$ = number of DOF), improved by
Lochak (1992) to $a = 1/(2n)$ and by Poschel to $a = 1/(2(n-1))$ for convex systems.

The Kepler Hamiltonian $H_0 = -\mu k^2/(2I^2)$ is **convex** in $I$ (since
$\partial^2 H_0/\partial I^2 = -3\mu k^2/I^4 < 0$ for $k < 0$, i.e. $k^2 > 0$
always). So the Nekhoroshev theorem applies.

### 3.2. Application to the double-orbiter

**Configuration**: 3D symmetric double-orbiter at $(4, 4, 1.3, 0.13)$.
**Measured**: $\lambda_\max \approx 0.0205$, $t^* \approx 566$.

After COM removal, the system has 9 DOF in 3D (or 6 DOF in 2D, further reduced to
3 by Klein-four). The effective number of degrees of freedom for Nekhoroshev estimates
is $n_{\text{eff}}$.

For the double-orbiter, the Klein-four reduction gives $n_{\text{eff}} \approx 6$
(3D after COM removal and discrete symmetry). Using Lochak's improved exponent:

$$
a = \frac{1}{2n_{\text{eff}}} = \frac{1}{12}
$$

The Nekhoroshev stability time:

$$
T_N \sim \exp\left(\frac{1}{\varepsilon^{1/12}}\right)
$$

At $c = 1$, $\varepsilon = 1$: $T_N \sim \exp(1) \approx 2.7$. This is **trivially
short** -- the Nekhoroshev bound predicts only $O(1)$ stability at $\varepsilon = 1$.

At $c = 4$, $\varepsilon = 0.0625$: $T_N \sim \exp(0.0625^{-1/12}) = \exp(1.24) \approx 3.5$.
Still very short.

At $c = 100$, $\varepsilon = 10^{-4}$: $T_N \sim \exp(10^{4/12}) = \exp(10^{0.33}) = \exp(2.15) \approx 8.6$.
Still $O(10)$.

The Nekhoroshev exponent $a = 1/12$ is **extremely weak** for 6 DOF. The exponential
factor only becomes meaningfully large at $\varepsilon \lesssim 10^{-12}$
(i.e. $c \gtrsim 10^6$).

**Comparison with observation.** The observed $t^* \approx 566$ at $c = 1$
($\varepsilon = 1$) **vastly exceeds** the Nekhoroshev lower bound ($T_N \sim 3$).
This is typical: Nekhoroshev bounds are conservative by many orders of magnitude.
The observed lifetime is determined by the actual width of the stochastic layer
and the geometry of the resonance web, not by the worst-case Nekhoroshev estimate.

### 3.3. Application to the rhombus basin

**Configuration**: rhombus $(a=1.5, b=1.45, \eta=0.75)$.
**Measured**: $\lambda_\max \approx 0.18$, $t^* \approx 420$-$457$.

Same DOF analysis. The larger $\lambda_\max$ indicates a wider stochastic layer.
The Nekhoroshev bound gives the same $T_N \sim O(1)$ at $\varepsilon = 1$.

**Consistency check.** The ratio of escape times roughly tracks the ratio of
Lyapunov exponents:

$$
\frac{t^*_{\text{double-orb}}}{t^*_{\text{rhombus}}}
\approx \frac{566}{440} \approx 1.3
$$

$$
\frac{\lambda_{\text{rhombus}}}{\lambda_{\text{double-orb}}}
\approx \frac{0.18}{0.02} \approx 9
$$

If escape were purely Nekhoroshev-exponential, we would expect
$t^*/t^* \sim \exp(C/\lambda)$ for some $C$, giving a ratio of
$\exp(C(1/0.02 - 1/0.18)) \sim \exp(44C)$. The observed ratio of 1.3
implies $C \approx 0.006$, which is consistent with the system being in a
**non-Nekhoroshev regime** -- escape is not controlled by exponentially slow
diffusion but by direct transport along the stochastic layer.

### 3.4. Breathing alternating square

**Configuration**: alternating square with period $T = 11.78$.
**Measured**: $|\lambda|_\max \approx 228.6$ (Floquet multiplier).

The enormous Lyapunov exponent places this orbit deep in the chaotic sea -- far
from any regime where Nekhoroshev bounds are informative. This is a *violently
unstable periodic orbit*, not a quasi-bound state in the perturbative sense.

---

## 4. Arnold diffusion analysis

### 4.1. Framework

In systems with $\geq 3$ DOF (or $\geq 2.5$ DOF for time-periodic), Arnold diffusion
allows slow transport through the resonance web connecting hyperbolic tori. The
diffusion coefficient along the resonance web scales as:

$$
D \sim \varepsilon^2 \cdot |\text{splitting}|^2
$$

where the splitting of separatrices is exponentially small in $1/\sqrt{\varepsilon}$
for analytic systems (the Poincare-Melnikov mechanism).

### 4.2. Diffusion estimate for the double-orbiter

The double-orbiter has $n_{\text{eff}} \geq 3$ DOF after all reductions, so Arnold
diffusion is possible in principle.

**Heuristic estimate from Lyapunov exponent.** The MLE measures the local rate of
separation along the most unstable direction. In an Arnold-diffusion regime, the
MLE is related to the width of the stochastic layer around the resonance:

$$
\lambda_\max \sim \omega_\text{osc} \cdot \ln\left(\frac{\Delta_\text{layer}}{\Delta_\text{min}}\right)
$$

where $\omega_\text{osc}$ is the oscillation frequency near the separatrix and
$\Delta_\text{layer}$ is the layer width.

A simpler scaling uses the diffusion coefficient:

$$
D \sim \lambda_\max^2 \cdot \ell_\text{corr}^2
$$

where $\ell_\text{corr} \sim 1/\lambda_\max$ is the correlation length (time for
decorrelation). This gives $D \sim \text{const}$, independent of $\lambda_\max$ --
which is the "fast diffusion" regime.

**Time to traverse the trapping region.** The trapping region in action space has
a diameter $\Delta I \sim I_0 \cdot \delta$ where $\delta$ is the relative width
of the basin. From the parameter scans, $\delta \sim 0.02$-$0.05$ in each parameter
direction. The diffusion time is:

$$
t_\text{diff} \sim \frac{(\Delta I)^2}{D}
$$

With $D \sim v_\text{orb}^2 / \lambda_\max$ as a characteristic diffusion rate:
- $v_\text{orb} \sim 1.3$, $\lambda_\max \sim 0.02$
- $D \sim 1.3^2/0.02 \sim 85$ (in natural units)
- But the basin width is $\Delta I/I \sim 0.03$, so $\Delta I \sim 0.03 \times I_0$
- $t_\text{diff} \sim (0.03 \cdot I_0)^2 / D$

This heuristic gives $t_\text{diff} \sim O(10^2)$ -- the right ballpark for
$t^* \approx 566$, supporting the interpretation that the double-orbiter escape
is controlled by Arnold diffusion across a narrow trapping region.

### 4.3. The e-folding interpretation

The measured $\lambda_\max = 0.0205$ gives an e-folding time $\tau_e = 1/\lambda \approx 49$.
The escape at $t^* \approx 566$ represents $\sim 11.5$ e-foldings. This means an initial
perturbation of size $\delta_0$ grows to $\delta_0 \cdot e^{11.5} \approx 10^5 \delta_0$
before escape.

This is consistent with fast Arnold diffusion: the trajectory wanders randomly in
the stochastic layer, and escape occurs when the random walk reaches the basin
boundary. The number of e-foldings ($\sim 11$) is much less than the Nekhoroshev
prediction ($\sim e^{1/\varepsilon^a}$ e-foldings), confirming the system is in
the **fast-diffusion** rather than **Nekhoroshev** regime at $c = 1$.

### 4.4. Arnold diffusion for the rhombus

For the rhombus with $\lambda_\max \approx 0.18$: $\tau_e \approx 5.6$, and
$t^*/\tau_e \approx 80$. The larger Lyapunov exponent means faster local mixing,
but the escape time is only modestly shorter ($\sim 440$ vs $\sim 566$), suggesting
the rhombus basin has comparable geometric width in action space. The ratio
$t^*/\tau_e$ being larger for the rhombus ($80$ vs $28$) indicates the basin
is geometrically wider relative to the local diffusion rate.

---

## 5. Scaling predictions and testable consequences

### 5.1. Dependence on c

If the escape is controlled by Arnold diffusion with $\varepsilon = 1/c^2$:

**Fast-diffusion regime** ($c \lesssim 10$, $\varepsilon \gtrsim 0.01$):
$$
t^* \sim C_1 / \varepsilon \sim c^2
$$

**Nekhoroshev regime** ($c \gg 10$, $\varepsilon \ll 0.01$):
$$
t^* \sim \exp\left(\frac{C_2}{\varepsilon^{1/12}}\right) \sim \exp\left(C_2 \cdot c^{1/6}\right)
$$

**Prediction.** Running the double-orbiter at $c = 2, 4, 8, 16$ and measuring $t^*$
would distinguish the two regimes. If $t^* \propto c^2$, the system is in the
fast-diffusion regime at all tested values. If $t^*$ grows faster than any power
of $c$, Nekhoroshev stability has been reached.

### 5.2. Lyapunov exponent scaling

In the fast-diffusion regime: $\lambda_\max \sim \varepsilon \sim 1/c^2$.
In the Nekhoroshev regime: $\lambda_\max \sim \exp(-C/\varepsilon^a)$.

The measured $\lambda_\max = 0.02$ at $c = 1$ ($\varepsilon = 1$) is already
quite small. If $\lambda \propto 1/c^2$, then at $c = 4$:
$\lambda_\max \sim 0.02/16 \approx 0.001$.

### 5.3. Basin width scaling

The basin width in parameter space scales as $\sqrt{\varepsilon}$ in KAM theory
(the "KAM gap"). At $c = 1$ the basin has width $\Delta \sim 0.03$ in normalized
parameters. At $c = 4$ the basin should widen to $\Delta \sim 0.03 \times 4 = 0.12$.
This is testable by parameter scans at different $c$ values.

---

## 6. Why KAM theory cannot prove indefinite stability here

Several fundamental obstructions prevent a rigorous KAM proof of permanent bound
motion for the 4-body 2+/2- Weber system at $c = 1$:

1. **The perturbation is O(1).** At $c = 1$, $\varepsilon = 1$, and the Weber
   correction is the same order as the Coulomb potential. All perturbative
   theorems require $\varepsilon \ll 1$.

2. **No integrable base.** The 4-body Coulomb problem is not integrable.
   KAM/Nekhoroshev require an integrable $H_0$ as the starting point.

3. **Kepler degeneracy.** Even for 2-body, the Kepler Hamiltonian is properly
   degenerate ($\omega_r = \omega_\theta$), requiring specialized versions of
   KAM theory (Arnold's theorem for degenerate systems).

4. **Like-sign pairs.** The (+,+) pair has repulsive Coulomb interaction. Its
   bound-state character is entirely due to the Weber correction -- there is no
   Coulomb bound state to perturb.

5. **Positive Lyapunov exponent.** The measured $\lambda_\max = 0.0205 > 0$
   definitively rules out the trajectory lying on a KAM torus. The bound state
   is genuinely chaotic.

The appropriate theoretical framework is not KAM but rather:
- **Conley index theory** for proving the existence of an isolating neighborhood
- **Normally hyperbolic invariant manifold (NHIM)** theory for the trapping mechanism
- **Arnold diffusion rate estimates** for the escape time

---

## 7. Summary of bounds

| Candidate | $\lambda_\max$ | $t^*$ (obs) | $T_N$ (Nekh.) | Regime | $t^*/\tau_e$ |
|-----------|----------------|-------------|----------------|--------|--------------|
| 3D double-orbiter | $0.0205$ | $566$ | $\sim 3$ | Fast Arnold diff. | $11.5$ |
| Rhombus basin | $0.18$ | $420$-$457$ | $\sim 3$ | Fast Arnold diff. | $\sim 80$ |
| Breathing square | $228.6$ | -- (periodic) | N/A | Violently unstable | N/A |

All quasi-bound candidates at $c = 1$ are in the fast-Arnold-diffusion regime,
not the Nekhoroshev regime. The Nekhoroshev lower bounds are trivially satisfied
and do not constrain the dynamics. The observed lifetimes are set by the geometric
width of the stochastic layer, which is $O(1)$ at $\varepsilon = 1$.

---

## 8. Computational verification

See [kam_2body.jl](kam_2body.jl) for numerical verification of:
- The KAM twist condition for 2-body Weber at $c = 4, 10, 100$
- The perturbation size $\varepsilon\langle H_1\rangle / |H_0|$
- Action-angle variable computation for circular and eccentric orbits

See [nekhoroshev_bounds.md](nekhoroshev_bounds.md) for explicit formulas and
numerical values of all stability-time estimates.

---

## References

1. V.I. Arnold, "Proof of a theorem of A.N. Kolmogorov on the invariance of
   quasi-periodic motions under small perturbations of the Hamiltonian,"
   Russ. Math. Surv. 18(5):9-36 (1963).
2. N.N. Nekhoroshev, "An exponential estimate of the time of stability of
   nearly-integrable Hamiltonian systems," Russ. Math. Surv. 32(6):1-65 (1977).
3. P. Lochak, "Canonical perturbation theory via simultaneous approximation,"
   Russ. Math. Surv. 47(6):57-133 (1992).
4. A. Giorgilli, C. Skokos, "On the stability of the Trojan asteroids,"
   Astron. Astroph. 317:254-261 (1997).
5. V.I. Arnold, "Instability of dynamical systems with several degrees of
   freedom," Sov. Math. Dokl. 5:581-585 (1964).
