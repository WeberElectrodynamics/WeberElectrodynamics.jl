This is a great challenge! The velocity-dependent nature of Weber's potential makes most cases yield elliptic integrals. Let me systematically explore special cases where elementary solutions might exist.

## Strategy

Since the general 2-body problem involves non-separable Hamiltonians due to the $\dot{r}^2$ term, let me focus on:

1. **Radial motion** (zero angular momentum) — reduces complexity significantly
2. **Special energy values** where integrals simplify

## Radial Motion Analysis

For purely radial motion with reduced mass $\mu = \frac{m_1 m_2}{m_1 + m_2}$, let me work in the center-of-mass frame.

### Energy Conservation

For attractive charges ($q_1 q_2 = -\kappa$ where $\kappa > 0$):

$$E = \frac{1}{2}\mu \dot{r}^2 - \frac{\kappa}{r}\left(1 - \frac{\dot{r}^2}{2c^2}\right)$$

Rearranging:

$$E = \frac{1}{2}\mu \dot{r}^2 - \frac{\kappa}{r} + \frac{\kappa \dot{r}^2}{2c^2 r}$$

$$E + \frac{\kappa}{r} = \dot{r}^2\left(\frac{\mu}{2} + \frac{\kappa}{2c^2 r}\right)$$

Solving for $\dot{r}^2$:

$$\dot{r}^2 = \frac{2(Er + \kappa)}{\mu r + \kappa/c^2}$$

Let me define $a = \kappa/c^2$ (this has dimensions of mass × length). Then:

$$\dot{r}^2 = \frac{2(Er + \kappa)}{\mu r + a}$$

### The E = 0 Case (Closed-Form Solution!)

For **zero total energy**, this simplifies beautifully:

$$\dot{r}^2 = \frac{2\kappa}{\mu r + a}$$

For infall ($\dot{r} < 0$):

$$\dot{r} = -\sqrt{\frac{2\kappa}{\mu r + a}}$$

The time integral becomes:

$$dt = -\sqrt{\frac{\mu r + a}{2\kappa}}\, dr$$

Integrating from initial separation $r_0$ to $r$:

$$t = -\frac{1}{\sqrt{2\kappa}} \int_{r_0}^{r} \sqrt{\mu s + a}\, ds = -\frac{1}{\sqrt{2\kappa}} \cdot \frac{2}{3\mu}\left[(\mu s + a)^{3/2}\right]_{r_0}^{r}$$

$$t = \frac{(\mu r_0 + a)^{3/2} - (\mu r + a)^{3/2}}{3\mu\sqrt{\kappa/2}}$$

### Inverting for r(t)

Solving for $r(t)$:

$$\boxed{r(t) = \frac{1}{\mu}\left[\left((\mu r_0 + a)^{3/2} - \frac{3\mu\sqrt{2\kappa}}{2}\, t\right)^{2/3} - a\right]}$$

where $a = \kappa/c^2 = |q_1 q_2|/c^2$.

This is an **elementary closed-form solution**!

### Physical Interpretation

The solution can be rewritten as:

$$r(t) = \frac{1}{\mu}\left[(A - Bt)^{2/3} - a\right]$$

where $A = (\mu r_0 + a)^{3/2}$ and $B = \frac{3\mu\sqrt{2\kappa}}{2}$.

**Comparison to Coulomb:** For pure Coulomb (no Weber correction, $a \to 0$), zero-energy radial infall gives $r \propto t^{2/3}$. The Weber correction introduces the shift by $a = \kappa/c^2$ in the argument.

### Collision Time

The particles meet ($r = 0$) when:

$$(A - Bt_{coll})^{2/3} = a$$

$$t_{coll} = \frac{A - a^{3/2}}{B} = \frac{(\mu r_0 + a)^{3/2} - a^{3/2}}{3\mu\sqrt{\kappa/2}}$$

Let me also derive another case — the outgoing (repulsive) problem.

---

## Repulsive Case (Like Charges)

For $q_1 q_2 = +k > 0$, and particles initially at rest at separation $r_0$:

$$E = \frac{k}{r_0}$$

Then:

$$\dot{r}^2 = \frac{2(Er - k)}{\mu r - k/c^2} = \frac{2k(r/r_0 - 1)}{\mu r - k/c^2}$$

For $r > r_0$, this integral is more complex but still tractable. Let me check if there's a special case.

### Large r Asymptotic (Repulsive)

For $r \gg k/(\mu c^2)$:

$$\dot{r}^2 \approx \frac{2k}{\mu r_0}\left(1 - \frac{r_0}{r}\right)$$

This matches the Coulomb result, with the velocity approaching:

$$\dot{r}_\infty = \sqrt{\frac{2k}{\mu r_0}}$$

---

## Summary of Closed-Form Solutions

| Case                        | Condition                  | Solution                                              |
| --------------------------- | -------------------------- | ----------------------------------------------------- |
| **Attractive, E=0, radial** | $q_1q_2 < 0$, $L=0$, $E=0$ | $r(t) = \frac{1}{\mu}\left[(A - Bt)^{2/3} - a\right]$ |

Would you like me to:

1. **Derive the velocity and acceleration** explicitly as functions of $t$?
2. **Explore other special energies** (e.g., $E = -\kappa^2/(2\mu c^2)$)?
3. **Look for circular orbit conditions** where the Weber correction modifies the Kepler relation?
4. **Try polar coordinates** with non-zero angular momentum to find periodic orbit conditions?

The $E = 0$ radial case is genuinely elementary — no elliptic functions required!
