# Hypergeometric Structure of the N-Body Weber Hamiltonian

**Reference**: A=B — Petkovšek, Wilf & Zeilberger (1996)  
**Tools**: Python / SymPy · Gosper's algorithm · Zeilberger / creative telescoping  
**Scripts**: `01_ratio_tests.py` · `02_closed_forms.py` · `03_recurrence.py` · `04_weber_analysis.py` · `05_ring_and_summary.py`

---

## The Equations Under Study

Hamilton's equations for the N-body Weber Hamiltonian (paper Eqs. 211, 217):

$$\dot x_i = \frac{1}{m_i}\!\left(p_{x_i} - \sum_{j\neq i}\frac{q_i q_j}{c^2}\,\frac{\dot r_{ij}(x_i-x_j)}{r_{ij}^2}\right)$$

$$\dot p_{x_i} = \sum_{j\neq i}\frac{q_i q_j}{r_{ij}^2}\!\left[\frac{x_i-x_j}{r_{ij}}\!\left(1-\frac{3\dot r_{ij}^2}{2c^2}\right) + \frac{\dot r_{ij}(\dot x_i-\dot x_j)}{c^2}\right]$$

**Three questions answered below:**
1. Are the summation expressions hypergeometric in the A=B sense?
2. Can Gosper/Zeilberger find closed-form N-dependence?
3. What do those closed forms look like?

---

## Part 1 — Hypergeometric Ratio Test (A=B §3.2)

A sequence $f(k)$ is **hypergeometric** iff $f(k+1)/f(k)$ is rational in $k$.

**Symmetric configuration**: N identical particles (mass $m$, charge $q$) in a 1D chain, spacing $d$, pair distance $r_{ij} = |i-j| \cdot d$, reindexed by $k = |i-j|$.

### Results

| Term | $f(k)$ | $f(k+1)/f(k)$ | Hypergeometric? |
|------|--------|---------------|-----------------|
| Coulomb energy | $(N-k)/k$ | $\dfrac{k(N-k-1)}{(k+1)(N-k)}$ | **YES** ✓ |
| Force ($\dot p$, static) | $1/k^2$ | $k^2/(k+1)^2$ | **YES** ✓ |
| Weber correction | $(N-k)/k \cdot (1-v^2/2c^2)$ | same as Coulomb | **YES** ✓ |
| Ring csc | $\csc(\pi k/N)$ | $\sin(\pi k/N)/\sin(\pi(k+1)/N)$ | **NO** ✗ |

The 1D-chain pair sums in Hamilton's equations **are hypergeometric** in $k$.  
The ring configuration is **not** — Gosper/Zeilberger do not apply there.

---

## Part 2 — Closed Forms via Gosper / Zeilberger (A=B §5–6)

### Important subtlety: Gosper on diagonal sums

When $N$ appears in **both** the summand and the upper limit ("diagonal sum"),
`gosper_sum` from SymPy gives incorrect results:

```
gosper_sum((N−k)/k, k=1..N−1) = −(N²−3N+1)/(N−1)

Verification:
  N=2: direct=1.0000,  gosper formula= 1.0000  ✓
  N=3: direct=2.5000,  gosper formula=−0.5000  ✗ WRONG
  N=4: direct=4.3333,  gosper formula=−1.6667  ✗ WRONG
```

`summation()` handles this correctly via its Zeilberger-based backend.

---

### Result 1: Total Coulomb energy

$$U_C(N) = \frac{q^2}{d} \sum_{k=1}^{N-1} \frac{N-k}{k} = \frac{q^2}{d}\bigl[N\,H_{N-1} - (N-1)\bigr]$$

where $H_n = \sum_{j=1}^n 1/j$ is the $n$-th harmonic number.

**Validation** (`summation()` result = direct pair sum):

| N  | Direct sum | Closed form $N H_{N-1}-(N-1)$ | Match |
|----|-----------|-------------------------------|-------|
| 2  | 1.000000  | 1.000000                      | ✓ |
| 3  | 2.500000  | 2.500000                      | ✓ |
| 4  | 4.333333  | 4.333333                      | ✓ |
| 5  | 6.416667  | 6.416667                      | ✓ |
| 10 | 19.289683 | 19.289683                     | ✓ |
| 20 | 51.954793 | 51.954793                     | ✓ |

Growth: $U_C(N)/N \to \ln N + \gamma$ as $N\to\infty$ (Euler–Mascheroni $\gamma \approx 0.5772$).

---

### Result 2: Weber-corrected energy

For a symmetric configuration with uniform relative radial speed $v$ (same for all pairs):

$$U_W(N) = \frac{q^2}{d}\left(1 - \frac{v^2}{2c^2}\right)\bigl[N\,H_{N-1} - (N-1)\bigr]$$

`summation()` result:
```
(2c² − v²)(N·harmonic(N−1) − N + 1) / (2c²)
```

**Key finding**: $S_W(N) / S_C(N) = 1 - v^2/(2c^2)$

The Weber velocity correction is a **constant scalar** — it does not alter the $N$-dependence at all.

Weber correction fraction $\Delta U / U_C = v^2/2c^2$:

| $v/c$ | $\Delta U/U_C$ | $U_W/U_C$ |
|--------|----------------|-----------|
| 0.001  | 0.000001       | 0.999999  |
| 0.010  | 0.000050       | 0.999950  |
| 0.100  | 0.005000       | 0.995000  |
| 0.200  | 0.020000       | 0.980000  |
| 0.300  | 0.045000       | 0.955000  |
| 0.500  | 0.125000       | 0.875000  |

---

### Result 3: Force on end particle

Gosper's algorithm **returns None** for $1/k^2$ — no rational anti-difference exists.  
`summation()` recognises it as a generalised harmonic sum:

$$F_\text{end}(N) = \frac{q^2}{d^2}\sum_{k=1}^{N-1}\frac{1}{k^2} = \frac{q^2}{d^2}\,H_{N-1}^{(2)}$$

where $H_n^{(2)} = \sum_{k=1}^n 1/k^2$ is the 2nd-order generalised harmonic number.

**Basel limit**: $F_\text{end} \to (q^2/d^2) \cdot \pi^2/6 \approx 1.6449\, q^2/d^2$ as $N\to\infty$.

**Validation**:

| N  | Direct sum  | $H_{N-1}^{(2)}$  | Match |
|----|-------------|-----------------|-------|
| 2  | 1.00000000  | 1.00000000      | ✓ |
| 3  | 1.25000000  | 1.25000000      | ✓ |
| 5  | 1.42361111  | 1.42361111      | ✓ |
| 10 | 1.53976773  | 1.53976773      | ✓ |
| 20 | 1.59366324  | 1.59366324      | ✓ |

---

### Result 4: Force on interior particle $i$

$$F_i(N) = \frac{q^2}{d^2}\Bigl[H_i^{(2)} - H_{N-1-i}^{(2)}\Bigr]$$

(positive = rightward; left neighbours push right, right neighbours push left)

Special cases:
- $i = 0$ (leftmost): $F = -H_{N-1}^{(2)}$ (pushed left)
- $i = (N-1)/2$ (middle): $F = 0$ (symmetry)
- $i = N-1$ (rightmost): $F = +H_{N-1}^{(2)}$ (pushed right)

**Validation for $N=5$**:

| $i$ | Direct force | $H_i^{(2)} - H_{N-1-i}^{(2)}$ | Match |
|-----|-------------|-------------------------------|-------|
| 0   | −1.423611   | −1.423611                     | ✓ |
| 1   | −0.361111   | −0.361111                     | ✓ |
| 2   |  0.000000   |  0.000000                     | ✓ |
| 3   | +0.361111   | +0.361111                     | ✓ |
| 4   | +1.423611   | +1.423611                     | ✓ |

---

## Part 3 — Zeilberger Recurrence (A=B §6)

**Theorem**: $E(N) := N H_{N-1} - (N-1)$ satisfies the recurrence

$$E(N+1) - E(N) = H_N, \qquad E(2) = 1$$

**Proof** (creative telescoping):

$$E(N+1) - E(N) = (N+1)H_N - N - N H_{N-1} + (N-1)$$
$$= H_{N-1} + \frac{N+1}{N} - 1 \quad\bigl[\text{using }H_N = H_{N-1} + 1/N\bigr]$$
$$= H_{N-1} + \frac{1}{N} = H_N \qquad \square$$

This is the **Zeilberger certificate**: a recurrence that proves the closed form without evaluating the sum from scratch.

**Numerical verification** (all 11 cases exact to 12 decimal places):

| $N$ | $E(N+1)$    | $E(N)$      | $E(N+1)-E(N)$ | $H_N$       | OK? |
|-----|------------|------------|--------------|------------|-----|
| 2   | 2.50000000 | 1.00000000 | 1.50000000   | 1.50000000 | ✓ |
| 3   | 4.33333333 | 2.50000000 | 1.83333333   | 1.83333333 | ✓ |
| 4   | 6.41666667 | 4.33333333 | 2.08333333   | 2.08333333 | ✓ |
| 5   | 8.70000000 | 6.41666667 | 2.28333333   | 2.28333333 | ✓ |
| 6   | 11.1500000 | 8.70000000 | 2.45000000   | 2.45000000 | ✓ |
| 7   | 13.7428571 | 11.1500000 | 2.59285714   | 2.59285714 | ✓ |
| 8   | 16.4607143 | 13.7428571 | 2.71785714   | 2.71785714 | ✓ |
| 10  | 22.2186508 | 19.2896825 | 2.92896825   | 2.92896825 | ✓ |
| 12  | 28.3417388 | 25.2385281 | 3.10321068   | 3.10321068 | ✓ |

**Physical meaning**: each time a particle is added to the chain, the total Coulomb energy increases by exactly $H_N$ (in units of $q^2/d$).

### Growth table ($N = 2\ldots20$, $q=d=1$)

| $N$ | $U_C(N)$   | $U_C/N$  | $\Delta E = H_N$ |
|-----|-----------|---------|-----------------|
| 2   | 1.000000  | 0.5000  | 1.500000 |
| 3   | 2.500000  | 0.8333  | 1.833333 |
| 4   | 4.333333  | 1.0833  | 2.083333 |
| 5   | 6.416667  | 1.2833  | 2.283333 |
| 6   | 8.700000  | 1.4500  | 2.450000 |
| 8   | 13.742857 | 1.7179  | 2.717857 |
| 10  | 19.289683 | 1.9290  | 2.928968 |
| 15  | 34.773435 | 2.3182  | 3.318229 |
| 20  | 51.954793 | 2.5977  | 3.597740 |

---

## Part 4 — Weber Velocity Analysis

### Case A: Collinear uniform motion

All particles move at the same velocity $v_0$ along the chain axis.  
Relative radial velocity $\dot r_{ij} = 0$ for all pairs.  
**Weber factor = 1** → $U_W = U_C$ (pure Coulomb).

Position equation (Eq. 211): $\dot x_i = p_{x_i}/m$ (velocity-coupling sum vanishes).  
Momentum equation (Eq. 217): reduces to static Coulomb force.

### Case B: Head-on approach at speed $v$

Particle 0 moves at $+v$, particle $N-1$ at $-v$, rest stationary.  
For particle 0 vs particle $k$: $\dot r_{0k} = v$.

Force on particle 0 from the Weber momentum equation:

$$\dot p_{x_0} = \frac{q^2}{d^2}\left(1 - \frac{3v^2}{2c^2}\right) H_{N-1}^{(2)}$$

Note the **different** Weber factor $(1 - 3v^2/2c^2)$ compared to the energy term $(1 - v^2/2c^2)$.

Numerical values ($v/c = 0.1$, $q=d=1$):

| $N$ | $\dot p_{x_0}$ |
|-----|---------------|
| 3   | 1.231250 |
| 5   | 1.402257 |
| 10  | 1.516671 |

### Case C: Uniform relative radial speed (symmetric)

`summation()` result:

```
S_W(N) = (2c² − v²)(N·H_{N-1} − N + 1) / (2c²)
S_W(N) / S_C(N) = 1 − v²/(2c²)
```

The Weber correction is a **constant independent of $N$** — the harmonic-number structure is fully preserved.

---

## Part 5 — Ring Configuration (Non-Hypergeometric)

For $N$ charges on a circle of radius $R$, pair distance $r_k = 2R\sin(\pi k/N)$:

$$U_\text{ring}(N) = \frac{Nq^2}{4R}\sum_{k=1}^{N-1}\csc\!\left(\frac{\pi k}{N}\right)$$

The csc ratio is **not rational** in $k$: Gosper/Zeilberger **do not apply**.

Numerical values of $T(N) = \sum_{k=1}^{N-1}\csc(\pi k/N)$:

| $N$  | $T(N)$      | $T(N) \div (N\ln N/\pi)$ |
|------|------------|--------------------------|
| 2    | 1.000000   | 2.2662 |
| 4    | 3.828427   | 2.1690 |
| 6    | 7.309401   | 2.1360 |
| 10   | 15.449800  | 2.1079 |
| 20   | 39.738094  | 2.0836 |
| 50   | 128.520836 | 2.0642 |

**Large-$N$ asymptotic**: $T(N) \sim (N/\pi)\ln N$ (ratio converges to ≈ 2.07 as $N\to\infty$; the true limit involves the Euler–Mascheroni constant).

**Exact special cases**:
- $N=2$: $T = 1$
- $N=4$: $T = 1 + 2\sqrt{2} \approx 3.8284$
- $N=6$: $T = 5 + 4/\sqrt{3} \approx 7.3094$

No general elementary closed form exists for arbitrary $N$.

---

## Complete Summary of Closed Forms

> **Configuration**: $N$ identical particles, mass $m$, charge $q$, spacing $d$ (1D chain).  
> All results are **exact** and verified against direct pair sums for $N = 2\ldots20$.

| # | Quantity | Closed Form |
|---|----------|-------------|
| 1 | Total Coulomb energy | $U_C(N) = \dfrac{q^2}{d}\bigl[N H_{N-1} - (N-1)\bigr]$ |
| 2 | Weber-corrected energy | $U_W(N) = \Bigl(1 - \dfrac{v^2}{2c^2}\Bigr) U_C(N)$ |
| 3 | Force on end particle (magnitude) | $\lvert F_\text{end}(N)\rvert = \dfrac{q^2}{d^2} H_{N-1}^{(2)} \;\to\; \dfrac{q^2\pi^2}{6d^2}$ |
| 4 | Force on particle $i$ (signed) | $F_i(N) = \dfrac{q^2}{d^2}\bigl[H_i^{(2)} - H_{N-1-i}^{(2)}\bigr]$ |
| 5 | Energy recurrence | $E(N+1) = E(N) + H_N,\quad E(2) = 1$ |
| 6 | Ring energy | $\propto T(N) \sim (N/\pi)\ln N$ (no closed form) |

**Notation**:
- $H_n = \text{harmonic}(n) = \sum_{k=1}^n 1/k$ (harmonic number)
- $H_n^{(2)} = \text{harmonic}(n,2) = \sum_{k=1}^n 1/k^2$ (generalised harmonic number)
- $\gamma \approx 0.5772$ (Euler–Mascheroni constant)
- $H_n \sim \ln n + \gamma$ for large $n$

---

## Scope and Limitations

These are **exact N-scaling formulas for the potential energy and forces at fixed symmetric configurations** — not solutions to the full nonlinear dynamics.

- The general N-body problem has no closed-form solution (Poincaré–Bruns theorem).
- A=B techniques give the exact $N$-dependence of the Hamiltonian at a snapshot in time.
- The Weber velocity correction enters only as a scalar prefactor; the harmonic-number structure of the $N$-body sums is preserved in all cases examined.
- The ring configuration falls outside the hypergeometric framework — only asymptotic results are available there.

---

## Algorithm Notes

| Algorithm | When it applies | SymPy call |
|-----------|----------------|------------|
| **Gosper** (indefinite sum) | $f(k+1)/f(k)$ rational in $k$, $N$ not in upper limit | `gosper_sum(f, k)` |
| **Zeilberger** (definite sum) | $f(k+1)/f(k)$ rational in $k$, $N$ may appear in limit | `summation(f, (k, a, b))` |
| **Neither** | $f(k+1)/f(k)$ not rational (e.g. csc term) | numerical only |

`gosper_sum` **fails silently** on diagonal sums (N in both summand and limit) — always validate against direct computation. `summation()` is the reliable choice.
