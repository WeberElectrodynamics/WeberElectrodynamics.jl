# Weber Electrodynamics: Two-Body Problem Formulas

Reference: Clemente & Assis (1991)

## Core Definitions

**Reduced mass:**
```
μ = m₁m₂/(m₁ + m₂)
```

**Interaction constant:**
```
U₀ = q₁q₂  (for electromagnetic interaction)
K = U₀/(μc²)
```

**Position variables:**
```
r = |r₁ - r₂|
r̂ = (r₁ - r₂)/r
```

## Weber Force Law

```
F₁,₂ = -U₀(r̂/r²)(1 + rṙ/c² - ṙ²/2c²)
```

## Conservation Laws

**Angular momentum:**
```
L = μr²θ̇
```

**Energy:**
```
W = (μ/2)(ṙ² + r²θ̇²) + (U₀/r)(1 - ṙ²/2c²)
```

## Trajectory Equations

**Variable transformation:**
```
x² = 1 - K/r
```

**Turning points:**
```
x₁,₂² = 1 + (μKU₀/L²)[1 ± (1 + 2WL²/(μU₀²))^(1/2)]
```

**Differential equation:**
```
dx/dθ = ±(1/2x²)[(x₁² - x²)(x² - x₂²)]^(1/2)
```

## Trajectory Solutions

**Attractive force (U₀ < 0):**
```
θᴬ = ±2|x₁|E(φ, k)
```

**Repulsive force (U₀ > 0):**
```
θᴿ = ±2|x₁|[E(k) - E(φ, k)]
```

**Where:**
```
φ = arcsin[(x₁² - x²)/(x₁² - x₂²)]^(1/2)
k = [(x₁² - x₂²)/x₁²]^(1/2)
```

E(φ, k) = incomplete elliptic integral of the second kind
E(k) = complete elliptic integral of the second kind

## Perihelion Precession (Limited Trajectories)

**Angle per cycle:**
```
Δθ = 4|x₁|E(k)  (always > 2π)
```

**Shift (small |K| approximation):**
```
δθ ≈ π|K|/[a(1 - e²)]
```

where a = semimajor axis, e = eccentricity
