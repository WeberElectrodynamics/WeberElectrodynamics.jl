# Zöllner Extension

> **Research / Experimental.** The Zöllner extension implements a 19th-century
> electrogravitational hypothesis. It is not part of standard Weber
> electrodynamics and is intended for comparative research, not general
> simulation. It is **off by default**.

## Background

Johann Karl Friedrich Zöllner (1837–1882) proposed that matter consists of
binary electrical dyads with a slight *mismatch*: the unlike-sign attraction is
fractionally stronger than the like-sign repulsion. The residual attraction
between neutral dyads produces an effective gravitational interaction.

The mismatch is parameterized by `a > 0`:

- Unlike-sign pairs: κᵢⱼ = 1 + a  (stronger attraction)
- Like-sign pairs:   κᵢⱼ = 1.0    (standard repulsion)

Each pair potential becomes:

```
U_ij = κ_ij · qᵢqⱼ / r · (1 - ṙ²/(2c²))
```

When `a = 0` (or Zöllner disabled) all κ = 1 and the standard Weber
potential is recovered exactly.

See [Zöllner Electrogravitational Theory](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/ZollnerElectrogravitationalTheory.md)
for the full theoretical derivation.

## Enabling Zöllner

```julia
prob = HamiltonianProblem(
    sys, tspan, q0, p0;
    masses  = [1.0, 1.0],
    charges = [1.0, -1.0],
    c       = 10.0,
    dt      = 0.01,
    zollner = ZollnerOptions(enabled = true, a = 0.01),   # mismatch parameter; must be > 0
)
```

`a` must be strictly positive when `enabled = true`.

## Combining with regularization

Zöllner is fully compatible with regularization:

```julia
prob = HamiltonianProblem(sys, tspan, q0, p0; ...
    zollner        = ZollnerOptions(enabled = true, a = 0.01),
    regularization = RegularizationOptions(enabled = true),
)
```

κ values are automatically included in the parameter vector passed to regularization
sub-steps, so the Zöllner coupling is respected even during close-encounter
regularization. No extra configuration is needed.

## Zöllner-specific plots

The Plots.jl extension provides four functions for visualizing the Zöllner
contribution:

```julia
using Plots

energy = compute_energy_timeseries(sol)
forces = compute_pair_force_timeseries(sol, 1, 2)

plot_zollner_energy(energy)               # pair potentials + Zöllner extra
plot_zollner_force_residual(forces)       # force magnitude vs Zöllner fraction
plot_weber_vs_zollner(sol_weber, sol_z)   # overlay two solutions
plot_zollner_phase_space(f1, f2)          # (r, ṙ) portrait comparison
```

## API reference

```@docs
ZollnerOptions
```
