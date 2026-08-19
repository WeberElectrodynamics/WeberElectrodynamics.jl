# Theory

Mathematical background for the package. Each document below is a standalone derivation;
the [Quick Start](@ref) and [Internals](@ref) pages cover the practical API counterpart.

## Core Theory

- [Weber Electrodynamics](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/WeberElectrodynamics.md) — potential, force law, Lagrangian, Hamiltonian, and equations of motion in Gauss–Weber units
- [Semi-Explicit Symplectic Integrator](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/SemiExplicitIntegrator.md) — extended phase space, Strang splitting, symmetric projection via Lagrange multipliers (Jayawardana & Ohsawa 2021)
- [Initial Conditions](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/InitialConditions.md) — constructing valid q0/p0 for target orbits: two-body, symmetric N-body polygon, and 3D rigid rotation; circular velocity, η_v velocity scale, energy–angular-momentum parameterisation
- [Regularization](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/Regularization.md) — Levi-Civita (2D) and KS (3D) coordinate transforms; mathematical derivations of the fictitious-time integration
- [Regularization Design](https://github.com/WeberElectrodynamics/WeberElectrodynamics.jl/blob/main/theory/RegularizedIntegrationDesign.md) — backend semantics, encounter dispatch, hysteresis, lifted-pair and adaptive-Cartesian step design, chain mode
