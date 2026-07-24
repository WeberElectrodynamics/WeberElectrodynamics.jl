# Solver

## One-shot interface

```@docs
solve
```

## Step-by-step interface

```@docs
init
step!
solve!
```

## Types

```@docs
HamiltonianIntegrator
HamiltonianSolution
SymmetricProjectionIntegrator
RegularizedIntegrator
```

## Archives

JLD2 archive helpers are provided by an optional extension. They currently
archive default Weber systems by storing primitive solution arrays and
problem metadata, then reconstructing the system on load. Load `JLD2` before
calling them:

```julia
using WeberElectrodynamics, JLD2
save_solution("run.jld2", sol; metadata = (case = "two-body",))
sol2 = load_solution("run.jld2")
```

Archive format v2 is Weber-only. Archives created by the 0.5.x series use
format v1 and must be exported or regenerated before upgrading.

```@docs
save_solution
load_solution
```
