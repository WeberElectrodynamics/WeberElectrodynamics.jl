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
archive default Weber/Zöllner systems by storing primitive solution arrays and
problem metadata, then reconstructing the system on load. Load `JLD2` before
calling them:

```julia
using WeberElectrodynamics, JLD2
save_solution("run.jld2", sol; metadata = (case = "two-body",))
sol2 = load_solution("run.jld2")
```

```@docs
save_solution
load_solution
```
