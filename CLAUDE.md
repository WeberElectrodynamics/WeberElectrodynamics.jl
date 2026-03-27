# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Julia package (v0.2.0) for n-body Weber electrodynamics simulation with Zöllner electrogravitational extension. Implements a symplectic Strang-splitting symmetric-projection integrator with Levi-Civita/KS regularization for close encounters and collision bounce for head-on singularities.

## Commands

```bash
# Run full test suite (from project root)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file (e.g., test_physics.jl)
julia --project=. -e 'using Test; using WeberElectrodynamics; using WeberElectrodynamics: SymmetricProjectionIntegrator; using LinearAlgebra; using Symbolics; @testset "single" begin include("test/test_utils.jl"); include("test/test_physics.jl") end'

# Format all Julia files
julia -e 'using JuliaFormatter; format(".")'

# Release a new version (patch / minor / major)
./release.sh patch
```

## Repository Structure

```
WeberElectrodynamics/
├── src/                    # Package source
├── ext/                    # Weak-dependency extensions (Plots, Makie)
├── test/                   # Test suite (20,297 tests)
├── examples/               # Jupyter notebooks (run from default Julia env)
├── docs/                   # Documenter.jl scaffold (make.jl, src/, build/)
├── research/
│   ├── theory/             # Mathematical derivations and theory documents
│   ├── exploratory/        # Research notes and lessons learned
│   └── sub_critical_weber_research/  # Sub-critical Weber research
├── papers/
│   └── Computational-Weber-Electrodynamics/   # LaTeX paper with own Project.toml
├── CHANGELOG.md            # Versioned changelog (semver)
└── Project.toml            # Package metadata and compat bounds
```

## Architecture

**Pipeline**: Symbolic Hamiltonian → compiled equations of motion → symplectic integration → statistics/plotting/animation

### Source files (`src/`)

- `WeberElectrodynamics.jl` — Module definition, exports, extension stubs (`plot_*`, `animate_weber`)
- `weber_system.jl` — `WeberSystem`: uses Symbolics.jl to build the Weber Hamiltonian symbolically, then compiles `dq_dt`, `dp_dt`, and `hamiltonian` functions via `build_function`
- `types.jl` — All core structs: `WeberProblem`, `WeberSolution`, `WeberIntegrator`, `SymmetricProjectionIntegrator`, `RegularizationOptions`, `ZollnerOptions`, buffer/diagnostics types
- `regularization.jl` — Internal helpers: pair distance detection, adjacency graph (BFS), Levi-Civita 2D projection, KS quaternion helpers
- `solve.jl` — Main integrator: Strang splitting flow, symmetric projection via fixed-point iteration on Lagrange multipliers, regularization dispatch, collision bounce, CommonSolve interface (`init`/`step!`/`solve!`/`solve`)
- `statistics/` — `energy.jl`, `forces.jl`, `momentum.jl`, `trajectories.jl` — post-solution analysis producing typed data structs

### Extensions (`ext/`)

- `WeberElectrodynamicsPlotsExt.jl` — Plots.jl weak dependency; provides `plot_trajectories`, `plot_energy`, `plot_pair_energy`, `plot_energy_errors`, `plot_pair_forces`, `plot_phase_space`, `plot_momentum`, and Zöllner-specific plot functions.
- `WeberElectrodynamicsMakieExt.jl` — Makie weak dependency (any backend: GLMakie, CairoMakie, WGLMakie); provides `animate_weber` for real-time streaming or solution replay with rolling trajectory/energy/momentum/phase-space dashboard.

### Tests (`test/`)

- `test_utils.jl` — Problem builders (`make_weber_problem()`, `make_coulomb_like_problem()`) and reference energy functions; **must be included before other test files**
- `runtests.jl` — Entry point, includes all test files in order
- Test files: `test_types.jl`, `test_weber_system.jl`, `test_solve.jl`, `test_statistics.jl`, `test_integration.jl`, `test_physics.jl`, `test_regularization.jl`, `test_zollner.jl`

### Examples (`examples/`)

Jupyter notebooks run via IJulia from the default Julia environment (where this package is `dev`'d). They use `Plots` for static figures and optionally `GLMakie` for animation.

### Docs (`docs/`)

Documenter.jl scaffold: `make.jl`, `Project.toml`, `src/` (page sources), `build/` (generated output).

### Research (`research/`)

- `theory/` — Weber electrodynamics, semi-explicit integrator, regularization, critical radius, initial conditions, Zöllner theory
- `exploratory/` — Collision bounce lessons learned, three-body bound states
- `sub_critical_weber_research/` — Sub-critical Weber exploration and literature searches

## Critical Conventions

### Params vector layout

```
params = [m₁, ..., mₙ, q₁, ..., qₙ, c, κ₁₂, κ₁₃, ..., κ_{N-1,N}]
length = 2N + 1 + N*(N-1)/2
```

Any code calling `sys.dq_dt_compiled(out, q, p, params)` or `sys.dp_dt_compiled(out, q, p, params)` directly **must** include the κ (kappa) entries. When Zöllner is disabled, all κ values are 1.0.

Pair index: `_pair_index(i, j, n) = (i-1)*(2n-i)÷2 + (j-i)` (1-based, i < j)

### Regularization backends

Only two valid values for `regularization_backend`:
- `:adaptive_cartesian` — KS-style, works for 2D and 3D
- `:lifted_pair` — Levi-Civita, **2D only** (auto-falls back to `:adaptive_cartesian` for 3D)

Neither backend regularizes Weber's velocity-dependent force — only the Coulomb/Kepler singularity.

### Collision bounce

- Enabled via `regularization_collision_bounce_radius` kwarg on `WeberProblem` (default 0.0 = off)
- Only valid for ℓ=0 (head-on) collisions
- Works best with the **unregularized** integrator (symplectic error stays bounded)

### Zöllner extension

- `ZollnerOptions(enabled, a)` — mismatch parameter `a`
- κ_ij = 1+a for unlike-sign charge pairs, 1.0 for like-sign
- Stored in `WeberProblem.kappas` and appended to the params vector automatically

### Makie animation extension

- Weak dependency is `Makie` (not `GLMakie`) — any backend triggers the extension
- `animate_weber(prob)` for live streaming, `animate_weber(sol)` for replay
- Compat: `Makie = "0.21, 0.22, 0.23, 0.24"`

### Immutable options pattern

`RegularizationOptions`, `ZollnerOptions` are immutable structs created once per problem. Pass configuration through `WeberProblem` keyword arguments rather than mutating options.

## EnergyStatistics fields

`en.statistics` has: `local_error_max`, `local_error_min`, `local_error_avg`, `global_error_ratio_max/min/avg`, `global_error_percent_max/min/avg`. There is **no** `local_error_percent_max`.

## Environment

- Julia version managed via `juliaup` (default: `release` channel)
- Package is `dev`'d in default Julia environment (`~/.julia/environments/v1.12/`)
- Notebooks run via IJulia `julia-1.12` kernel from the same default environment
- Tests run via `julia --project=. -e 'using Pkg; Pkg.test()'`

## Releasing a New Version

**Claude Code handles all CHANGELOG entries for this project.** See conventions below.

### Steps

1. Write release notes under `## [Unreleased]` in `CHANGELOG.md` (see conventions below)
2. Update the version number in the `CLAUDE.md` Project description
3. Commit all changed source/doc/test files
4. Run `./release.sh patch|minor|major`

`release.sh` handles everything else: bumps `Project.toml`, moves `[Unreleased]` entries to the
new version heading, commits, pushes, and posts the `@JuliaRegistrator register` comment with
the changelog content embedded.

### Full automated pipeline (after `release.sh`)

```
release.sh posts @JuliaRegistrator register comment
  → JuliaRegistrator opens PR to JuliaRegistries/General (~15 min)
  → AutoMerge validates and merges (see requirements below)
  → JuliaTagBot creates the git tag automatically
  → TagBot.yml creates a GitHub Release using CHANGELOG.md (changelog: true)
  → Docs.yml deploys /stable/ on the new tag
```

No manual `git tag` or GitHub Release creation needed — TagBot and Docs deploy automatically.

### CHANGELOG conventions

**Format**: Keep a Changelog (`## [Unreleased]` at top, sections like `### Added`, `### Fixed`,
`### Changed`, `### Breaking changes`).

**AutoMerge requirement**: JuliaRegistries AutoMerge requires the words **"breaking"** or
**"changelog"** to appear in the release notes for any breaking version bump. The release script
prepends `"See CHANGELOG.md for full details."` to every registrator comment, so the word
"changelog" is always present. Still, use `### Breaking changes` as the section header whenever
there are actual breaking changes — it makes the notes self-documenting and clear to users.

**What counts as breaking (Julia semver)**:
- Package is pre-1.0 (`0.x.x`): every **minor** bump (`0.2 → 0.3`) is treated as breaking by
  the registry. Patch bumps (`0.2.0 → 0.2.1`) are not.
- Post-1.0: only **major** bumps are breaking.

**Example CHANGELOG entry for a breaking minor release**:

```markdown
## [Unreleased]

### Breaking changes
- `foo` now returns `Bar` instead of `Baz`. Upgrade: replace `x::Baz` with `x::Bar`.

### Added
- New `bar` function for ...
```

**Guard**: `release.sh` will exit with an error (before touching git) if the new version section
is empty. Always populate `[Unreleased]` before running the release script.

### GitHub account

Always use the **WeberElectrodynamics** account for `gh` CLI operations on this repo.
Only `WeberElectrodynamics` has collaborator status — JuliaRegistrator will reject comments
from any other account. Before running `release.sh` or any `gh api` command, verify:

```bash
gh auth switch --user WeberElectrodynamics
gh auth status   # confirm Active account: WeberElectrodynamics
```
