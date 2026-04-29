# AGENTS.md

This file provides guidance to AI coding agents working in this repository. Follow the closest
`AGENTS.md` file in the tree; explicit user instructions take precedence.

## Project

Julia package for n-body Weber electrodynamics simulation with Zöllner electrogravitational extension. Version is tracked in `Project.toml` — do not duplicate it here.

Implements a symplectic Strang-splitting symmetric-projection integrator. Levi-Civita/KS regularization handles close encounters; collision bounce handles head-on singularities.

## Commands

```bash
# Run full test suite (from project root)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file (e.g., test_physics.jl)
julia --project=. -e 'using Test; using WeberElectrodynamics; using WeberElectrodynamics: SymmetricProjectionIntegrator; using LinearAlgebra; using Symbolics; @testset "single" begin include("test/test_utils.jl"); include("test/test_physics.jl") end'

# Format all Julia files (requires JuliaFormatter in the active environment)
julia -e 'using JuliaFormatter; format(".")'

# Build documentation locally
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl

# Release a new version — dry-run (safe, no changes)
./release.sh patch
# Release a new version — actually execute
./release.sh patch --execute
```

## Repository Structure

`src/`, `ext/`, `test/`, `examples/` — package source, weak-dep extensions,
tests, reference notebook. `docs/` — Documenter scaffold. `theory/` — math
derivations. `_research/` — **exploratory sandbox**, not authoritative.
`papers/<name>/` — LaTeX papers with their own `VERSION` files. See
[docs/src/internals.md](docs/src/internals.md) for per-file descriptions.

## Architecture

**Pipeline**: Symbolic Hamiltonian → compiled equations of motion → symplectic integration → statistics/plotting/animation

See [docs/src/internals.md](docs/src/internals.md) for per-file descriptions of `src/`, extensions, tests, and the `_research/` sandbox.

## Critical Conventions

Quick-reference mirror — source of truth: [docs/src/internals.md](docs/src/internals.md).

### Params and κ vector layout

`params = [m₁…mₙ, q₁…qₙ, c]` — length `2N + 1`.
`kappas = [κ₁₂, κ₁₃, …, κ_{N-1,N}]` — length `N*(N-1)/2`, indexed by `_pair_index(i, j, n)`.
Compiled EOM signature: `dq_dt_compiled(out, q, p, t, params, kappas)` (same for `dp_dt` and `hamiltonian_compiled`). Direct callers **must** pass both; when Zöllner is disabled, `kappas = ones(N*(N-1)÷2)`.
Per-pair accessor: `kappa(prob, i, j)`.

### Algorithms and callbacks

Regularization and collision bounce are composed outside the problem:

- **Algorithm wrapper**: `RegularizedIntegrator(SymmetricProjectionIntegrator(); backend=:lifted_pair|:adaptive_cartesian, r_on_factor=..., r_off_factor=..., max_substeps=..., ...)`.
  Pass to `solve(prob, alg)`.
- **Callback**: `CollisionBounce(radius)` passed via `solve(prob, alg; callbacks=CollisionBounce(r))`.
- Regularization backends: `:adaptive_cartesian` (2D+3D) and `:lifted_pair` (2D only). Neither regularizes Weber's velocity-dependent force.

### Accessor API

`HamiltonianProblem` stores `params` and `kappas` as fields; `masses`, `charges`,
and `c` are O(1) views into `params`. `ZollnerOptions` and `RegularizationOptions`
are construction-time only and not retained on the problem. Use the exported
accessors: `masses(prob)`, `charges(prob)`, `speed_of_light(prob)`, `kappas(prob)`,
`kappa(prob, i, j)`, `params(prob)`, plus `n_particles(sys)`, `dims(sys)`,
`degrees_of_freedom(sys)` on `HamiltonianSystem`.

### EnergyStatistics gotcha

No `local_error_percent_max` field — use `local_error_max`. Full field list in
[internals.md](docs/src/internals.md).

## Environment

- Package supports Julia `1.9+`; CI tests Julia `1.9` and latest stable on Ubuntu, plus
  latest stable on macOS and Windows.
- Prefer `julia --project=.` from the repository root for package work.
- Notebooks should run from a Julia environment where this package is available, usually via
  `Pkg.develop(path=pwd())`; source, docs, and tests remain authoritative over notebook state.
- Tests run via `julia --project=. -e 'using Pkg; Pkg.test()'`.

## Git Workflow

### Branch strategy

`main` is always stable and releasable. All non-trivial changes go through a branch and PR.

**When to push directly to `main`:** only for truly trivial changes — a typo fix, a single-word
doc tweak, a one-liner correction. Everything else gets a branch.

**Branch naming:**

| Prefix | Use for |
|--------|---------|
| `feature/<name>` | New functionality |
| `fix/<name>` | Bug fixes |
| `docs/<name>` | Documentation only |
| `refactor/<name>` | Internal restructuring, no behaviour change |
| `experiment/<name>` | Exploratory / research branches — may not merge |

Examples: `feature/3d-regularization`, `fix/energy-drift`, `docs/zollner-theory`,
`experiment/three-body-bound-states`

### Pull requests

- Open a PR to `main` for all non-trivial branches — **open it immediately after pushing, by default**
- Title should match the eventual commit message (conventional commit format, see below)
- Tests must pass (CI runs automatically on push)
- Merge strategy: **regular merge commit** (not squash)
- Never force-push `main`

### Commit message conventions

This repo uses [Conventional Commits](https://www.conventionalcommits.org/):

```
feat:      new user-facing functionality
fix:       bug fix
docs:      documentation only
refactor:  internal restructuring, no behaviour change
test:      adding or fixing tests
chore:     tooling, CI, dependencies
release:   version bump commits made by release.sh
```

Examples:
```
feat: add 3D Levi-Civita regularization backend
fix: correct sign flip in 1D regularization lift
docs: add Zöllner theory page
refactor: extract pair-distance helpers into regularization.jl
test: add collision bounce smoke test for 3D case
chore: bump Symbolics compat to 7
release: v0.2.3
```

## Releasing a New Version

Agents write `RELEASENOTES.md`; `release.sh` promotes it into `CHANGELOG.md`. See conventions below.

### Pre-release checklist

1. Commit all changed source/doc/test files
2. Write release notes in `RELEASENOTES.md` (the developer-facing input point)
3. If the paper changed, bump `papers/<name>/VERSION` (see [Paper versioning](#paper-versioning))
4. Dry-run: `./release.sh patch` (or `minor`/`major`) — safe, no changes made
5. Execute: `./release.sh patch --execute`

`release.sh --execute` handles everything else: bumps `Project.toml`, moves notes from
`RELEASENOTES.md` into `CHANGELOG.md`, resets `RELEASENOTES.md`, commits, pushes, and posts
the `@JuliaRegistrator register` comment.

### release.sh usage

```bash
./release.sh [patch|minor|major]            # dry-run (default, safe)
./release.sh [patch|minor|major] --execute  # actually release
```

**Dry-run** prints a preview of the version bump, the release notes, and every action that would
be taken — no files are modified.

**Gate**: both modes fail immediately if `RELEASENOTES.md` is empty (blank lines and
`<!-- comments -->` don't count).

### Full automated pipeline (after `release.sh --execute`)

```
release.sh --execute
  → commits Project.toml + CHANGELOG.md + RELEASENOTES.md, pushes
  → posts @JuliaRegistrator register comment
  → JuliaRegistrator opens PR to JuliaRegistries/General (~15 min)
  → AutoMerge validates and merges
  → JuliaTagBot creates the git tag automatically
  → TagBot.yml creates a GitHub Release (description = release notes)
  → Papers.yml fires on release published → compiles PDFs → uploads as assets
  → Docs.yml deploys /stable/ on the new tag
```

No manual `git tag`, GitHub Release creation, or PDF upload needed — all automated.

### CHANGELOG conventions

`CHANGELOG.md` is populated **automatically** from `RELEASENOTES.md` by `release.sh`.
Agents still write `RELEASENOTES.md` content when preparing a release.

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

**Example `RELEASENOTES.md` for a breaking minor release**:

```markdown
### Breaking changes
- `foo` now returns `Bar` instead of `Baz`. Upgrade: replace `x::Baz` with `x::Bar`.

### Added
- New `bar` function for ...
```

### GitHub account

Always use the **WeberElectrodynamics** account for `gh` CLI operations on this repo.
Only `WeberElectrodynamics` has collaborator status — JuliaRegistrator will reject comments
from any other account. Before running `release.sh` or any `gh api` command, verify:

```bash
gh auth switch --user WeberElectrodynamics
gh auth status   # confirm Active account: WeberElectrodynamics
```

## Paper versioning

Papers live in `papers/<name>/` and each has its own independent version in `papers/<name>/VERSION`.

**Convention**: `MAJOR.MINOR` (e.g. `1.0`, `1.1`, `2.0`) — independent of the Julia package version.

| Bump | When |
|------|------|
| `MINOR` (1.0 → 1.1) | Corrections, new sections, new references |
| `MAJOR` (1.0 → 2.0) | Substantial restructuring or a new paper version |

The VERSION file is committed as part of the normal pre-release commit. The PDF filename on the
Releases page is `<paper-name>-v<version>.pdf` (e.g. `Computational-Weber-Electrodynamics-v1.0.pdf`).

If the paper has **not** changed since the last release, leave VERSION as-is — the same filename
will simply be re-uploaded to the new release.

### Adding a new paper

1. Create `papers/<new-paper>/` with `.tex`, `.bib`, and a `VERSION` file (start at `1.0`)
2. Add `- <new-paper>` to the matrix in `.github/workflows/Papers.yml`
3. That's it — the paper compiles and uploads on every future release automatically

### Zenodo & arXiv notes

- **Zenodo**: `.zenodo.json` is already present at the repo root → Zenodo auto-archives every
  GitHub Release. The PDF asset is included automatically. Each paper will eventually get its
  own separate Zenodo record.
- **arXiv**: arXiv requires source submission (`.tex` + `.bib` + figures), not the compiled PDF.
  This is a deliberate manual step when ready — no automation needed.
