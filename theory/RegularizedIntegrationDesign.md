# Regularized Integration Design for Weber Hamiltonian

## Scope

This document specifies the regularized integration design implemented for the Weber Hamiltonian.

- Scope is limited to the Weber Hamiltonian produced by `HamiltonianSystem`.
- Public solver signatures (`solve/init/step!/solve!`) remain unchanged.
- Regularization is configured by wrapping the base algorithm in `RegularizedIntegrator` and is disabled by default.
- Far-field evolution remains on the existing Cartesian projected kernel.

## Backend Semantics

Regularization backend is configured through `RegularizedIntegrator(base_alg; backend = ...)`:

- `:lifted_pair`
- `:adaptive_cartesian`

Effective backend is resolved at init. For supported dimensions (`1`, `2`, and
`3`), both configured backends are available for binary pair mode. Fallback
behavior is reserved for unsupported future dimensions:

- one init-time warning when `warn_on_fallback=true`
- diagnostics counter increment on fallback pair steps

Support matrix:

- 1D pair mode: lifted square-root or adaptive Cartesian backend
- 2D pair mode: lifted Levi-Civita or adaptive Cartesian backend
- 3D pair mode: lifted KS or adaptive Cartesian backend
- chain mode (all dims): adaptive Cartesian backend

## Encounter Dispatch and Hysteresis

Each outer step performs:

1. Compute all pair distances.
2. Build encounter graph for pairs with `r <= r_on`.
3. Select the connected component containing the minimum-distance pair.
4. Apply hysteresis:
   - activate at `r_on`
   - remain active while active-component max distance `<= r_off`
   - deactivate otherwise

   Note: the anchor pair `(i, j)` that triggered activation is held fixed for
   the duration of the active period. During the deactivation check the
   component is always rooted at the original anchor (using `r_off` adjacency),
   not re-derived from the current global minimum-distance pair. This prevents
   mode-hopping when two different pairs trade the minimum-distance position
   across consecutive steps.

5. Choose mode:
   - component size 2: pair mode
   - component size > 2 with chain enabled: chain mode
   - otherwise: Cartesian mode

## Pair Mode: `:adaptive_cartesian`

The adaptive Cartesian pair backend keeps the previous robust path:

1. Detect active pair.
2. Compute adaptive substep count: `substeps = clamp(⌈r_on / max(r, g_floor)⌉, 1, max_substeps)`.
3. For each substep, run the existing projected Cartesian kernel.
4. In 3D, apply KS constraint projection in lifted diagnostics path.

## Pair Mode: `:lifted_pair`

Lifted pair mode uses a split method:

- `A`: external perturbation half-step (physical time, midpoint)
- `B`: lifted pair full-step (1D square-root, 2D LC, or 3D KS coordinates)
- `A`: external perturbation half-step

### Derivative decomposition

For active pair `(i,j)`, build `params_pair`:

- masses unchanged
- charges zero for all particles except `i,j`
- speed of light unchanged

Then evaluate with existing compiled Weber RHS:

- full derivatives: `dq_full, dp_full` with full params
- pair derivatives: `dq_pair, dp_pair` with `params_pair`
- external perturbation: `dq_ext = dq_full - dq_pair`, `dp_ext = dp_full - dp_pair`

### External half-step (`A`)

Use explicit midpoint for perturbation over `dt/2`:

1. evaluate ext derivatives at start
2. build midpoint state
3. re-evaluate at midpoint
4. update physical `q,p`

### Multi-substep composition

Multiple A-B-A substeps are applied back-to-back within one macro-step. The
substep size is adaptive: `dt_sub = min(r_current * dtau_target, t_remaining)`
where `dtau_target = dt / (1.5 * r_on)`. Adjacent A half-steps between
substeps are not merged, so the structure is
`A(h₁/2)–B(h₁)–A(h₁/2)–A(h₂/2)–B(h₂)–A(h₂/2)–…`. If `max_substeps` is
exhausted before `t_remaining` reaches zero, the remainder is taken as a
single final substep.

### Lifted pair step (`B`)

For each pair substep:

1. Convert pair to COM + relative variables.
2. Lift relative state with Levi-Civita map.
3. Enforce LC sheet continuity (`u,U` sign) against previous substep state.
4. Use monitor `g = max(r, g_floor)` at substep start and freeze it within that substep.
5. Integrate `(u,U)` with explicit midpoint in fictitious time (`dt = g dτ`).
6. Project LC state back to relative Cartesian and reconstruct pair particles in-place.

The frozen-per-substep monitor avoids distortion from re-scaling midpoint stages with different `g` values.

## Chain Mode

Chain mode uses adaptive Cartesian integration:

1. Build deterministic chain ordering from closest-link traversal.
2. Compute active-component monitor
   `Ω = Σ 1/r_ij` over active component pairs.
3. Use `g = max(1/Ω, g_floor)` for substep sizing.
4. Advance each substep with projected Cartesian kernel.

## KS Constraint Handling

3D regularized paths have two KS-related modes:

- `:lifted_pair` uses a KS pair chart with fictitious time for binary close
  encounters.
- `:adaptive_cartesian` still lifts to KS variables for active-pair diagnostics,
  projects momentum to satisfy the bilinear KS constraint, and then advances the
  Cartesian projected kernel.
- Diagnostics track max constraint violation plus KS projection/iteration counts.

True chain-coordinate KS regularization for multi-particle clusters remains
deferred; chain mode is adaptive Cartesian over the active component.

## Diagnostics

`RegularizationDiagnostics` includes:

- `requested_backend`, `used_backend`
- activation/deactivation counters
- `active_steps`: steps taken while regularization was active (`pair_steps + chain_steps`)
- `pair_steps`, `adaptive_pair_steps`, `lifted_pair_steps`, `chain_steps`, `unregularized_steps`
- `backend_fallback_steps`
- `total_substeps`, `max_substeps_used`
- `min_encounter_distance`, `max_constraint_violation`
- `ks_constraint_projection_count`, `ks_constraint_iteration_count`
- per-step mode history (`0` Cartesian, `1` pair, `2` chain)

`used_backend` semantics:

- `:adaptive_cartesian` or `:lifted_pair` when only one regularized backend is used
- `:mixed` when both are used in one run
- `:disabled` when regularization is disabled, or when regularization is enabled
  but no encounter was detected during the run (zero regularized steps taken)

## Collision Bounce

`RegularizationOptions` accepts `collision_bounce_radius::Float64 = 0.0`
(disabled by default). When positive, a pre-step reflection is applied to any
pair closer than this radius: `q_rel → -q_rel` with momenta unchanged. This is
the C⁰-continuation of a head-on (ℓ=0) collision and preserves energy exactly.
The feature is independent of the regularization backend and may be used with
`enabled = false`.

## Memory Model

Regularization workspaces are preallocated and reused:

- encounter graph and pair-distance buffers
- backend state and chain ordering buffers
- full/pair/external derivative buffers
- midpoint and split-step scratch buffers
- LC/KS temporary state buffers

No per-step heap allocation is required in regularization helper kernels; unregularized allocation profile is unchanged.

See also: [Regularization.md](Regularization.md), [SemiExplicitIntegrator.md](SemiExplicitIntegrator.md), [WeberElectrodynamics.md](WeberElectrodynamics.md).
