# Regularized Integration Design for Weber Hamiltonian

## Scope

This document specifies the regularized integration design implemented for WeberElectrodynamics v3 continuation.

- Scope is limited to the Weber Hamiltonian produced by `WeberSystem`.
- Public solver signatures (`solve/init/step!/solve!`) remain unchanged.
- Regularization is problem-level configurable and enabled by default.
- Far-field evolution remains on the existing Cartesian projected kernel.

## Backend Semantics

Regularization backend is configured through `RegularizationOptions.backend` (or `WeberProblem(...; regularization_backend=...)`):

- `:lifted_pair`
- `:adaptive_cartesian`

Effective backend is resolved at init:

- 2D: requested backend is honored.
- 1D/3D with requested `:lifted_pair`: fallback to `:adaptive_cartesian`.

Fallback behavior:

- one init-time warning when `warn_on_fallback=true`
- diagnostics counter increment on fallback pair steps

Support matrix in this milestone:

- 2D pair mode: true lifted pair backend available
- 1D pair mode: adaptive Cartesian backend
- 3D pair mode: adaptive Cartesian backend
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
5. Choose mode:
   - component size 2: pair mode
   - component size > 2 with chain enabled: chain mode
   - otherwise: Cartesian mode

## Pair Mode: `:adaptive_cartesian`

The adaptive Cartesian pair backend keeps the previous robust path:

1. Detect active pair.
2. Compute adaptive substep count from encounter scale.
3. For each substep, run the existing projected Cartesian kernel.
4. In 3D, apply KS constraint projection in lifted diagnostics path.

## Pair Mode: `:lifted_pair` (2D)

2D lifted pair mode uses a split method:

- `A`: external perturbation half-step (physical time, midpoint)
- `B`: lifted pair full-step (LC coordinates)
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

Chain mode remains on adaptive Cartesian integration in this milestone:

1. Build deterministic chain ordering from closest-link traversal.
2. Compute active-component monitor
   `Ω = Σ 1/r_ij` over active component pairs.
3. Use `g = max(1/Ω, g_floor)` for substep sizing.
4. Advance each substep with projected Cartesian kernel.

## KS Constraint Handling

3D regularized paths still use KS diagnostics support in adaptive Cartesian mode:

- lift to KS variables for active pair checks
- project momentum to satisfy bilinear KS constraint
- track max constraint violation in diagnostics

True 3D lifted KS stepping is deferred.

## Diagnostics

`RegularizationDiagnostics` includes:

- `requested_backend`, `used_backend`
- activation/deactivation counters
- `pair_steps`, `adaptive_pair_steps`, `lifted_pair_steps`, `chain_steps`, `unregularized_steps`
- `backend_fallback_steps`
- `total_substeps`, `max_substeps_used`
- `min_encounter_distance`, `max_constraint_violation`
- per-step mode history (`0` Cartesian, `1` pair, `2` chain)

`used_backend` semantics:

- `:adaptive_cartesian` or `:lifted_pair` when only one regularized backend is used
- `:mixed` when both are used in one run
- `:disabled` when regularization is disabled

## Memory Model

Regularization workspaces are preallocated and reused:

- encounter graph and pair-distance buffers
- backend state and chain ordering buffers
- full/pair/external derivative buffers
- midpoint and split-step scratch buffers
- LC/KS temporary state buffers

No per-step heap allocation is required in regularization helper kernels; unregularized allocation profile is unchanged.

See also: [Regularization.md](Regularization.md), [SemiExplicitIntegrator.md](SemiExplicitIntegrator.md), [WeberElectrodynamics.md](WeberElectrodynamics.md).
