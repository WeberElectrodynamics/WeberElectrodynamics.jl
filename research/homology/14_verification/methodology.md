# Verification Protocol for Bound Weber Orbits

## Purpose

Classify trajectories of the Weber electrodynamics n-body system into one of four
categories based on phase-space topology and dynamical stability.

## Classification categories

| Category       | Definition |
|----------------|------------|
| PERIODIC       | Phase-space return: relative distance < 0.01 at detected period T, energy drift < 1e-4 over integration |
| QUASI-PERIODIC | Bounded trajectory with no exact phase-space return; energy conserved; perturbation-stable |
| CHAOTIC-BOUND  | Bounded trajectory with positive Lyapunov exponent (amplification > 10^4) |
| ESCAPE         | At least one pair distance exceeds escape_ratio times the initial maximum separation |

## 5-step verification pipeline

### Step 1: Base integration and period detection

Integrate the initial condition (q0, p0) using the symplectic Strang-splitting
symmetric-projection integrator at timestep dt for duration tmax.

**Period detection** uses two methods:
- **Autocorrelation**: Compute the autocorrelation of the primary pair-distance signal.
  The first significant peak (height > 0.3) gives the candidate period.
- **Circular orbit detection**: If the pair-distance variance is negligible
  (var / (mean^2 * N) < 1e-10), scan phase-space return distance and take the first
  local minimum below 0.01 relative return.
- **Phase-space scan fallback**: If autocorrelation finds no period, scan for the first
  local minimum of ||z(t) - z(0)|| / ||z(0)|| below 0.05.

**Escape check**: If max pair distance > escape_ratio * max(r0_max, 1.0), classify as ESCAPE
immediately.

**Energy drift**: Compute max |H(t) - H(0)| / |H(0)| over all saved timesteps using
the full Weber Hamiltonian (kinetic + velocity-dependent potential).

### Step 2: Convergence check

Re-integrate with dt/2 for min(2T, tmax). Compare the phase-space state at the end:

    convergence_discrepancy = ||z_base(t) - z_half(t)||

This validates that the base timestep is sufficient. A large discrepancy indicates the
orbit requires finer resolution.

### Step 3: Extended integration

If a period T was detected and n_periods_check * T > tmax, re-integrate for
n_periods_check * T (default: 5 periods). Check:
- Energy drift does not grow
- No escape occurs in the extended window

### Step 4: Perturbation test (Lyapunov-type)

Apply a small deterministic perturbation (eps = 1e-4, direction seeded by orbit ID hash)
to the initial condition. Integrate both reference and perturbed trajectories for
min(2T, tmax). Measure:

    lyapunov_separation = ||z_ref(t_end) - z_pert(t_end)||
    amplification = lyapunov_separation / eps

Amplification >> 1 indicates sensitive dependence on initial conditions.

### Step 5: Classification

Apply the decision tree:

1. If max pair distance exceeds escape threshold: **ESCAPE**
2. If phase-space return < periodic_tol (relative) AND extended energy drift < 1e-8: **PERIODIC**
3. If phase-space return < 0.01 (relative) AND energy drift < energy_tol: **PERIODIC**
4. If Lyapunov amplification > 10^4: **CHAOTIC-BOUND**
5. If bounded with phase return < 0.5 AND good energy conservation: **QUASI-PERIODIC**
6. If bounded, no period detected, good energy: **QUASI-PERIODIC**
7. Otherwise: **Unclassified**

## Parameters

| Parameter         | Default | Description |
|-------------------|---------|-------------|
| tmax              | 100.0   | Base integration time |
| dt                | 1e-3    | Base timestep |
| n_periods_check   | 5       | Number of periods for extended check |
| energy_tol        | 1e-4    | Energy drift tolerance for classification |
| periodic_tol      | 1e-6    | Relative phase-space return for strict periodic |
| escape_ratio      | 10.0    | Max separation / initial separation threshold |
| perturbation_eps  | 1e-4    | Perturbation magnitude for Lyapunov test |
| bounce_r          | 0.0     | Collision bounce radius (for like-charge encounters) |

## Known limitations

1. **Unlike-charge collisions**: The Weber integrator cannot handle head-on collisions
   between oppositely charged particles. Neither collision bounce (designed for
   like-charge repulsion) nor Levi-Civita regularization (does not regularize
   velocity-dependent Weber force) resolves this. Avoid initial conditions where
   opposite charges can collide.

2. **Weber precession**: Elliptical orbits in Weber electrodynamics precess due to the
   velocity-dependent correction. This means Kepler-periodic orbits become quasi-periodic
   under Weber dynamics. This is physical, not numerical.

3. **Period detection for chaotic orbits**: Autocorrelation-based detection fails for
   chaotic trajectories (no periodic signal). The pipeline correctly returns NaN for
   the period in such cases.

4. **Lyapunov estimate is single-shot**: The perturbation test computes separation at
   one time only, not a proper maximum Lyapunov exponent. It distinguishes exponential
   from polynomial sensitivity but does not give a quantitative exponent.

## API

```julia
verify_orbit(id, q0, p0, masses, charges, c;
             n_particles, dims, tmax=100.0, dt=1e-3, ...) -> VerificationResult

# VerificationResult fields:
#   .classification :: Symbol    (:Periodic, :QuasiPeriodic, :ChaoticBound, :Escape, :Failure)
#   .period         :: Float64   (detected period, NaN if none)
#   .energy_drift_base      :: Float64
#   .energy_drift_extended  :: Float64
#   .lyapunov_separation    :: Float64
#   .convergence_discrepancy :: Float64
#   .phase_return           :: Float64
```
