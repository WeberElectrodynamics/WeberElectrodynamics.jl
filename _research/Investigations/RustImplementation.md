# High-Performance Rust Implementation of the Semi-Explicit Symplectic Integrator

Research notes on porting the semi-explicit symmetric-projection integrator to Rust with multi-threading and Newton's third law optimization.

## 1. Problem Statement

### 1.1 The N-Squared Bottleneck

The Julia implementation uses Symbolics.jl to build the Weber Hamiltonian symbolically and auto-differentiate it into compiled equations of motion (`dq_dt_compiled`, `dp_dt_compiled`). These compiled functions are **particle-centric**: each call computes all N*d derivatives in a single monolithic evaluation. Internally, every particle i sums contributions from all N-1 partners, yielding N*(N-1) pairwise force evaluations per RHS call.

But Newton's third law guarantees $\mathbf{F}_{ij} = -\mathbf{F}_{ji}$, meaning only N*(N-1)/2 unique pair interactions exist. The symbolic compilation obscures this symmetry because `Symbolics.derivative(H, q_vars[k])` differentiates the full Hamiltonian with respect to each coordinate independently, producing code that re-evaluates the same pair interaction from both sides.

### 1.2 Amplification by the Integrator

The Strang splitting flow (see `solve.jl:4-50`) calls `dq_dt_compiled` and `dp_dt_compiled` **6 times** per step:

```
Flow A(dt/2):  dq_dt(buf, Q, Y, params)  +  dp_dt(buf, Q, Y, params)
Flow B(dt):    dq_dt(buf, X, P, params)  +  dp_dt(buf, X, P, params)
Flow A(dt/2):  dq_dt(buf, Q, Y, params)  +  dp_dt(buf, Q, Y, params)
```

The symmetric projection iteration (`solve.jl:92-233`) re-evaluates the full Strang splitting at every iteration (typically 2-3 iterations). Including the convergence-check evaluation, the total is:

| Projection iterations | RHS evaluations per timestep |
|----------------------|------------------------------|
| 2                    | 18                           |
| 3                    | 24                           |

Each of these RHS evaluations performs N*(N-1) pair computations instead of N*(N-1)/2.

### 1.3 Why Rust

| Concern | Julia (current) | Rust (proposed) |
|---------|-----------------|-----------------|
| Force symmetry | N*(N-1) per RHS call (symbolic compilation obscures symmetry) | N*(N-1)/2 via pair-centric loop with explicit accumulation |
| Parallelism | Single-threaded RHS evaluation | Rayon work-stealing over pair loop |
| Memory control | GC-managed, pre-allocated buffers mitigate but don't eliminate | Zero-cost ownership, no GC pauses in hot path |
| SIMD | LLVM auto-vectorization (hit-or-miss) | Explicit `std::simd` or `std::arch` intrinsics |
| First-call latency | Symbolics.jl compilation: seconds to minutes for large N | AOT compilation, zero startup cost |

### 1.4 Target Speedup

Three multiplicative factors:

1. **2x** from Newton's 3rd law: N*(N-1)/2 vs N*(N-1)
2. **~P** from multi-threading on P cores (sub-linear due to Amdahl's law on serial sections)
3. **1.5-3x** from SoA memory layout, explicit SIMD, and elimination of GC

Conservative estimate for N=100 on 8 cores: **20-40x** per timestep.

## 2. Architecture

### 2.1 Design Inversion: Pair-Centric vs Particle-Centric

The Julia implementation is **particle-centric**: `dq_dt_compiled(out, q, p, params)` writes all N*d derivatives at once. The Rust implementation inverts this to be **pair-centric**: iterate over unique pairs (i,j) with i < j, compute the pair contribution once, and accumulate to both particles.

The Julia code in `weber_system.jl:77-106` already iterates `for i=1:n, j=(i+1):n` when building the Hamiltonian. But Symbolics.jl then differentiates per-component, losing the pair structure. In Rust, we differentiate analytically at the pair level and keep it.

### 2.2 Crate Structure

```
weber-integrator/
  src/
    lib.rs              -- Public API, re-exports
    types.rs            -- Particle, System, Problem, Solution, IntegratorState
    force.rs            -- PairwiseForce trait, WeberForce implementation
    accumulator.rs      -- Thread-safe force accumulation strategies
    pair_loop.rs        -- N(N-1)/2 pair iteration, parallel dispatch
    integrator/
      mod.rs            -- init/step/solve interface
      strang.rs         -- Strang splitting flow in extended phase space
      projection.rs     -- Symmetric projection fixed-point iteration
      composition.rs    -- Higher-order: triple jump, Suzuki, Yoshida
    params.rs           -- Params vector layout (Julia-compatible)
  benches/
    pair_loop.rs        -- Criterion benchmarks for pair loop
    integration.rs      -- End-to-end timestep benchmarks
```

Regularization is out of scope for the initial Rust port (it operates on isolated pair subsystems and can be added later without architectural changes).

## 3. Force Computation Interface

### 3.1 Splitting Kinetic and Interaction Contributions

The Hamiltonian $H = T(p) + V(q, p)$ has two types of contributions to the equations of motion:

**Kinetic (per-particle, trivial):**
$$\frac{\partial T}{\partial p_i} = \frac{p_i}{m_i}$$

This contributes to $\dot{q}_i$ and is computed in a simple O(N) per-particle loop with no pairwise interaction.

**Interaction (pairwise, expensive):**

The Weber potential couples position and momentum:
$$U_{ij} = \frac{\kappa_{ij} \, q_i \, q_j}{r_{ij}} \left(1 - \frac{\dot{r}_{ij}^2}{2c^2}\right)$$

where $\dot{r}_{ij} = \hat{r}_{ij} \cdot (v_i - v_j)$ and $v_k = p_k / m_k$.

This contributes to both $\dot{q}$ (via $\partial U / \partial p$, the Weber velocity correction) and $\dot{p}$ (via $-\partial U / \partial q$, the force). The pairwise structure guarantees:

$$\frac{\partial U_{ij}}{\partial q_i} = -\frac{\partial U_{ij}}{\partial q_j}, \qquad \frac{\partial U_{ij}}{\partial p_i} = -\frac{m_j}{m_i} \frac{\partial U_{ij}}{\partial p_j} \quad \text{(up to mass factors)}$$

More precisely, the momentum force satisfies exact Newton's third law: $\dot{p}_i^{(ij)} = -\dot{p}_j^{(ij)}$.

### 3.2 The PairwiseForce Trait

```rust
/// Contribution from a single pair (i, j) to the equations of motion.
/// All arrays are d-dimensional (1, 2, or 3).
pub struct PairContribution<const D: usize> {
    /// Additive correction to dq/dt for particle i (Weber velocity term)
    pub dq_i: [f64; D],
    /// Additive correction to dq/dt for particle j
    pub dq_j: [f64; D],
    /// Force on particle i from this pair: dp_i/dt contribution
    pub dp_i: [f64; D],
    /// Force on particle j (must equal -dp_i for momentum conservation)
    pub dp_j: [f64; D],
}

/// Trait for pairwise velocity-dependent interactions.
/// Implementors compute the contribution of a single pair (i, j) with i < j.
pub trait PairwiseForce<const D: usize>: Send + Sync {
    fn evaluate(
        &self,
        q_i: &[f64; D], q_j: &[f64; D],
        p_i: &[f64; D], p_j: &[f64; D],
        mass_i: f64, mass_j: f64,
        charge_i: f64, charge_j: f64,
        kappa_ij: f64,
    ) -> PairContribution<D>;
}
```

The trait is generic over dimension `D` (const generic) and requires `Send + Sync` for multi-threaded pair loops. The integrator framework calls this once per unique pair and accumulates the results.

### 3.3 Weber Force Implementation

The Weber pair force, analytically derived from $U_{ij}$, computes shared intermediates once:

```rust
pub struct WeberForce {
    pub c: f64,
}

impl PairwiseForce<3> for WeberForce {
    fn evaluate(
        &self,
        q_i: &[f64; 3], q_j: &[f64; 3],
        p_i: &[f64; 3], p_j: &[f64; 3],
        mass_i: f64, mass_j: f64,
        charge_i: f64, charge_j: f64,
        kappa_ij: f64,
    ) -> PairContribution<3> {
        let c2 = self.c * self.c;

        // Relative position and velocity
        let dr = [q_i[0] - q_j[0], q_i[1] - q_j[1], q_i[2] - q_j[2]];
        let r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        let r = r2.sqrt();
        let r_inv = 1.0 / r;

        let vi = [p_i[0]/mass_i, p_i[1]/mass_i, p_i[2]/mass_i];
        let vj = [p_j[0]/mass_j, p_j[1]/mass_j, p_j[2]/mass_j];
        let dv = [vi[0] - vj[0], vi[1] - vj[1], vi[2] - vj[2]];

        // Radial velocity: rdot = (dr . dv) / r
        let dr_dot_dv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
        let rdot = dr_dot_dv * r_inv;

        let k = kappa_ij * charge_i * charge_j;

        // --- dp/dt contributions (force) ---
        // F_ij on particle i from the Weber potential:
        //   dp_i/dt = -dU/dq_i
        //   dp_j/dt = -dU/dq_j = +dU/dq_i = -dp_i/dt  (Newton's 3rd law)
        //
        // The full Weber force expression (from analytical differentiation):
        //   F_ij = k/r^2 * r_hat * [ -(1 - rdot^2/(2c^2)) + rdot^2/c^2 ]
        //        + k/(r * c^2) * rdot * dv_perp_terms
        //        ... (see theory/WeberElectrodynamics.md for full expression)
        //
        // Implementation computes dp_i; dp_j = -dp_i exactly.

        let coulomb = k * r_inv * r_inv; // k / r^2
        let weber_r = rdot * rdot / c2;
        let r_hat = [dr[0] * r_inv, dr[1] * r_inv, dr[2] * r_inv];

        // Radial force component
        let f_radial = -coulomb * (1.0 - weber_r / 2.0 + weber_r);

        // Velocity-dependent transverse component
        let f_vel_coeff = k * r_inv * rdot / c2;

        let mut dp_i = [0.0f64; 3];
        for d in 0..3 {
            dp_i[d] = f_radial * r_hat[d] + f_vel_coeff * dv[d];
        }

        // Newton's third law: dp_j = -dp_i
        let dp_j = [-dp_i[0], -dp_i[1], -dp_i[2]];

        // --- dq/dt contributions (Weber velocity correction) ---
        // dq_i/dt has a correction: d(U_ij)/dp_i
        // This arises because U_ij depends on p through rdot.
        // The correction is: (k / (r * c^2)) * rdot * dr / m_i
        // with opposite sign for particle j.
        let vel_correction_coeff = k * r_inv * rdot / c2;
        let dq_i = [
            vel_correction_coeff * dr[0] / mass_i,
            vel_correction_coeff * dr[1] / mass_i,
            vel_correction_coeff * dr[2] / mass_i,
        ];
        let dq_j = [
            -vel_correction_coeff * dr[0] / mass_j,
            -vel_correction_coeff * dr[1] / mass_j,
            -vel_correction_coeff * dr[2] / mass_j,
        ];

        PairContribution { dq_i, dq_j, dp_i, dp_j }
    }
}
```

> **Note:** The force expression above is schematic. The exact analytical derivatives of the Weber potential must be carefully verified against the Julia symbolic output (or against the formulas in `theory/WeberElectrodynamics.md`) before use. The key structural point is that `dp_j = -dp_i` holds exactly.

### 3.4 The Full RHS Evaluation

```rust
fn evaluate_rhs<const D: usize>(
    q: &SoaVec<D>,           // positions for N particles
    p: &SoaVec<D>,           // momenta for N particles
    force: &impl PairwiseForce<D>,
    system: &System,
    dq_out: &mut SoaVec<D>,  // output: dq/dt
    dp_out: &mut SoaVec<D>,  // output: dp/dt
) {
    let n = system.n_particles;

    // Step 1: Kinetic contribution dq_i = p_i / m_i (O(N), per-particle)
    for i in 0..n {
        for d in 0..D {
            dq_out[i][d] = p[i][d] / system.masses[i];
        }
    }

    // Step 2: Clear force accumulator
    dp_out.fill_zero();

    // Step 3: Pair loop — N(N-1)/2 unique pairs
    for i in 0..n {
        for j in (i+1)..n {
            let kappa = system.kappas[pair_index(i, j, n)];
            let contrib = force.evaluate(
                &q[i], &q[j], &p[i], &p[j],
                system.masses[i], system.masses[j],
                system.charges[i], system.charges[j],
                kappa,
            );
            // Accumulate to both particles
            for d in 0..D {
                dq_out[i][d] += contrib.dq_i[d];
                dq_out[j][d] += contrib.dq_j[d];
                dp_out[i][d] += contrib.dp_i[d];
                dp_out[j][d] += contrib.dp_j[d];
            }
        }
    }
}
```

This performs exactly N*(N-1)/2 pair evaluations — half the Julia implementation's count.

## 4. Parallelism Strategy

### 4.1 The Data Race Problem

When pair (i,j) writes to accumulators for particles i and j, and pair (i,k) concurrently writes to particle i, we have a write-write conflict. Three strategies address this.

### 4.2 Option A: Atomic Accumulators

Each particle's force accumulator slot is an `AtomicU64` (reinterpreted as `f64` via `fetch_add` on the bit pattern). Pairs write via atomic addition with `Ordering::Relaxed`.

```rust
// Pseudocode
dp_atomic[i].fetch_add(contrib.dp_i, Relaxed);
dp_atomic[j].fetch_add(contrib.dp_j, Relaxed);
```

**Pros:** Simple, scales well for large N where contention is rare.
**Cons:** Atomic `f64` requires bit-level tricks (`f64::to_bits` / `from_bits` CAS loop or the `atomic_float` crate). Non-deterministic summation order means floating-point results vary across runs at the ULP level.

### 4.3 Option B: Per-Thread Local Accumulators (Recommended)

Each Rayon thread maintains a private accumulator array of size N*d. After the parallel pair loop, reduce by summing across threads.

```rust
use rayon::prelude::*;

let result = pairs
    .par_chunks(chunk_size)
    .fold(
        || Accumulator::zeros(n_particles, D),
        |mut acc, chunk| {
            for &(i, j) in chunk {
                let contrib = force.evaluate(/*...*/);
                acc.add_dq(i, &contrib.dq_i);
                acc.add_dq(j, &contrib.dq_j);
                acc.add_dp(i, &contrib.dp_i);
                acc.add_dp(j, &contrib.dp_j);
            }
            acc
        },
    )
    .reduce(
        || Accumulator::zeros(n_particles, D),
        |a, b| a.merge(&b),
    );
```

**Pros:** Fully safe Rust. Deterministic if chunk assignment is fixed. No atomic overhead.
**Cons:** O(T * N * d) extra memory for T threads. The reduction pass is O(T * N * d).

**Memory cost:** For N=1000, d=3, T=8 threads: 8 * 1000 * 3 * 8 bytes = 192 KB — fits comfortably in L2 cache.

### 4.4 Option C: Tiled Pair Decomposition

Partition the upper-triangular pair matrix into B*B tiles. Tiles that don't share row or column indices can execute concurrently without conflicts. This is a graph-coloring problem on the tile dependency structure.

**Pros:** No extra memory, no atomics, fully deterministic.
**Cons:** Complex scheduling logic. Tiles near the diagonal are smaller (load imbalance). Only worthwhile for very large N (> 10,000).

### 4.5 Recommendation

| N range | Strategy | Rationale |
|---------|----------|-----------|
| N < 20 | Sequential | Pair count too small for threading overhead |
| 20 ≤ N < 5000 | Option B (per-thread accumulators) | Best balance of simplicity, safety, and performance |
| N ≥ 5000 | Option A (atomics) or Option C (tiled) | Accumulator memory becomes significant; atomics have low contention at high N |

For the Weber electrodynamics use case (typically N = 2-100), **Option B is the clear winner**, with a fallback to sequential for very small N.

### 4.6 What Can Be Parallelized

The projection iteration itself is inherently sequential ($\mu^{(k+1)}$ depends on $\mu^{(k)}$). Within each iteration, the expensive work is the Strang splitting flow, which calls `evaluate_rhs` three times. Each `evaluate_rhs` call contains the parallelizable pair loop.

For higher-order compositions (triple jump with 9 stages, Suzuki with 15 stages), the stages are sequential but each stage's pair loop is independently parallelizable.

The overall parallelism structure is:

```
for each projection iteration (sequential, 2-3):
    for each Strang sub-step (sequential, 3):
        evaluate_rhs:
            kinetic loop (sequential, O(N))
            pair loop (PARALLEL, O(N^2/2))  <-- this is where threading lives
```

## 5. Memory Layout

### 5.1 Structure of Arrays (SoA)

The Julia implementation stores positions as a flat `Vector{Float64}` with interleaved coordinates: `[x1, y1, z1, x2, y2, z2, ...]`. This is an AoS (Array of Structs) layout that hinders SIMD because loading particle i's coordinates requires a stride-3 gather.

The Rust implementation should use SoA:

```rust
/// SoA storage for N particles in D dimensions.
/// Each coordinate axis is a contiguous array.
pub struct SoaVec<const D: usize> {
    /// data[d][i] = d-th coordinate of particle i
    pub data: [[f64; MAX_PARTICLES]; D],
    pub len: usize,
}
```

Or more flexibly with heap allocation:

```rust
pub struct ParticleCoords {
    pub x: Vec<f64>,  // [x_0, x_1, ..., x_{N-1}]
    pub y: Vec<f64>,  // [y_0, y_1, ..., y_{N-1}]
    pub z: Vec<f64>,  // [z_0, z_1, ..., z_{N-1}]
}
```

**SIMD benefit:** Computing distances between particles i and j loads `x[i], x[j]` from adjacent cache lines. With AVX-256, process 4 pairs simultaneously by loading `x[i0..i3]` and `x[j0..j3]` into SIMD registers.

### 5.2 Extended Phase Space Buffers

The extended state $Z = (Q, X, P, Y)$ from the theory document (SemiExplicitIntegrator.md) is stored as four separate `ParticleCoords` structs:

```rust
pub struct ExtendedState {
    pub q: ParticleCoords,  // Original position (Q)
    pub x: ParticleCoords,  // Auxiliary position (X)
    pub p: ParticleCoords,  // Original momentum (P)
    pub y: ParticleCoords,  // Auxiliary momentum (Y)
}
```

This replaces the Julia `Z::Vector{Float64}` of length 4d with structured access. Total memory: 4 * N * d * 8 bytes = 96N bytes for 3D.

### 5.3 Implicit A and A^T Operations

The Julia implementation constructs explicit dense matrices `A` (2d x 4d) and `A^T` (4d x 2d) and uses `mul!` for matrix-vector products (`solve.jl:72,87`). But the structure of A is trivially sparse — it's just identity-block arithmetic:

$$A \cdot Z = \begin{bmatrix} Q - X \\ P - Y \end{bmatrix}, \qquad A^T \cdot \mu = \begin{bmatrix} \mu_q \\ -\mu_q \\ \mu_p \\ -\mu_p \end{bmatrix}$$

In Rust, implement these as inline element-wise operations:

```rust
/// Compute f(mu) = A * Z_result, where Z_result = (Q, X, P, Y)
#[inline(always)]
fn apply_A(q: &[f64], x: &[f64], p: &[f64], y: &[f64], out: &mut [f64]) {
    let d = q.len();
    for i in 0..d {
        out[i] = q[i] - x[i];         // Q - X
        out[d + i] = p[i] - y[i];     // P - Y
    }
}

/// Compute A^T * mu and add to extended state
#[inline(always)]
fn apply_AT(mu: &[f64], d: usize, q: &mut [f64], x: &mut [f64],
            p: &mut [f64], y: &mut [f64]) {
    for i in 0..d {
        q[i] += mu[i];          // +mu_q
        x[i] -= mu[i];          // -mu_q
        p[i] += mu[d + i];      // +mu_p
        y[i] -= mu[d + i];      // -mu_p
    }
}
```

This eliminates the O(d^2) matrix storage and the O(d^2) matrix-vector multiply, replacing both with O(d) inline arithmetic. For N=100 in 3D (d=300), this saves a 600x1200 matrix allocation and a 600x1200 GEMV per projection iteration.

### 5.4 Zero-Allocation Hot Path

All buffers must be pre-allocated at `init` time. The `step` function must never call `Vec::push`, `Box::new`, or any allocating operation. The Rust ownership model makes this easy to enforce: all scratch space lives in the `IntegratorState` struct, and `step(&mut self)` borrows it mutably.

```rust
pub struct IntegratorState<const D: usize, F: PairwiseForce<D>> {
    // Problem definition (immutable after init)
    system: System,
    force: F,
    dt: f64,
    relaxation: f64,
    tol: f64,
    max_iter: usize,

    // Current state
    t: f64,
    q: ParticleCoords,
    p: ParticleCoords,

    // Pre-allocated scratch buffers
    extended: ExtendedState,          // Q, X, P, Y
    extended_hat: ExtendedState,      // Z_hat for Strang flow
    mu: Vec<f64>,                     // Lagrange multiplier (2d)
    mu_prev: Vec<f64>,
    f_mu: Vec<f64>,                   // Constraint residual (2d)
    dq_buf: ParticleCoords,           // RHS evaluation buffers
    dp_buf: ParticleCoords,
    dq_buf2: ParticleCoords,
    dp_buf2: ParticleCoords,

    // Thread-local accumulators (for parallel pair loop)
    accumulators: Vec<Accumulator>,   // One per Rayon thread

    // Solution history
    t_history: Vec<f64>,
    q_history: Vec<ParticleCoords>,
    p_history: Vec<ParticleCoords>,
}
```

## 6. Integrator Structure

### 6.1 Strang Splitting with Pair-Centric RHS

The Strang splitting flow from the theory document (`SemiExplicitIntegrator.md:55-59`) composes flows A and B:

$$\hat{\Phi}_{\Delta t} = \hat{\Phi}^A_{\Delta t/2} \circ \hat{\Phi}^B_{\Delta t} \circ \hat{\Phi}^A_{\Delta t/2}$$

In the extended phase space $(Q, X, P, Y)$:
- **Flow A** uses Hamiltonian $H(Q, Y)$: freezes Q and Y, evolves X and P
- **Flow B** uses Hamiltonian $H(X, P)$: freezes X and P, evolves Q and Y

```rust
fn strang_splitting_flow<const D: usize>(
    ext: &mut ExtendedState,
    dt: f64,
    force: &impl PairwiseForce<D>,
    system: &System,
    dq_buf: &mut ParticleCoords,
    dp_buf: &mut ParticleCoords,
    accumulators: &mut [Accumulator],
) {
    // Flow A(dt/2): evaluate RHS at (Q, Y), evolve X and P
    evaluate_rhs_parallel(&ext.q, &ext.y, force, system, dq_buf, dp_buf, accumulators);
    ext.x.axpy(dt / 2.0, dq_buf);
    ext.p.axpy(dt / 2.0, dp_buf);

    // Flow B(dt): evaluate RHS at (X, P), evolve Q and Y
    evaluate_rhs_parallel(&ext.x, &ext.p, force, system, dq_buf, dp_buf, accumulators);
    ext.q.axpy(dt, dq_buf);
    ext.y.axpy(dt, dp_buf);

    // Flow A(dt/2): evaluate RHS at (Q, Y), evolve X and P
    evaluate_rhs_parallel(&ext.q, &ext.y, force, system, dq_buf, dp_buf, accumulators);
    ext.x.axpy(dt / 2.0, dq_buf);
    ext.p.axpy(dt / 2.0, dp_buf);
}
```

Each `evaluate_rhs_parallel` call does the N(N-1)/2 pair loop with optional Rayon parallelism (gated on `n_particles > PARALLEL_THRESHOLD`).

### 6.2 Symmetric Projection Iteration

The fixed-point iteration solves for the Lagrange multiplier $\mu \in \mathbb{R}^{2d}$ such that the Strang-evolved state lands on the constraint manifold $\mathcal{N}$ where $Q = X$ and $P = Y$.

```rust
fn projected_step<const D: usize>(
    state: &mut IntegratorState<D, impl PairwiseForce<D>>,
) -> Result<(), ConvergenceError> {
    let d = state.system.dof();

    // 1. Embed: Z = (q, q, p, p)
    state.extended.embed(&state.q, &state.p);

    // 2. Fixed-point iteration
    state.mu.fill(0.0);  // Or warm-start from previous step

    for _ in 0..state.max_iter {
        // Copy mu for convergence check
        state.mu_prev.copy_from_slice(&state.mu);

        // Compute Z_hat = Z + A^T mu
        state.extended_hat.copy_from(&state.extended);
        apply_AT(&state.mu, d, &mut state.extended_hat);

        // Apply Strang splitting flow
        strang_splitting_flow(
            &mut state.extended_hat, state.dt, &state.force,
            &state.system, &mut state.dq_buf, &mut state.dp_buf,
            &mut state.accumulators,
        );

        // Compute f(mu) = A * (Z_hat + A^T mu)
        //   Z_result = Z_hat + A^T mu  (shift endpoint)
        apply_AT(&state.mu, d, &mut state.extended_hat);
        //   f_mu = A * Z_result = [Q-X, P-Y]
        apply_A(&state.extended_hat, &mut state.f_mu);

        // Relaxed Newton update: mu -= 0.25 * f(mu)
        for i in 0..(2 * d) {
            state.mu[i] -= state.relaxation * state.f_mu[i];
        }

        // Check convergence
        let step_norm = norm_diff(&state.mu, &state.mu_prev);
        if step_norm < state.tol {
            let residual_norm = norm(&state.f_mu);
            if residual_norm < state.tol {
                break;
            }
        }
    }

    // 3. Extract result: q = Q, p = P from the converged Z_result
    state.q.copy_from(&state.extended_hat.q);
    state.p.copy_from(&state.extended_hat.p);

    Ok(())
}
```

**Warm-starting:** The converged $\mu$ from the previous timestep can be reused as the initial guess for the next timestep. This typically reduces iterations from 2-3 to 1-2 (see `SemiExplicitIntegrator.md:140`).

### 6.3 Higher-Order Compositions

Higher-order methods are sequences of 2nd-order Strang splitting calls with pre-computed gamma coefficients:

```rust
/// Pre-computed composition coefficients.
pub struct CompositionCoeffs {
    pub gammas: Vec<f64>,  // Sub-step fractions
    pub order: usize,
}

impl CompositionCoeffs {
    /// Triple jump: 3 stages, order n from order n-2
    pub fn triple_jump(order: usize) -> Self {
        let gamma1 = 1.0 / (2.0 - 2.0_f64.powf(1.0 / (order as f64 - 1.0)));
        let gamma2 = -2.0_f64.powf(1.0 / (order as f64 - 1.0))
                    / (2.0 - 2.0_f64.powf(1.0 / (order as f64 - 1.0)));
        CompositionCoeffs {
            gammas: vec![gamma1, gamma2, gamma1],
            order,
        }
    }

    /// Suzuki fractal: 5 stages
    pub fn suzuki(order: usize) -> Self {
        let gamma_outer = 1.0 / (4.0 - 4.0_f64.powf(1.0 / (order as f64 - 1.0)));
        let gamma_center = -4.0_f64.powf(1.0 / (order as f64 - 1.0))
                         / (4.0 - 4.0_f64.powf(1.0 / (order as f64 - 1.0)));
        CompositionCoeffs {
            gammas: vec![gamma_outer, gamma_outer, gamma_center,
                         gamma_outer, gamma_outer],
            order,
        }
    }
}
```

The composed integrator replaces `strang_splitting_flow` with a loop over gamma-scaled sub-steps. The projection iteration wraps around the full composed flow.

## 7. Performance Analysis

### 7.1 Operation Count Per RHS Evaluation

| Operation | Julia (particle-centric) | Rust (pair-centric) | Savings |
|-----------|------------------------|--------------------|---------|
| Pair force evaluations | N*(N-1) | N*(N-1)/2 | 2x |
| `sqrt` calls | N*(N-1) | N*(N-1)/2 | 2x |
| Division by r | N*(N-1) | N*(N-1)/2 | 2x |
| Kinetic dq/dt | N (same) | N (same) | 1x |
| A/A^T matrix ops | O(d^2) GEMV | O(d) inline | ~d/4 x |

### 7.2 Total Pair Evaluations Per Timestep

For a 2nd-order integrator with $k$ projection iterations:

| | Julia | Rust | Speedup |
|-|-------|------|---------|
| RHS calls per Strang step | 6 | 6 | 1x |
| Strang steps per timestep | k+1 | k+1 | 1x |
| Pairs per RHS call | N(N-1) | N(N-1)/2 | 2x |
| **Total pair evals** | **6(k+1)N(N-1)** | **3(k+1)N(N-1)** | **2x** |

With k=2 (typical), N=50 (25 particles in 2D):
- Julia: 6 * 3 * 2450 = 44,100 pair evaluations per timestep
- Rust: 3 * 3 * 2450 = 22,050 pair evaluations per timestep

### 7.3 Threading Speedup

The pair loop is embarrassingly parallel. With Rayon's work-stealing scheduler:

| N (particles) | Pairs | 1 thread | 4 threads | 8 threads | 16 threads |
|---------------|-------|----------|-----------|-----------|------------|
| 10 | 45 | 1x | ~1x (overhead dominates) | ~1x | ~1x |
| 50 | 1,225 | 1x | ~3.5x | ~6x | ~8x |
| 100 | 4,950 | 1x | ~3.8x | ~7x | ~12x |
| 500 | 124,750 | 1x | ~3.9x | ~7.5x | ~14x |
| 1000 | 499,500 | 1x | ~4.0x | ~7.8x | ~15x |

(Estimates assume ~30 FLOPs per pair and include reduction overhead.)

### 7.4 SIMD Potential

The pair force inner loop performs:
- 3 subtractions (dr), 3 multiplies + 2 adds (r^2), 1 sqrt, 1 reciprocal
- ~15 more FLOPs for the Weber-specific terms
- Total: ~25-30 FLOPs per pair

With AVX-256 processing 4 `f64` values simultaneously:
- **Horizontal SIMD** (one pair, vectorize across dimensions): limited benefit for d=3 (only 3 lanes used)
- **Vertical SIMD** (4 pairs simultaneously, same dimension): 4x throughput on the arithmetic, but requires gather loads for non-contiguous particle data

SoA layout enables vertical SIMD: `x[i0..i3]` and `x[j0..j3]` are contiguous if pairs are sorted by particle index. Expected SIMD speedup: **1.5-2.5x** for the pair loop body.

### 7.5 Projected Composite Speedup (N=100, 8 cores)

| Factor | Multiplier |
|--------|------------|
| Newton's 3rd law | 2.0x |
| Multi-threading (8 cores) | ~7.0x |
| SIMD (conservative) | ~1.5x |
| No GC + cache optimization | ~1.3x |
| Implicit A/A^T (vs dense GEMV) | ~1.1x |
| **Composite** | **~30x** |

## 8. Rust-Specific Considerations

### 8.1 Inlining and Bounds-Check Elision

The pair force evaluation is called N(N-1)/2 times per RHS — it must be fully inlined:

```rust
#[inline(always)]
fn evaluate(&self, ...) -> PairContribution<D> { ... }
```

Bounds checking on SoA arrays adds branch overhead in the inner loop. Two safe approaches:

1. **Iterator-based access:** Use `.iter().zip()` patterns that the compiler can prove are in-bounds
2. **Assert once, index freely:** Place a single `assert!(n <= data.len())` before the loop, then use `unsafe { *data.get_unchecked(i) }` inside

The second approach is more natural for the pair loop where both `i` and `j` index into the same arrays.

### 8.2 Portable SIMD (Nightly)

Rust's `std::simd` (unstable as of early 2026) provides portable SIMD types:

```rust
#![feature(portable_simd)]
use std::simd::f64x4;

// Process 4 pairs simultaneously
let dx = f64x4::from_array([
    x[pairs[0].0] - x[pairs[0].1],
    x[pairs[1].0] - x[pairs[1].1],
    x[pairs[2].0] - x[pairs[2].1],
    x[pairs[3].0] - x[pairs[3].1],
]);
let r2 = dx * dx + dy * dy + dz * dz;
let r = r2.sqrt();  // 4 sqrts in one instruction
```

For stable Rust, use the `packed_simd2` or `wide` crates, or rely on LLVM auto-vectorization with `#[target_feature(enable = "avx2")]`.

### 8.3 Borrow Checker Patterns for the Integrator

The Strang splitting flow needs mutable access to the extended state while simultaneously reading system parameters. Rust's borrow checker requires careful structuring:

```rust
// Problem: strang_splitting_flow needs &mut ext AND &system,
// but both live in IntegratorState.
// Solution: destructure the state before calling.

fn step(state: &mut IntegratorState) {
    let IntegratorState {
        ref system, ref force, ref mut extended, ref mut dq_buf, ref mut dp_buf,
        ref mut accumulators, ..
    } = *state;

    strang_splitting_flow(extended, dt, force, system, dq_buf, dp_buf, accumulators);
}
```

This pattern (destructuring borrows) is idiomatic Rust for complex mutable state machines.

### 8.4 FFI for Julia Interop

A `#[no_mangle] extern "C"` API allows the Rust integrator to be called from Julia via `ccall`, enabling incremental adoption:

```rust
#[no_mangle]
pub extern "C" fn weber_step(
    state: *mut WeberIntegratorOpaque,
) -> i32 {
    let state = unsafe { &mut *state };
    match state.step() {
        Ok(true) => 1,   // continue
        Ok(false) => 0,  // done
        Err(_) => -1,    // convergence failure
    }
}
```

Julia side:
```julia
retcode = ccall((:weber_step, libweber), Cint, (Ptr{Cvoid},), state_ptr)
```

This allows keeping the Julia front-end (symbolic construction, Plots.jl, Makie animation) while offloading the compute kernel to Rust. The params vector layout is already identical (`[m1..mN, q1..qN, c, kappa_12..]`).

### 8.5 Testing Against the Julia Reference

Port the key physics tests from the Julia test suite:

1. **Energy conservation:** Run N-body Weber problem, verify $|H(t) - H(0)| / |H(0)| < \varepsilon$ with $\varepsilon$ matching Julia's tolerance
2. **Newton's third law:** For each pair evaluation, assert `dp_i + dp_j == 0` to machine precision
3. **Momentum conservation:** Verify $\sum_i p_i(t) = \sum_i p_i(0)$ (Weber conserves linear momentum)
4. **Trajectory matching:** Run identical ICs in Julia and Rust, compare $(q, p)$ at each timestep to within roundoff

### 8.6 Benchmarking

Use [Criterion.rs](https://github.com/bheisler/criterion.rs) for statistical benchmarks:

- **Micro:** Single pair force evaluation latency (~5-10 ns expected)
- **Pair loop:** N(N-1)/2 force accumulation for N = 10, 50, 100, 500
- **Timestep:** Full `projected_step` including projection iterations
- **Scaling:** Vary thread count from 1 to `num_cpus` for the pair loop
- **Comparison:** Julia `@benchmark` on the same problems for apples-to-apples comparison

## 9. Summary of Recommendations

1. **Start with the pair-centric `PairwiseForce` trait** — this is the most impactful architectural decision, halving computation regardless of threading
2. **Use Rayon fold/reduce** for the parallel pair loop — safe, deterministic, good scaling
3. **SoA memory layout** for particle data — enables SIMD and improves cache utilization
4. **Implicit A/A^T** as inline arithmetic — eliminates unnecessary matrix storage and GEMV
5. **Pre-allocate everything** at init time — zero allocations in the hot path
6. **Warm-start mu** across timesteps — reduces projection iterations
7. **FFI bridge to Julia** for incremental adoption — keep the front-end, replace the kernel
8. **Benchmark against Julia** at every stage — quantify actual vs projected speedups

See also: [SemiExplicitIntegrator.md](../theory/SemiExplicitIntegrator.md), [WeberElectrodynamics.md](../theory/WeberElectrodynamics.md).
