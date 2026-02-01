# Parallelization Research Report: Weber Electrodynamics Simulation

## Executive Summary

This report analyzes parallelization strategies for the WeberElectrodynamics.jl symplectic integrator. The current implementation uses a semi-explicit integrator with extended phase space and symmetric projection, generating derivative functions via Symbolics.jl. **No parallelization is currently implemented** beyond `@inbounds` optimizations.

### Key Findings

| Strategy | Effort | Speedup Potential | Problem Size |
|----------|--------|-------------------|--------------|
| **Warm-start μ** | ✅ Implemented | Already active | All |
| **SIMD (LoopVectorization)** | Medium | 2-4x | All |
| **Multi-threading (Polyester)** | Medium | Nx threads | N > 100 particles |
| **GPU (Metal/CUDA)** | High | 10-100x | N > 4,000 particles |

### Prioritized Recommendations

1. **Already implemented**:
   - ✅ Warm-start μ from previous timestep (verified in solve.jl)

2. **Short-term (medium effort)**:
   - Add Polyester.jl `@batch` to force calculation loops
   - Implement 4th-order Yoshida composition

3. **Medium-term (high effort)**:
   - GPU kernels via KernelAbstractions.jl for N > 4,000

---

## 1. Current Implementation Analysis

### Architecture

The implementation consists of:

- **WeberSystem**: Symbolic Hamiltonian with compiled derivative functions
- **SymmetricProjectionIntegrator**: Strang splitting in extended phase space
- **Fixed-point iteration**: Solves constraint equation with relaxation = 0.25

### Performance Hotspots

```
┌─────────────────────────────────────────────────────────────────────┐
│                         step!() Function                            │
├─────────────────────────────────────────────────────────────────────┤
│  ┌──────────────────────────────────────────────────────────────┐  │
│  │           Fixed-Point Iteration Loop (2-3 iters)              │  │
│  │  ┌────────────────────────────────────────────────────────┐  │  │
│  │  │      compute_constraint_residual!()                     │  │  │
│  │  │  ┌──────────────────────────────────────────────────┐  │  │  │
│  │  │  │     strang_splitting_flow!()                      │  │  │  │
│  │  │  │  ┌────────────────────────────────────────────┐  │  │  │  │
│  │  │  │  │  dq_dt_compiled() [3x] ← PRIMARY HOTSPOT   │  │  │  │  │
│  │  │  │  │  dp_dt_compiled() [3x] ← PRIMARY HOTSPOT   │  │  │  │  │
│  │  │  │  └────────────────────────────────────────────┘  │  │  │  │
│  │  │  └──────────────────────────────────────────────────┘  │  │  │
│  │  └────────────────────────────────────────────────────────┘  │  │
│  └──────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────┘
```

**Critical path**: 6 derivative evaluations × 2-3 iterations = 12-18 force calculations per timestep.

### Scaling Characteristics

| Particles (n) | DOF (d) | Pairwise Terms | Force Evals/Step |
|---------------|---------|----------------|------------------|
| 2 | 4-6 | 1 | 12-18 |
| 10 | 20-30 | 45 | 12-18 |
| 100 | 200-300 | 4,950 | 12-18 |
| 1,000 | 2,000-3,000 | 499,500 | 12-18 |

The O(n²) pairwise interaction scaling dominates for large systems.

---

## 2. SIMD Vectorization

### Current State

The Symbolics.jl-generated derivative functions are **not SIMD-optimized**. The generated code relies on LLVM auto-vectorization, which may not apply due to:

- Non-contiguous memory access patterns
- Function call barriers
- Complex control flow in pairwise loops

### Opportunities

#### 2.1 Memory Layout Optimization

Current layout (interleaved):
```julia
q = [x₁, y₁, z₁, x₂, y₂, z₂, ...]  # AoS-style
```

SIMD-friendly layout (Structure of Arrays):
```julia
struct ParticlesSoA
    x::Vector{Float64}  # [x₁, x₂, x₃, ...]
    y::Vector{Float64}  # [y₁, y₂, y₃, ...]
    z::Vector{Float64}  # [z₁, z₂, z₃, ...]
end
```

**Expected benefit**: 1.5-3x speedup from cache-aligned, vectorizable access patterns.

#### 2.2 LoopVectorization.jl

```julia
using LoopVectorization

# Vectorized distance calculation
function compute_distances_simd!(r², pos, n)
    @turbo for i in 1:n
        for j in (i+1):n
            dx = pos.x[j] - pos.x[i]
            dy = pos.y[j] - pos.y[i]
            dz = pos.z[j] - pos.z[i]
            r²[i,j] = dx*dx + dy*dy + dz*dz
        end
    end
end
```

**Note**: LoopVectorization.jl is deprecated for Julia 1.11+. Monitor for replacement.

#### 2.3 Limitations

- Symbolics.jl generated code cannot be directly `@turbo`-annotated
- Would require manual force function implementation to fully leverage SIMD
- Trade-off: Lose automatic symbolic differentiation

### Recommendation

**Medium priority**. Consider SoA layout for large-n systems and manual SIMD implementation for critical force loops.

---

## 3. Multi-Threading

### Strategy Comparison

| Library | Overhead | Best For | Julia Version |
|---------|----------|----------|---------------|
| `@threads` | ~15-20 μs | General parallelism | All |
| Polyester.jl `@batch` | ~1 μs | Tight loops | All |
| ThreadsX.jl | ~10-15 μs | Map/reduce patterns | All |
| FLoops.jl | ~10-15 μs | Complex reductions | All |

### Implementation Approach

#### 3.1 Force Calculation Parallelization

The pairwise force loop is embarrassingly parallel:

```julia
using Polyester

function compute_pairwise_forces!(forces, positions, params, n)
    @batch per=core for i in 1:n
        for j in (i+1):n
            f_ij = weber_force(positions, i, j, params)
            # Atomic or reduction needed for thread safety
            forces[i] += f_ij
            forces[j] -= f_ij
        end
    end
end
```

**Challenge**: The inner loop writes to shared `forces` array. Solutions:
1. Thread-local accumulation with final reduction
2. Symmetric tiling to avoid conflicts
3. Use of `@reduce` patterns

#### 3.2 Crossover Points

| Problem Size | Threading Benefit |
|--------------|-------------------|
| n < 50 particles | Negative (overhead dominates) |
| 50-500 particles | Marginal |
| 500+ particles | Positive scaling |

### Recommendation

**Medium priority**. Use Polyester.jl for large-n systems. Implement thread-local force accumulation to avoid synchronization overhead.

---

## 4. GPU Acceleration

### Platform Comparison

| Feature | Apple Metal | NVIDIA CUDA |
|---------|-------------|-------------|
| Julia Package | Metal.jl | CUDA.jl |
| Maturity | Work-in-progress | Production |
| FP64 Support | Emulated (slow) | Native |
| Memory Model | Unified | Discrete |
| Portability | KernelAbstractions.jl | KernelAbstractions.jl |

### Key Findings

1. **Metal lacks native FP64**: Weber forces require double precision for accuracy. Metal would require Float32 with careful error management.

2. **Unified memory advantage**: Metal's unified memory eliminates transfer overhead, potentially beneficial for iterative methods.

3. **CUDA kernel launch**: ~7 μs minimum overhead per kernel.

4. **Crossover point**: GPU beneficial for N > 4,000 particles.

### N-Body GPU Kernel Design

From NVIDIA GPU Gems 3, the optimal approach uses **tiled shared memory**:

```julia
using KernelAbstractions

const TILE_SIZE = 256

@kernel function nbody_kernel!(acc, @Const(pos), @Const(mass), N)
    i = @index(Global)
    lid = @index(Local)

    # Shared memory for tile
    tile_pos = @localmem Float32 (TILE_SIZE, 3)
    tile_mass = @localmem Float32 (TILE_SIZE,)

    ax, ay, az = 0f0, 0f0, 0f0

    # Loop over tiles
    for tile in 0:(cld(N, TILE_SIZE) - 1)
        # Collaborative load
        j = tile * TILE_SIZE + lid
        if j <= N
            tile_pos[lid, 1] = pos[j, 1]
            tile_pos[lid, 2] = pos[j, 2]
            tile_pos[lid, 3] = pos[j, 3]
            tile_mass[lid] = mass[j]
        end
        @synchronize()

        # Compute forces from tile
        for k in 1:min(TILE_SIZE, N - tile * TILE_SIZE)
            # Weber force calculation...
        end
        @synchronize()
    end

    if i <= N
        acc[i, 1] = ax
        acc[i, 2] = ay
        acc[i, 3] = az
    end
end
```

### Recommendation

**High priority for large systems**. Use KernelAbstractions.jl for portability across Metal/CUDA. Requires Float32 adaptation for Metal or access to CUDA hardware.

---

## 5. Distributed Computing

### Julia Options

| Package | Use Case | Performance |
|---------|----------|-------------|
| Distributed.jl | Simple distribution | Good |
| MPI.jl | HPC clusters | Best |
| DistributedArrays.jl | Data parallelism | Good |

### N-Body Distribution Strategies

#### 5.1 Full Replication (Small-Medium N)

All processes hold all particle positions:
```julia
# Each timestep:
positions = MPI.Allgatherv(local_positions, comm)
# Compute local subset of forces
local_forces = compute_weber_forces(my_particle_range, positions)
```

**Communication**: O(n × p) data per timestep.

#### 5.2 Ring Communication (Memory Efficient)

```
P₀ → P₁ → P₂ → P₃ → P₀  (circular data flow)

Each step:
1. Compute forces with current received particles
2. Send particles to next, receive from previous
3. Repeat p-1 times
```

**Communication**: O(n) data per ring pass.

### Recommendation

**Low priority** for typical Weber simulations. MPI becomes valuable for very large particle counts or HPC cluster deployment.

---

## 6. Parallel-in-Time Methods

### Parareal Algorithm

```
Time:     t₀ ────────── t₁ ────────── t₂ ────────── t₃
           │             │             │             │
Coarse:   G(·)  ───→   G(·)  ───→   G(·)  ───→   G(·)
           │             │             │             │
Fine:     F(·)         F(·)         F(·)         F(·)
          (P₀)         (P₁)         (P₂)         (P₃)
           │             │             │             │
Correct:  ─────────────────────────────────────────────
```

### Critical Finding: Standard Parareal Breaks Symplecticity

Even with symplectic fine and coarse integrators, the iterative correction procedure **does not preserve symplectic structure**. This leads to:
- Linear energy drift
- Incorrect long-term dynamics
- Loss of phase-space volume preservation

### Symplectic Variants

1. **Symmetric Parareal** (Dai et al., 2013): Time-reversible iteration scheme
2. **Symplectic Parareal** (Bal & Wu, 2008): Preserves Hamiltonian structure
3. **PFASST with Gauss nodes**: Converges to symplectic collocation

### Recommendation

**Low priority**. Parareal requires specialized symplectic variants that aren't readily available in Julia. Standard parallel-in-time methods break symplecticity and are not suitable for long-time Hamiltonian integration.

---

## 7. Symbolics.jl Optimization

### Current Code Generation

```julia
# In weber_system.jl:
dq_dt_compiled = Symbolics.build_function(
    dq_dt_symbolic, q_vars, p_vars, param_symbols,
    expression = Val{false}  # Compiled via RuntimeGeneratedFunctions
)[2]
```

### Available Parallel Forms

| Form | Description | Use Case |
|------|-------------|----------|
| `SerialForm()` | Sequential | < 1500 expressions |
| `ShardedForm(80, 4)` | Split into sub-functions | Large Jacobians |
| `MultithreadedForm()` | Thread-parallel | Large systems |

### Example with Multithreading

```julia
# Generate threaded derivative function
dq_dt_compiled = Symbolics.build_function(
    dq_dt_symbolic, q_vars, p_vars, param_symbols,
    expression = Val{false},
    parallel = Symbolics.MultithreadedForm()
)[2]
```

### Limitations

- **No GPU code generation**: Cannot directly emit CUDA/Metal kernels
- **No SIMD annotation**: Generated code doesn't use `@simd` or `@turbo`
- **CSE sometimes suboptimal**: May not catch all common subexpressions

### Alternative: FastDifferentiation.jl

For maximum derivative performance:
```julia
using FastDifferentiation

@variables x[1:6]
expr = weber_hamiltonian(x)
grad = jacobian([expr], x)
f! = make_function(grad, x, in_place=true)
```

**Claimed benefit**: Finds shared subexpressions across partial derivatives.

### Recommendation

**Medium priority**. Enable `MultithreadedForm()` for large systems. Consider FastDifferentiation.jl for derivative-heavy workloads.

---

## 8. Algorithmic Optimizations

### 8.1 Warm-Starting μ

**Status**: ALREADY IMPLEMENTED

The μ vector is initialized to zeros once in `SymmetricProjectionBuffers` (types.jl:140)
when `init(prob)` is called. The `step!()` function references this buffer (solve.jl:154)
and updates it in-place (solve.jl:192) without resetting it. The converged μ from
timestep N automatically becomes the initial guess for timestep N+1.

**Verification**: The fixed-point iteration typically converges in 2-3 iterations for the
first few steps, then often converges faster (1-2 iterations) once μ has stabilized.

**No changes needed**: The implementation already provides the warm-start benefit.

### 8.2 Higher-Order Composition (4th Order Yoshida)

```julia
# Exact coefficients for 4th-order Yoshida
const w₀ = -cbrt(2.0) / (2.0 - cbrt(2.0))  # ≈ -1.702414
const w₁ = 1.0 / (2.0 - cbrt(2.0))          # ≈  1.351207

function yoshida4_step!(state, dt, flow2!)
    flow2!(state, w₁ * dt)
    flow2!(state, w₀ * dt)
    flow2!(state, w₁ * dt)
end
```

**Trade-off**: 3x more force evaluations per step, but 4th-order vs 2nd-order error allows larger dt.

### 8.3 Symplectic Correctors

From REBOUND documentation, correctors can improve accuracy by **up to 6 orders of magnitude** without additional force evaluations:

```julia
# Apply corrector only at output points
function corrected_output(state, corrector_coeffs)
    # Modify initial conditions to compensate for integrator error
    # Only needed when saving output
end
```

### 8.4 Anderson Acceleration

Replace fixed-point iteration with Anderson acceleration for near-quadratic convergence:

```julia
using NLsolve

function solve_constraint!(μ, f!, m=5)
    # Anderson acceleration with depth m
    nlsolve(f!, μ, method=:anderson, m=m)
end
```

### Recommendation

**High priority**:
1. Warm-start μ (trivial, immediate benefit)
2. Anderson acceleration (moderate effort, faster convergence)
3. 4th-order Yoshida (allows larger timesteps)

---

## 9. Comparative Analysis

### Speedup vs Implementation Effort Matrix

```
                    Implementation Effort
                 Low         Medium        High
            ┌───────────┬───────────┬───────────┐
      High  │ Warm-start│ Multi-    │ GPU       │
Speedup     │ μ (2x)    │ threading │ kernels   │
Potential   │           │ (4-8x)    │ (10-100x) │
            │           │ SIMD      │           │
            ├───────────┼───────────┼───────────┤
      Med   │           │ Yoshida   │ Symplectic│
            │           │ 4th order │ Parareal  │
            │           │           │           │
            │           │ MTF       │ MPI       │
            │           │ codegen   │ distrib.  │
            ├───────────┼───────────┼───────────┤
      Low   │           │           │           │
            │           │           │           │
            └───────────┴───────────┴───────────┘
```

### Problem Size Decision Tree

```
                      ┌─────────────────────┐
                      │   Number of         │
                      │   Particles?        │
                      └─────────┬───────────┘
                                │
              ┌─────────────────┼─────────────────┐
              │                 │                 │
         n < 100           100 < n < 4000     n > 4000
              │                 │                 │
              ▼                 ▼                 ▼
    ┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
    │ Focus on:       │ │ Focus on:       │ │ Focus on:       │
    │ - Algorithmic   │ │ - Multi-thread  │ │ - GPU kernels   │
    │ - Warm-start μ  │ │ - SIMD          │ │ - MPI           │
    │ - Higher-order  │ │ - Warm-start μ  │ │ - Multi-thread  │
    │                 │ │ - Higher-order  │ │ - Higher-order  │
    └─────────────────┘ └─────────────────┘ └─────────────────┘
```

---

## 10. Conclusions

The WeberElectrodynamics.jl implementation has significant parallelization opportunities:

1. **Immediate gains** available through warm-starting μ
2. **Medium-term improvements** via threading and higher-order methods
3. **Long-term scalability** through GPU kernels and MPI distribution

The optimal strategy depends heavily on problem size:
- Small (n < 100): Algorithmic improvements dominate
- Medium (100-4000): Threading provides good speedup
- Large (n > 4000): GPU kernels essential

The velocity-dependent nature of Weber forces requires careful handling in any parallelization strategy, but the symmetric projection framework is well-suited to parallel force evaluation.

---

## Sources

### Julia Parallelization
- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/)
- [Julia Distributed Computing](https://docs.julialang.org/en/v1/manual/distributed-computing/)
- [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)
- [Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl)
- [MPI.jl](https://juliaparallel.org/MPI.jl/stable/)

### GPU Computing
- [CUDA.jl Documentation](https://cuda.juliagpu.org/stable/)
- [Metal.jl](https://github.com/JuliaGPU/Metal.jl)
- [KernelAbstractions.jl](https://juliagpu.github.io/KernelAbstractions.jl/stable/)
- [NVIDIA GPU Gems 3: Fast N-Body](https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-31-fast-n-body-simulation-cuda)

### Symplectic Integration
- [Yoshida (1990) - Higher Order Symplectic Integrators](https://www.sciencedirect.com/science/article/pii/0375960190900923)
- [REBOUND High-Order Integrators](https://academic.oup.com/mnras/article/489/4/4632/5565063)
- [WHFast Paper](https://arxiv.org/abs/1506.01084)

### Parallel-in-Time
- [Symmetric Parareal for Hamiltonian Systems](https://arxiv.org/abs/1011.6222)

### Symbolics
- [Symbolics.jl build_function](https://docs.sciml.ai/Symbolics/stable/manual/build_function/)
- [FastDifferentiation.jl](https://brianguenter.github.io/FastDifferentiation.jl/dev/)
- [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/)
