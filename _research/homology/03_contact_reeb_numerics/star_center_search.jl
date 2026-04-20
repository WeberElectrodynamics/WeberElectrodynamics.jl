"""
Star-center optimizer for the Weber Hill region {V ≤ E}.

Given a Weber system, energy E, and candidate star-center q*, shoots N random
rays and checks:
  (a) monotonicity of V along ray (V non-decreasing outward)
  (b) single-crossing of {V = E}

Then optimizes over q* to maximize the fraction of monotone rays using a
coordinate-descent / random-perturbation scheme.

Tested on: 2-body (2D), 4-body alternating square, 4-body rhombus.
"""

using Random, LinearAlgebra, Printf
using WeberElectrodynamics

# ============================================================================
# Core ray-shooting analysis
# ============================================================================

"""
    shoot_rays(V, q_star, E, n_rays, t_max, n_samples; seed=nothing)

Shoot `n_rays` random rays from `q_star` and return:
  - frac_monotone: fraction of rays where V is non-decreasing
  - frac_single_cross: fraction of rays crossing {V=E} at most once
  - details: per-ray (monotone::Bool, n_crossings::Int)
"""
function shoot_rays(V, q_star::Vector{Float64}, E::Float64, n_rays::Int;
                    t_max::Float64=50.0, n_samples::Int=400, seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end
    d = length(q_star)
    n_monotone = 0
    n_single = 0
    details = Tuple{Bool,Int}[]

    for _ in 1:n_rays
        v = randn(d)
        v ./= norm(v)
        ts = range(1e-4, t_max; length=n_samples)
        Vs = [V(q_star .+ t .* v) for t in ts]

        # Monotonicity check
        diffs = diff(Vs)
        is_monotone = all(d -> d >= -1e-10, diffs)

        # Crossing count
        crossings = 0
        for i in 1:(length(Vs)-1)
            if (Vs[i] - E) * (Vs[i+1] - E) < 0
                crossings += 1
            end
        end

        if is_monotone
            n_monotone += 1
        end
        if crossings <= 1
            n_single += 1
        end
        push!(details, (is_monotone, crossings))
    end

    return (frac_monotone = n_monotone / n_rays,
            frac_single_cross = n_single / n_rays,
            details = details)
end

"""
    optimize_star_center(V, E, q_init, n_rays; kwargs...)

Random-perturbation coordinate descent to find q* maximizing monotone-ray fraction.
Returns (best_q, best_frac_monotone, best_frac_single_cross).
"""
function optimize_star_center(V, E::Float64, q_init::Vector{Float64}, n_rays::Int;
                              n_iters::Int=200, step_size::Float64=0.5,
                              t_max::Float64=50.0, n_samples::Int=200,
                              cooling::Float64=0.995, seed::Int=42)
    Random.seed!(seed)
    d = length(q_init)
    q_best = copy(q_init)
    res = shoot_rays(V, q_best, E, n_rays; t_max, n_samples, seed=seed+1)
    best_score = res.frac_monotone + 0.5 * res.frac_single_cross
    best_mono = res.frac_monotone
    best_single = res.frac_single_cross

    σ = step_size
    for iter in 1:n_iters
        q_trial = q_best .+ σ .* randn(d)
        # Ensure trial point is inside the Hill region
        if V(q_trial) > E
            σ *= cooling
            continue
        end
        res_trial = shoot_rays(V, q_trial, E, n_rays; t_max, n_samples,
                               seed=seed+iter+1)
        score = res_trial.frac_monotone + 0.5 * res_trial.frac_single_cross
        if score > best_score
            q_best = q_trial
            best_score = score
            best_mono = res_trial.frac_monotone
            best_single = res_trial.frac_single_cross
        end
        σ *= cooling
    end

    return (q_star = q_best, frac_monotone = best_mono,
            frac_single_cross = best_single)
end

# ============================================================================
# 2-body potential helper (for standalone 2-body tests without 4-body IC)
# ============================================================================

function make_2body_V(charges, masses, c)
    sys = HamiltonianSystem(2, 2)
    kappas = [1.0]
    params = vcat(masses, charges, [c], kappas)
    zero_p = zeros(4)
    return q -> sys.hamiltonian_compiled(q, zero_p, params)
end

# ============================================================================
# Run tests
# ============================================================================

function run_2body_test()
    println("="^60)
    println("2-BODY UNLIKE CHARGES (q1=+1, q2=-1) in 2D")
    println("="^60)

    V = make_2body_V([1.0, -1.0], [1.0, 1.0], 1.0)

    # Place particles at r=1 apart along x in COM frame
    q_eq = [0.5, 0.0, -0.5, 0.0]  # separation = 1.0
    E_eq = V(q_eq)
    @printf "V at r=1 separation: %.6f\n" E_eq

    # Test several energies
    for E in [-5.0, -2.0, -1.0, -0.5, -0.1]
        # Hill region for unlike charges: V(r) = -1/r <= E (E<0)
        # => -1/r <= E => 1/r >= |E| => r <= 1/|E|.
        # Hill region is the COMPACT set (0, 1/|E|]. V → -∞ as r→0.
        # Place center at r = 1/(2|E|) where V = -2|E| = 2E < E (deep inside).
        r_boundary = 1.0 / abs(E)
        r_center = r_boundary / 2.0   # V = -2|E| < E (inside Hill region)
        q_star = [r_center/2, 0.0, -r_center/2, 0.0]

        # Verify q_star is inside Hill region
        V_star = V(q_star)
        if V_star > E
            @printf "E=%.2f: q* outside Hill region (V(q*)=%.4f > E), skipping\n" E V_star
            continue
        end

        res = shoot_rays(V, q_star, E, 200; t_max=max(30.0, 5*r_center), n_samples=400, seed=42)
        @printf "E=%.2f: monotone=%.1f%%, single-cross=%.1f%% (r_center=%.3f, V*=%.4f)\n" E 100*res.frac_monotone 100*res.frac_single_cross r_center V_star

        # Try optimization
        opt = optimize_star_center(V, E, q_star, 100; n_iters=150, step_size=r_center*0.3,
                                   t_max=max(30.0, 5*r_center), n_samples=200)
        @printf "  optimized: monotone=%.1f%%, single-cross=%.1f%%\n" 100*opt.frac_monotone 100*opt.frac_single_cross
    end

    println()
    println("="^60)
    println("2-BODY LIKE CHARGES (q1=+1, q2=+1) in 2D — supercritical")
    println("="^60)

    V_like = make_2body_V([1.0, 1.0], [1.0, 1.0], 1.0)

    # Critical radius: ρ = q1*q2/(μ*c²) = 1/(0.5*1) = 2.0
    ρ = 2.0
    @printf "Critical radius ρ = %.2f\n" ρ

    # For like charges V(r) = +1/r (at p=0). The Hill region {V <= E} is {r >= 1/E}.
    # The supercritical zone is r > ρ = 2. So we need 1/E > ρ, i.e. E < 1/ρ = 0.5.
    # For E >= 0.5, the entire supercritical zone satisfies V <= E.
    for E in [0.1, 0.2, 0.3, 0.45, 0.5, 1.0, 2.0]
        r_boundary = 1.0 / E  # V(r) = E at r = 1/E
        # Hill region in supercritical zone: r >= max(ρ, 1/E)
        r_inner = max(ρ, r_boundary)
        # Place center at r_inner + generous offset
        r_center = r_inner + 3.0
        q_star = [r_center/2, 0.0, -r_center/2, 0.0]

        V_star = V_like(q_star)
        if V_star > E
            @printf "E=%.2f: q* outside Hill region (V(q*)=%.4f), skipping\n" E V_star
            continue
        end

        res = shoot_rays(V_like, q_star, E, 200; t_max=max(20.0, 3*r_center), n_samples=300, seed=42)
        @printf "E=%.2f: monotone=%.1f%%, single-cross=%.1f%% (center at r=%.3f, V*=%.4f)\n" E 100*res.frac_monotone 100*res.frac_single_cross r_center V_star
    end
end

function run_4body_test()
    println()
    println("="^60)
    println("4-BODY ALTERNATING SQUARE")
    println("="^60)

    # Build system
    sys = HamiltonianSystem(4, 2)

    # Alternating square configuration
    s = 0.5  # side/2 = 0.5
    q0 = [s, s, -s, -s, -s, s, s, -s]
    masses = [1.0, 1.0, 1.0, 1.0]
    charges = [1.0, 1.0, -1.0, -1.0]
    c = 1.0
    kappas = ones(6)
    params = vcat(masses, charges, [c], kappas)
    zero_p = zeros(8)

    V(q) = sys.hamiltonian_compiled(q, zero_p, params)
    E_sq = V(q0)
    @printf "V(square) = %.6f\n" E_sq

    res = shoot_rays(V, q0, E_sq, 200; t_max=50.0, n_samples=400, seed=42)
    @printf "Square as star-center: monotone=%.1f%%, single-cross=%.1f%%\n" 100*res.frac_monotone 100*res.frac_single_cross

    # Try inflated square (larger separation, deeper in Hill region)
    for scale in [1.5, 2.0, 3.0, 5.0]
        q_inflated = q0 .* scale
        V_infl = V(q_inflated)
        if V_infl > E_sq
            @printf "  scale=%.1f: outside Hill region (V=%.4f > E=%.4f)\n" scale V_infl E_sq
            continue
        end
        res_infl = shoot_rays(V, q_inflated, E_sq, 200; t_max=50.0, n_samples=400, seed=42)
        @printf "  scale=%.1f: V=%.4f, monotone=%.1f%%, single-cross=%.1f%%\n" scale V_infl 100*res_infl.frac_monotone 100*res_infl.frac_single_cross
    end

    # Optimization from square
    opt = optimize_star_center(V, E_sq, q0, 100; n_iters=200, step_size=1.0,
                               t_max=50.0, n_samples=200)
    @printf "Optimized: monotone=%.1f%%, single-cross=%.1f%%\n" 100*opt.frac_monotone 100*opt.frac_single_cross

    println()
    println("="^60)
    println("4-BODY RHOMBUS (a=1.0, b=0.6)")
    println("="^60)

    q_rh = [1.0, 0.0, -1.0, 0.0, 0.0, 0.6, 0.0, -0.6]
    E_rh = V(q_rh)
    @printf "V(rhombus) = %.6f\n" E_rh

    res_rh = shoot_rays(V, q_rh, E_rh, 200; t_max=50.0, n_samples=400, seed=42)
    @printf "Rhombus as star-center: monotone=%.1f%%, single-cross=%.1f%%\n" 100*res_rh.frac_monotone 100*res_rh.frac_single_cross

    opt_rh = optimize_star_center(V, E_rh, q_rh, 100; n_iters=200, step_size=1.0,
                                  t_max=50.0, n_samples=200)
    @printf "Optimized: monotone=%.1f%%, single-cross=%.1f%%\n" 100*opt_rh.frac_monotone 100*opt_rh.frac_single_cross
end

# Main
if abspath(PROGRAM_FILE) == @__FILE__
    run_2body_test()
    run_4body_test()
end
