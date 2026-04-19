"""
Agent 7 — Maximal Lyapunov Exponent (MLE) via two-trajectory shadow method
for the 4-body 2+/2− Weber Hamiltonian.

Method
------
For each IC (q0, p0):
  1. Build state x = [q; p] ∈ ℝ^16 (4 particles × 2 dims × 2 phase-space copies).
  2. Pick random unit vector u ∈ ℝ^16; perturbed state x̃₀ = x₀ + ε u, ε = 1e-8.
  3. Integrate reference and shadow for Δt_renorm = 0.5 with dt = 1e-3.
  4. Measure δ = ||x̃ − x||; accumulate γ = log(δ/ε).
  5. Renormalize x̃ ← x + ε (x̃ − x)/δ, repeat for N intervals.
  6. λ_max ≈ (1/(N Δt)) Σ γ_k.

Each leg uses `SharedSurvey.run` (tmax = Δt_renorm) so we reuse the canonical
integrator and respect Weber/Zöllner/regularization options.

Outputs
-------
- `lyapunov_results.csv` — config × η rows with λ_max, retcode, δ_final.
- `lyapunov_heatmap.png`  — 3 × 4 heatmap of λ_max.
"""

using LinearAlgebra
using Random
using Printf
using Statistics
using DelimitedFiles
using Plots

using WeberElectrodynamics

const SHARED = joinpath(@__DIR__, "..", "shared")
include(joinpath(SHARED, "ic_generators.jl"))
include(joinpath(SHARED, "metrics.jl"))
include(joinpath(SHARED, "run_survey.jl"))
using .SharedSurvey

# --------------------------------------------------------------------------
# Single-shot integration of [q,p] across one renorm interval
# --------------------------------------------------------------------------
function integrate_leg(q0, p0, masses, charges; dt_renorm::Float64, dt::Float64)
    sol = SharedSurvey.run(
        q0, p0, masses, charges;
        tmax = dt_renorm, dt = dt, c = 1.0, dims = 2,
    )
    return sol
end

state(q, p) = vcat(q, p)
splitstate(x, n) = (x[1:n], x[n+1:end])

# --------------------------------------------------------------------------
# Compute MLE for one IC
# --------------------------------------------------------------------------
function compute_mle(
    q0::Vector{Float64}, p0::Vector{Float64},
    masses::Vector{Float64}, charges::Vector{Float64};
    tmax::Float64 = 20.0,
    dt_renorm::Float64 = 0.5,
    dt::Float64 = 1e-3,
    epsilon::Float64 = 1e-8,
    rng::AbstractRNG = MersenneTwister(42),
)
    n_q = length(q0)
    dim_state = 2 * n_q  # 16 for 2D 4-body
    N = Int(round(tmax / dt_renorm))

    # initial perturbation
    u = randn(rng, dim_state)
    u ./= norm(u)

    x_ref = state(q0, p0)
    x_sh  = x_ref .+ epsilon .* u

    log_sum = 0.0
    delta_final = NaN
    for k = 1:N
        qr, pr = splitstate(x_ref, n_q)
        qs, ps = splitstate(x_sh,  n_q)

        sol_r = integrate_leg(qr, pr, masses, charges; dt_renorm = dt_renorm, dt = dt)
        sol_s = integrate_leg(qs, ps, masses, charges; dt_renorm = dt_renorm, dt = dt)

        if sol_r.retcode != :Success || sol_s.retcode != :Success
            # Partial estimate from successfully completed intervals
            if k > 1
                lam = log_sum / ((k - 1) * dt_renorm)
                return (lambda = lam, retcode = :partial,
                        completed = k - 1, delta_final = delta_final)
            else
                return (lambda = NaN, retcode = :failed_immediately,
                        completed = 0, delta_final = delta_final)
            end
        end

        x_ref = state(sol_r.q[end], sol_r.p[end])
        x_sh  = state(sol_s.q[end], sol_s.p[end])

        delta = norm(x_sh .- x_ref)
        if !(delta > 0) || !isfinite(delta)
            return (lambda = NaN, retcode = :degenerate, completed = k, delta_final = delta)
        end
        log_sum += log(delta / epsilon)
        delta_final = delta

        # renormalize
        x_sh = x_ref .+ epsilon .* (x_sh .- x_ref) ./ delta
    end

    lambda = log_sum / (N * dt_renorm)
    return (lambda = lambda, retcode = :Success, completed = N, delta_final = delta_final)
end

# --------------------------------------------------------------------------
# IC grid: 3 configurations × 4 energy fractions (mirrors Agent 6).
# --------------------------------------------------------------------------
const CONFIGS = [
    :alternating_square,
    :rhombus,
    :two_dimers,
]
const ETAS = [0.10, 0.25, 0.45, 0.70]

function build_ic(config::Symbol, eta::Float64)
    if config == :alternating_square
        return alternating_square(side = 1.0, energy_fraction = eta, velocity_mode = :rotating)
    elseif config == :rhombus
        return rhombus(a = 1.0, b = 0.6, energy_fraction = eta, velocity_mode = :rotating)
    elseif config == :two_dimers
        # map η to intra_fraction so we sweep energy in a meaningful way
        # Wider dyad to avoid immediate fixed-point failure; intra_fraction
        # scales energy — dimers with too little KE collapse to head-on collision.
        return two_dimers(
            dyad_length = 0.8, separation = 3.0,
            intra_fraction = clamp(0.3 + 0.6*eta, 0.2, 0.9),
            inter_velocity = 0.0,
        )
    else
        error("unknown config $config")
    end
end

# --------------------------------------------------------------------------
# Run the grid
# --------------------------------------------------------------------------
function run_grid(; tmax::Float64 = 20.0, dt_renorm::Float64 = 0.5, dt::Float64 = 1e-3)
    results = Matrix{Float64}(undef, length(CONFIGS), length(ETAS))
    fill!(results, NaN)
    rows = String[]
    push!(rows, "config,energy_fraction,lambda_max,retcode,completed_intervals,delta_final")

    for (i, cfg) in enumerate(CONFIGS), (j, eta) in enumerate(ETAS)
        q0, p0, m, qch = build_ic(cfg, eta)
        @printf("[%-22s η=%.2f] computing MLE … ", String(cfg), eta)
        t0 = time()
        res = compute_mle(
            collect(q0), collect(p0), collect(m), collect(qch);
            tmax = tmax, dt_renorm = dt_renorm, dt = dt,
            rng = MersenneTwister(1000 + 10*i + j),
        )
        dt_wall = time() - t0
        @printf("λ=%s retcode=%s (%.1fs)\n",
                isnan(res.lambda) ? "NaN" : @sprintf("%.4f", res.lambda),
                String(res.retcode), dt_wall)
        results[i, j] = res.lambda
        push!(rows,
            join([String(cfg), @sprintf("%.3f", eta),
                  isnan(res.lambda) ? "NaN" : @sprintf("%.6f", res.lambda),
                  String(res.retcode), string(res.completed),
                  isnan(res.delta_final) ? "NaN" : @sprintf("%.3e", res.delta_final)],
                 ","))
    end
    return results, rows
end

function plot_heatmap(results::Matrix{Float64}, outfile::String)
    # Replace NaN with a sentinel for display only
    display_mat = copy(results)
    label_mat = similar(results, String)
    for i in eachindex(results)
        if isnan(results[i])
            display_mat[i] = -1.0
            label_mat[i] = "NaN"
        else
            label_mat[i] = @sprintf("%.2f", results[i])
        end
    end
    cfg_labels = ["alt. square", "rhombus", "two dimers"]
    eta_labels = [@sprintf("η=%.2f", e) for e in ETAS]
    h = heatmap(
        eta_labels, cfg_labels, display_mat,
        c = :viridis, clims = (-0.5, max(1.0, maximum(filter(!isnan, results); init = 1.0))),
        xlabel = "energy fraction", ylabel = "configuration",
        title  = "Finite-time MLE λ_max  (4-body 2+/2−, tmax=20, Δt=0.5)",
        size = (700, 420),
        right_margin = 5Plots.mm,
    )
    for i = 1:size(results, 1), j = 1:size(results, 2)
        annotate!(h, j, i, text(label_mat[i, j], 9, :white))
    end
    savefig(h, outfile)
    return outfile
end

# --------------------------------------------------------------------------
# Main entry
# --------------------------------------------------------------------------
function main()
    outdir = @__DIR__
    results, rows = run_grid()
    csvpath = joinpath(outdir, "lyapunov_results.csv")
    open(csvpath, "w") do io
        for r in rows
            println(io, r)
        end
    end
    @info "wrote $csvpath"
    pngpath = plot_heatmap(results, joinpath(outdir, "lyapunov_heatmap.png"))
    @info "wrote $pngpath"
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
