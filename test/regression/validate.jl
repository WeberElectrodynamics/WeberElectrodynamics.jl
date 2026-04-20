#!/usr/bin/env julia
# Numerical-equivalence validation against captured regression fixtures.
#
# For each .jld2 fixture under test/regression/fixtures/, this script:
#   1. Reads the captured setup (masses, charges, c, ICs, dt, tspan,
#      regularization/Zöllner options).
#   2. Rebuilds the problem using the CURRENT API — whatever `main` looks like
#      at the time of running.
#   3. Runs `solve(prob, SymmetricProjectionIntegrator())`.
#   4. Asserts the re-run trajectory matches the captured trajectory
#      pointwise to 1e-12 on q and p.
#
# This is the phase-gate invariant for the architectural refactor: any change
# that drifts the integrator output past 1e-12 blocks the phase merge.
#
# Run standalone (reports result + exits non-zero on failure):
#   julia test/regression/validate.jl
#
# Or include from a testset:
#   include("test/regression/validate.jl"); validate_all()

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using JLD2
using Printf

const FIXTURE_DIR = joinpath(@__DIR__, "fixtures")
const FIXTURES = [
    "twobody_ellipse",
    "threebody_mixed",
    "close_approach_lifted",
    "zollner_offmatch",
]

# Hard threshold — any drift above this blocks a phase merge. The integrator
# is deterministic at full Float64 precision, so any nonzero drift signals a
# bug, not a tolerance issue.
const TOLERANCE = 1e-12

# Reconstruct RegularizationOptions from a fixture setup dict.
function _rebuild_reg_options(d::Dict{String,Any})
    backend_sym = Symbol(d["backend"])
    RegularizationOptions(
        enabled = d["enabled"],
        r_on = d["r_on"],
        r_off = d["r_off"],
        r_on_factor = d["r_on_factor"],
        r_off_factor = d["r_off_factor"],
        max_substeps = d["max_substeps"],
        constraint_tolerance = d["constraint_tolerance"],
        g_floor = d["g_floor"],
        chain_enabled = d["chain_enabled"],
        backend = backend_sym,
        warn_on_fallback = d["warn_on_fallback"],
        collision_bounce_radius = d["collision_bounce_radius"],
    )
end

function _rebuild_zollner_options(d::Dict{String,Any})
    ZollnerOptions(enabled = d["enabled"], a = d["a"])
end

# Rebuild the HamiltonianProblem corresponding to a fixture.
#
# This function is the single point that must be updated each time the public
# problem-building API changes (e.g. Phase 1 of the refactor will rewrite it
# to use HamiltonianSystem + weber_term + HamiltonianProblem).
function rebuild_problem(setup::Dict{String,Any})
    n_particles = setup["n_particles"]::Int
    dims = setup["dims"]::Int
    system = HamiltonianSystem(n_particles, dims)

    reg = _rebuild_reg_options(setup["regularization"])
    zol = _rebuild_zollner_options(setup["zollner"])

    tspan = Tuple(setup["tspan"]::Vector{Float64})

    return HamiltonianProblem(
        system,
        tspan,
        setup["q_initial"]::Vector{Float64},
        setup["p_initial"]::Vector{Float64};
        masses = setup["masses"]::Vector{Float64},
        charges = setup["charges"]::Vector{Float64},
        c = setup["c"]::Float64,
        dt = setup["dt"]::Float64,
        convergence_tolerance = setup["convergence_tolerance"]::Float64,
        maximum_iterations = setup["maximum_iterations"]::Int,
        regularization = reg,
        zollner = zol,
    )
end

# Compare a fresh solution against a captured trajectory.
# Returns (max_t_err, max_q_err, max_p_err, retcode_match).
function compare_trajectory(sol, captured::Dict{String,Any})
    n_new = length(sol.t)
    n_ref = captured["n_steps"]::Int
    @assert n_new == n_ref "time-grid length mismatch: new=$n_new captured=$n_ref"

    t_ref = captured["t"]::Vector{Float64}
    q_ref = captured["q"]::Matrix{Float64}
    p_ref = captured["p"]::Matrix{Float64}

    max_t_err = maximum(abs.(sol.t .- t_ref))
    max_q_err = 0.0
    max_p_err = 0.0
    @inbounds for i in eachindex(sol.t)
        for k in axes(q_ref, 1)
            max_q_err = max(max_q_err, abs(sol.q[i][k] - q_ref[k, i]))
            max_p_err = max(max_p_err, abs(sol.p[i][k] - p_ref[k, i]))
        end
    end
    retcode_match = String(sol.retcode) == captured["retcode"]::String
    return max_t_err, max_q_err, max_p_err, retcode_match
end

function validate_fixture(name::String)
    path = joinpath(FIXTURE_DIR, "$(name).jld2")
    isfile(path) || error("fixture missing: $path (run capture.jl first)")

    fixture = jldopen(path, "r") do file
        Dict{String,Any}(
            "metadata" => read(file, "metadata"),
            "setup" => read(file, "setup"),
            "trajectory" => read(file, "trajectory"),
            "diagnostics" => read(file, "diagnostics"),
        )
    end

    prob = rebuild_problem(fixture["setup"])
    sol = solve(prob, SymmetricProjectionIntegrator())
    max_t, max_q, max_p, ok_ret = compare_trajectory(sol, fixture["trajectory"])

    pass = max_q ≤ TOLERANCE && max_p ≤ TOLERANCE && max_t ≤ TOLERANCE && ok_ret
    return pass, max_t, max_q, max_p, ok_ret
end

function validate_all()
    println("# Regression validation")
    println("  tolerance: $(TOLERANCE) on |Δt|, |Δq|, |Δp|")
    println()

    all_pass = true
    for name in FIXTURES
        pass, mt, mq, mp, okret = validate_fixture(name)
        status = pass ? "PASS" : "FAIL"
        @printf("  %s  %-26s  |Δt|=%.2e  |Δq|=%.2e  |Δp|=%.2e  retcode=%s\n",
                status, name, mt, mq, mp, okret ? "ok" : "MISMATCH")
        all_pass &= pass
    end

    println()
    if all_pass
        println("All fixtures agree to tolerance $(TOLERANCE).")
    else
        println("FAILURE: one or more fixtures drifted beyond tolerance.")
    end
    return all_pass
end

# Run standalone.
if abspath(PROGRAM_FILE) == @__FILE__
    ok = validate_all()
    exit(ok ? 0 : 1)
end
