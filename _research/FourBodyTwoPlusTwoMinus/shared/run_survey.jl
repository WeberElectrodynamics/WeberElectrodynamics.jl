"""
Unified runner for 4-body 2+/2− Weber surveys. Builds a HamiltonianProblem from
IC arrays, solves, and returns solution + basic statistics.

Usage (from any agent subdir):
    include(joinpath(@__DIR__, "..", "shared", "run_survey.jl"))
    using .SharedSurvey
    sol = SharedSurvey.run(q0, p0, masses, charges; tmax=50.0, dt=1e-3, c=1.0)
"""
module SharedSurvey

using WeberElectrodynamics
using LinearAlgebra

const _SYSTEMS = Dict{Tuple{Int,Int},Any}()

function _system(n::Int, d::Int)
    key = (n, d)
    haskey(_SYSTEMS, key) && return _SYSTEMS[key]
    sys = HamiltonianSystem(n, d)
    _SYSTEMS[key] = sys
    return sys
end

"""
    run(q0, p0, masses, charges; tmax, dt, c, dims=2, zollner_enabled=false,
        zollner_a=0.0, bounce_r=0.0, alg=nothing)
Returns (sol, energy, momentum).
"""
function run(
    q0::AbstractVector,
    p0::AbstractVector,
    masses::AbstractVector,
    charges::AbstractVector;
    tmax::Float64 = 50.0,
    dt::Float64 = 1e-3,
    c::Float64 = 1.0,
    dims::Int = 2,
    zollner_enabled::Bool = false,
    zollner_a::Float64 = 0.0,
    bounce_r::Float64 = 0.0,
)
    n = length(masses)
    sys = _system(n, dims)
    reg = RegularizationOptions(collision_bounce_radius = bounce_r)
    zol = ZollnerOptions(enabled = zollner_enabled, a = zollner_a)
    prob = HamiltonianProblem(
        sys,
        (0.0, tmax),
        q0,
        p0;
        masses = masses,
        charges = charges,
        c = c,
        dt = dt,
        regularization = reg,
        zollner = zol,
    )
    sol = solve(prob, SymmetricProjectionIntegrator())
    return sol
end

"""
    summarize(sol) -> NamedTuple
"""
function summarize(sol)
    en = compute_energy_timeseries(sol)
    mom = compute_momentum_timeseries(sol)
    return (
        retcode = sol.retcode,
        n_steps = length(sol.t),
        t_final = sol.t[end],
        E_drift_pct = en.statistics.global_error_percent_max,
        L_series_available = true,
        energy = en,
        momentum = mom,
    )
end

end  # module
