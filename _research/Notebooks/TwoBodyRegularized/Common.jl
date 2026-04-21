# Common.jl — shared setup for the TwoBodyRegularized notebook tour.
#
# Each notebook in this folder includes this file so the narrative cells stay
# focused on the physics and regularization being demonstrated rather than on
# boilerplate. Nothing here reimplements package functionality — all physics,
# statistics, and plotting come directly from WeberElectrodynamics.

using WeberElectrodynamics
using Plots
using Printf
using LinearAlgebra

# Canonical physical constants — identical to examples/two_body_reference.ipynb
const M_CANON            = (1.0, 1.0)
const Q_CANON_UNLIKE     = (1.0, -1.0)
const Q_CANON_LIKE       = (1.0, 1.0)
const C_CANON            = 4.0
const DT_CANON           = 0.001

# Canonical Kepler orbit (e = 3/4, a = 8/7, μ = 1/2, |k| = 1)
const A_SEMI_CANON       = 8 / 7
const MU_CANON           = 0.5
const T_ORBIT_CANON      = 2π * sqrt(A_SEMI_CANON^3 * MU_CANON / 1.0)  # ≈ 5.4282

"""
    canonical_2d_problem(; n_orbits=5, dt=DT_CANON, zollner=ZollnerOptions())

Build the canonical 2D two-body Kepler-like Weber problem from
`examples/two_body_reference.ipynb`: equal unit masses, unlike unit charges,
`c = 4`, apoapsis initial conditions giving `e = 3/4`, `a = 8/7`, period
`T ≈ 5.4282`.

Regularization is no longer a problem-level option; see
[`canonical_2d_algorithm`](@ref) for the matching algorithm wrapper.
"""
function canonical_2d_problem(;
    n_orbits::Real = 5,
    dt::Real = DT_CANON,
    zollner::ZollnerOptions = ZollnerOptions(),
)
    system = HamiltonianSystem(2, 2)
    q0 = [-1.0, 0.0, 1.0, 0.0]
    p0 = [ 0.0, -0.25, 0.0, 0.25]
    tspan = (0.0, n_orbits * T_ORBIT_CANON)
    return HamiltonianProblem(system, tspan, q0, p0;
        masses  = collect(M_CANON),
        charges = collect(Q_CANON_UNLIKE),
        c       = C_CANON,
        dt      = dt,
        zollner = zollner,
    )
end

"""
    canonical_2d_algorithm(; regularization=nothing, collision_bounce_radius=0.0)

Build the canonical algorithm for the two-body regularized tour.

- `regularization = nothing` and `collision_bounce_radius = 0.0` →
  bare `SymmetricProjectionIntegrator()`.
- `regularization::NamedTuple` → wraps the base algorithm in
  `RegularizedIntegrator(...; regularization..., collision_bounce_radius=...)`.
  NamedTuple kwargs are forwarded verbatim (`backend`, `r_on`, `r_off`,
  `r_on_factor`, `r_off_factor`, `chain_enabled`, …).

Pass the result as the second argument to `solve(prob, alg)`.
"""
function canonical_2d_algorithm(;
    regularization::Union{Nothing,NamedTuple} = nothing,
    collision_bounce_radius::Real = 0.0,
)
    if regularization === nothing && collision_bounce_radius == 0.0
        return SymmetricProjectionIntegrator()
    end
    reg_kwargs = regularization === nothing ? NamedTuple() : regularization
    return RegularizedIntegrator(
        SymmetricProjectionIntegrator();
        reg_kwargs...,
        collision_bounce_radius = collision_bounce_radius,
    )
end

"""
    figures_dir()

Return (creating if necessary) the `figures/` directory next to this file.
"""
function figures_dir()
    dir = joinpath(@__DIR__, "figures")
    isdir(dir) || mkpath(dir)
    return dir
end

"""
    print_reg_diagnostics(sol; label="")

Pretty-print every field of `sol.regularization` plus the retcode. Used by
every notebook in the tour to show what the backend actually did.
"""
function print_reg_diagnostics(sol::HamiltonianSolution; label::AbstractString = "")
    d = sol.regularization
    prefix = isempty(label) ? "" : "[$label] "
    @printf("%sretcode                  = %s\n", prefix, sol.retcode)
    @printf("%srequested_backend        = %s\n", prefix, d.requested_backend)
    @printf("%sused_backend             = %s\n", prefix, d.used_backend)
    @printf("%sactive_steps             = %d\n", prefix, d.active_steps)
    @printf("%spair_steps               = %d\n", prefix, d.pair_steps)
    @printf("%sadaptive_pair_steps      = %d\n", prefix, d.adaptive_pair_steps)
    @printf("%slifted_pair_steps        = %d\n", prefix, d.lifted_pair_steps)
    @printf("%schain_steps              = %d\n", prefix, d.chain_steps)
    @printf("%sunregularized_steps      = %d\n", prefix, d.unregularized_steps)
    @printf("%sbackend_fallback_steps   = %d\n", prefix, d.backend_fallback_steps)
    @printf("%stotal_substeps           = %d\n", prefix, d.total_substeps)
    @printf("%smax_substeps_used        = %d\n", prefix, d.max_substeps_used)
    @printf("%sactivation_count         = %d\n", prefix, d.activation_count)
    @printf("%sdeactivation_count       = %d\n", prefix, d.deactivation_count)
    @printf("%smin_encounter_distance   = %.6g\n", prefix, d.min_encounter_distance)
    @printf("%smax_constraint_violation = %.3e\n", prefix, d.max_constraint_violation)
    return nothing
end

"""
    energy_drift_percent(energy)

Shorthand for the max-global-energy-error percentage (for printed tables).
"""
energy_drift_percent(energy) = energy.statistics.global_error_percent_max

"""
    final_position_delta(sol_a, sol_b)

L2 norm of the position difference at the final step between two solutions
that started from identical initial conditions. Used to quantify how close
two regularized integrations track each other.
"""
function final_position_delta(sol_a::HamiltonianSolution, sol_b::HamiltonianSolution)
    return norm(sol_a.q[end] .- sol_b.q[end])
end

nothing
