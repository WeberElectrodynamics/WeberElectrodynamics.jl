#!/usr/bin/env julia
# Capture regression fixtures against the current API.
#
# This script is run ONCE against `main` before Phase 1 of the refactor begins.
# It produces four JLD2 fixture files under test/regression/fixtures/ covering:
#
#   1. twobody_ellipse.jld2        — 2-body unregularized elliptic orbit
#   2. threebody_mixed.jld2        — 3-body unregularized mixed-charge system
#   3. close_approach_lifted.jld2  — 2-body close approach, :lifted_pair backend
#   4. zollner_offmatch.jld2       — 2-body Zöllner with a ≠ 0, regularized
#
# After the refactor, test/regression/validate.jl rebuilds each problem with
# the current (new-API) code, runs it, and asserts pointwise agreement with the
# captured trajectory to 1e-12. This is the single most important invariant of
# the refactor.
#
# Run from project root using the default Julia env (where the package is dev'd
# and JLD2 is available):
#   julia test/regression/capture.jl
#
# JLD2 is listed in [extras] for the test target so validate.jl can load these
# fixtures during `Pkg.test()` without polluting the package's runtime deps.

using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator, RegularizedIntegrator
using JLD2
using Dates
using Symbolics: pkgversion

const FIXTURE_DIR = joinpath(@__DIR__, "fixtures")
const GIT_HASH = try
    strip(read(`git -C $(@__DIR__) rev-parse HEAD`, String))
catch
    "unknown"
end

# ---------------------------------------------------------------------------
# Fixture schema — stored top-level as a Dict{String,Any} for JLD2 robustness.
# ---------------------------------------------------------------------------

function reg_opts_to_dict(r::RegularizationOptions)
    Dict{String,Any}(
        "enabled" => r.enabled,
        "r_on" => r.r_on,
        "r_off" => r.r_off,
        "r_on_factor" => r.r_on_factor,
        "r_off_factor" => r.r_off_factor,
        "max_substeps" => r.max_substeps,
        "constraint_tolerance" => r.constraint_tolerance,
        "g_floor" => r.g_floor,
        "chain_enabled" => r.chain_enabled,
        "backend" => String(r.backend),
        "warn_on_fallback" => r.warn_on_fallback,
        "collision_bounce_radius" => r.collision_bounce_radius,
    )
end

function zollner_opts_to_dict(z::ZollnerOptions)
    Dict{String,Any}("enabled" => z.enabled, "a" => z.a)
end

function diagnostics_to_dict(d::RegularizationDiagnostics)
    Dict{String,Any}(
        "enabled" => d.enabled,
        "requested_backend" => String(d.requested_backend),
        "used_backend" => String(d.used_backend),
        "activation_count" => d.activation_count,
        "deactivation_count" => d.deactivation_count,
        "active_steps" => d.active_steps,
        "pair_steps" => d.pair_steps,
        "adaptive_pair_steps" => d.adaptive_pair_steps,
        "lifted_pair_steps" => d.lifted_pair_steps,
        "chain_steps" => d.chain_steps,
        "unregularized_steps" => d.unregularized_steps,
        "backend_fallback_steps" => d.backend_fallback_steps,
        "total_substeps" => d.total_substeps,
        "max_substeps_used" => d.max_substeps_used,
        "max_constraint_violation" => d.max_constraint_violation,
        "min_encounter_distance" => d.min_encounter_distance,
    )
end

function trajectory_to_dict(sol)
    dof = length(sol.q[1])
    nsteps = length(sol.t)
    q_mat = Matrix{Float64}(undef, dof, nsteps)
    p_mat = Matrix{Float64}(undef, dof, nsteps)
    @inbounds for i = 1:nsteps
        q_mat[:, i] = sol.q[i]
        p_mat[:, i] = sol.p[i]
    end
    Dict{String,Any}(
        "t" => copy(sol.t),
        "q" => q_mat,
        "p" => p_mat,
        "retcode" => String(sol.retcode),
        "n_steps" => nsteps,
        "dof" => dof,
    )
end

function metadata_dict(name::String, description::String)
    Dict{String,Any}(
        "name" => name,
        "description" => description,
        "git_commit" => GIT_HASH,
        "captured_at" => string(Dates.now(Dates.UTC)),
        "julia_version" => string(VERSION),
        "package_version" => "0.4.3",
    )
end

# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------

function fixture_twobody_ellipse()
    # 2-body bound elliptic orbit, unregularized, finite c.
    m1, m2 = 1.0, 0.1
    q1, q2 = sqrt(0.1), -sqrt(0.1)
    c = 4.0
    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(abs(q1 * q2) * M / (m1 * m2 * r0))
    v_scale = 0.9  # sub-circular → mild eccentricity
    q_initial = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p_initial =
        [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, m2 * (m1 / M * v_circ * v_scale)]
    tspan = (0.0, 2.0)
    dt = 5e-4

    system = HamiltonianSystem(2, 2)
    prob = HamiltonianProblem(
        system,
        tspan,
        q_initial,
        p_initial;
        masses = [m1, m2],
        charges = [q1, q2],
        c = c,
        dt = dt,
    )
    alg = SymmetricProjectionIntegrator()
    sol = solve(prob, alg)
    return prob,
    alg,
    sol,
    ZollnerOptions(),
    "twobody_ellipse",
    "2-body unregularized bound elliptic orbit (finite c, mild eccentricity)"
end

function fixture_threebody_mixed()
    # 3-body mixed-charge, unregularized. Starts from rest; well-separated so no
    # close encounters trigger. Weber velocity corrections are small.
    system = HamiltonianSystem(3, 2)
    q_initial = [1.0, 0.0, -0.5, 0.0, 0.2, 1.5]
    p_initial = zeros(6)
    masses = [1.0, 1.0, 0.5]
    charges = [0.1, -0.1, 0.05]
    c = 10.0
    tspan = (0.0, 0.5)
    dt = 5e-4
    prob = HamiltonianProblem(
        system,
        tspan,
        q_initial,
        p_initial;
        masses = masses,
        charges = charges,
        c = c,
        dt = dt,
    )
    alg = SymmetricProjectionIntegrator()
    sol = solve(prob, alg)
    return prob,
    alg,
    sol,
    ZollnerOptions(),
    "threebody_mixed",
    "3-body unregularized mixed-charge system starting from rest"
end

function fixture_close_approach_lifted()
    # 2-body close-approach exercising :lifted_pair Levi-Civita (2D).
    # Apocenter-start orbit with v_tan = 0.5·v_circ → eccentricity 0.75,
    # pericenter ≈ 0.286, period ≈ 7.3. tspan=8.0 guarantees ≥1 pericenter pass.
    m1, m2 = 1.0, 0.1
    q1, q2 = sqrt(0.1), -sqrt(0.1)
    c = 4.0
    r0 = 2.0
    M = m1 + m2
    v_circ = sqrt(abs(q1 * q2) * M / (m1 * m2 * r0))
    v_scale = 0.5
    q_initial = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p_initial =
        [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, m2 * (m1 / M * v_circ * v_scale)]
    tspan = (0.0, 8.0)
    dt = 2e-3

    system = HamiltonianSystem(2, 2)
    prob = HamiltonianProblem(
        system,
        tspan,
        q_initial,
        p_initial;
        masses = [m1, m2],
        charges = [q1, q2],
        c = c,
        dt = dt,
    )
    alg = RegularizedIntegrator(
        SymmetricProjectionIntegrator();
        backend = :lifted_pair,
        r_on_factor = 0.3,
        r_off_factor = 0.45,
        max_substeps = 512,
        constraint_tolerance = 1e-12,
        warn_on_fallback = false,
    )
    sol = solve(prob, alg)
    return prob,
    alg,
    sol,
    ZollnerOptions(),
    "close_approach_lifted",
    "2-body close approach with :lifted_pair Levi-Civita regularization"
end

function fixture_zollner_offmatch()
    # 2-body attractive pair with Zöllner mismatch; adaptive-Cartesian regularized.
    # Apocenter-start eccentric orbit that plunges through pericenter, so κ=1+a
    # propagates through the regularization substep code. Short period so
    # tspan=5.0 spans a full orbit.
    m1, m2 = 1.0, 1.0
    q1, q2 = 0.1, -0.1
    c = 10.0
    r0 = 0.5
    M = m1 + m2
    v_circ = sqrt(abs(q1 * q2) * M / (m1 * m2 * r0))
    v_scale = 0.5
    q_initial = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
    p_initial =
        [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, m2 * (m1 / M * v_circ * v_scale)]
    tspan = (0.0, 5.0)
    dt = 1e-3

    zol = ZollnerOptions(enabled = true, a = 0.05)

    system = HamiltonianSystem(2, 2)
    prob = HamiltonianProblem(
        system,
        tspan,
        q_initial,
        p_initial;
        masses = [m1, m2],
        charges = [q1, q2],
        c = c,
        dt = dt,
        zollner = zol,
    )
    alg = RegularizedIntegrator(
        SymmetricProjectionIntegrator();
        backend = :adaptive_cartesian,
        r_on_factor = 0.4,
        r_off_factor = 0.6,
        max_substeps = 512,
        warn_on_fallback = false,
    )
    sol = solve(prob, alg)
    return prob,
    alg,
    sol,
    zol,
    "zollner_offmatch",
    "2-body Zöllner off-match (a≠0) with adaptive-Cartesian regularization"
end

# ---------------------------------------------------------------------------
# Write JLD2 fixture file
# ---------------------------------------------------------------------------

function _alg_reg_opts(::SymmetricProjectionIntegrator)
    RegularizationOptions()
end
function _alg_reg_opts(alg::RegularizedIntegrator)
    alg.options
end

function save_fixture(
    prob::HamiltonianProblem,
    alg,
    sol::HamiltonianSolution,
    zol::ZollnerOptions,
    name::String,
    desc::String,
)
    setup = Dict{String,Any}(
        "n_particles" => prob.system.n_particles,
        "dims" => prob.system.dims,
        "q_initial" => copy(prob.q_initial),
        "p_initial" => copy(prob.p_initial),
        "masses" => collect(masses(prob)),
        "charges" => collect(charges(prob)),
        "c" => speed_of_light(prob),
        "kappas" => collect(kappas(prob)),
        "params" => copy(prob.params),
        "tspan" => collect(prob.tspan),
        "dt" => prob.dt,
        "convergence_tolerance" => prob.convergence_tolerance,
        "maximum_iterations" => prob.maximum_iterations,
        "regularization" => reg_opts_to_dict(_alg_reg_opts(alg)),
        "zollner" => zollner_opts_to_dict(zol),
    )

    fixture = Dict{String,Any}(
        "metadata" => metadata_dict(name, desc),
        "setup" => setup,
        "trajectory" => trajectory_to_dict(sol),
        "diagnostics" => diagnostics_to_dict(sol.regularization),
    )

    path = joinpath(FIXTURE_DIR, "$(name).jld2")
    jldopen(path, "w") do file
        for (k, v) in fixture
            write(file, k, v)
        end
    end
    sz = filesize(path)
    println(
        "  wrote $path  ($(round(sz / 1024; digits=1)) KiB)  " *
        "retcode=$(sol.retcode)  n_steps=$(length(sol.t))",
    )
    return path
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    mkpath(FIXTURE_DIR)
    println("# Regression fixture capture")
    println("Julia $(VERSION)  |  git $(GIT_HASH)  |  $(Dates.now(Dates.UTC))")
    println()

    for builder in (
        fixture_twobody_ellipse,
        fixture_threebody_mixed,
        fixture_close_approach_lifted,
        fixture_zollner_offmatch,
    )
        name_guess = string(nameof(builder))
        print("Building $name_guess ... ")
        t0 = time()
        prob, alg, sol, zol, name, desc = builder()
        elapsed = time() - t0
        println("solve=$(round(elapsed; digits=2))s")
        save_fixture(prob, alg, sol, zol, name, desc)
        if sol.retcode != :Success
            error("Fixture $name produced retcode=$(sol.retcode); aborting.")
        end
    end

    println()
    println("All fixtures captured successfully.")
end

main()
