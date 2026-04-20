using JLD2

@testset "RegularizedIntegrator wrapper" begin
    @testset "subtyping + unwrap" begin
        base = SymmetricProjectionIntegrator(relaxation = 0.25)
        alg = WeberElectrodynamics.RegularizedIntegrator(base; backend = :lifted_pair)
        @test alg isa WeberElectrodynamics.HamiltonianAlgorithm
        @test WeberElectrodynamics.base_algorithm(alg) === base
        @test WeberElectrodynamics.base_algorithm(base) === base
        @test alg.options.enabled === true
        @test alg.options.backend === :lifted_pair
        @test alg.options.r_on_factor == 0.15
    end

    @testset "options kwargs forwarded" begin
        alg = WeberElectrodynamics.RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            r_on_factor = 0.2,
            r_off_factor = 0.5,
            max_substeps = 64,
            backend = :adaptive_cartesian,
            chain_enabled = false,
            warn_on_fallback = false,
            collision_bounce_radius = 0.01,
        )
        @test alg.options.r_on_factor == 0.2
        @test alg.options.r_off_factor == 0.5
        @test alg.options.max_substeps == 64
        @test alg.options.backend === :adaptive_cartesian
        @test alg.options.chain_enabled === false
        @test alg.options.warn_on_fallback === false
        @test alg.options.collision_bounce_radius == 0.01
    end

    # The wrapper must produce trajectories numerically identical to those from
    # the equivalent prob-level `regularization = ...` kwarg. This verifies the
    # Phase 3c.1 shim — rewriting `prob.regularization` at init time — preserves
    # bit-exact output against the existing regularization dispatch code.
    @testset "equivalence with prob-level regularization" begin
        fixture_path =
            joinpath(@__DIR__, "regression", "fixtures", "close_approach_lifted.jld2")
        fixture = jldopen(fixture_path, "r") do file
            Dict{String,Any}(
                "setup" => read(file, "setup"),
                "trajectory" => read(file, "trajectory"),
            )
        end
        setup = fixture["setup"]
        reg_dict = setup["regularization"]
        reg_opts = RegularizationOptions(
            enabled = reg_dict["enabled"],
            r_on = reg_dict["r_on"],
            r_off = reg_dict["r_off"],
            r_on_factor = reg_dict["r_on_factor"],
            r_off_factor = reg_dict["r_off_factor"],
            max_substeps = reg_dict["max_substeps"],
            constraint_tolerance = reg_dict["constraint_tolerance"],
            g_floor = reg_dict["g_floor"],
            chain_enabled = reg_dict["chain_enabled"],
            backend = Symbol(reg_dict["backend"]),
            warn_on_fallback = reg_dict["warn_on_fallback"],
            collision_bounce_radius = reg_dict["collision_bounce_radius"],
        )

        sys = HamiltonianSystem(setup["n_particles"]::Int, setup["dims"]::Int)
        base_kwargs = (
            masses = setup["masses"]::Vector{Float64},
            charges = setup["charges"]::Vector{Float64},
            c = setup["c"]::Float64,
            dt = setup["dt"]::Float64,
            convergence_tolerance = setup["convergence_tolerance"]::Float64,
            maximum_iterations = setup["maximum_iterations"]::Int,
        )

        # Problem A: regularization baked into the problem (existing API).
        prob_a = HamiltonianProblem(
            sys,
            Tuple(setup["tspan"]::Vector{Float64}),
            setup["q_initial"]::Vector{Float64},
            setup["p_initial"]::Vector{Float64};
            base_kwargs...,
            regularization = reg_opts,
        )
        sol_a = solve(prob_a, SymmetricProjectionIntegrator())

        # Problem B: same setup but regularization via the algorithm wrapper.
        prob_b = HamiltonianProblem(
            sys,
            Tuple(setup["tspan"]::Vector{Float64}),
            setup["q_initial"]::Vector{Float64},
            setup["p_initial"]::Vector{Float64};
            base_kwargs...,
        )
        alg = WeberElectrodynamics.RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            r_on = reg_opts.r_on,
            r_off = reg_opts.r_off,
            r_on_factor = reg_opts.r_on_factor,
            r_off_factor = reg_opts.r_off_factor,
            max_substeps = reg_opts.max_substeps,
            constraint_tolerance = reg_opts.constraint_tolerance,
            g_floor = reg_opts.g_floor,
            chain_enabled = reg_opts.chain_enabled,
            backend = reg_opts.backend,
            warn_on_fallback = reg_opts.warn_on_fallback,
            collision_bounce_radius = reg_opts.collision_bounce_radius,
        )
        sol_b = solve(prob_b, alg)

        @test sol_a.retcode === sol_b.retcode
        @test length(sol_a.t) == length(sol_b.t)
        max_q = maximum(maximum(abs, sol_a.q[i] .- sol_b.q[i]) for i in eachindex(sol_a.q))
        max_p = maximum(maximum(abs, sol_a.p[i] .- sol_b.p[i]) for i in eachindex(sol_a.p))
        @test max_q == 0.0
        @test max_p == 0.0
    end
end
