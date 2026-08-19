@testset "CommonSolve Interface" begin
    @testset "solve basic" begin
        prob = make_weber_problem(tspan = (0.0, 0.5), dt = 0.001)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.t) > 1
        @test sol.t[1] == 0.0
        @test sol.t[end] ≈ 0.5 atol = 0.01
    end

    @testset "solve with custom algorithm" begin
        prob = make_weber_problem(tspan = (0.0, 0.1), dt = 0.001)
        alg = SymmetricProjectionIntegrator(relaxation = 0.3)
        sol = solve(prob, alg)

        @test sol.retcode == :Success
    end

    @testset "init" begin
        prob = make_weber_problem()
        integrator = init(prob)

        @test integrator.t == prob.tspan[1]
        @test integrator.t_end == prob.tspan[2]
        @test integrator.q == prob.q_initial
        @test integrator.p == prob.p_initial
        @test integrator.step_count == 0
        @test length(integrator.t_history) > 0
    end

    @testset "step!" begin
        prob = make_weber_problem(tspan = (0.0, 1.0), dt = 0.01)
        integrator = init(prob)

        # First step
        result = step!(integrator)
        @test result == true  # More steps to go
        @test integrator.step_count == 1
        @test integrator.t ≈ 0.01
        @test integrator.t_history[2] ≈ 0.01

        # State changed from initial
        @test integrator.q != prob.q_initial || integrator.p != prob.p_initial

        # Continue stepping
        for _ = 1:5
            step!(integrator)
        end
        @test integrator.step_count == 6
        @test integrator.t ≈ 0.06
    end

    @testset "step! returns false at end" begin
        prob = make_weber_problem(tspan = (0.0, 0.05), dt = 0.01)
        integrator = init(prob)

        # Step until done
        steps = 0
        while step!(integrator)
            steps += 1
        end

        @test integrator.t >= prob.tspan[2] - eps(prob.tspan[2])

        # Calling step! again returns false
        @test step!(integrator) == false
    end

    @testset "step! clamps final step to tspan end" begin
        prob = make_weber_problem(tspan = (0.0, 1.0), dt = 0.3)
        integrator = init(prob)
        while step!(integrator)
        end

        @test integrator.t ≈ prob.tspan[2]
        @test integrator.t <= prob.tspan[2] + eps(prob.tspan[2])

        sol = solve(prob)
        @test sol.t[end] ≈ prob.tspan[2]
        @test sol.t[end] <= prob.tspan[2] + eps(prob.tspan[2])
    end

    @testset "solve!" begin
        prob = make_weber_problem(tspan = (0.0, 0.5), dt = 0.001)
        integrator = init(prob)

        # Step a few times first
        step!(integrator)
        step!(integrator)
        @test integrator.step_count == 2

        # Complete integration
        sol = solve!(integrator)

        @test sol.retcode == :Success
        @test sol.t[1] == 0.0
        @test sol.t[end] ≈ 0.5 atol = 0.01
        @test length(sol.t) == integrator.step_count + 1
    end

    @testset "solve vs init+solve! equivalence" begin
        prob = make_weber_problem(tspan = (0.0, 0.5), dt = 0.001)

        sol1 = solve(prob)

        integrator = init(prob)
        sol2 = solve!(integrator)

        @test length(sol1.t) == length(sol2.t)
        @test sol1.t ≈ sol2.t
        @test all(sol1.q .≈ sol2.q)
        @test all(sol1.p .≈ sol2.p)
    end

    @testset "History pre-allocation" begin
        prob = make_weber_problem(tspan = (0.0, 1.0), dt = 0.001)
        integrator = init(prob)

        # Check history vectors are pre-allocated
        expected_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
        @test length(integrator.t_history) >= expected_steps
        @test length(integrator.q_history) >= expected_steps
        @test length(integrator.p_history) >= expected_steps
    end

    @testset "save_stride stores sparse history and final state" begin
        prob = make_weber_problem(tspan = (0.0, 0.1), dt = 0.01)
        sol = solve(prob; save_stride = 4)

        @test sol.retcode == :Success
        @test sol.t ≈ [0.0, 0.04, 0.08, 0.1]
        @test length(sol.q) == 4
        @test length(sol.p) == 4

        sol_alias = solve(prob; save_every = 5)
        @test sol_alias.t ≈ [0.0, 0.05, 0.1]
        @test_throws ArgumentError solve(prob; save_stride = 2, save_every = 5)
        @test_throws ArgumentError solve(prob; save_stride = 0)
    end

    @testset "stream_sink receives every macro-step" begin
        prob = make_weber_problem(tspan = (0.0, 0.03), dt = 0.01)
        seen = Tuple{Float64,Int}[]
        sink = (t, q, p, step) -> push!(seen, (t, step))

        sol = solve(prob; save_stride = 10, stream_sink = sink)

        @test sol.retcode == :Success
        @test sol.t ≈ [0.0, 0.03]
        @test first.(seen) ≈ [0.0, 0.01, 0.02, 0.03]
        @test last.(seen) == [0, 1, 2, 3]
    end

    @testset "stream_sink fires per macro-step under regularized integrator" begin
        # Regression guard: when the regularized integrator subdivides a
        # macro-step into many lifted substeps, stream_sink must still be
        # invoked once per macro-step (not once per substep). Use a tight
        # 2D close-encounter setup with a coarse macro-dt so :lifted_pair
        # actually subdivides each macro-step into many substeps.
        sys = HamiltonianSystem(2, 2)
        q0 = [0.4, 0.0, -0.4, 0.0]
        p0 = [0.0, 0.15, 0.0, -0.15]
        prob = HamiltonianProblem(
            sys,
            (0.0, 1.0),
            q0,
            p0;
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.05,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :lifted_pair,
            warn_on_fallback = false,
            r_on = 0.5,
            r_off = 0.6,
            chain_enabled = false,
        )

        seen = Tuple{Float64,Int}[]
        sink = (t, q, p, step) -> push!(seen, (t, step))
        sol = solve(prob, alg; stream_sink = sink)

        @test sol.retcode == :Success
        # One emission per macro-step plus one for the initial state.
        @test length(seen) == length(sol.t)
        @test last.(seen) == collect(0:(length(sol.t)-1))
        # Lifted-pair was actually used and subdivided into multiple substeps
        # (otherwise the regression-guard premise is moot).
        @test sol.regularization.lifted_pair_steps > 0
        @test sol.regularization.total_substeps > sol.regularization.lifted_pair_steps
    end

    @testset "JLD2 solution archive helpers" begin
        prob = make_weber_problem(tspan = (0.0, 0.02), dt = 0.01)
        sol = solve(prob)
        path = tempname() * ".jld2"

        @test save_solution(path, sol; metadata = (case = "archive",)) == path
        archive = JLD2.load(path)["archive"]
        @test archive.format_version == 2
        # Pin the exact v2 problem schema. This is stronger than asserting the
        # absence of the removed fields, and it keeps the removed identifiers
        # out of the tree entirely.
        @test propertynames(archive.problem) == (
            :n_particles,
            :dims,
            :tspan,
            :q_initial,
            :p_initial,
            :params,
            :dt,
            :convergence_tolerance,
            :maximum_iterations,
        )
        loaded = load_solution(path)

        @test loaded isa HamiltonianSolution
        @test loaded.t == sol.t
        @test loaded.q == sol.q
        @test loaded.p == sol.p
        @test loaded.retcode == sol.retcode
        rm(path; force = true)

        old_path = tempname() * ".jld2"
        JLD2.jldsave(old_path; archive = (format_version = 1,))
        @test_throws ArgumentError load_solution(old_path)
        rm(old_path; force = true)
    end

    @testset "Two-body system solve" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.q[1]) == 4  # 2 particles × 2 dims
        @test length(sol.p[1]) == 4
    end

    @testset "Weber system solve" begin
        prob = make_weber_problem(tspan = (0.0, 0.5), dt = 0.001)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.t) > 1
    end

    @testset "Custom tolerance" begin
        system = HamiltonianSystem(2, 2)
        q0 = [1.0, 0.0, -1.0, 0.0]
        p0 = [0.0, 0.05, 0.0, -0.05]

        prob_tight = HamiltonianProblem(
            system,
            (0.0, 0.1),
            q0,
            p0;
            masses = [1.0, 0.1],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
            convergence_tolerance = 1e-14,
        )
        sol = solve(prob_tight)
        @test sol.retcode == :Success
    end
end
