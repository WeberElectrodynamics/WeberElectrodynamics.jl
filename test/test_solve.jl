@testset "CommonSolve Interface" begin
    @testset "solve basic" begin
        prob = make_harmonic_problem(tspan=(0.0, 0.5), dt=0.01)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.t) > 1
        @test sol.t[1] == 0.0
        @test sol.t[end] ≈ 0.5 atol = 0.01
    end

    @testset "solve with custom algorithm" begin
        # Use a short problem with small dt for reliable convergence
        prob = make_harmonic_problem(tspan=(0.0, 0.1), dt=0.001)
        alg = SymmetricProjectionIntegrator(solver=RelaxedFixedPointSolver(relaxation=0.3))
        sol = solve(prob, alg)

        @test sol.retcode == :Success
    end

    @testset "init" begin
        prob = make_harmonic_problem()
        integrator = init(prob)

        @test integrator.t == prob.tspan[1]
        @test integrator.t_end == prob.tspan[2]
        @test integrator.q == prob.q_initial
        @test integrator.p == prob.p_initial
        @test integrator.step_count == 0
        @test length(integrator.t_history) > 0
    end

    @testset "step!" begin
        prob = make_harmonic_problem(tspan=(0.0, 1.0), dt=0.1)
        integrator = init(prob)

        # First step
        result = step!(integrator)
        @test result == true  # More steps to go
        @test integrator.step_count == 1
        @test integrator.t ≈ 0.1
        @test integrator.t_history[2] ≈ 0.1

        # State changed from initial
        @test integrator.q != prob.q_initial || integrator.p != prob.p_initial

        # Continue stepping
        for _ in 1:5
            step!(integrator)
        end
        @test integrator.step_count == 6
        @test integrator.t ≈ 0.6
    end

    @testset "step! returns false at end" begin
        prob = make_harmonic_problem(tspan=(0.0, 0.05), dt=0.01)
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

    @testset "solve!" begin
        prob = make_harmonic_problem(tspan=(0.0, 0.5), dt=0.01)
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
        prob = make_harmonic_problem(tspan=(0.0, 0.5), dt=0.01)

        sol1 = solve(prob)

        integrator = init(prob)
        sol2 = solve!(integrator)

        @test length(sol1.t) == length(sol2.t)
        @test sol1.t ≈ sol2.t
        @test all(sol1.q .≈ sol2.q)
        @test all(sol1.p .≈ sol2.p)
    end

    @testset "params as Vector" begin
        prob = make_harmonic_problem()
        sol = solve(prob)
        @test sol.retcode == :Success
    end

    @testset "params as NamedTuple" begin
        H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])
        prob = WeberProblem(H, (0.0, 1.0), [1.0], [0.0]; params=(m=1.0, k=1.0), dt=0.01)
        sol = solve(prob)
        @test sol.retcode == :Success
    end

    @testset "History pre-allocation" begin
        prob = make_harmonic_problem(tspan=(0.0, 1.0), dt=0.01)
        integrator = init(prob)

        # Check history vectors are pre-allocated
        expected_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / prob.dt))
        @test length(integrator.t_history) >= expected_steps
        @test length(integrator.q_history) >= expected_steps
        @test length(integrator.p_history) >= expected_steps
    end

    @testset "Two-body system solve" begin
        prob = make_coulomb_problem(tspan=(0.0, 1.0), dt=0.01)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.q[1]) == 4  # 2 particles × 2 dims
        @test length(sol.p[1]) == 4
    end

    @testset "Weber system solve" begin
        prob = make_weber_problem(tspan=(0.0, 0.5), dt=0.001)
        sol = solve(prob)

        @test sol.retcode == :Success
        @test length(sol.t) > 1
    end

    @testset "Custom tolerance" begin
        prob = make_harmonic_problem(tspan=(0.0, 0.1), dt=0.01)
        # Tighter tolerance
        H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])
        prob_tight = WeberProblem(H, (0.0, 0.1), [1.0], [0.0];
            params=[1.0, 1.0], dt=0.01, convergence_tolerance=1e-14)
        sol = solve(prob_tight)
        @test sol.retcode == :Success
    end
end
