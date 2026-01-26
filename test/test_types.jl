@testset "Types" begin
    @testset "WeberAlgorithm & SymmetricProjectionIntegrator" begin
        # Default construction
        alg = SymmetricProjectionIntegrator()
        @test alg isa WeberAlgorithm
        @test alg.relaxation == 0.25

        # Custom relaxation
        alg2 = SymmetricProjectionIntegrator(relaxation=0.5)
        @test alg2.relaxation == 0.5

        # Edge case: relaxation = 1.0 is valid
        @test SymmetricProjectionIntegrator(relaxation=1.0).relaxation == 1.0

        # Validation: relaxation must be in (0, 1]
        @test_throws AssertionError SymmetricProjectionIntegrator(relaxation=0.0)
        @test_throws AssertionError SymmetricProjectionIntegrator(relaxation=-0.1)
        @test_throws AssertionError SymmetricProjectionIntegrator(relaxation=1.5)
    end

    @testset "WeberHamiltonian" begin
        H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])

        @test H isa WeberHamiltonian
        @test H.degrees_of_freedom == 1
        @test H.parameter_names == [:m, :k]
        @test !isnothing(H.hamiltonian_symbolic)
        @test !isnothing(H.dq_dt_symbolic)
        @test !isnothing(H.dp_dt_symbolic)
        @test !isnothing(H.dq_dt_compiled)
        @test !isnothing(H.dp_dt_compiled)
    end

    @testset "WeberProblem" begin
        H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])

        # Valid construction
        prob = WeberProblem(H, (0.0, 1.0), [1.0], [0.0]; params=[1.0, 1.0], dt=0.01)
        @test prob.tspan == (0.0, 1.0)
        @test prob.q_initial == [1.0]
        @test prob.p_initial == [0.0]
        @test prob.dt == 0.01
        @test prob.convergence_tolerance == 1e-13  # default
        @test prob.maximum_iterations == 100  # default

        # Custom convergence_tolerance and maximum_iterations
        prob2 = WeberProblem(H, (0.0, 1.0), [1.0], [0.0];
            params=[1.0, 1.0], dt=0.01,
            convergence_tolerance=1e-10, maximum_iterations=50)
        @test prob2.convergence_tolerance == 1e-10
        @test prob2.maximum_iterations == 50

        # NamedTuple params
        prob_nt = WeberProblem(H, (0.0, 1.0), [1.0], [0.0]; params=(m=1.0, k=1.0), dt=0.01)
        @test prob_nt isa WeberProblem

        # Validation errors
        @test_throws AssertionError WeberProblem(H, (1.0, 0.0), [1.0], [0.0]; params=[1.0, 1.0], dt=0.01)  # tspan reversed
        @test_throws AssertionError WeberProblem(H, (0.0, 1.0), [1.0], [0.0]; params=[1.0, 1.0], dt=-0.01)  # negative dt
        @test_throws AssertionError WeberProblem(H, (0.0, 1.0), [1.0, 2.0], [0.0]; params=[1.0, 1.0], dt=0.01)  # q/p mismatch
        @test_throws AssertionError WeberProblem(H, (0.0, 1.0), [1.0, 2.0], [0.0, 0.0]; params=[1.0, 1.0], dt=0.01)  # DOF mismatch
    end

    @testset "WeberSolution" begin
        prob = make_harmonic_problem(tspan=(0.0, 0.1))
        sol = solve(prob)

        @test sol isa WeberSolution
        @test sol.retcode == :Success
        @test length(sol) > 1
        @test length(sol.t) == length(sol.q) == length(sol.p)

        # Indexing
        t, q, p = sol[1]
        @test t == sol.t[1]
        @test q == sol.q[1]
        @test p == sol.p[1]

        # firstindex/lastindex
        @test firstindex(sol) == 1
        @test lastindex(sol) == length(sol)

        # Iteration
        count = 0
        for (t_i, q_i, p_i) in sol
            count += 1
            @test t_i isa Float64
            @test q_i isa Vector{Float64}
            @test p_i isa Vector{Float64}
        end
        @test count == length(sol)

        # show methods
        io = IOBuffer()
        show(io, sol)
        @test occursin("WeberSolution", String(take!(io)))

        io = IOBuffer()
        show(io, MIME"text/plain"(), sol)
        str = String(take!(io))
        @test occursin("WeberSolution", str)
        @test occursin("retcode", str)
    end

    @testset "WeberIntegrator" begin
        prob = make_harmonic_problem()
        integrator = init(prob)

        @test integrator isa WeberIntegrator
        @test integrator.t == prob.tspan[1]
        @test integrator.q == prob.q_initial
        @test integrator.p == prob.p_initial
        @test integrator.step_count == 0

        # show method
        io = IOBuffer()
        show(io, integrator)
        @test occursin("WeberIntegrator", String(take!(io)))
    end

    @testset "NonlinearSolveError" begin
        err = NonlinearSolveError(50, 1e-12, 1e-8, 100, 0.5, "RelaxedFixedPointSolver")

        @test err.iterations == 50
        @test err.convergence_tolerance == 1e-12
        @test err.final_residual == 1e-8
        @test err.step == 100
        @test err.time == 0.5
        @test err.solver_name == "RelaxedFixedPointSolver"

        # showerror
        io = IOBuffer()
        showerror(io, err)
        msg = String(take!(io))
        @test occursin("NonlinearSolveError", msg)
        @test occursin("RelaxedFixedPointSolver", msg)
        @test occursin("step 100", msg)
    end
end
