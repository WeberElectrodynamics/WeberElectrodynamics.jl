@testset "Types" begin
    @testset "WeberAlgorithm & SymmetricProjection" begin
        # Default construction
        alg = SymmetricProjection()
        @test alg isa WeberAlgorithm
        @test alg.solver isa RelaxedFixedPoint
        @test alg.solver.relaxation == 0.25

        # Custom solver
        custom_solver = RelaxedFixedPoint(relaxation=0.5)
        alg2 = SymmetricProjection(solver=custom_solver)
        @test alg2.solver.relaxation == 0.5
    end

    @testset "RelaxedFixedPoint" begin
        # Default relaxation
        @test RelaxedFixedPoint().relaxation == 0.25

        # Custom relaxation
        @test RelaxedFixedPoint(relaxation=0.5).relaxation == 0.5

        # Edge case: relaxation = 1.0 is valid
        @test RelaxedFixedPoint(relaxation=1.0).relaxation == 1.0

        # Validation: relaxation must be in (0, 1]
        @test_throws AssertionError RelaxedFixedPoint(relaxation=0.0)
        @test_throws AssertionError RelaxedFixedPoint(relaxation=-0.1)
        @test_throws AssertionError RelaxedFixedPoint(relaxation=1.5)

        # show method
        io = IOBuffer()
        show(io, RelaxedFixedPoint(relaxation=0.3))
        @test occursin("0.3", String(take!(io)))
    end

    @testset "WeberHamiltonian" begin
        H = build_hamiltonian(harmonic_oscillator_H, 1, 1; param_names=[:m, :k])

        @test H isa WeberHamiltonian
        @test H.n_dof == 1
        @test H.param_names == [:m, :k]
        @test !isnothing(H.H_sym)
        @test !isnothing(H.qdot_sym)
        @test !isnothing(H.pdot_sym)
        @test !isnothing(H.qdot_func)
        @test !isnothing(H.pdot_func)
    end

    @testset "WeberProblem" begin
        H = build_hamiltonian(harmonic_oscillator_H, 1, 1; param_names=[:m, :k])

        # Valid construction
        prob = WeberProblem(H, (0.0, 1.0), [1.0], [0.0]; params=[1.0, 1.0], dt=0.01)
        @test prob.tspan == (0.0, 1.0)
        @test prob.q₀ == [1.0]
        @test prob.p₀ == [0.0]
        @test prob.dt == 0.01
        @test prob.tolerance == 1e-13  # default
        @test prob.max_iterations == 100  # default

        # Custom tolerance and max_iterations
        prob2 = WeberProblem(H, (0.0, 1.0), [1.0], [0.0];
            params=[1.0, 1.0], dt=0.01,
            tolerance=1e-10, max_iterations=50)
        @test prob2.tolerance == 1e-10
        @test prob2.max_iterations == 50

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
        @test integrator.q == prob.q₀
        @test integrator.p == prob.p₀
        @test integrator.step_count == 0

        # show method
        io = IOBuffer()
        show(io, integrator)
        @test occursin("WeberIntegrator", String(take!(io)))
    end

    @testset "IntegratorBuffers" begin
        params_vec = [1.0, 2.0]
        d = 3
        buffers = IntegratorBuffers(d, params_vec)

        @test buffers.d == 3
        @test size(buffers.A) == (2d, 4d)  # 6 × 12
        @test length(buffers.Z_current) == 4d  # 12
        @test length(buffers.Z_post_phi) == 4d
        @test length(buffers.Z_result) == 4d
        @test length(buffers.qs_buf) == d
        @test length(buffers.xs_buf) == d
        @test length(buffers.ps_buf) == d
        @test length(buffers.ys_buf) == d
        @test length(buffers.ATμ) == 4d
        @test length(buffers.μ) == 2d
        @test length(buffers.μ_old) == 2d
        @test length(buffers.f_val) == 2d
        @test buffers.params_vec === params_vec
    end

    @testset "NonlinearSolveError" begin
        err = NonlinearSolveError(50, 1e-12, 1e-8, 100, 0.5, "RelaxedFixedPoint")

        @test err.iterations == 50
        @test err.tolerance == 1e-12
        @test err.final_residual == 1e-8
        @test err.step == 100
        @test err.time == 0.5
        @test err.solver_name == "RelaxedFixedPoint"

        # showerror
        io = IOBuffer()
        showerror(io, err)
        msg = String(take!(io))
        @test occursin("NonlinearSolveError", msg)
        @test occursin("RelaxedFixedPoint", msg)
        @test occursin("step 100", msg)

        # Deprecated alias
        @test NewtonConvergenceError === NonlinearSolveError
    end
end
