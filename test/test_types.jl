@testset "Types" begin
    @testset "SymmetricProjectionIntegrator" begin
        # Default construction
        alg = SymmetricProjectionIntegrator()
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

    @testset "WeberSystem" begin
        # Basic construction
        system = WeberSystem(2, 2; masses=[1.0, 0.5], charges=[1.0, -1.0], c=4.0)

        @test system isa WeberSystem
        @test system.n_particles == 2
        @test system.dims == 2
        @test system.degrees_of_freedom == 4
        @test system.masses == [1.0, 0.5]
        @test system.charges == [1.0, -1.0]
        @test system.c == 4.0

        # Symbolic fields exist
        @test !isnothing(system.hamiltonian_symbolic)
        @test !isnothing(system.dq_dt_symbolic)
        @test !isnothing(system.dp_dt_symbolic)

        # Compiled functions exist
        @test !isnothing(system.dq_dt_compiled)
        @test !isnothing(system.dp_dt_compiled)
    end

    @testset "WeberSystem validation" begin
        # Valid cases
        @test WeberSystem(1, 1; masses=[1.0], charges=[1.0], c=1.0) isa WeberSystem
        @test WeberSystem(3, 3; masses=[1.0, 2.0, 3.0], charges=[1.0, -1.0, 0.5], c=299792458.0) isa WeberSystem

        # Invalid: dims must be 1, 2, or 3
        @test_throws AssertionError WeberSystem(2, 0; masses=[1.0, 1.0], charges=[1.0, -1.0], c=1.0)
        @test_throws AssertionError WeberSystem(2, 4; masses=[1.0, 1.0], charges=[1.0, -1.0], c=1.0)

        # Invalid: masses and charges must match n_particles
        @test_throws AssertionError WeberSystem(2, 2; masses=[1.0], charges=[1.0, -1.0], c=1.0)
        @test_throws AssertionError WeberSystem(2, 2; masses=[1.0, 1.0], charges=[1.0], c=1.0)

        # Invalid: masses must be positive
        @test_throws AssertionError WeberSystem(2, 2; masses=[1.0, -1.0], charges=[1.0, -1.0], c=1.0)
        @test_throws AssertionError WeberSystem(2, 2; masses=[1.0, 0.0], charges=[1.0, -1.0], c=1.0)

        # Invalid: c must be positive
        @test_throws AssertionError WeberSystem(2, 2; masses=[1.0, 1.0], charges=[1.0, -1.0], c=0.0)
        @test_throws AssertionError WeberSystem(2, 2; masses=[1.0, 1.0], charges=[1.0, -1.0], c=-1.0)

        # Valid: charges can be zero, negative, or positive
        @test WeberSystem(2, 2; masses=[1.0, 1.0], charges=[0.0, 0.0], c=1.0) isa WeberSystem
        @test WeberSystem(2, 2; masses=[1.0, 1.0], charges=[-1.0, -1.0], c=1.0) isa WeberSystem
    end

    @testset "WeberSystem different configurations" begin
        # 1D, 1 particle (minimal)
        sys1 = WeberSystem(1, 1; masses=[1.0], charges=[1.0], c=1.0)
        @test sys1.degrees_of_freedom == 1

        # 3D, 3 particles
        sys2 = WeberSystem(3, 3; masses=[1.0, 2.0, 3.0], charges=[1.0, -1.0, 0.5], c=1.0)
        @test sys2.degrees_of_freedom == 9
        @test length(sys2.dq_dt_symbolic) == 9
        @test length(sys2.dp_dt_symbolic) == 9

        # 2D, 2 particles (common case)
        sys3 = WeberSystem(2, 2; masses=[1.0, 0.1], charges=[-0.1, 0.1], c=4.0)
        @test sys3.degrees_of_freedom == 4
    end

    @testset "WeberProblem" begin
        system = WeberSystem(2, 2; masses=[1.0, 0.5], charges=[1.0, -1.0], c=4.0)

        # Valid construction
        prob = WeberProblem(system, (0.0, 1.0), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1]; dt=0.01)
        @test prob.tspan == (0.0, 1.0)
        @test prob.q_initial == [1.0, 0.0, -1.0, 0.0]
        @test prob.p_initial == [0.0, 0.1, 0.0, -0.1]
        @test prob.dt == 0.01
        @test prob.convergence_tolerance == 1e-13  # default
        @test prob.maximum_iterations == 100  # default

        # Custom convergence_tolerance and maximum_iterations
        prob2 = WeberProblem(system, (0.0, 1.0), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            dt=0.01, convergence_tolerance=1e-10, maximum_iterations=50)
        @test prob2.convergence_tolerance == 1e-10
        @test prob2.maximum_iterations == 50

        # Validation errors
        @test_throws AssertionError WeberProblem(system, (1.0, 0.0), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1]; dt=0.01)  # tspan reversed
        @test_throws AssertionError WeberProblem(system, (0.0, 1.0), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1]; dt=-0.01)  # negative dt
        @test_throws AssertionError WeberProblem(system, (0.0, 1.0), [1.0, 2.0], [0.0, 0.1, 0.0, -0.1]; dt=0.01)  # q length mismatch
        @test_throws AssertionError WeberProblem(system, (0.0, 1.0), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1]; dt=0.01)  # p length mismatch
    end

    @testset "WeberSolution" begin
        prob = make_weber_problem(tspan=(0.0, 0.1))
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
        prob = make_weber_problem()
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

    @testset "WeberSystem display" begin
        system = WeberSystem(2, 2; masses=[1.0, 0.5], charges=[1.0, -1.0], c=4.0)

        # show(io, system)
        io = IOBuffer()
        show(io, system)
        str = String(take!(io))
        @test occursin("WeberSystem", str)
        @test occursin("2 particles", str)
        @test occursin("2D", str)
        @test occursin("4 DOF", str)

        # show(io, MIME"text/plain", system)
        io = IOBuffer()
        show(io, MIME"text/plain"(), system)
        str = String(take!(io))
        @test occursin("WeberSystem", str)
        @test occursin("Particles: 2", str)
        @test occursin("Dimensions: 2", str)
        @test occursin("Masses:", str)
        @test occursin("Charges:", str)
    end
end
