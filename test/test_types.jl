@testset "Types" begin
    @testset "SymmetricProjectionIntegrator" begin
        # Default construction
        alg = SymmetricProjectionIntegrator()
        @test alg.relaxation == 0.25

        # Custom relaxation
        alg2 = SymmetricProjectionIntegrator(relaxation = 0.5)
        @test alg2.relaxation == 0.5

        # Edge case: relaxation = 1.0 is valid
        @test SymmetricProjectionIntegrator(relaxation = 1.0).relaxation == 1.0

        # Validation: relaxation must be in (0, 1]
        @test_throws AssertionError SymmetricProjectionIntegrator(relaxation = 0.0)
        @test_throws AssertionError SymmetricProjectionIntegrator(relaxation = -0.1)
        @test_throws AssertionError SymmetricProjectionIntegrator(relaxation = 1.5)
    end

    @testset "HamiltonianSystem" begin
        # Basic construction (now purely symbolic)
        system = HamiltonianSystem(2, 2)

        @test system isa HamiltonianSystem
        @test system.n_particles == 2
        @test system.dims == 2
        @test system.degrees_of_freedom == 4

        # Symbolic fields exist
        @test !isnothing(system.hamiltonian_symbolic)
        @test !isnothing(system.dq_dt_symbolic)
        @test !isnothing(system.dp_dt_symbolic)

        # Compiled functions exist
        @test !isnothing(system.dq_dt_compiled)
        @test !isnothing(system.dp_dt_compiled)

        # Symbolic parameters exist: [m1, m2, q1, q2, c, kappa_1_2]
        @test length(system.param_symbols) == 6  # m1, m2, q1, q2, c, κ₁₂
        @test length(system.q_symbols) == 4  # x1, y1, x2, y2
        @test length(system.p_symbols) == 4  # px1, py1, px2, py2
    end

    @testset "HamiltonianSystem validation" begin
        # Valid cases
        @test HamiltonianSystem(1, 1) isa HamiltonianSystem
        @test HamiltonianSystem(3, 3) isa HamiltonianSystem

        # Invalid: dims must be 1, 2, or 3
        @test_throws AssertionError HamiltonianSystem(2, 0)
        @test_throws AssertionError HamiltonianSystem(2, 4)
    end

    @testset "HamiltonianSystem different configurations" begin
        # 1D, 1 particle (minimal)
        sys1 = HamiltonianSystem(1, 1)
        @test sys1.degrees_of_freedom == 1

        # 3D, 3 particles
        sys2 = HamiltonianSystem(3, 3)
        @test sys2.degrees_of_freedom == 9
        @test length(sys2.dq_dt_symbolic) == 9
        @test length(sys2.dp_dt_symbolic) == 9

        # 2D, 2 particles (common case)
        sys3 = HamiltonianSystem(2, 2)
        @test sys3.degrees_of_freedom == 4
    end

    @testset "HamiltonianProblem" begin
        system = HamiltonianSystem(2, 2)

        # Valid construction
        prob = HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
        )
        @test prob.tspan == (0.0, 1.0)
        @test prob.q_initial == [1.0, 0.0, -1.0, 0.0]
        @test prob.p_initial == [0.0, 0.1, 0.0, -0.1]
        @test prob.masses == [1.0, 0.5]
        @test prob.charges == [1.0, -1.0]
        @test prob.c == 4.0
        @test prob.dt == 0.01
        @test prob.convergence_tolerance == 1e-13  # default
        @test prob.maximum_iterations == 100  # default
        # params = [m1, m2, q1, q2, c, κ₁₂]; charges +1/-1 are unlike so κ=1 (Zöllner disabled)
        @test prob.params == [1.0, 0.5, 1.0, -1.0, 4.0, 1.0]
        @test prob.kappas == [1.0]  # unlike charges but Zöllner disabled → κ=1
        @test prob.zollner.enabled == false
        @test prob.zollner.a == 0.0
        @test prob.regularization.enabled == false
        @test prob.regularization.r_on === nothing
        @test prob.regularization.r_off === nothing
        @test prob.regularization.max_substeps == 512

        # Custom convergence_tolerance and maximum_iterations
        prob2 = HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
            convergence_tolerance = 1e-10,
            maximum_iterations = 50,
            regularization = RegularizationOptions(
                enabled = false,
                r_on = 0.2,
                r_off = 0.3,
                max_substeps = 64,
            ),
        )
        @test prob2.convergence_tolerance == 1e-10
        @test prob2.maximum_iterations == 50
        @test prob2.regularization.enabled == false
        @test prob2.regularization.r_on == 0.2
        @test prob2.regularization.r_off == 0.3
        @test prob2.regularization.max_substeps == 64

        # Validation errors
        @test_throws AssertionError HamiltonianProblem(
            system,
            (1.0, 0.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
        )  # tspan reversed
        @test_throws AssertionError HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = -0.01,
        )  # negative dt
        @test_throws AssertionError HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 2.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
        )  # q length mismatch
        @test_throws AssertionError HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
        )  # p length mismatch
        @test_throws AssertionError HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0],  # wrong length
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
        )
        @test_throws AssertionError HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, -0.5],  # negative mass
            charges = [1.0, -1.0],
            c = 4.0,
            dt = 0.01,
        )
        @test_throws AssertionError HamiltonianProblem(
            system,
            (0.0, 1.0),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 0.5],
            charges = [1.0, -1.0],
            c = -1.0,  # negative c
            dt = 0.01,
        )
    end

    @testset "HamiltonianSolution" begin
        prob = make_weber_problem(tspan = (0.0, 0.1))
        sol = solve(prob)

        @test sol isa HamiltonianSolution
        @test sol.retcode == :Success
        @test length(sol) > 1
        @test length(sol.t) == length(sol.q) == length(sol.p)
        @test sol.regularization isa RegularizationDiagnostics

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
        @test occursin("HamiltonianSolution", String(take!(io)))

        io = IOBuffer()
        show(io, MIME"text/plain"(), sol)
        str = String(take!(io))
        @test occursin("HamiltonianSolution", str)
        @test occursin("retcode", str)
    end

    @testset "HamiltonianIntegrator" begin
        prob = make_weber_problem()
        integrator = init(prob)

        @test integrator isa HamiltonianIntegrator
        @test integrator.t == prob.tspan[1]
        @test integrator.q == prob.q_initial
        @test integrator.p == prob.p_initial
        @test integrator.step_count == 0

        # show method
        io = IOBuffer()
        show(io, integrator)
        @test occursin("HamiltonianIntegrator", String(take!(io)))
    end

    @testset "HamiltonianSystem display" begin
        system = HamiltonianSystem(2, 2)

        # show(io, system)
        io = IOBuffer()
        show(io, system)
        str = String(take!(io))
        @test occursin("HamiltonianSystem", str)
        @test occursin("2 particles", str)
        @test occursin("2D", str)
        @test occursin("4 DOF", str)

        # show(io, MIME"text/plain", system)
        io = IOBuffer()
        show(io, MIME"text/plain"(), system)
        str = String(take!(io))
        @test occursin("HamiltonianSystem", str)
        @test occursin("Particles: 2", str)
        @test occursin("Dimensions: 2", str)
    end
end
