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

    @testset "RegularizationOptions validation" begin
        opts = RegularizationOptions()
        @test opts.enabled == false
        @test opts.backend == WeberElectrodynamics.REG_BACKEND_LIFTED

        @test_throws AssertionError RegularizationOptions(r_on = 0.0)
        @test_throws AssertionError RegularizationOptions(r_off = -1.0)
        @test_throws AssertionError RegularizationOptions(r_on = 0.2, r_off = 0.1)
        @test_throws AssertionError RegularizationOptions(r_on_factor = 0.0)
        @test_throws AssertionError RegularizationOptions(r_off_factor = 0.0)
        @test_throws AssertionError RegularizationOptions(max_substeps = 0)
        @test_throws AssertionError RegularizationOptions(constraint_tolerance = 0.0)
        @test_throws AssertionError RegularizationOptions(g_floor = 0.0)
        @test_throws AssertionError RegularizationOptions(backend = :unknown)
        @test_throws AssertionError RegularizationOptions(collision_bounce_radius = -0.1)
    end

    @testset "HamiltonianSystem" begin
        # Basic construction (now purely symbolic)
        system = HamiltonianSystem(2, 2)

        @test system isa HamiltonianSystem
        @test system.n_particles == 2
        @test system.dims == 2
        @test system.degrees_of_freedom == 4

        # The default Weber system is analytic: its exact canonical Hamiltonian
        # needs a per-evaluation pair solve and carries no symbolic expression.
        @test !has_symbolic_hamiltonian(system)
        @test isnothing(system.hamiltonian_symbolic)
        @test isnothing(system.dq_dt_symbolic)
        @test isnothing(system.dp_dt_symbolic)

        # Compiled functions exist
        @test !isnothing(system.dq_dt_compiled)
        @test !isnothing(system.dp_dt_compiled)

        # Symbolic parameters exist: [m1, m2, q1, q2, c].
        @test length(system.param_symbols) == 5  # m1, m2, q1, q2, c
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
        @test length(sys2.q_symbols) == 9
        @test length(sys2.p_symbols) == 9

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
        @test masses(prob) == [1.0, 0.5]
        @test charges(prob) == [1.0, -1.0]
        @test speed_of_light(prob) == 4.0
        @test prob.dt == 0.01
        @test prob.convergence_tolerance == 1e-13  # default
        @test prob.maximum_iterations == 100  # default
        @test fieldtype(typeof(prob), :system) == typeof(system)
        # params = [m1, m2, q1, q2, c]
        @test params(prob) == [1.0, 0.5, 1.0, -1.0, 4.0]

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
        )
        @test prob2.convergence_tolerance == 1e-10
        @test prob2.maximum_iterations == 50

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
        @test fieldtype(typeof(integrator), :prob) == typeof(prob)

        # show method
        io = IOBuffer()
        show(io, integrator)
        @test occursin("HamiltonianIntegrator", String(take!(io)))
    end

    @testset "_with_tspan clones a problem with a new tspan" begin
        prob = HamiltonianProblem(
            HamiltonianSystem(2, 2),
            (0.0, 1.5),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 2.0],
            charges = [1.0, -1.0],
            c = 10.0,
            dt = 0.01,
        )
        new_prob = WeberElectrodynamics._with_tspan(prob, (0.5, 4.0))

        # tspan replaced.
        @test new_prob.tspan == (0.5, 4.0)

        # All other fields are equal-by-value to the original.
        @test new_prob.system === prob.system
        @test new_prob.q_initial == prob.q_initial
        @test new_prob.p_initial == prob.p_initial
        @test params(new_prob) == params(prob)
        @test new_prob.dt == prob.dt
        @test new_prob.convergence_tolerance == prob.convergence_tolerance
        @test new_prob.maximum_iterations == prob.maximum_iterations

        # ICs and params are *copies*, not aliases — mutating the
        # clone must not propagate back to the original.
        @test new_prob.q_initial !== prob.q_initial
        @test new_prob.p_initial !== prob.p_initial
        @test params(new_prob) !== params(prob)
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
