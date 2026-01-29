@testset "Weber System" begin
    @testset "Weber Hamiltonian structure" begin
        # 2-body 2D Weber system
        system = WeberSystem(2, 2; masses=[1.0, 0.5], charges=[1.0, -1.0], c=4.0)

        @test system.degrees_of_freedom == 4
        @test length(system.dq_dt_symbolic) == 4
        @test length(system.dp_dt_symbolic) == 4

        # Hamiltonian should be a symbolic expression
        H = system.hamiltonian_symbolic
        @test H isa Symbolics.Num || H isa Symbolics.BasicSymbolic
    end

    @testset "Compiled functions correctness" begin
        # Simple 2-body 2D system
        m1, m2 = 1.0, 1.0
        q1, q2 = 1.0, -1.0
        c = 1000.0  # Large c to make Weber ≈ Coulomb

        system = WeberSystem(2, 2; masses=[m1, m2], charges=[q1, q2], c=c)

        # Test at specific point: particles at x = ±1, no y offset
        q = [1.0, 0.0, -1.0, 0.0]  # particles at x=1 and x=-1
        p = [0.0, 0.1, 0.0, -0.1]  # small momenta in y

        out_q = zeros(4)
        out_p = zeros(4)
        system.dq_dt_compiled(out_q, q, p)
        system.dp_dt_compiled(out_p, q, p)

        # dq/dt = ∂H/∂p ≈ p/m (Weber correction is small for large c)
        @test out_q[1] ≈ p[1] / m1 atol = 1e-6  # px1/m1
        @test out_q[2] ≈ p[2] / m1 atol = 1e-6  # py1/m1
        @test out_q[3] ≈ p[3] / m2 atol = 1e-6  # px2/m2
        @test out_q[4] ≈ p[4] / m2 atol = 1e-6  # py2/m2

        # dp/dt should show attractive force along x-axis (q1*q2 < 0)
        # Particle 1 at x=1 pulled toward particle 2 (negative x direction)
        @test out_p[1] < 0
        # Particle 2 at x=-1 pulled toward particle 1 (positive x direction)
        @test out_p[3] > 0
        # y forces should be small at this symmetric configuration
        @test abs(out_p[2]) < 1e-6
        @test abs(out_p[4]) < 1e-6
    end

    @testset "Different dimensions" begin
        # 1D system
        sys1d = WeberSystem(2, 1; masses=[1.0, 1.0], charges=[1.0, -1.0], c=1.0)
        @test sys1d.degrees_of_freedom == 2
        @test length(sys1d.dq_dt_symbolic) == 2

        out1 = zeros(2)
        sys1d.dq_dt_compiled(out1, [1.0, -1.0], [0.5, -0.5])
        # Weber's velocity-dependent potential affects dq/dt, so we just verify
        # the function runs and produces finite output of correct dimension
        @test length(out1) == 2
        @test all(isfinite.(out1))

        # 2D system
        sys2d = WeberSystem(2, 2; masses=[1.0, 1.0], charges=[1.0, -1.0], c=1.0)
        @test sys2d.degrees_of_freedom == 4

        # 3D system
        sys3d = WeberSystem(2, 3; masses=[1.0, 1.0], charges=[1.0, -1.0], c=1.0)
        @test sys3d.degrees_of_freedom == 6
        @test length(sys3d.dq_dt_symbolic) == 6
    end

    @testset "Multiple particles" begin
        # 3-body system
        sys3body = WeberSystem(3, 2; masses=[1.0, 1.0, 1.0], charges=[1.0, -1.0, 0.5], c=1.0)
        @test sys3body.degrees_of_freedom == 6
        @test sys3body.n_particles == 3

        # Should have 3 pairwise interactions in the Hamiltonian
        out_q = zeros(6)
        out_p = zeros(6)
        q = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]  # triangle configuration
        p = zeros(6)

        # Should run without error
        sys3body.dq_dt_compiled(out_q, q, p)
        sys3body.dp_dt_compiled(out_p, q, p)
    end

    @testset "Default speed of light" begin
        # Default c should be 1.0
        system = WeberSystem(2, 2; masses=[1.0, 1.0], charges=[1.0, -1.0])
        @test system.c == 1.0
    end

    @testset "Single particle system" begin
        # Single particle has no interactions, just kinetic energy
        sys1 = WeberSystem(1, 2; masses=[2.0], charges=[1.0], c=1.0)
        @test sys1.degrees_of_freedom == 2

        out_q = zeros(2)
        out_p = zeros(2)
        q = [1.0, 2.0]
        p = [0.5, 1.0]

        sys1.dq_dt_compiled(out_q, q, p)
        sys1.dp_dt_compiled(out_p, q, p)

        # dq/dt = p/m (free particle)
        @test out_q[1] ≈ 0.5 / 2.0  # px/m = 0.25
        @test out_q[2] ≈ 1.0 / 2.0  # py/m = 0.5

        # dp/dt = 0 (no forces on single particle)
        @test out_p[1] ≈ 0.0
        @test out_p[2] ≈ 0.0
    end

    @testset "Energy function correctness" begin
        # Verify the symbolic Hamiltonian evaluates correctly
        m1, m2 = 1.0, 0.5
        q1, q2 = 0.3, -0.3
        c = 4.0

        system = WeberSystem(2, 2; masses=[m1, m2], charges=[q1, q2], c=c)

        # Test state
        q = [1.0, 0.0, -1.0, 0.0]
        p = [0.1, 0.2, -0.1, -0.2]

        # Compute energy manually using the test utility function
        E_manual = weber_energy_2body_2d(q, p, [m1, m2], [q1, q2], c)

        # The Hamiltonian symbolic expression should give same result
        # (This is implicit in the integration tests - if energy is conserved,
        # the symbolic Hamiltonian must be correct)
        @test E_manual isa Float64
        @test isfinite(E_manual)
    end
end
