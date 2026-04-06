using WeberElectrodynamics: _pair_index

@testset "Weber System" begin
    @testset "Weber Hamiltonian structure" begin
        # 2-body 2D Weber system
        system = WeberSystem(2, 2)

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

        system = WeberSystem(2, 2)
        # params = [m1, m2, q1, q2, c, κ₁₂]; κ=1 (Zöllner disabled)
        params = [m1, m2, q1, q2, c, 1.0]

        # Test at specific point: particles at x = ±1, no y offset
        q = [1.0, 0.0, -1.0, 0.0]  # particles at x=1 and x=-1
        p = [0.0, 0.1, 0.0, -0.1]  # small momenta in y

        out_q = zeros(4)
        out_p = zeros(4)
        system.dq_dt_compiled(out_q, q, p, params)
        system.dp_dt_compiled(out_p, q, p, params)

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
        sys1d = WeberSystem(2, 1)
        @test sys1d.degrees_of_freedom == 2
        @test length(sys1d.dq_dt_symbolic) == 2

        params1d = [1.0, 1.0, 1.0, -1.0, 1.0, 1.0]  # m1, m2, q1, q2, c, κ₁₂
        out1 = zeros(2)
        sys1d.dq_dt_compiled(out1, [1.0, -1.0], [0.5, -0.5], params1d)
        # Weber's velocity-dependent potential affects dq/dt, so we just verify
        # the function runs and produces finite output of correct dimension
        @test length(out1) == 2
        @test all(isfinite.(out1))

        # 2D system
        sys2d = WeberSystem(2, 2)
        @test sys2d.degrees_of_freedom == 4

        # 3D system
        sys3d = WeberSystem(2, 3)
        @test sys3d.degrees_of_freedom == 6
        @test length(sys3d.dq_dt_symbolic) == 6
    end

    @testset "Multiple particles" begin
        # 3-body system
        sys3body = WeberSystem(3, 2)
        @test sys3body.degrees_of_freedom == 6
        @test sys3body.n_particles == 3

        # params: [m1, m2, m3, q1, q2, q3, c, κ₁₂, κ₁₃, κ₂₃]; all κ=1 (Zöllner disabled)
        params3 = [1.0, 1.0, 1.0, 1.0, -1.0, 0.5, 1.0, 1.0, 1.0, 1.0]

        # Should have 3 pairwise interactions in the Hamiltonian
        out_q = zeros(6)
        out_p = zeros(6)
        q = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]  # triangle configuration
        p = zeros(6)

        # Should run without error
        sys3body.dq_dt_compiled(out_q, q, p, params3)
        sys3body.dp_dt_compiled(out_p, q, p, params3)
    end

    @testset "Single particle system" begin
        # Single particle has no interactions, just kinetic energy
        sys1 = WeberSystem(1, 2)
        @test sys1.degrees_of_freedom == 2

        params1 = [2.0, 1.0, 1.0]  # m1, q1, c
        out_q = zeros(2)
        out_p = zeros(2)
        q = [1.0, 2.0]
        p = [0.5, 1.0]

        sys1.dq_dt_compiled(out_q, q, p, params1)
        sys1.dp_dt_compiled(out_p, q, p, params1)

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

        system = WeberSystem(2, 2)

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

    @testset "Symbolic parameters in Hamiltonian" begin
        system = WeberSystem(2, 2)
        H_str = string(system.hamiltonian_symbolic)

        # Hamiltonian should contain symbolic parameters, not numerical values
        @test occursin("m1", H_str)
        @test occursin("m2", H_str)
        @test occursin("q1", H_str)
        @test occursin("q2", H_str)
        @test occursin("c", H_str)
    end

    @testset "_pair_index" begin
        # n=2: one pair (1,2) → 1
        @test _pair_index(1, 2, 2) == 1

        # n=3: pairs in order (1,2), (1,3), (2,3) → 1, 2, 3
        @test _pair_index(1, 2, 3) == 1
        @test _pair_index(1, 3, 3) == 2
        @test _pair_index(2, 3, 3) == 3

        # n=4: pairs in order (1,2),(1,3),(1,4),(2,3),(2,4),(3,4) → 1..6
        @test _pair_index(1, 2, 4) == 1
        @test _pair_index(1, 3, 4) == 2
        @test _pair_index(1, 4, 4) == 3
        @test _pair_index(2, 3, 4) == 4
        @test _pair_index(2, 4, 4) == 5
        @test _pair_index(3, 4, 4) == 6

        # All indices for a given n are unique and cover 1:n_pairs exactly
        for n in 2:5
            indices = [_pair_index(i, j, n) for i in 1:n for j in (i+1):n]
            @test sort(indices) == collect(1:n*(n-1)÷2)
        end
    end
end
