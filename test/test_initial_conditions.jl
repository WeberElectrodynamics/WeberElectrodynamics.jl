@testset "Initial-condition helpers" begin
    function total_com(q, masses, dims)
        out = zeros(dims)
        for i in eachindex(masses)
            for d = 1:dims
                out[d] += masses[i] * q[(i-1)*dims+d]
            end
        end
        return out
    end

    function total_momentum(p, n, dims)
        out = zeros(dims)
        for i = 1:n
            for d = 1:dims
                out[d] += p[(i-1)*dims+d]
            end
        end
        return out
    end

    function pair_rdot(q, p, masses, dims, i, j)
        rel_q = [q[(i-1)*dims+d] - q[(j-1)*dims+d] for d = 1:dims]
        rel_v = [p[(i-1)*dims+d] / masses[i] - p[(j-1)*dims+d] / masses[j] for d = 1:dims]
        return dot(rel_q, rel_v) / norm(rel_q)
    end

    @testset "center_of_mass_frame" begin
        framed = center_of_mass_frame(
            [2.0, 0.0, -1.0, 0.0],
            [3.0, 0.0, -1.0, 0.0],
            [1.0, 2.0],
            2,
        )

        @test total_com(framed.q, [1.0, 2.0], 2) ≈ zeros(2)
        @test total_momentum(framed.p, 2, 2) ≈ zeros(2) atol = 1e-14
    end

    @testset "two_body_initial_conditions" begin
        ic = two_body_initial_conditions(
            [1.0, 1.0],
            [1.0, -1.0];
            separation = 2.0,
            dims = 2,
            velocity_scale = 1.0,
        )

        @test ic.q ≈ [1.0, 0.0, -1.0, 0.0]
        @test ic.p ≈ [0.0, 0.5, 0.0, -0.5]
        @test total_com(ic.q, ic.masses, 2) ≈ zeros(2)
        @test total_momentum(ic.p, 2, 2) ≈ zeros(2)
        @test pair_rdot(ic.q, ic.p, ic.masses, 2, 1, 2) ≈ 0.0 atol = 1e-14
    end

    @testset "polygon_initial_conditions" begin
        ic = polygon_initial_conditions(5; radius = 2.0, mass = 1.5)

        @test length(ic.q) == 10
        @test length(ic.p) == 10
        @test ic.charges == [1.0, -1.0, 1.0, -1.0, 1.0]
        @test total_com(ic.q, ic.masses, 2) ≈ zeros(2) atol = 1e-12
        @test total_momentum(ic.p, 5, 2) ≈ zeros(2) atol = 1e-12
        for i = 1:5, j = (i+1):5
            @test pair_rdot(ic.q, ic.p, ic.masses, 2, i, j) ≈ 0.0 atol = 1e-12
        end
    end

    @testset "rigid_rotation_initial_conditions" begin
        q = [1.0, 0.0, 0.0, -0.5, sqrt(3) / 2, 0.2, -0.5, -sqrt(3) / 2, -0.4]
        masses = [1.0, 1.0, 1.0]
        ic = rigid_rotation_initial_conditions(
            q,
            masses;
            angular_velocity = [0.0, 0.0, 0.2],
            dims = 3,
            charges = [1.0, -1.0, 1.0],
        )

        @test total_com(ic.q, ic.masses, 3) ≈ zeros(3) atol = 1e-12
        @test total_momentum(ic.p, 3, 3) ≈ zeros(3) atol = 1e-12
        for i = 1:3, j = (i+1):3
            @test pair_rdot(ic.q, ic.p, ic.masses, 3, i, j) ≈ 0.0 atol = 1e-12
        end
    end
end
