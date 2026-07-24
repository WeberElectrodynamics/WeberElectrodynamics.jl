@testset "Accessor API" begin
    sys = HamiltonianSystem(3, 2)

    @testset "HamiltonianSystem shape accessors" begin
        @test n_particles(sys) == 3
        @test dims(sys) == 2
        @test degrees_of_freedom(sys) == 6
        @test n_pairs(sys) == 3
        @test pair_indices(sys) == [(1, 2), (1, 3), (2, 3)]
    end

    @testset "HamiltonianProblem parameter accessors" begin
        prob = HamiltonianProblem(
            sys,
            (0.0, 1.0),
            zeros(6),
            zeros(6);
            masses = [1.0, 2.0, 3.0],
            charges = [1.0, -1.0, 0.5],
            c = 10.0,
            dt = 0.01,
        )

        @test n_particles(prob) == 3
        @test dims(prob) == 2
        @test masses(prob) == [1.0, 2.0, 3.0]
        @test charges(prob) == [1.0, -1.0, 0.5]
        @test speed_of_light(prob) == 10.0
        @test n_pairs(prob) == 3
        @test pair_indices(prob) == [(1, 2), (1, 3), (2, 3)]

        # params layout: [m₁…mₙ, q₁…qₙ, c]
        @test length(params(prob)) == 2*3 + 1
        @test params(prob)[1:3] == masses(prob)
        @test params(prob)[4:6] == charges(prob)
        @test params(prob)[7] == speed_of_light(prob)
    end

    @testset "accessor view aliasing" begin
        prob = HamiltonianProblem(
            sys,
            (0.0, 1.0),
            zeros(6),
            zeros(6);
            masses = [1.0, 2.0, 3.0],
            charges = [1.0, -1.0, 0.5],
            c = 10.0,
            dt = 0.01,
        )
        @test parent(masses(prob)) === params(prob)
        @test parent(charges(prob)) === params(prob)
    end
end
