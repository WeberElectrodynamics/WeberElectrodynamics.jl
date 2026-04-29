@testset "Accessor API" begin
    sys = HamiltonianSystem(3, 2)

    @testset "HamiltonianSystem shape accessors" begin
        @test n_particles(sys) == 3
        @test dims(sys) == 2
        @test degrees_of_freedom(sys) == 6
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
        @test kappas(prob) == ones(3)

        # params layout: [m₁…mₙ, q₁…qₙ, c]; kappas stored separately
        @test length(params(prob)) == 2*3 + 1
        @test params(prob)[1:3] == masses(prob)
        @test params(prob)[4:6] == charges(prob)
        @test params(prob)[7] == speed_of_light(prob)
        @test length(kappas(prob)) == 3

        # kappa(prob, i, j) per-pair accessor
        @test kappa(prob, 1, 2) == 1.0
        @test kappa(prob, 1, 3) == 1.0
        @test kappa(prob, 2, 3) == 1.0
    end

    @testset "zollner injects κ accessor" begin
        prob = HamiltonianProblem(
            sys,
            (0.0, 1.0),
            zeros(6),
            zeros(6);
            masses = [1.0, 1.0, 1.0],
            charges = [1.0, -1.0, 1.0],
            c = 10.0,
            dt = 0.01,
            zollner = ZollnerOptions(enabled = true, a = 0.1),
        )
        # pairs (1,2) and (2,3) are unlike-sign → κ = 1.1; (1,3) like-sign → κ = 1
        κ = kappas(prob)
        @test κ ≈ [1.1, 1.0, 1.1]
    end

    @testset "kappa argument validation" begin
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
        # Valid index ordering — no throw.
        @test kappa(prob, 1, 2) == 1.0
        @test kappa(prob, 2, 3) == 1.0
        # Reversed pairs are canonicalized to the same κ.
        @test kappa(prob, 2, 1) == kappa(prob, 1, 2)
        @test kappa(prob, 3, 2) == kappa(prob, 2, 3)
        @test_throws AssertionError kappa(prob, 1, 1)
        # Out-of-range indices.
        @test_throws AssertionError kappa(prob, 0, 1)
        @test_throws AssertionError kappa(prob, 1, 4)
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
        # masses(prob) and charges(prob) are O(1) views into prob.params —
        # parent identity proves they share memory rather than copy.
        @test parent(masses(prob)) === params(prob)
        @test parent(charges(prob)) === params(prob)
        # kappas lives on its own field, NOT a slice of params.
        @test kappas(prob) === prob.kappas
        @test parent(kappas(prob)) !== params(prob)
    end

end
