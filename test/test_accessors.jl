@testset "Accessor API" begin
    sys = HamiltonianSystem(3, 2)

    @testset "HamiltonianSystem shape accessors" begin
        @test n_particles(sys) == 3
        @test dims(sys) == 2
        @test degrees_of_freedom(sys) == 6
    end

    @testset "HamiltonianProblem parameter accessors" begin
        prob = HamiltonianProblem(
            sys, (0.0, 1.0), zeros(6), zeros(6);
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

        # params layout: [m₁…mₙ, q₁…qₙ, c, κ₁₂, κ₁₃, κ₂₃]
        @test length(params(prob)) == 2*3 + 1 + 3
        @test params(prob)[1:3] == masses(prob)
        @test params(prob)[4:6] == charges(prob)
        @test params(prob)[7] == speed_of_light(prob)
        @test params(prob)[8:10] == kappas(prob)
    end

    @testset "zollner injects κ accessor" begin
        prob = HamiltonianProblem(
            sys, (0.0, 1.0), zeros(6), zeros(6);
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

    @testset "regularization/zollner options accessors" begin
        reg = RegularizationOptions(enabled = true, r_on_factor = 0.2, r_off_factor = 0.3)
        zopts = ZollnerOptions(enabled = true, a = 0.05)
        prob = HamiltonianProblem(
            sys, (0.0, 1.0), zeros(6), zeros(6);
            masses = [1.0, 1.0, 1.0],
            charges = [1.0, -1.0, 1.0],
            c = 10.0,
            dt = 0.01,
            regularization = reg,
            zollner = zopts,
        )
        @test regularization(prob) === reg
        @test zollner(prob) === zopts
    end
end
