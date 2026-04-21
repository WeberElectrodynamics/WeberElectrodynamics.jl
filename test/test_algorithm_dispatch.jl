@testset "Algorithm dispatch hook" begin
    # Dummy algorithm with no _allocate_cache / _step_core! methods —
    # should hit the default error path at init time so unsupported
    # algorithms fail loudly.
    struct _UnsupportedAlgorithm <: WeberElectrodynamics.HamiltonianAlgorithm end

    @testset "SymmetricProjectionIntegrator <: HamiltonianAlgorithm" begin
        @test SymmetricProjectionIntegrator <: WeberElectrodynamics.HamiltonianAlgorithm
    end

    @testset "dispatch methods exist for the built-in algorithm" begin
        @test hasmethod(
            WeberElectrodynamics._step_core!,
            Tuple{
                WeberElectrodynamics.HamiltonianIntegrator,
                SymmetricProjectionIntegrator,
                Float64,
            },
        )
        @test hasmethod(
            WeberElectrodynamics._allocate_cache,
            Tuple{HamiltonianProblem,SymmetricProjectionIntegrator},
        )
    end

    sys = HamiltonianSystem(2, 2)
    prob = HamiltonianProblem(
        sys,
        (0.0, 0.05),
        [1.0, 0.0, -1.0, 0.0],
        [0.0, 0.1, 0.0, -0.1];
        masses = [1.0, 1.0],
        charges = [1.0, -1.0],
        c = 100.0,
        dt = 0.01,
    )

    @testset "supported algorithm solves" begin
        sol = solve(prob, SymmetricProjectionIntegrator())
        @test sol.retcode === :Success
    end

    @testset "unsupported algorithm raises at init" begin
        @test_throws ArgumentError solve(prob, _UnsupportedAlgorithm())
    end
end
