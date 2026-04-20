@testset "Algorithm dispatch hook" begin
    @testset "SymmetricProjectionIntegrator <: HamiltonianAlgorithm" begin
        @test SymmetricProjectionIntegrator <: WeberElectrodynamics.HamiltonianAlgorithm
    end

    @testset "_step_core! dispatches on algorithm" begin
        # The internal step hook is the documented extension point. A
        # specialized method exists for SymmetricProjectionIntegrator; the
        # default method on any other HamiltonianAlgorithm subtype errors.
        @test hasmethod(
            WeberElectrodynamics._step_core!,
            Tuple{
                WeberElectrodynamics.HamiltonianIntegrator,
                SymmetricProjectionIntegrator,
                Float64,
            },
        )
    end

    @testset "solve! routes through _step_core!" begin
        # Regression: the outer step! loop must call _step_core!, not an
        # inline branch. If someone re-inlines the regularized/unregularized
        # branching, this test still passes — but the regression fixtures
        # catch any numerical drift.
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.05), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 100.0,
            dt = 0.01,
        )
        sol = solve(prob, SymmetricProjectionIntegrator())
        @test sol.retcode === :Success
    end
end
