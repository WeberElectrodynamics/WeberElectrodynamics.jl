using Aqua

@testset "Aqua quality checks" begin
    Aqua.test_all(
        WeberElectrodynamics;
        # Symbolics / Latexify can trigger method ambiguities outside our control.
        ambiguities = false,
        # Piracies from overriding Base.show for our types are fine.
        piracies = false,
        # Known-flaky precompilation check inside Pkg.test sandbox.
        persistent_tasks = false,
    )
end
