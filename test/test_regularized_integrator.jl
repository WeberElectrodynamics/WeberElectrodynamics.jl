@testset "RegularizedIntegrator wrapper" begin
    @testset "subtyping + unwrap" begin
        base = SymmetricProjectionIntegrator(relaxation = 0.25)
        alg = WeberElectrodynamics.RegularizedIntegrator(base; backend = :lifted_pair)
        @test alg isa WeberElectrodynamics.HamiltonianAlgorithm
        @test WeberElectrodynamics.base_algorithm(alg) === base
        @test WeberElectrodynamics.base_algorithm(base) === base
        @test alg.options.enabled === true
        @test alg.options.backend === :lifted_pair
        @test alg.options.r_on_factor == 0.15
    end

    @testset "options kwargs forwarded" begin
        alg = WeberElectrodynamics.RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            r_on_factor = 0.2,
            r_off_factor = 0.5,
            max_substeps = 64,
            backend = :adaptive_cartesian,
            chain_enabled = false,
            warn_on_fallback = false,
            collision_bounce_radius = 0.01,
        )
        @test alg.options.r_on_factor == 0.2
        @test alg.options.r_off_factor == 0.5
        @test alg.options.max_substeps == 64
        @test alg.options.backend === :adaptive_cartesian
        @test alg.options.chain_enabled === false
        @test alg.options.warn_on_fallback === false
        @test alg.options.collision_bounce_radius == 0.01
    end

    @testset "default base algorithm shorthand" begin
        alg = WeberElectrodynamics.RegularizedIntegrator(; backend = :adaptive_cartesian)
        @test alg.base_alg isa SymmetricProjectionIntegrator
        @test alg.options.enabled === true
        @test alg.options.backend === :adaptive_cartesian
    end
end
