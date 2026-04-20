using Test
using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using LinearAlgebra
using Random
using Symbolics

@testset "WeberElectrodynamics.jl" begin
    include("test_utils.jl")
    include("test_types.jl")
    include("test_hamiltonian_system.jl")
    include("test_builders.jl")
    include("test_named_term.jl")
    include("test_accessors.jl")
    include("test_algorithm_dispatch.jl")
    include("test_solve.jl")
    include("test_statistics.jl")
    include("test_integration.jl")
    include("test_physics.jl")
    include("test_regularization.jl")
    include("test_zollner.jl")
    include("test_aqua.jl")

    @testset "regression fixtures" begin
        include("regression/validate.jl")
        @test validate_all()
    end
end
