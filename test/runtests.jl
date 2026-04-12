using Test
using WeberElectrodynamics
using WeberElectrodynamics: SymmetricProjectionIntegrator
using LinearAlgebra
using Random
using Symbolics

@testset "WeberElectrodynamics.jl" begin
    include("test_utils.jl")
    include("test_types.jl")
    include("test_weber_system.jl")
    include("test_solve.jl")
    include("test_statistics.jl")
    include("test_integration.jl")
    include("test_physics.jl")
    include("test_regularization.jl")
    include("test_zollner.jl")
    include("test_aqua.jl")
end
