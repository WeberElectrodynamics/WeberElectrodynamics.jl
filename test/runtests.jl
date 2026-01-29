using Test
using WeberElectrodynamics
using LinearAlgebra

@testset "WeberElectrodynamics.jl" begin
    include("test_utils.jl")
    include("test_types.jl")
    include("test_hamiltonian_systems.jl")
    include("test_solve.jl")
    include("test_statistics.jl")
    include("test_integration.jl")
    include("test_physics.jl")
end
