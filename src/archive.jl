"""
    save_solution(path, sol; metadata=NamedTuple())

Save a [`HamiltonianSolution`](@ref) to a JLD2 archive.

The archive stores primitive solution arrays plus problem metadata for default
Weber/Zöllner systems. Custom symbolic Hamiltonians are not reconstructable by
this helper yet.

This function is implemented by the optional JLD2 extension. Load it with
`using JLD2` before calling `save_solution`.
"""
function save_solution end

"""
    load_solution(path) -> HamiltonianSolution

Load a [`HamiltonianSolution`](@ref) from a JLD2 archive created by
[`save_solution`](@ref). The default Weber/Zöllner `HamiltonianSystem` is
reconstructed from the archived metadata.

This function is implemented by the optional JLD2 extension. Load it with
`using JLD2` before calling `load_solution`.
"""
function load_solution end
