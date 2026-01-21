using Symbolics
using Latexify: latexify

# =============================================================================
# Phase Space Variable Creation
# =============================================================================

"""
    generate_phase_space_symbols(n_particles, dims) -> (coordinate_names, momentum_names)

Generate symbolic variable names for phase space coordinates.
Returns tuple of Symbol vectors for positions (x,y,z) and momenta (px,py,pz).
"""
function generate_phase_space_symbols(n_particles::Int, dims::Int)::Tuple{Vector{Symbol}, Vector{Symbol}}
    coord_names = [:x, :y, :z]
    momentum_names = [:px, :py, :pz]

    # Pre-size arrays (avoids reallocations from push!)
    n_total = n_particles * dims
    coordinate_symbols = Vector{Symbol}(undef, n_total)
    momentum_symbols = Vector{Symbol}(undef, n_total)

    idx = 1
    for i in 1:n_particles
        for d in 1:dims
            coordinate_symbols[idx] = Symbol(string(coord_names[d]) * string(i))
            momentum_symbols[idx] = Symbol(string(momentum_names[d]) * string(i))
            idx += 1
        end
    end

    return (coordinate_symbols, momentum_symbols)
end

# =============================================================================
# Hamiltonian Compilation
# =============================================================================

"""
    @hamiltonian n_particles dims parameter_names H_expr -> WeberHamiltonian

Create a compiled Hamiltonian. `H_expr` is `(q, p, params) -> H`.

    H = @hamiltonian 2 2 [:m, :k] (q, p, params) -> p'p/(2params[1]) + params[2]*q'q/2
"""
macro hamiltonian(n_particles, dims, parameter_names, H_expr)
    quote
        _compile_hamiltonian_internal($(esc(n_particles)), $(esc(dims)), $(esc(parameter_names)), $(esc(H_expr)))
    end
end

"""
    _compile_hamiltonian_internal(n_particles, dims, parameter_names, H_func)

Internal function to compile WeberHamiltonian from a Hamiltonian function.
"""
function _compile_hamiltonian_internal(n_particles::Int, dims::Int, parameter_names::Vector{Symbol}, H_func::Function)
    coordinate_symbols, momentum_symbols = generate_phase_space_symbols(n_particles, dims)

    q_vars = [Symbolics.variable(sym) for sym in coordinate_symbols]
    p_vars = [Symbolics.variable(sym) for sym in momentum_symbols]
    param_vars = [Symbolics.variable(param_name) for param_name in parameter_names]

    hamiltonian_symbolic = H_func(q_vars, p_vars, param_vars)

    # Hamilton's equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
    dq_dt_symbolic = [Symbolics.derivative(hamiltonian_symbolic, p_vars[i]) for i in eachindex(p_vars)]
    dp_dt_symbolic = [-Symbolics.derivative(hamiltonian_symbolic, q_vars[i]) for i in eachindex(q_vars)]

    # Compile to fast numeric functions
    # Signature: func(out, q, p, params) - GI-compatible
    dq_dt_compiled = build_function(dq_dt_symbolic, q_vars, p_vars, param_vars, expression=Val{false})[2]
    dp_dt_compiled = build_function(dp_dt_symbolic, q_vars, p_vars, param_vars, expression=Val{false})[2]

    WeberHamiltonian(
        hamiltonian_symbolic,
        dq_dt_symbolic,
        dp_dt_symbolic,
        dq_dt_compiled,
        dp_dt_compiled,
        length(q_vars),
        parameter_names
    )
end

# =============================================================================
# Alternative: Function-based API (without macro)
# =============================================================================

"""
    compile_hamiltonian(H_func, n_particles, dims; parameter_names=Symbol[]) -> WeberHamiltonian

Compile a WeberHamiltonian from a function `(q, p, params) -> H` without using the macro.
"""
function compile_hamiltonian(H_func::Function, n_particles::Int, dims::Int; parameter_names::Vector{Symbol}=Symbol[])
    _compile_hamiltonian_internal(n_particles, dims, parameter_names, H_func)
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, H::WeberHamiltonian)
    print(io, "WeberHamiltonian($(H.degrees_of_freedom) DOF, params=$(H.parameter_names))")
end

function Base.show(io::IO, ::MIME"text/plain", H::WeberHamiltonian)
    println(io, "WeberHamiltonian")
    println(io, "  DOF: $(H.degrees_of_freedom)")
    println(io, "  Parameters: $(H.parameter_names)")
    println(io, "  H = $(H.hamiltonian_symbolic)")
end

# LaTeX display for Jupyter notebooks
function Base.show(io::IO, ::MIME"text/latex", H::WeberHamiltonian)
    print(io, latexify(H.hamiltonian_symbolic))
end
