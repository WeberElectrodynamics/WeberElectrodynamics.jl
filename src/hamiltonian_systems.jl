using Symbolics
using Latexify: latexify

# =============================================================================
# Phase Space Variable Creation
# =============================================================================

"""
    create_phase_space_variables(n_particles, dims) -> (q_syms, p_syms)

Create symbolic variable names for phase space coordinates.
Returns tuple of Symbol vectors for positions (x,y,z) and momenta (px,py,pz).
"""
function create_phase_space_variables(n_particles::Int, dims::Int)::Tuple{Vector{Symbol}, Vector{Symbol}}
    coord_names = [:x, :y, :z]
    momentum_names = [:px, :py, :pz]

    # Pre-size arrays (avoids reallocations from push!)
    n_total = n_particles * dims
    q_syms = Vector{Symbol}(undef, n_total)
    p_syms = Vector{Symbol}(undef, n_total)

    idx = 1
    for i in 1:n_particles
        for d in 1:dims
            q_syms[idx] = Symbol(string(coord_names[d]) * string(i))
            p_syms[idx] = Symbol(string(momentum_names[d]) * string(i))
            idx += 1
        end
    end

    return (q_syms, p_syms)
end

# =============================================================================
# Hamiltonian Compilation
# =============================================================================

"""
    @hamiltonian n_particles dims param_names H_expr -> WeberHamiltonian

Create a compiled Hamiltonian. `H_expr` is `(q, p, params) -> H`.

    H = @hamiltonian 2 2 [:m, :k] (q, p, params) -> p'p/(2params[1]) + params[2]*q'q/2
"""
macro hamiltonian(n_particles, dims, param_names, H_expr)
    quote
        _build_hamiltonian($(esc(n_particles)), $(esc(dims)), $(esc(param_names)), $(esc(H_expr)))
    end
end

"""
    _build_hamiltonian(n_particles, dims, param_names, H_func)

Internal function to build WeberHamiltonian from a Hamiltonian function.
"""
function _build_hamiltonian(n_particles::Int, dims::Int, param_names::Vector{Symbol}, H_func::Function)
    q_syms, p_syms = create_phase_space_variables(n_particles, dims)

    q_vars = [Symbolics.variable(q_sym) for q_sym in q_syms]
    p_vars = [Symbolics.variable(p_sym) for p_sym in p_syms]
    param_vars = [Symbolics.variable(param_name) for param_name in param_names]

    H_sym = H_func(q_vars, p_vars, param_vars)

    qdot_sym = [Symbolics.derivative(H_sym, p_vars[i]) for i in eachindex(p_vars)]
    pdot_sym = [-Symbolics.derivative(H_sym, q_vars[i]) for i in eachindex(q_vars)]

    # Compile to fast numeric functions
    # Signature: func(out, q, p, params) - GI-compatible
    qdot_func = build_function(qdot_sym, q_vars, p_vars, param_vars, expression=Val{false})[2]
    pdot_func = build_function(pdot_sym, q_vars, p_vars, param_vars, expression=Val{false})[2]

    WeberHamiltonian(
        H_sym,
        qdot_sym,
        pdot_sym,
        qdot_func,
        pdot_func,
        length(q_vars),
        param_names
    )
end

# =============================================================================
# Alternative: Function-based API (without macro)
# =============================================================================

"""
    build_hamiltonian(H_func, n_particles, dims; param_names=Symbol[]) -> WeberHamiltonian

Build a WeberHamiltonian from a function `(q, p, params) -> H` without using the macro.
"""
function build_hamiltonian(H_func::Function, n_particles::Int, dims::Int; param_names::Vector{Symbol}=Symbol[])
    _build_hamiltonian(n_particles, dims, param_names, H_func)
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, H::WeberHamiltonian)
    print(io, "WeberHamiltonian($(H.n_dof) DOF, params=$(H.param_names))")
end

function Base.show(io::IO, ::MIME"text/plain", H::WeberHamiltonian)
    println(io, "WeberHamiltonian")
    println(io, "  DOF: $(H.n_dof)")
    println(io, "  Parameters: $(H.param_names)")
    println(io, "  H = $(H.H_sym)")
end

# LaTeX display for Jupyter notebooks
function Base.show(io::IO, ::MIME"text/latex", H::WeberHamiltonian)
    print(io, latexify(H.H_sym))
end
