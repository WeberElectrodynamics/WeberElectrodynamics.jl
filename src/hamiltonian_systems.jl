using Symbolics
using Latexify: latexify

# =============================================================================
# Phase Space Variable Creation
# =============================================================================

"""
    create_phase_space_variables(n_particles::Int, dims::Int)

Create symbolic variable names for phase space coordinates.

# Arguments
- `n_particles`: Number of particles
- `dims`: Spatial dimensions (1, 2, or 3)

# Returns
Tuple `(q_syms, p_syms)` of Symbol vectors.

# Example
```julia
q_syms, p_syms = create_phase_space_variables(2, 3)
# q_syms = [:x1, :y1, :z1, :x2, :y2, :z2]
# p_syms = [:px1, :py1, :pz1, :px2, :py2, :pz2]
```
"""
function create_phase_space_variables(n_particles::Int, dims::Int)::Tuple{Vector{Symbol}, Vector{Symbol}}
    coord_names = [:x, :y, :z]
    momentum_names = [:px, :py, :pz]

    q_syms = Symbol[]
    p_syms = Symbol[]

    for i in 1:n_particles
        for d in 1:dims
            push!(q_syms, Symbol(string(coord_names[d]) * string(i)))
            push!(p_syms, Symbol(string(momentum_names[d]) * string(i)))
        end
    end

    return (q_syms, p_syms)
end

# =============================================================================
# Hamiltonian Compilation
# =============================================================================

"""
    @hamiltonian(n_particles, dims, param_names, H_expr)

Create a compiled Hamiltonian from a symbolic expression.

# Arguments
- `n_particles`: Number of particles
- `dims`: Spatial dimensions (1, 2, or 3)
- `param_names`: Vector of parameter names as symbols
- `H_expr`: Function `(q, p, params) -> H` returning Hamiltonian expression

# Returns
`WeberHamiltonian` with compiled vector field functions.

# Example
```julia
H = @hamiltonian 2 2 [:m1, :m2, :k, :c] (q, p, params) -> begin
    m1, m2, k, c = params
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    KE = (px1^2 + py1^2)/(2m1) + (px2^2 + py2^2)/(2m2)
    # ... potential energy
    return KE + PE
end
```
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
    build_hamiltonian(H_func, n_particles, dims; param_names=Symbol[])

Build a WeberHamiltonian from a function without using the macro.

# Arguments
- `H_func`: Function `(q, p, params) -> H` returning symbolic Hamiltonian
- `n_particles`: Number of particles
- `dims`: Spatial dimensions
- `param_names`: Parameter names (optional, for documentation)

# Example
```julia
function my_hamiltonian(q, p, params)
    m1, m2, k, c = params
    # ... define H
    return KE + PE
end

H = build_hamiltonian(my_hamiltonian, 2, 2; param_names=[:m1, :m2, :k, :c])
```
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
