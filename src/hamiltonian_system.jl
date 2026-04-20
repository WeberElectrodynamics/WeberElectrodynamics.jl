using Symbolics
using Latexify: latexify

"""
    HamiltonianSystem

Symbolic and compiled representation of the n-body Weber Hamiltonian.

Constructed once for a given `(n_particles, dims)` pair via `HamiltonianSystem(n_particles, dims)`.
The compiled equations of motion are reused across many `HamiltonianProblem` instances
with different physical parameters. Construction involves symbolic differentiation
and code generation via Symbolics.jl; expect a few seconds for the first call.

# Fields
- `n_particles::Int`: Number of charged particles.
- `dims::Int`: Spatial dimension (1, 2, or 3).
- `q_symbols::Vector{Num}`: Symbolic coordinate variables `[x1, y1, ..., xN, yN, ...]`.
- `p_symbols::Vector{Num}`: Symbolic momentum variables `[px1, py1, ...]`.
- `t_symbol::Num`: Symbolic time variable. Reserved for time-dependent terms;
  the current Weber Hamiltonian is autonomous and does not use it.
- `param_symbols`: Symbolic parameter vector `[m1…mN, q1…qN, c, κ12, κ13, …]`.
- `hamiltonian_symbolic`: Full symbolic Weber Hamiltonian expression.
- `dq_dt_symbolic`, `dp_dt_symbolic`: Symbolic Hamilton's equations.
- `dq_dt_compiled(out, q, p, t, params)`, `dp_dt_compiled(out, q, p, t, params)`:
  In-place compiled equations of motion. `t` is currently unused.
- `hamiltonian_compiled(q, p, t, params)`: Compiled scalar Hamiltonian function.
- `degrees_of_freedom::Int`: Total DOF = `n_particles × dims`.
"""
struct HamiltonianSystem{H,QD,PD,QF,PF,HF,PS}
    n_particles::Int
    dims::Int

    q_symbols::Vector{Num}
    p_symbols::Vector{Num}
    t_symbol::Num
    param_symbols::PS

    hamiltonian_symbolic::H
    dq_dt_symbolic::QD
    dp_dt_symbolic::PD

    dq_dt_compiled::QF
    dp_dt_compiled::PF
    hamiltonian_compiled::HF

    degrees_of_freedom::Int
end

function _generate_phase_space_symbols(n_particles::Int, dims::Int)
    coord_names = [:x, :y, :z]
    momentum_names = [:px, :py, :pz]

    n_total = n_particles * dims
    coordinate_symbols = Vector{Symbol}(undef, n_total)
    momentum_symbols = Vector{Symbol}(undef, n_total)

    idx = 1
    @inbounds for i = 1:n_particles
        for d = 1:dims
            coordinate_symbols[idx] = Symbol(string(coord_names[d]) * string(i))
            momentum_symbols[idx] = Symbol(string(momentum_names[d]) * string(i))
            idx += 1
        end
    end

    return (coordinate_symbols, momentum_symbols)
end

function _generate_param_symbols(n_particles::Int)
    mass_symbols = [Symbol("m$i") for i = 1:n_particles]
    charge_symbols = [Symbol("q$i") for i = 1:n_particles]
    kappa_symbols = [
        Symbol("kappa_$(i)_$(j)") for i = 1:n_particles for j = (i+1):n_particles
    ]
    return (mass_symbols, charge_symbols, :c, kappa_symbols)
end

# Pair (i,j) with i<j maps to this 1-based index in the kappas vector.
@inline function _pair_index(i::Int, j::Int, n::Int)::Int
    return (i - 1) * (2n - i) ÷ 2 + (j - i)
end

include("hamiltonian/builders/weber.jl")
include("hamiltonian/builders/zollner.jl")

"""
    HamiltonianSystem(H, q_vars, p_vars;
                      param_symbols, t, n_particles, dims) -> HamiltonianSystem

Generic constructor from a pre-built symbolic Hamiltonian `H`.

Derives Hamilton's equations analytically via `Symbolics.derivative`, then
compiles them to in-place Julia functions with
`Symbolics.build_function(…, q_vars, p_vars, t, param_symbols)`. The signature
of the compiled EOMs is `(out, q, p, t, params)`; `t` is presently unused but
reserved for time-dependent terms.

Use this overload to build a custom Hamiltonian (e.g. `H = weber_term(…) +
zollner_term(…)`). For the default pure Weber case, the convenience constructor
`HamiltonianSystem(n_particles, dims)` wraps this path.

# Arguments
- `H`: Symbolic Hamiltonian expression.
- `q_vars`, `p_vars`: Phase-space symbolic variables, each length `n_particles*dims`.

# Keywords
- `param_symbols`: Symbolic parameter vector passed to `build_function`.
- `t`: Symbolic time variable (reserved; may be unused).
- `n_particles::Int`, `dims::Int`: Problem shape.
"""
function HamiltonianSystem(
    H,
    q_vars::AbstractVector,
    p_vars::AbstractVector;
    param_symbols::AbstractVector,
    t,
    n_particles::Int,
    dims::Int,
)
    dq_dt_symbolic = [Symbolics.derivative(H, p_vars[i]) for i in eachindex(p_vars)]
    dp_dt_symbolic = [-Symbolics.derivative(H, q_vars[i]) for i in eachindex(q_vars)]

    dq_dt_compiled = Symbolics.build_function(
        dq_dt_symbolic, q_vars, p_vars, t, param_symbols, expression = Val{false},
    )[2]
    dp_dt_compiled = Symbolics.build_function(
        dp_dt_symbolic, q_vars, p_vars, t, param_symbols, expression = Val{false},
    )[2]
    hamiltonian_compiled = Symbolics.build_function(
        H, q_vars, p_vars, t, param_symbols, expression = Val{false},
    )

    degrees_of_freedom = n_particles * dims

    HamiltonianSystem(
        n_particles,
        dims,
        q_vars,
        p_vars,
        t,
        param_symbols,
        H,
        dq_dt_symbolic,
        dp_dt_symbolic,
        dq_dt_compiled,
        dp_dt_compiled,
        hamiltonian_compiled,
        degrees_of_freedom,
    )
end

"""
    HamiltonianSystem(n_particles::Int, dims::Int) -> HamiltonianSystem

Convenience constructor for the default n-body Weber Hamiltonian in `dims`
spatial dimensions. Equivalent to building `weber_term(…)` with all κ=1 symbols
and passing it to the generic `HamiltonianSystem(H, q, p; …)` constructor.

# Arguments
- `n_particles`: Number of particles (≥ 1).
- `dims`: Spatial dimension; must be 1, 2, or 3.
"""
function HamiltonianSystem(n_particles::Int, dims::Int)
    @assert n_particles >= 1 "Must have at least 1 particle"
    @assert dims in (1, 2, 3) "Dimensions must be 1, 2, or 3, got $dims"

    coordinate_symbols, momentum_symbols = _generate_phase_space_symbols(n_particles, dims)
    q_vars = [Symbolics.variable(sym) for sym in coordinate_symbols]
    p_vars = [Symbolics.variable(sym) for sym in momentum_symbols]

    mass_symbols, charge_symbols, c_symbol, kappa_syms = _generate_param_symbols(n_particles)
    m_vars = [Symbolics.variable(sym) for sym in mass_symbols]
    charge_vars = [Symbolics.variable(sym) for sym in charge_symbols]
    c_var = Symbolics.variable(c_symbol)
    kappa_vars = [Symbolics.variable(sym) for sym in kappa_syms]

    t_var = Symbolics.variable(:t)

    param_symbols = vcat(m_vars, charge_vars, [c_var], kappa_vars)

    H = weber_term(
        q_vars,
        p_vars;
        masses = m_vars,
        charges = charge_vars,
        c = c_var,
        kappas = kappa_vars,
        n_particles = n_particles,
        dims = dims,
    )

    return HamiltonianSystem(
        H,
        q_vars,
        p_vars;
        param_symbols = param_symbols,
        t = t_var,
        n_particles = n_particles,
        dims = dims,
    )
end

function Base.show(io::IO, sys::HamiltonianSystem)
    print(
        io,
        "HamiltonianSystem($(sys.n_particles) particles, $(sys.dims)D, $(sys.degrees_of_freedom) DOF)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", sys::HamiltonianSystem)
    println(io, "HamiltonianSystem")
    println(io, "  Particles: $(sys.n_particles)")
    println(io, "  Dimensions: $(sys.dims)")
    println(io, "  DOF: $(sys.degrees_of_freedom)")
    println(io, "  H = $(sys.hamiltonian_symbolic)")
end

function Base.show(io::IO, ::MIME"text/latex", sys::HamiltonianSystem)
    print(io, latexify(sys.hamiltonian_symbolic))
end
