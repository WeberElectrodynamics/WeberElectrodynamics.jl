using Symbolics
using Latexify: latexify

"""
    WeberSystem

Symbolic and compiled representation of the n-body Weber Hamiltonian.

Constructed once for a given `(n_particles, dims)` pair via `WeberSystem(n_particles, dims)`.
The compiled equations of motion are reused across many `WeberProblem` instances
with different physical parameters. Construction involves symbolic differentiation
and code generation via Symbolics.jl; expect a few seconds for the first call.

# Fields
- `n_particles::Int`: Number of charged particles.
- `dims::Int`: Spatial dimension (1, 2, or 3).
- `q_symbols::Vector{Num}`: Symbolic coordinate variables `[x1, y1, ..., xN, yN, ...]`.
- `p_symbols::Vector{Num}`: Symbolic momentum variables `[px1, py1, ...]`.
- `param_symbols`: Symbolic parameter vector `[m1…mN, q1…qN, c, κ12, κ13, …]`.
- `hamiltonian_symbolic`: Full symbolic Weber Hamiltonian expression.
- `dq_dt_symbolic`, `dp_dt_symbolic`: Symbolic Hamilton's equations.
- `dq_dt_compiled`, `dp_dt_compiled`: In-place compiled equations of motion.
- `hamiltonian_compiled`: Compiled scalar Hamiltonian function.
- `degrees_of_freedom::Int`: Total DOF = `n_particles × dims`.
"""
struct WeberSystem{H,QD,PD,QF,PF,HF,PS}
    n_particles::Int
    dims::Int

    q_symbols::Vector{Num}
    p_symbols::Vector{Num}
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

function _build_weber_hamiltonian(
    q_vars::Vector,
    p_vars::Vector,
    m_vars::Vector,
    charge_vars::Vector,
    c_var,
    kappa_vars::Vector,
    n_particles::Int,
    dims::Int,
)
    H = zero(eltype(q_vars))

    @inbounds for i = 1:n_particles
        p_start = (i - 1) * dims + 1
        p_end = i * dims
        p_squared = sum(p_vars[p_start:p_end] .^ 2)
        H = H + p_squared / (2 * m_vars[i])
    end

    c_squared = c_var^2
    @inbounds for i = 1:n_particles
        for j = (i+1):n_particles
            qi_start = (i - 1) * dims + 1
            qj_start = (j - 1) * dims + 1

            pi_start = (i - 1) * dims + 1
            pj_start = (j - 1) * dims + 1

            r_squared = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                r_squared = r_squared + dq^2
            end
            r = sqrt(r_squared)

            r_dot_v = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                dv = p_vars[pi_start+d-1] / m_vars[i] - p_vars[pj_start+d-1] / m_vars[j]
                r_dot_v = r_dot_v + dq * dv
            end
            r_dot = r_dot_v / r

            pair_idx = _pair_index(i, j, n_particles)
            kappa = kappa_vars[pair_idx]
            k = kappa * charge_vars[i] * charge_vars[j]
            U_ij = k / r * (1 - r_dot^2 / (2 * c_squared))
            H = H + U_ij
        end
    end

    return H
end

"""
    WeberSystem(n_particles::Int, dims::Int) -> WeberSystem

Build and compile the symbolic Weber Hamiltonian for `n_particles` particles in
`dims` dimensions.

Uses Symbolics.jl to derive Hamilton's equations analytically, then compiles
them to efficient in-place Julia functions via `build_function`. The resulting
`WeberSystem` can be reused across multiple `WeberProblem` instances.

# Arguments
- `n_particles`: Number of particles (≥ 1).
- `dims`: Spatial dimension; must be 1, 2, or 3.

# Returns
- `WeberSystem` with compiled equations of motion ready for `WeberProblem`.
"""
function WeberSystem(n_particles::Int, dims::Int)
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

    param_symbols = vcat(m_vars, charge_vars, [c_var], kappa_vars)

    hamiltonian_symbolic = _build_weber_hamiltonian(
        q_vars,
        p_vars,
        m_vars,
        charge_vars,
        c_var,
        kappa_vars,
        n_particles,
        dims,
    )

    dq_dt_symbolic =
        [Symbolics.derivative(hamiltonian_symbolic, p_vars[i]) for i in eachindex(p_vars)]
    dp_dt_symbolic =
        [-Symbolics.derivative(hamiltonian_symbolic, q_vars[i]) for i in eachindex(q_vars)]

    dq_dt_compiled =
        Symbolics.build_function(dq_dt_symbolic, q_vars, p_vars, param_symbols, expression = Val{false})[2]
    dp_dt_compiled =
        Symbolics.build_function(dp_dt_symbolic, q_vars, p_vars, param_symbols, expression = Val{false})[2]
    hamiltonian_compiled =
        Symbolics.build_function(hamiltonian_symbolic, q_vars, p_vars, param_symbols, expression = Val{false})

    degrees_of_freedom = n_particles * dims

    WeberSystem(
        n_particles,
        dims,
        q_vars,
        p_vars,
        param_symbols,
        hamiltonian_symbolic,
        dq_dt_symbolic,
        dp_dt_symbolic,
        dq_dt_compiled,
        dp_dt_compiled,
        hamiltonian_compiled,
        degrees_of_freedom,
    )
end

function Base.show(io::IO, sys::WeberSystem)
    print(
        io,
        "WeberSystem($(sys.n_particles) particles, $(sys.dims)D, $(sys.degrees_of_freedom) DOF)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", sys::WeberSystem)
    println(io, "WeberSystem")
    println(io, "  Particles: $(sys.n_particles)")
    println(io, "  Dimensions: $(sys.dims)")
    println(io, "  DOF: $(sys.degrees_of_freedom)")
    println(io, "  H = $(sys.hamiltonian_symbolic)")
end

function Base.show(io::IO, ::MIME"text/latex", sys::WeberSystem)
    print(io, latexify(sys.hamiltonian_symbolic))
end
