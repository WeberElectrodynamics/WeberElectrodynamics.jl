using Symbolics
using Latexify: latexify

include("hamiltonian/terms.jl")

"""
    HamiltonianSystem

Compiled representation of an n-body Hamiltonian, with an optional symbolic
backing.

There are two construction paths:

- **Analytic** — `HamiltonianSystem(n_particles, dims)` builds the default
  *exact canonical Weber* system. Its equations of motion are hand-derived
  closed forms evaluated numerically, because recovering physical velocities
  from canonical momenta requires solving a coupled linear system that has no
  practical closed symbolic form for general `n`. The three symbolic fields are
  `nothing` for this path.
- **Symbolic** — `HamiltonianSystem(H, q_vars, p_vars; …)` takes any symbolic
  Hamiltonian, differentiates it with Symbolics.jl, and compiles the result.
  Use this for custom Hamiltonians.

Both paths expose the same compiled-function signatures, so everything
downstream (solver, statistics, plotting, animation) is agnostic to which was
used. Construct a system once for a given `(n_particles, dims)` and reuse it
across many `HamiltonianProblem` instances.

# Fields
- `n_particles::Int`: Number of charged particles.
- `dims::Int`: Spatial dimension (1, 2, or 3).
- `q_symbols::Vector{Num}`: Symbolic coordinate variables `[x1, y1, ..., xN, yN, ...]`.
- `p_symbols::Vector{Num}`: Symbolic momentum variables `[px1, py1, ...]`.
- `t_symbol::Num`: Symbolic time variable. Reserved for time-dependent terms;
  the Weber Hamiltonian is autonomous and does not use it.
- `param_symbols`: Symbolic parameter vector `[m1…mN, q1…qN, c]`.
- `hamiltonian_symbolic`: Symbolic Hamiltonian expression, or `nothing` for an
  analytic system.
- `dq_dt_symbolic`, `dp_dt_symbolic`: Symbolic Hamilton's equations, or
  `nothing` for an analytic system.
- `dq_dt_compiled(out, q, p, t, params)`, `dp_dt_compiled(out, q, p, t, params)`:
  In-place compiled equations of motion. `t` is currently unused.
- `hamiltonian_compiled(q, p, t, params)`: Compiled scalar Hamiltonian function.
- `degrees_of_freedom::Int`: Total DOF = `n_particles × dims`.
- `terms::Vector{NamedTerm}`: Named components of the Hamiltonian preserving
  the decomposition for per-term statistics and plotting. The generic
  constructor assigns a single `:hamiltonian` term by default.

See also [`has_symbolic_hamiltonian`](@ref).
"""
struct HamiltonianSystem{H,QD,PD,QF,PF,HF,PS,TS<:AbstractVector{<:NamedTerm}}
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

    terms::TS
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
    return (mass_symbols, charge_symbols, :c)
end

include("hamiltonian/builders/weber.jl")
include("hamiltonian/builders/basic.jl")

"""
    HamiltonianSystem(H, q_vars, p_vars;
                      param_symbols, t, n_particles, dims) -> HamiltonianSystem

Generic constructor from a pre-built symbolic Hamiltonian `H`.

Derives Hamilton's equations analytically via `Symbolics.derivative`, then
compiles them to in-place Julia functions with
`Symbolics.build_function(…, q_vars, p_vars, t, param_symbols)`.
The signature of the compiled EOMs is `(out, q, p, t, params)`;
`t` is presently unused but reserved for time-dependent terms.

Use this overload to build a custom Hamiltonian. For the default Weber case,
the convenience constructor `HamiltonianSystem(n_particles, dims)` wraps this
path.

# Arguments
- `H`: Symbolic Hamiltonian expression.
- `q_vars`, `p_vars`: Phase-space symbolic variables, each length `n_particles*dims`.

# Keywords
- `param_symbols`: Symbolic parameter vector `[m1…mN, q1…qN, c]`.
- `t`: Symbolic time variable (reserved; may be unused).
- `n_particles::Int`, `dims::Int`: Problem shape.
- `terms`: Optional `Vector{NamedTerm}` naming the components of `H`. If
  omitted, a single `NamedTerm(:hamiltonian, H)` is stored so queries like
  `get_term(sys, :hamiltonian)` still work.
"""
function HamiltonianSystem(
    H,
    q_vars::AbstractVector,
    p_vars::AbstractVector;
    param_symbols::AbstractVector,
    t,
    n_particles::Int,
    dims::Int,
    terms::Union{Nothing,AbstractVector{<:NamedTerm}} = nothing,
)
    dq_dt_symbolic = [Symbolics.derivative(H, p_vars[i]) for i in eachindex(p_vars)]
    dp_dt_symbolic = [-Symbolics.derivative(H, q_vars[i]) for i in eachindex(q_vars)]

    dq_dt_compiled = Symbolics.build_function(
        dq_dt_symbolic,
        q_vars,
        p_vars,
        t,
        param_symbols,
        expression = Val{false},
    )[2]
    dp_dt_compiled = Symbolics.build_function(
        dp_dt_symbolic,
        q_vars,
        p_vars,
        t,
        param_symbols,
        expression = Val{false},
    )[2]
    hamiltonian_compiled = Symbolics.build_function(
        H,
        q_vars,
        p_vars,
        t,
        param_symbols,
        expression = Val{false},
    )

    degrees_of_freedom = n_particles * dims

    resolved_terms = terms === nothing ? [NamedTerm(:hamiltonian, H)] : collect(terms)

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
        resolved_terms,
    )
end

"""
    HamiltonianSystem(n_particles::Int, dims::Int) -> HamiltonianSystem

Convenience constructor for the default n-body **exact canonical Weber** system
in `dims` spatial dimensions.

The Weber Lagrangian is velocity dependent, so its canonical momentum is
`p_i = ∂L/∂v_i = m_i v_i − Σ_{j≠i} (q_i q_j/c²) ṙ_ij (r_i − r_j)/r_ij²`, not
`m_i v_i`. Recovering the physical velocities therefore requires solving a
coupled linear system over the `n(n−1)/2` particle pairs at every evaluation.
That solve has no practical closed symbolic form for general `n`, so this
constructor takes the **analytic** path: the exact Hamiltonian and both
canonical equations are evaluated by hand-derived closed-form routines rather
than by Symbolics.jl differentiation. `hamiltonian_symbolic`, `dq_dt_symbolic`,
and `dp_dt_symbolic` are `nothing` for the resulting system.

The compiled-function signatures are identical to the symbolic path, so
`HamiltonianProblem`, `solve`, statistics, plotting, and animation are
unaffected.

Construction is cheap — there is no symbolic differentiation or code generation.

# Arguments
- `n_particles`: Number of particles (≥ 1).
- `dims`: Spatial dimension; must be 1, 2, or 3.

# Throws
- [`WeberCriticalRadiusError`](@ref) during evaluation if the canonical mass
  matrix becomes singular (a like-charge pair exactly at Weber's critical
  radius `ρ = q₁q₂/(μc²)`).
"""
function HamiltonianSystem(n_particles::Int, dims::Int)
    @assert n_particles >= 1 "Must have at least 1 particle"
    @assert dims in (1, 2, 3) "Dimensions must be 1, 2, or 3, got $dims"

    coordinate_symbols, momentum_symbols = _generate_phase_space_symbols(n_particles, dims)
    q_vars = [Symbolics.variable(sym) for sym in coordinate_symbols]
    p_vars = [Symbolics.variable(sym) for sym in momentum_symbols]

    mass_symbols, charge_symbols, c_symbol = _generate_param_symbols(n_particles)
    m_vars = [Symbolics.variable(sym) for sym in mass_symbols]
    charge_vars = [Symbolics.variable(sym) for sym in charge_symbols]
    c_var = Symbolics.variable(c_symbol)

    t_var = Symbolics.variable(:t)

    param_symbols = vcat(m_vars, charge_vars, [c_var])

    # Each compiled entry point owns a workspace; a workspace is not reentrant.
    ws_dq = WeberWorkspace(n_particles, dims)
    ws_dp = WeberWorkspace(n_particles, dims)
    ws_H = WeberWorkspace(n_particles, dims)
    ws_decomp = WeberWorkspace(n_particles, dims)
    ws_ke = WeberWorkspace(n_particles, dims)

    dq_dt_compiled = (out, q, p, t, params) -> weber_dq_dt!(out, ws_dq, q, p, params)
    dp_dt_compiled = (out, q, p, t, params) -> weber_dp_dt!(out, ws_dp, q, p, params)
    hamiltonian_compiled = (q, p, t, params) -> weber_hamiltonian(ws_H, q, p, params)

    weber_decomp = (q, p, params) -> _weber_pair_decomposition(ws_decomp, q, p, params)
    weber_kinetic = (q, p, params) -> _weber_kinetic_energy(ws_ke, q, p, params)

    return HamiltonianSystem(
        n_particles,
        dims,
        q_vars,
        p_vars,
        t_var,
        param_symbols,
        nothing,
        nothing,
        nothing,
        dq_dt_compiled,
        dp_dt_compiled,
        hamiltonian_compiled,
        n_particles * dims,
        [
            NamedTerm(
                :weber,
                nothing;
                pair_decomposition = weber_decomp,
                kinetic_energy = weber_kinetic,
            ),
        ],
    )
end

"""
    has_symbolic_hamiltonian(sys::HamiltonianSystem) -> Bool

Return `true` when `sys` carries a Symbolics.jl expression for `H` and its
derivatives, i.e. when it was built through the generic symbolic constructor.

The default Weber system is analytic and returns `false`: its exact canonical
Hamiltonian requires a per-evaluation linear solve and has no closed symbolic
form for general `n`. Guard any code that reads `hamiltonian_symbolic`,
`dq_dt_symbolic`, or `dp_dt_symbolic` with this predicate.
"""
has_symbolic_hamiltonian(sys::HamiltonianSystem) = sys.hamiltonian_symbolic !== nothing

"""
    term_names(sys::HamiltonianSystem) -> Vector{Symbol}
    has_term(sys::HamiltonianSystem, name::Symbol) -> Bool
    get_term(sys::HamiltonianSystem, name::Symbol) -> NamedTerm

Query the named Hamiltonian components stored on `sys`.
"""
term_names(sys::HamiltonianSystem) = term_names(sys.terms)
has_term(sys::HamiltonianSystem, name::Symbol) = has_term(sys.terms, name)
get_term(sys::HamiltonianSystem, name::Symbol) = get_term(sys.terms, name)

"""
    n_particles(sys::HamiltonianSystem) -> Int
    dims(sys::HamiltonianSystem) -> Int
    degrees_of_freedom(sys::HamiltonianSystem) -> Int

Read-only accessors for the problem shape. Prefer these over direct field
access so extensions stay stable across future Hamiltonian shapes.
"""
n_particles(sys::HamiltonianSystem) = sys.n_particles
dims(sys::HamiltonianSystem) = sys.dims
degrees_of_freedom(sys::HamiltonianSystem) = sys.degrees_of_freedom

function Base.show(io::IO, sys::HamiltonianSystem)
    print(
        io,
        "HamiltonianSystem($(sys.n_particles) particles, $(sys.dims)D, $(sys.degrees_of_freedom) DOF)",
    )
end

# Plain-text and LaTeX renderings of the exact canonical Weber Hamiltonian.
# Used when the system is analytic and therefore carries no symbolic expression.
const _WEBER_H_TEXT =
    "H = Σᵢ |pᵢ|²/(2mᵢ) + ½ Σᵢ<ⱼ kᵢⱼ ṙᵢⱼ sᵢⱼ + Σᵢ<ⱼ qᵢqⱼ/rᵢⱼ" *
    "   [kᵢⱼ = qᵢqⱼ/(c²rᵢⱼ), sᵢⱼ = r̂ᵢⱼ·(pᵢ/mᵢ − pⱼ/mⱼ), " *
    "ṙᵢⱼ from (I − GK)ṙ = s]"

const _WEBER_H_LATEX =
    raw"$$H = \sum_i \frac{\lVert \vec p_i\rVert^2}{2 m_i}" *
    raw" + \frac{1}{2}\sum_{i<j} k_{ij}\,\dot r_{ij}\,s_{ij}" *
    raw" + \sum_{i<j} \frac{q_i q_j}{r_{ij}},\qquad" *
    raw" k_{ij} = \frac{q_i q_j}{c^2 r_{ij}},\quad" *
    raw" s_{ij} = \hat r_{ij}\cdot\left(\frac{\vec p_i}{m_i}" *
    raw" - \frac{\vec p_j}{m_j}\right)$$"

function Base.show(io::IO, ::MIME"text/plain", sys::HamiltonianSystem)
    println(io, "HamiltonianSystem")
    println(io, "  Particles: $(sys.n_particles)")
    println(io, "  Dimensions: $(sys.dims)")
    println(io, "  DOF: $(sys.degrees_of_freedom)")
    if has_symbolic_hamiltonian(sys)
        println(io, "  Backing: symbolic (Symbolics.jl)")
        println(io, "  H = $(sys.hamiltonian_symbolic)")
    else
        println(io, "  Backing: analytic (exact canonical Weber)")
        println(io, "  $(_WEBER_H_TEXT)")
    end
end

function Base.show(io::IO, ::MIME"text/latex", sys::HamiltonianSystem)
    if has_symbolic_hamiltonian(sys)
        print(io, latexify(sys.hamiltonian_symbolic))
    else
        print(io, _WEBER_H_LATEX)
    end
end
