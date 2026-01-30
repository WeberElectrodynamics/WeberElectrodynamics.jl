using Symbolics
using Latexify: latexify

struct WeberSystem{H,QD,PD,QF,PF}
    # Physics configuration
    n_particles::Int
    dims::Int
    masses::Vector{Float64}
    charges::Vector{Float64}
    c::Float64

    # Symbolic representation
    hamiltonian_symbolic::H
    dq_dt_symbolic::QD
    dp_dt_symbolic::PD

    # Compiled functions: func(out, q, p) -> Nothing
    dq_dt_compiled::QF
    dp_dt_compiled::PF

    # Derived
    degrees_of_freedom::Int
end

# =============================================================================
# Internal: Phase Space Symbol Generation
# =============================================================================

"""Generate symbolic variable names for positions and momenta."""
function _generate_phase_space_symbols(n_particles::Int, dims::Int)
    coord_names = [:x, :y, :z]
    momentum_names = [:px, :py, :pz]

    n_total = n_particles * dims
    coordinate_symbols = Vector{Symbol}(undef, n_total)
    momentum_symbols = Vector{Symbol}(undef, n_total)

    idx = 1
    for i = 1:n_particles
        for d = 1:dims
            coordinate_symbols[idx] = Symbol(string(coord_names[d]) * string(i))
            momentum_symbols[idx] = Symbol(string(momentum_names[d]) * string(i))
            idx += 1
        end
    end

    return (coordinate_symbols, momentum_symbols)
end

function _build_weber_hamiltonian(
    q_vars::Vector,      # Symbolic position variables
    p_vars::Vector,      # Symbolic momentum variables
    n_particles::Int,
    dims::Int,
    masses::Vector{Float64},
    charges::Vector{Float64},
    c::Float64,
)
    H = zero(eltype(q_vars))  # Start with symbolic zero

    # Kinetic energy: Σᵢ |pᵢ|² / (2mᵢ)
    for i = 1:n_particles
        p_start = (i - 1) * dims + 1
        p_end = i * dims
        p_squared = sum(p_vars[p_start:p_end] .^ 2)
        H = H + p_squared / (2 * masses[i])
    end

    # Weber potential energy: Σᵢ<ⱼ Uᵢⱼ
    c_squared = c^2
    for i = 1:n_particles
        for j = (i+1):n_particles
            # Position indices for particles i and j
            qi_start = (i - 1) * dims + 1
            qj_start = (j - 1) * dims + 1

            # Momentum indices for particles i and j
            pi_start = (i - 1) * dims + 1
            pj_start = (j - 1) * dims + 1

            # Relative position squared: Σ_d (qi_d - qj_d)²
            r_squared = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                r_squared = r_squared + dq^2
            end
            r = sqrt(r_squared)

            r_dot_v = zero(eltype(q_vars))
            for d = 1:dims
                dq = q_vars[qi_start+d-1] - q_vars[qj_start+d-1]
                dv = p_vars[pi_start+d-1] / masses[i] - p_vars[pj_start+d-1] / masses[j]
                r_dot_v = r_dot_v + dq * dv
            end
            r_dot = r_dot_v / r

            # Weber potential: U = (qi*qj/r)(1 - ṙ²/(2c²))
            k = charges[i] * charges[j]
            U_ij = k / r * (1 - r_dot^2 / (2 * c_squared))
            H = H + U_ij
        end
    end

    return H
end

# =============================================================================
# Constructor
# =============================================================================

function WeberSystem(
    n_particles::Int,
    dims::Int;
    masses::AbstractVector{<:Real},
    charges::AbstractVector{<:Real},
    c::Real = 1.0,
)
    # Validation
    @assert n_particles >= 1 "Must have at least 1 particle"
    @assert dims in (1, 2, 3) "Dimensions must be 1, 2, or 3, got $dims"
    @assert length(masses) == n_particles "masses must have length $n_particles, got $(length(masses))"
    @assert length(charges) == n_particles "charges must have length $n_particles, got $(length(charges))"
    @assert all(m -> m > 0, masses) "All masses must be positive"
    @assert c > 0 "Speed of light must be positive, got $c"

    # Convert to Float64
    masses_f64 = Vector{Float64}(masses)
    charges_f64 = Vector{Float64}(charges)
    c_f64 = Float64(c)

    # Generate symbolic variables
    coordinate_symbols, momentum_symbols = _generate_phase_space_symbols(n_particles, dims)
    q_vars = [Symbolics.variable(sym) for sym in coordinate_symbols]
    p_vars = [Symbolics.variable(sym) for sym in momentum_symbols]

    # Build Weber Hamiltonian symbolically
    hamiltonian_symbolic = _build_weber_hamiltonian(
        q_vars,
        p_vars,
        n_particles,
        dims,
        masses_f64,
        charges_f64,
        c_f64,
    )

    # Derive Hamilton's equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
    dq_dt_symbolic =
        [Symbolics.derivative(hamiltonian_symbolic, p_vars[i]) for i in eachindex(p_vars)]
    dp_dt_symbolic =
        [-Symbolics.derivative(hamiltonian_symbolic, q_vars[i]) for i in eachindex(q_vars)]

    # Compile to fast numeric functions
    # Signature: func(out, q, p) - no params needed, they're baked in
    dq_dt_compiled =
        Symbolics.build_function(dq_dt_symbolic, q_vars, p_vars, expression = Val{false})[2]
    dp_dt_compiled =
        Symbolics.build_function(dp_dt_symbolic, q_vars, p_vars, expression = Val{false})[2]

    degrees_of_freedom = n_particles * dims

    WeberSystem(
        n_particles,
        dims,
        masses_f64,
        charges_f64,
        c_f64,
        hamiltonian_symbolic,
        dq_dt_symbolic,
        dp_dt_symbolic,
        dq_dt_compiled,
        dp_dt_compiled,
        degrees_of_freedom,
    )
end

# =============================================================================
# Display
# =============================================================================

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
    println(io, "  Masses: $(sys.masses)")
    println(io, "  Charges: $(sys.charges)")
    println(io, "  c: $(sys.c)")
    println(io, "  H = $(sys.hamiltonian_symbolic)")
end

# LaTeX display for Jupyter notebooks
function Base.show(io::IO, ::MIME"text/latex", sys::WeberSystem)
    print(io, latexify(sys.hamiltonian_symbolic))
end
