"""
    center_of_mass_frame(q, p, masses, dims) -> NamedTuple

Return `(q, p)` shifted into the centre-of-mass frame.

Positions are translated so `sum(m_i * r_i) == 0`. Momenta are shifted by the
centre-of-mass velocity so `sum(p_i) == 0`. Inputs are not mutated.
"""
function center_of_mass_frame(
    q::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    masses_in::AbstractVector{<:Real},
    dims_in::Integer,
)
    dims_value = Int(dims_in)
    dims_value >= 1 || throw(ArgumentError("dims must be positive, got $dims_value"))

    n = length(masses_in)
    dof = n * dims_value
    length(q) == dof ||
        throw(ArgumentError("q must have length $dof for $n particles in $dims_value D"))
    length(p) == dof ||
        throw(ArgumentError("p must have length $dof for $n particles in $dims_value D"))

    masses_vec = Vector{Float64}(masses_in)
    all(m -> m > 0, masses_vec) || throw(ArgumentError("all masses must be positive"))

    total_mass = sum(masses_vec)
    q_out = Vector{Float64}(q)
    p_out = Vector{Float64}(p)

    com = zeros(Float64, dims_value)
    total_p = zeros(Float64, dims_value)
    @inbounds for i = 1:n
        offset = (i - 1) * dims_value
        mi = masses_vec[i]
        for d = 1:dims_value
            idx = offset + d
            com[d] += mi * q_out[idx]
            total_p[d] += p_out[idx]
        end
    end
    @inbounds for d = 1:dims_value
        com[d] /= total_mass
        total_p[d] /= total_mass
    end

    @inbounds for i = 1:n
        offset = (i - 1) * dims_value
        mi = masses_vec[i]
        for d = 1:dims_value
            idx = offset + d
            q_out[idx] -= com[d]
            p_out[idx] -= mi * total_p[d]
        end
    end

    return (q = q_out, p = p_out)
end

"""
    two_body_initial_conditions(masses, charges; separation, dims=2,
                                velocity_scale=1.0,
                                radial_velocity=0.0, c=nothing) -> NamedTuple

Construct centre-of-mass two-body initial conditions in **canonical** variables.

The particles are placed on the x axis with separation `separation`. For
`dims >= 2`, the transverse relative momentum is
`velocity_scale * sqrt(mu * abs(q1*q2) / separation)`, where `mu` is the
reduced mass. Set `velocity_scale=1` for the circular Coulomb scale.

`radial_velocity` is a **physical** radial velocity ṙ. Because canonical
momentum in Weber electrodynamics is `p_i = ∂L/∂v_i`, the conjugate radial
momentum is

```
p_r = (mu - q1*q2/(r*c^2)) * rdot
```

not `mu * rdot`. Supplying a nonzero `radial_velocity` therefore requires `c`;
omitting it throws. Only the radial direction is affected — the transverse
momentum is unchanged, and for `radial_velocity = 0` canonical and kinetic
momenta coincide so `c` is irrelevant and may be omitted.

Throws if the pair sits exactly at Weber's critical radius
`rho = q1*q2/(mu*c^2)`, where `p_r` carries no information about ṙ.

Returns `(q, p, masses, charges)`.
"""
function two_body_initial_conditions(
    masses_in::AbstractVector{<:Real},
    charges_in::AbstractVector{<:Real};
    separation::Real,
    dims::Integer = 2,
    velocity_scale::Real = 1.0,
    radial_velocity::Real = 0.0,
    c::Union{Nothing,Real} = nothing,
)
    length(masses_in) == 2 || throw(ArgumentError("masses must have length 2"))
    length(charges_in) == 2 || throw(ArgumentError("charges must have length 2"))

    dims_value = Int(dims)
    1 <= dims_value <= 3 || throw(ArgumentError("dims must be 1, 2, or 3, got $dims_value"))
    separation_value = Float64(separation)
    separation_value > 0 ||
        throw(ArgumentError("separation must be positive, got $separation"))

    masses_vec = Vector{Float64}(masses_in)
    charges_vec = Vector{Float64}(charges_in)
    all(m -> m > 0, masses_vec) || throw(ArgumentError("all masses must be positive"))

    mi, mj = masses_vec
    total_mass = mi + mj
    reduced_mass = mi * mj / total_mass
    pair_coupling = charges_vec[1] * charges_vec[2]

    radial_velocity_value = Float64(radial_velocity)
    if radial_velocity_value != 0
        isnothing(c) && throw(
            ArgumentError(
                "c is required when radial_velocity != 0: canonical radial momentum " *
                "is p_r = (mu - q1*q2/(r*c^2)) * rdot, not mu * rdot",
            ),
        )
        Float64(c) > 0 || throw(ArgumentError("c must be positive, got $c"))
    end

    q = zeros(Float64, 2 * dims_value)
    p = zeros(Float64, 2 * dims_value)

    # q_rel = q1 - q2 points along +x and has norm `separation`.
    q[1] = mj / total_mass * separation_value
    q[dims_value+1] = -mi / total_mass * separation_value

    rel_p = zeros(Float64, dims_value)
    if radial_velocity_value != 0
        c_value = Float64(c)
        radial_inertia = reduced_mass - pair_coupling / (separation_value * c_value^2)
        if abs(radial_inertia) <= 1e-12 * reduced_mass
            rho = pair_coupling / (reduced_mass * c_value^2)
            throw(
                ArgumentError(
                    "separation = $separation_value is at Weber's critical radius " *
                    "rho = $rho, where canonical radial momentum does not determine " *
                    "a physical radial velocity",
                ),
            )
        end
        rel_p[1] = radial_inertia * radial_velocity_value
    end
    if dims_value >= 2
        circular_p = sqrt(reduced_mass * abs(pair_coupling) / separation_value)
        rel_p[2] = Float64(velocity_scale) * circular_p
    end

    @inbounds for d = 1:dims_value
        p[d] = rel_p[d]
        p[dims_value+d] = -rel_p[d]
    end

    framed = center_of_mass_frame(q, p, masses_vec, dims_value)
    return (q = framed.q, p = framed.p, masses = masses_vec, charges = charges_vec)
end

"""
    polygon_initial_conditions(n; radius, mass=1.0, charge_magnitude=1.0,
                               energy_ratio=0.5, speed=nothing,
                               clockwise=false) -> NamedTuple

Construct planar regular-polygon initial conditions with alternating charges.

Particles sit at regular `n`-gon vertices with circumradius `radius`. Momenta
are tangential. If `speed` is omitted, the speed is chosen from
`energy_ratio = T0 / abs(U0)` using the initial Coulomb potential.

Returns `(q, p, masses, charges)`.
"""
function polygon_initial_conditions(
    n_in::Integer;
    radius::Real,
    mass::Real = 1.0,
    charge_magnitude::Real = 1.0,
    energy_ratio::Real = 0.5,
    speed::Union{Nothing,Real} = nothing,
    clockwise::Bool = false,
)
    n = Int(n_in)
    n >= 2 || throw(ArgumentError("n must be at least 2, got $n"))
    radius_value = Float64(radius)
    radius_value > 0 || throw(ArgumentError("radius must be positive, got $radius"))
    mass_value = Float64(mass)
    mass_value > 0 || throw(ArgumentError("mass must be positive, got $mass"))
    Float64(energy_ratio) >= 0 ||
        throw(ArgumentError("energy_ratio must be non-negative, got $energy_ratio"))

    masses_vec = fill(mass_value, n)
    qmag = Float64(charge_magnitude)
    charges_vec = [isodd(i) ? qmag : -qmag for i = 1:n]

    q = Vector{Float64}(undef, 2n)
    angles = [2pi * (i - 1) / n for i = 1:n]
    @inbounds for i = 1:n
        theta = angles[i]
        q[2i-1] = radius_value * cos(theta)
        q[2i] = radius_value * sin(theta)
    end

    potential = 0.0
    @inbounds for i = 1:n
        qi = 2i - 1
        for j = (i+1):n
            qj = 2j - 1
            dx = q[qi] - q[qj]
            dy = q[qi+1] - q[qj+1]
            r = sqrt(dx * dx + dy * dy)
            potential += charges_vec[i] * charges_vec[j] / r
        end
    end

    speed_value =
        isnothing(speed) ?
        sqrt(2 * Float64(energy_ratio) * abs(potential) / (n * mass_value)) : Float64(speed)
    speed_value >= 0 || throw(ArgumentError("speed must be non-negative, got $speed"))

    sign = clockwise ? -1.0 : 1.0
    p = Vector{Float64}(undef, 2n)
    @inbounds for i = 1:n
        theta = angles[i]
        p[2i-1] = mass_value * speed_value * sign * (-sin(theta))
        p[2i] = mass_value * speed_value * sign * cos(theta)
    end

    framed = center_of_mass_frame(q, p, masses_vec, 2)
    return (q = framed.q, p = framed.p, masses = masses_vec, charges = charges_vec)
end

"""
    rigid_rotation_initial_conditions(q, masses; angular_velocity,
                                      dims=3, charges=zeros(n)) -> NamedTuple

Return centre-of-mass initial conditions for a rigidly rotating configuration.

The input `q` is flattened particle positions. In 3D, `angular_velocity` is a
3-vector and momenta are `p_i = m_i * (omega x r_i)`. In 2D,
`angular_velocity` may be a scalar z-angular velocity.

Returns `(q, p, masses, charges)`.
"""
function rigid_rotation_initial_conditions(
    q_in::AbstractVector{<:Real},
    masses_in::AbstractVector{<:Real};
    angular_velocity,
    dims::Integer = 3,
    charges::AbstractVector{<:Real} = zeros(Float64, length(masses_in)),
)
    dims_value = Int(dims)
    dims_value in (2, 3) ||
        throw(ArgumentError("rigid_rotation_initial_conditions supports dims=2 or dims=3"))

    n = length(masses_in)
    dof = n * dims_value
    length(q_in) == dof ||
        throw(ArgumentError("q must have length $dof for $n particles in $dims_value D"))
    length(charges) == n || throw(ArgumentError("charges must have length $n"))

    masses_vec = Vector{Float64}(masses_in)
    charges_vec = Vector{Float64}(charges)
    all(m -> m > 0, masses_vec) || throw(ArgumentError("all masses must be positive"))

    q_centered = center_of_mass_frame(q_in, zeros(Float64, dof), masses_vec, dims_value).q
    p = zeros(Float64, dof)

    if dims_value == 2
        omega =
            angular_velocity isa Real ? Float64(angular_velocity) :
            Float64(angular_velocity[end])
        @inbounds for i = 1:n
            offset = (i - 1) * 2
            x = q_centered[offset+1]
            y = q_centered[offset+2]
            mi = masses_vec[i]
            p[offset+1] = mi * (-omega * y)
            p[offset+2] = mi * (omega * x)
        end
    else
        length(angular_velocity) == 3 ||
            throw(ArgumentError("3D angular_velocity must have length 3"))
        wx = Float64(angular_velocity[1])
        wy = Float64(angular_velocity[2])
        wz = Float64(angular_velocity[3])
        @inbounds for i = 1:n
            offset = (i - 1) * 3
            x = q_centered[offset+1]
            y = q_centered[offset+2]
            z = q_centered[offset+3]
            mi = masses_vec[i]
            p[offset+1] = mi * (wy * z - wz * y)
            p[offset+2] = mi * (wz * x - wx * z)
            p[offset+3] = mi * (wx * y - wy * x)
        end
    end

    framed = center_of_mass_frame(q_centered, p, masses_vec, dims_value)
    return (q = framed.q, p = framed.p, masses = masses_vec, charges = charges_vec)
end
