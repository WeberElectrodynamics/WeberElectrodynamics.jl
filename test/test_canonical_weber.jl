# Independent checks of the exact canonical Weber system.
#
# Every expectation here is derived from the Weber Lagrangian by an independent
# route — the forward map p(q,v) = ∂L/∂v, finite differences of H, or the
# closed-form two-particle inverse — never from the routine under test. The same
# identities are verified symbolically in
# `papers/Computational-Weber-Electrodynamics/verify_formulas.py`.

using LinearAlgebra

_pairs(n) = [(i, j) for i = 1:n for j = (i+1):n]

# Forward map p(q, v) = ∂L/∂v, read straight off the Weber Lagrangian:
#   p_i = m_i v_i − Σ_{j≠i} (q_i q_j / c²) ṙ_ij (r_i − r_j)/r_ij²
function lagrangian_momenta(q, v, ms, chs, c, n, d)
    p = [ms[i] * v[(i-1)*d+t] for i = 1:n for t = 1:d]
    for (i, j) in _pairs(n)
        dq = [q[(i-1)*d+t] - q[(j-1)*d+t] for t = 1:d]
        r = norm(dq)
        rh = dq ./ r
        rd = dot(rh, [v[(i-1)*d+t] - v[(j-1)*d+t] for t = 1:d])
        k = chs[i] * chs[j] / (c^2 * r)
        for t = 1:d
            p[(i-1)*d+t] -= k * rd * rh[t]
            p[(j-1)*d+t] += k * rd * rh[t]
        end
    end
    return p
end

# Velocity-space energy E = Σ ½ m|v|² + Σ q_i q_j/r (1 − ṙ²/(2c²)).
function velocity_space_energy(q, v, ms, chs, c, n, d)
    E = sum(0.5 * ms[i] * sum(v[(i-1)*d+t]^2 for t = 1:d) for i = 1:n)
    for (i, j) in _pairs(n)
        dq = [q[(i-1)*d+t] - q[(j-1)*d+t] for t = 1:d]
        r = norm(dq)
        rd = dot(dq, [v[(i-1)*d+t] - v[(j-1)*d+t] for t = 1:d]) / r
        E += chs[i] * chs[j] / r * (1 - rd^2 / (2 * c^2))
    end
    return E
end

# Corrected canonical momentum rate, written out independently of the package:
#   ṗ_i = Σ_{j≠i} q_i q_j/r² [ r̂ (1 + 3ṙ²/(2c²)) − ṙ (v_i − v_j)/c² ]
function canonical_momentum_rate(q, v, ms, chs, c, n, d)
    out = zeros(n * d)
    for (i, j) in _pairs(n)
        dq = [q[(i-1)*d+t] - q[(j-1)*d+t] for t = 1:d]
        r = norm(dq)
        rh = dq ./ r
        dv = [v[(i-1)*d+t] - v[(j-1)*d+t] for t = 1:d]
        rd = dot(rh, dv)
        pref = chs[i] * chs[j] / r^2
        for t = 1:d
            term = pref * (rh[t] * (1 + 3 * rd^2 / (2 * c^2)) - rd * dv[t] / c^2)
            out[(i-1)*d+t] += term
            out[(j-1)*d+t] -= term
        end
    end
    return out
end

const CANONICAL_CASES = [
    (
        name = "two-body 2D",
        n = 2,
        d = 2,
        q = [1.0, 0.3, -1.2, 0.7],
        v = [0.21, -0.13, -0.09, 0.17],
        ms = [1.3, 0.7],
        chs = [0.9, -1.1],
        c = 5.0,
    ),
    (
        name = "two-body 2D like charges",
        n = 2,
        d = 2,
        q = [1.0, 0.3, -1.2, 0.7],
        v = [0.21, -0.13, -0.09, 0.17],
        ms = [1.0, 1.0],
        chs = [1.0, 1.0],
        c = 6.0,
    ),
    (
        name = "three-body 2D",
        n = 3,
        d = 2,
        q = [4.5, 2.4, 1.0, -1.0, -2.0, 2.6],
        v = [0.21, -0.13, -0.09, 0.17, 0.05, -0.11],
        ms = [1.3, 0.7, 1.9],
        chs = [0.9, -1.1, 0.6],
        c = 5.0,
    ),
    (
        name = "four-body 3D",
        n = 4,
        d = 3,
        q = [1.0, 0.0, 0.0, -1.0, 0.2, 0.1, 0.3, 1.4, -0.5, -0.4, -1.1, 0.9],
        v = [0.05, 0.11, -0.03, -0.07, 0.02, 0.08, 0.12, -0.06, 0.01, -0.02, 0.04, -0.09],
        ms = [1.0, 1.5, 0.8, 1.2],
        chs = [1.0, -1.0, 0.5, -0.5],
        c = 6.0,
    ),
    (
        name = "two-body 1D",
        n = 2,
        d = 1,
        q = [0.0, 2.5],
        v = [0.3, -0.2],
        ms = [1.0, 2.0],
        chs = [1.0, 1.0],
        c = 4.0,
    ),
]

@testset "Canonical Weber Hamiltonian" begin

    @testset "canonical momentum is ∂L/∂v, not m·v" begin
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)
            kinetic = [ms[i] * v[(i-1)*d+t] for i = 1:n for t = 1:d]
            # Every case here has nonzero pair radial velocity, so canonical and
            # kinetic momentum must genuinely differ.
            @test maximum(abs.(p .- kinetic)) > 1e-6
        end
    end

    @testset "physical_velocities inverts the Lagrangian forward map" begin
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            par = vcat(ms, chs, [c])
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)
            v_rec = physical_velocities(q, p, par; n_particles = n, dims = d)
            @test maximum(abs.(v_rec .- v)) < 1e-12
        end
    end

    @testset "H is the exact Legendre transform of the Weber Lagrangian" begin
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            sys = HamiltonianSystem(n, d)
            par = vcat(ms, chs, [c])
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)
            H = sys.hamiltonian_compiled(q, p, 0.0, par)
            E = velocity_space_energy(q, v, ms, chs, c, n, d)
            @test H ≈ E atol = 1e-12
        end
    end

    @testset "H differs from the naive p/m Hamiltonian" begin
        # Guards against regression to the pre-correction Hamiltonian.
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            sys = HamiltonianSystem(n, d)
            par = vcat(ms, chs, [c])
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)

            H_naive = sum(sum(p[(i-1)*d+t]^2 for t = 1:d) / (2 * ms[i]) for i = 1:n)
            for (i, j) in _pairs(n)
                dq = [q[(i-1)*d+t] - q[(j-1)*d+t] for t = 1:d]
                r = norm(dq)
                s = dot(dq, [p[(i-1)*d+t] / ms[i] - p[(j-1)*d+t] / ms[j] for t = 1:d]) / r
                H_naive += chs[i] * chs[j] / r * (1 - s^2 / (2 * c^2))
            end
            @test !isapprox(sys.hamiltonian_compiled(q, p, 0.0, par), H_naive; rtol = 1e-9)
        end
    end

    @testset "first canonical equation: q̇ = v" begin
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            sys = HamiltonianSystem(n, d)
            par = vcat(ms, chs, [c])
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)
            out = zeros(n * d)
            sys.dq_dt_compiled(out, q, p, 0.0, par)
            @test maximum(abs.(out .- v)) < 1e-12
        end
    end

    @testset "second canonical equation matches the corrected closed form" begin
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            sys = HamiltonianSystem(n, d)
            par = vcat(ms, chs, [c])
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)
            out = zeros(n * d)
            sys.dp_dt_compiled(out, q, p, 0.0, par)
            expected = canonical_momentum_rate(q, v, ms, chs, c, n, d)
            @test maximum(abs.(out .- expected)) < 1e-12
        end
    end

    @testset "finite differences of H reproduce both canonical equations" begin
        h = 1e-6
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs, c) = case
            sys = HamiltonianSystem(n, d)
            par = vcat(ms, chs, [c])
            p = lagrangian_momenta(q, v, ms, chs, c, n, d)

            dq = zeros(n * d)
            dp = zeros(n * d)
            sys.dq_dt_compiled(dq, q, p, 0.0, par)
            sys.dp_dt_compiled(dp, q, p, 0.0, par)

            for i = 1:(n*d)
                qp = copy(q);
                qp[i] += h
                qm = copy(q);
                qm[i] -= h
                num_pdot =
                    -(
                        sys.hamiltonian_compiled(qp, p, 0.0, par) -
                        sys.hamiltonian_compiled(qm, p, 0.0, par)
                    ) / (2h)
                @test num_pdot ≈ dp[i] atol = 1e-7

                pp = copy(p);
                pp[i] += h
                pm = copy(p);
                pm[i] -= h
                num_qdot =
                    (
                        sys.hamiltonian_compiled(q, pp, 0.0, par) -
                        sys.hamiltonian_compiled(q, pm, 0.0, par)
                    ) / (2h)
                @test num_qdot ≈ dq[i] atol = 1e-7
            end
        end
    end

    @testset "canonical equations reproduce the mechanical Weber force" begin
        # m_i a_i = ṗ_i + d(α_i)/dt must equal Weber's force law. The α time
        # derivative is taken by central differences along the exact flow.
        n, d = 2, 2
        ms = [1.3, 0.7]
        chs = [0.9, -1.1]
        c = 5.0
        par = vcat(ms, chs, [c])
        sys = HamiltonianSystem(n, d)

        q0 = [1.0, 0.3, -1.2, 0.7]
        v0 = [0.21, -0.13, -0.09, 0.17]
        p0 = lagrangian_momenta(q0, v0, ms, chs, c, n, d)

        # Advance the exact canonical flow with a tiny symmetric step to get dv/dt.
        h = 1e-6
        function flow(qs, ps)
            dq = zeros(n * d)
            dp = zeros(n * d)
            sys.dq_dt_compiled(dq, qs, ps, 0.0, par)
            sys.dp_dt_compiled(dp, qs, ps, 0.0, par)
            return dq, dp
        end
        dq0, dp0 = flow(q0, p0)
        q_plus = q0 .+ h .* dq0
        p_plus = p0 .+ h .* dp0
        q_minus = q0 .- h .* dq0
        p_minus = p0 .- h .* dp0
        v_plus = physical_velocities(q_plus, p_plus, par; n_particles = n, dims = d)
        v_minus = physical_velocities(q_minus, p_minus, par; n_particles = n, dims = d)
        accel = (v_plus .- v_minus) ./ (2h)

        # Weber's force on particle 1 from the vector form.
        r_vec = [q0[1] - q0[3], q0[2] - q0[4]]
        r = norm(r_vec)
        rh = r_vec ./ r
        v_rel = [v0[1] - v0[3], v0[2] - v0[4]]
        a_rel = [accel[1] - accel[3], accel[2] - accel[4]]
        weber_force =
            chs[1] * chs[2] / r^2 .* rh .*
            (1 + (dot(v_rel, v_rel) + dot(r_vec, a_rel) - 1.5 * dot(rh, v_rel)^2) / c^2)

        @test ms[1] * accel[1] ≈ weber_force[1] atol = 1e-6
        @test ms[1] * accel[2] ≈ weber_force[2] atol = 1e-6
        # Newton's third law in strong form.
        @test ms[1] * accel[1] ≈ -ms[2] * accel[3] atol = 1e-6
        @test ms[1] * accel[2] ≈ -ms[2] * accel[4] atol = 1e-6
    end

    @testset "Coulomb limit c → ∞" begin
        for case in CANONICAL_CASES
            (; name, n, d, q, v, ms, chs) = case
            sys = HamiltonianSystem(n, d)
            par = vcat(ms, chs, [1e9])
            p = lagrangian_momenta(q, v, ms, chs, 1e9, n, d)

            H_coulomb = sum(sum(p[(i-1)*d+t]^2 for t = 1:d) / (2 * ms[i]) for i = 1:n)
            for (i, j) in _pairs(n)
                H_coulomb +=
                    chs[i] * chs[j] / norm([q[(i-1)*d+t] - q[(j-1)*d+t] for t = 1:d])
            end
            @test sys.hamiltonian_compiled(q, p, 0.0, par) ≈ H_coulomb atol = 1e-9

            dq = zeros(n * d)
            sys.dq_dt_compiled(dq, q, p, 0.0, par)
            v_coulomb = [p[(i-1)*d+t] / ms[i] for i = 1:n for t = 1:d]
            @test maximum(abs.(dq .- v_coulomb)) < 1e-9
        end
    end

    @testset "zero radial velocity limit" begin
        # Rigid translation: every pair separation is constant, so all ṙ vanish
        # and canonical momentum coincides with kinetic momentum.
        n, d = 3, 2
        ms = [1.3, 0.7, 1.9]
        chs = [0.9, -1.1, 0.6]
        c = 5.0
        par = vcat(ms, chs, [c])
        q = [4.5, 2.4, 1.0, -1.0, -2.0, 2.6]
        v = repeat([0.2, -0.1], 3)
        sys = HamiltonianSystem(n, d)

        p = lagrangian_momenta(q, v, ms, chs, c, n, d)
        kinetic = [ms[i] * v[(i-1)*d+t] for i = 1:n for t = 1:d]
        @test maximum(abs.(p .- kinetic)) < 1e-14

        H = sys.hamiltonian_compiled(q, p, 0.0, par)
        H_coulomb = sum(sum(p[(i-1)*d+t]^2 for t = 1:d) / (2 * ms[i]) for i = 1:n)
        for (i, j) in _pairs(n)
            H_coulomb += chs[i] * chs[j] / norm([q[(i-1)*d+t] - q[(j-1)*d+t] for t = 1:d])
        end
        @test H ≈ H_coulomb atol = 1e-12
    end

    @testset "translation and rotation invariance" begin
        n, d = 3, 2
        ms = [1.3, 0.7, 1.9]
        chs = [0.9, -1.1, 0.6]
        c = 5.0
        par = vcat(ms, chs, [c])
        sys = HamiltonianSystem(n, d)
        q = [4.5, 2.4, 1.0, -1.0, -2.0, 2.6]
        p = [0.31, -0.22, -0.17, 0.28, -0.14, -0.06]
        H0 = sys.hamiltonian_compiled(q, p, 0.0, par)

        shift = [1.7, -0.9]
        q_shift = [q[(i-1)*d+t] + shift[t] for i = 1:n for t = 1:d]
        @test sys.hamiltonian_compiled(q_shift, p, 0.0, par) ≈ H0 atol = 1e-12

        θ = 0.7
        ct, st = cos(θ), sin(θ)
        rot(u) = [
            f(i) for i = 1:n for f in (
                i -> ct * u[(i-1)*d+1] - st * u[(i-1)*d+2],
                i -> st * u[(i-1)*d+1] + ct * u[(i-1)*d+2],
            )
        ]
        @test sys.hamiltonian_compiled(rot(q), rot(p), 0.0, par) ≈ H0 atol = 1e-12
    end

    @testset "two-particle scalar inverse p_r = (μ − q₁q₂/(rc²)) ṙ" begin
        ms = [1.3, 0.7]
        chs = [0.9, -1.1]
        c = 5.0
        q = [1.0, 0.3, -1.2, 0.7]
        v = [0.21, -0.13, -0.09, 0.17]
        p = lagrangian_momenta(q, v, ms, chs, c, 2, 2)

        r_vec = [q[1] - q[3], q[2] - q[4]]
        r = norm(r_vec)
        rh = r_vec ./ r
        rdot = dot(rh, [v[1] - v[3], v[2] - v[4]])
        mu = ms[1] * ms[2] / (ms[1] + ms[2])
        p_r = mu * dot(rh, [p[1] / ms[1] - p[3] / ms[2], p[2] / ms[1] - p[4] / ms[2]])

        @test p_r ≈ (mu - chs[1] * chs[2] / (r * c^2)) * rdot atol = 1e-12
    end

    @testset "behaviour at Weber's critical radius" begin
        # For a like-charge pair the canonical mass matrix is singular exactly at
        # ρ = q₁q₂/(μc²). Evaluation there must fail loudly, not silently return
        # garbage. Away from ρ (including below it) it must stay finite.
        ms = [1.0, 1.0]
        chs = [1.0, 1.0]
        c = 1.0
        mu = 0.5
        rho = chs[1] * chs[2] / (mu * c^2)          # = 2.0
        par = vcat(ms, chs, [c])
        sys = HamiltonianSystem(2, 1)

        p = [0.1, -0.1]
        @test_throws WeberCriticalRadiusError sys.hamiltonian_compiled(
            [0.0, rho],
            p,
            0.0,
            par,
        )

        err = try
            physical_velocities([0.0, rho], p, par; n_particles = 2, dims = 1)
            nothing
        catch e
            e
        end
        @test err isa WeberCriticalRadiusError
        @test err.pair == (1, 2)
        @test err.r ≈ rho
        @test err.rho ≈ rho
        @test occursin("critical radius", sprint(showerror, err))

        # Above ρ: ordinary positive radial inertia.
        v_above = physical_velocities([0.0, 2 * rho], p, par; n_particles = 2, dims = 1)
        @test all(isfinite, v_above)

        # Below ρ: the effective radial inertia is negative but finite, which is
        # the sub-critical "negative inertia" regime Weber described.
        v_below = physical_velocities([0.0, rho / 2], p, par; n_particles = 2, dims = 1)
        @test all(isfinite, v_below)
        rdot_below = v_below[1] - v_below[2]
        s_below = p[1] / ms[1] - p[2] / ms[2]
        @test sign(rdot_below) == -sign(s_below)
    end

    @testset "energy statistics reconstruct the compiled Hamiltonian" begin
        prob = make_weber_problem(tspan = (0.0, 0.5), dt = 0.001)
        sol = solve(prob)
        energy = compute_energy_timeseries(sol)

        # The manual velocity-space decomposition must match H(q,p) to machine
        # precision — this is the Legendre transform, checked per timestep.
        @test maximum(energy.hamiltonian_validation_error) < 1e-10

        # Reported kinetic energy is physical, so it differs from Σ|p|²/(2m).
        ms = masses(prob)
        d = dims(prob)
        canonical_KE = [
            sum(sum(p[(i-1)*d+t]^2 for t = 1:d) / (2 * ms[i]) for i = 1:n_particles(prob)) for p in sol.p
        ]
        @test !isapprox(energy.kinetic_energy[end], canonical_KE[end]; rtol = 1e-8)
    end
end
