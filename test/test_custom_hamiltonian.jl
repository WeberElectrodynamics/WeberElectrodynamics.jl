@testset "Custom Hamiltonian via term builders" begin
    @testset "kinetic_term + coulomb_term builds a pure-Coulomb system" begin
        # Build H = Σᵢ ‖pᵢ‖²/2mᵢ + Σᵢ<ⱼ qᵢqⱼ/rᵢⱼ without the
        # Weber velocity correction.
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 qq1 qq2 tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]

        H =
            kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2) +
            coulomb_term(q_syms; charges = [qq1, qq2], n_particles = 2, dims = 2)

        sys = HamiltonianSystem(
            H,
            q_syms,
            p_syms;
            param_symbols = [m1, m2, qq1, qq2],
            t = tt,
            n_particles = 2,
            dims = 2,
        )

        @test term_names(sys) == [:hamiltonian]
        @test n_particles(sys) == 2
        @test dims(sys) == 2

        # Compiled EOMs evaluate to the analytic 2-body Coulomb derivatives.
        out_q = zeros(4)
        out_p = zeros(4)
        q = [1.0, 0.0, -1.0, 0.0]
        p = [0.0, 0.5, 0.0, -0.5]
        params_vec = [1.0, 1.0, 1.0, -1.0]   # m1=m2=1, q1=+1, q2=-1

        sys.dq_dt_compiled(out_q, q, p, 0.0, params_vec)
        sys.dp_dt_compiled(out_p, q, p, 0.0, params_vec)

        # dq/dt = p/m
        @test out_q ≈ [0.0, 0.5, 0.0, -0.5]
        # Pair (q1=+1 at (1,0), q2=-1 at (-1,0)): attractive, r=2.
        # ∂U/∂x1 = -q1q2·(x1-x2)/r³ = -(-1)·2/8 = 1/4
        # ⇒ dp1/dt = -∂U/∂x1 = -1/4 (pulling toward q2 in -x direction).
        # By Newton's 3rd law, dp2/dt is +1/4 in x.
        @test out_p ≈ [-0.25, 0.0, 0.25, 0.0]
    end

    @testset "kinetic + coulomb integrates a periodic 2-body Kepler orbit" begin
        # `HamiltonianProblem` always packs `params = [m₁,m₂,q₁,q₂,c]` (length
        # `2N+1`), so the custom system must declare its `param_symbols` in the
        # same layout — here we include a `cc` symbol even though pure Coulomb
        # doesn't use it.
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 qq1 qq2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]
        H =
            kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2) +
            coulomb_term(q_syms; charges = [qq1, qq2], n_particles = 2, dims = 2)
        sys = HamiltonianSystem(
            H,
            q_syms,
            p_syms;
            param_symbols = [m1, m2, qq1, qq2, cc],
            t = tt,
            n_particles = 2,
            dims = 2,
        )

        m1v, m2v = 1.0, 1.0
        q1v, q2v = 1.0, -1.0
        r0 = 2.0
        M = m1v + m2v
        μ = m1v * m2v / M
        v_circ = sqrt(abs(q1v * q2v) / (μ * r0))
        v = 0.95 * v_circ

        q0 = [-m2v / M * r0, 0.0, m1v / M * r0, 0.0]
        p0 = [0.0, m1v * (-m2v / M * v), 0.0, m2v * (m1v / M * v)]

        # Period of the resulting attractive-Kepler orbit at this energy.
        # H = T + U with U = q1q2/r = -k/r where k = -q1q2 = +1 (attractive).
        # Bound state: E < 0; semi-major axis a = k / (-2E); period 2π√(μa³/k).
        E = 0.5 * (sum(p0[1:2] .^ 2) / m1v + sum(p0[3:4] .^ 2) / m2v) + q1v * q2v / r0
        k = -q1v * q2v
        @assert k > 0 && E < 0 "test setup must be attractive and bound"
        a = k / (-2 * E)
        T = 2π * sqrt(μ * a^3 / k)

        prob = HamiltonianProblem(
            sys,
            (0.0, T),
            q0,
            p0;
            masses = [m1v, m2v],
            charges = [q1v, q2v],
            c = 1.0,            # unused by the pure-Coulomb H, but required by params layout
            dt = T / 10_000,
        )
        sol = solve(prob)
        @test sol.retcode === :Success

        # After one Kepler period the orbit should be near its initial state.
        # Tolerance is loose (1e-3) because the Strang-projection integrator's
        # phase error accumulates over a full period at this dt.
        @test isapprox(sol.q[end], q0; atol = 1e-3)
        @test isapprox(sol.p[end], p0; atol = 1e-3)
    end

    @testset "user-defined NamedTerm with custom name + pair_decomposition" begin
        # Build a system tagged with a custom term whose pair_decomposition
        # returns an arbitrary NamedTuple shape — verify the round-trip.
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 qq1 qq2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]

        H =
            kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2) +
            coulomb_term(q_syms; charges = [qq1, qq2], n_particles = 2, dims = 2)

        # Custom decomposition: per-pair |Δq| only — exercises the closure
        # plumbing with the per-state signature (q, p, params).
        custom_pd = (q, p, params) -> begin
            r = sqrt((q[1] - q[3])^2 + (q[2] - q[4])^2)
            return (separation = [r], sentinel = :ok)
        end

        sys = HamiltonianSystem(
            H,
            q_syms,
            p_syms;
            param_symbols = [m1, m2, qq1, qq2, cc],
            t = tt,
            n_particles = 2,
            dims = 2,
            terms = [NamedTerm(:my_custom_term, H; pair_decomposition = custom_pd)],
        )

        @test term_names(sys) == [:my_custom_term]
        @test has_term(sys, :my_custom_term)
        @test !has_term(sys, :weber)
        custom = get_term(sys, :my_custom_term)
        @test custom.name === :my_custom_term
        @test custom.pair_decomposition !== nothing

        result = custom.pair_decomposition(
            [3.0, 0.0, 0.0, 4.0],
            zeros(4),
            [1.0, 1.0, 1.0, -1.0, 1.0],
        )
        @test result.separation[1] ≈ 5.0
        @test result.sentinel === :ok
        # A custom term may omit the kinetic_energy hook entirely.
        @test custom.kinetic_energy === nothing
    end

    @testset "energy statistics support generic Hamiltonian systems" begin
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 qq1 qq2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]

        H = kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2)
        sys = HamiltonianSystem(
            H,
            q_syms,
            p_syms;
            param_symbols = [m1, m2, qq1, qq2, cc],
            t = tt,
            n_particles = 2,
            dims = 2,
        )

        prob = HamiltonianProblem(
            sys,
            (0.0, 0.1),
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 0.5, 0.0, -0.5];
            masses = [1.0, 1.0],
            charges = [0.0, 0.0],
            c = 1.0,
            dt = 0.05,
        )
        sol = solve(prob)
        energy = compute_energy_timeseries(sol)

        @test energy.n_pairs == 1
        @test isempty(energy.pair_energies)
        @test energy.total_energy ≈ energy.kinetic_energy
        @test all(abs.(energy.total_potential_energy) .< 1e-14)
        @test all(energy.hamiltonian_validation_error .< 1e-14)
    end
end
