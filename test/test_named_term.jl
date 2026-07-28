@testset "NamedTerm queries" begin
    @testset "convenience ctor exposes the Weber term" begin
        sys = HamiltonianSystem(3, 2)
        @test term_names(sys) == [:weber]
        @test has_term(sys, :weber)
        weber = get_term(sys, :weber)
        @test weber.name === :weber
        @test isequal(weber.H_symbolic, sys.hamiltonian_symbolic)
        @test weber.pair_decomposition !== nothing
    end

    @testset "pair_decomposition matches Weber components" begin
        sys = HamiltonianSystem(2, 2)
        weber = get_term(sys, :weber)
        q = [1.0, 0.2, -0.5, 0.3]
        p = [0.1, 0.2, -0.15, -0.1]
        params = [1.0, 2.0, 1.0, -1.0, 5.0]

        # Reference built from the independent two-body closed-form inverse.
        v1, v2, rdot_ref, r_ref =
            weber_two_body_velocities_2d(q, p, params[1:2], params[3:4], params[5])
        v = [v1[1], v1[2], v2[1], v2[2]]
        coulomb, velocity, rdot = WeberElectrodynamics.compute_pair_weber_components(
            q,
            v,
            1,
            2,
            params[3:4],
            params[5],
            2,
        )

        result = weber.pair_decomposition(q, p, params)
        @test length(result.coulomb) == 1
        @test result.coulomb[1] ≈ coulomb
        @test result.velocity[1] ≈ velocity
        @test result.rdot[1] ≈ rdot ≈ rdot_ref
        @test result.r[1] ≈ r_ref

        # Physical rdot differs from the p/m surrogate once radial motion exists.
        s_naive =
            (
                (q[1] - q[3]) * (p[1] / params[1] - p[3] / params[2]) +
                (q[2] - q[4]) * (p[2] / params[1] - p[4] / params[2])
            ) / r_ref
        @test !isapprox(result.rdot[1], s_naive; rtol = 1e-8)
    end

    @testset "pair_decomposition and kinetic_energy reproduce H in 3D" begin
        sys = HamiltonianSystem(3, 3)
        weber = get_term(sys, :weber)
        q = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.5, 0.2]
        p = [0.1, 0.0, 0.0, -0.1, 0.05, 0.0, 0.0, -0.05, 0.02]
        params = [1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 8.0]

        result = weber.pair_decomposition(q, p, params)
        @test length(result.coulomb) == n_pairs(sys) == 3

        # The velocity-space decomposition must reconstruct the compiled
        # canonical Hamiltonian exactly (Legendre transform).
        KE = weber.kinetic_energy(q, p, params)
        E = KE + sum(result.coulomb .+ result.velocity)
        @test E ≈ sys.hamiltonian_compiled(q, p, 0.0, params) atol = 1e-12

        # Physical kinetic energy is not the canonical Σ|p|²/(2m).
        KE_canonical = sum(sum(p[(i-1)*3+t]^2 for t = 1:3) / (2 * params[i]) for i = 1:3)
        @test !isapprox(KE, KE_canonical; rtol = 1e-10)

        # Pair ordering matches pair_indices.
        for (a, (i, j)) in enumerate(pair_indices(sys))
            r_expected = sqrt(sum((q[(i-1)*3+t] - q[(j-1)*3+t])^2 for t = 1:3))
            @test result.r[a] ≈ r_expected
        end
    end

    @testset "generic ctor defaults to :hamiltonian" begin
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 tt
        q = [x1, y1, x2, y2]
        p = [px1, py1, px2, py2]
        H = sum(p .^ 2) / 2
        sys = HamiltonianSystem(
            H,
            q,
            p;
            param_symbols = [m1, m2],
            t = tt,
            n_particles = 2,
            dims = 2,
        )
        @test term_names(sys) == [:hamiltonian]
        @test_throws KeyError get_term(sys, :missing)
    end

    @testset "generic ctor accepts an explicit terms vector" begin
        @variables x1 y1 x2 y2 px1 py1 px2 py2 m1 m2 q1 q2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]
        kinetic = kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2)
        coulomb = coulomb_term(q_syms; charges = [q1, q2], n_particles = 2, dims = 2)
        sys = HamiltonianSystem(
            kinetic + coulomb,
            q_syms,
            p_syms;
            param_symbols = [m1, m2, q1, q2, cc],
            t = tt,
            n_particles = 2,
            dims = 2,
            terms = [NamedTerm(:kinetic, kinetic), NamedTerm(:coulomb, coulomb)],
        )

        @test term_names(sys) == [:kinetic, :coulomb]
        @test has_term(sys, :kinetic) && has_term(sys, :coulomb)
        @test isequal(get_term(sys, :kinetic).H_symbolic, kinetic)
        @test isequal(get_term(sys, :coulomb).H_symbolic, coulomb)
    end
end
