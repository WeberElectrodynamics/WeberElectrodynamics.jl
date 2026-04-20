@testset "Hamiltonian term builders" begin
    @testset "weber_term alone reproduces HamiltonianSystem(n, dims)" begin
        # Building via the generic ctor on weber_term output must yield the
        # same compiled EOMs as the convenience ctor (same symbolic expression).
        sys_conv = HamiltonianSystem(2, 2)
        params = [1.0, 1.0, 1.0, -1.0, 100.0, 1.0]  # m1, m2, q1, q2, c, κ12
        q0 = [1.0, 0.0, -1.0, 0.0]
        p0 = [0.0, 0.1, 0.0, -0.1]

        out_conv = zeros(4)
        sys_conv.dp_dt_compiled(out_conv, q0, p0, 0.0, params)

        # Sanity: attractive Coulomb, equal & opposite forces
        @test out_conv[1] < 0  # particle 1 pulled in −x
        @test out_conv[3] > 0  # particle 2 pulled in +x
        @test isapprox(out_conv[1], -out_conv[3]; atol = 1e-12)
    end

    @testset "weber_term(κ) ≡ weber_term(κ=1) + zollner_term(κ)" begin
        # Decomposition identity: splitting the full Weber H into a
        # pure-Weber term plus a Zöllner correction yields the same EOMs.
        @variables x1 y1 x2 y2 x3 y3 px1 py1 px2 py2 px3 py3
        @variables m1 m2 m3 q1 q2 q3 cc k12 k13 k23 tt
        q_syms = [x1, y1, x2, y2, x3, y3]
        p_syms = [px1, py1, px2, py2, px3, py3]
        m_syms = [m1, m2, m3]
        q_charge_syms = [q1, q2, q3]
        kappa_syms = [k12, k13, k23]
        param_syms = vcat(m_syms, q_charge_syms, [cc], kappa_syms)

        H_full = weber_term(
            q_syms, p_syms;
            masses = m_syms, charges = q_charge_syms, c = cc,
            kappas = kappa_syms, n_particles = 3, dims = 2,
        )
        H_pure = weber_term(
            q_syms, p_syms;
            masses = m_syms, charges = q_charge_syms, c = cc,
            kappas = [1.0, 1.0, 1.0], n_particles = 3, dims = 2,
        )
        H_corr = zollner_term(
            q_syms, p_syms;
            masses = m_syms, charges = q_charge_syms, c = cc,
            kappas = kappa_syms, n_particles = 3, dims = 2,
        )

        sys_full = HamiltonianSystem(
            H_full, q_syms, p_syms;
            param_symbols = param_syms, t = tt,
            n_particles = 3, dims = 2,
        )
        sys_comp = HamiltonianSystem(
            H_pure + H_corr, q_syms, p_syms;
            param_symbols = param_syms, t = tt,
            n_particles = 3, dims = 2,
        )

        q0 = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
        p0 = [0.1, 0.2, -0.1, 0.0, 0.0, -0.2]
        pvals = [1.0, 1.0, 0.5, 1.0, -1.0, 0.5, 10.0, 1.3, 1.3, 1.0]

        outA = zeros(6); outB = zeros(6)
        sys_full.dp_dt_compiled(outA, q0, p0, 0.0, pvals)
        sys_comp.dp_dt_compiled(outB, q0, p0, 0.0, pvals)
        @test maximum(abs.(outA .- outB)) < 1e-14

        outA .= 0; outB .= 0
        sys_full.dq_dt_compiled(outA, q0, p0, 0.0, pvals)
        sys_comp.dq_dt_compiled(outB, q0, p0, 0.0, pvals)
        @test maximum(abs.(outA .- outB)) < 1e-14
    end

    @testset "zollner_term vanishes when all κ = 1" begin
        # With all κ = 1 the Zöllner correction must be exactly zero.
        @variables x1 y1 x2 y2 px1 py1 px2 py2
        @variables m1 m2 q1 q2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]

        H_corr = zollner_term(
            q_syms, p_syms;
            masses = [m1, m2], charges = [q1, q2], c = cc,
            kappas = [1.0], n_particles = 2, dims = 2,
        )
        H_simplified = Symbolics.simplify(H_corr)
        @test iszero(H_simplified)
    end
end
