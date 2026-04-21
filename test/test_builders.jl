@testset "Hamiltonian term builders" begin
    @testset "weber_term alone reproduces HamiltonianSystem(n, dims)" begin
        # Building via the generic ctor on weber_term output must yield the
        # same compiled EOMs as the convenience ctor (same symbolic expression).
        sys_conv = HamiltonianSystem(2, 2)
        params = [1.0, 1.0, 1.0, -1.0, 100.0]  # m1, m2, q1, q2, c
        kappas = [1.0]                         # κ12
        q0 = [1.0, 0.0, -1.0, 0.0]
        p0 = [0.0, 0.1, 0.0, -0.1]

        out_conv = zeros(4)
        sys_conv.dp_dt_compiled(out_conv, q0, p0, 0.0, params, kappas)

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
        param_syms = vcat(m_syms, q_charge_syms, [cc])

        H_full = weber_term(
            q_syms,
            p_syms;
            masses = m_syms,
            charges = q_charge_syms,
            c = cc,
            kappas = kappa_syms,
            n_particles = 3,
            dims = 2,
        )
        H_pure = weber_term(
            q_syms,
            p_syms;
            masses = m_syms,
            charges = q_charge_syms,
            c = cc,
            kappas = [1.0, 1.0, 1.0],
            n_particles = 3,
            dims = 2,
        )
        H_corr = zollner_term(
            q_syms,
            p_syms;
            masses = m_syms,
            charges = q_charge_syms,
            c = cc,
            kappas = kappa_syms,
            n_particles = 3,
            dims = 2,
        )

        sys_full = HamiltonianSystem(
            H_full,
            q_syms,
            p_syms;
            param_symbols = param_syms,
            kappa_symbols = kappa_syms,
            t = tt,
            n_particles = 3,
            dims = 2,
        )
        sys_comp = HamiltonianSystem(
            H_pure + H_corr,
            q_syms,
            p_syms;
            param_symbols = param_syms,
            kappa_symbols = kappa_syms,
            t = tt,
            n_particles = 3,
            dims = 2,
        )

        q0 = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
        p0 = [0.1, 0.2, -0.1, 0.0, 0.0, -0.2]
        pvals = [1.0, 1.0, 0.5, 1.0, -1.0, 0.5, 10.0]  # 3 masses + 3 charges + c
        kvals = [1.3, 1.3, 1.0]                        # κ12, κ13, κ23

        outA = zeros(6);
        outB = zeros(6)
        sys_full.dp_dt_compiled(outA, q0, p0, 0.0, pvals, kvals)
        sys_comp.dp_dt_compiled(outB, q0, p0, 0.0, pvals, kvals)
        @test maximum(abs.(outA .- outB)) < 1e-14

        outA .= 0;
        outB .= 0
        sys_full.dq_dt_compiled(outA, q0, p0, 0.0, pvals, kvals)
        sys_comp.dq_dt_compiled(outB, q0, p0, 0.0, pvals, kvals)
        @test maximum(abs.(outA .- outB)) < 1e-14
    end

    @testset "zollner_term vanishes when all κ = 1" begin
        # With all κ = 1 the Zöllner correction must be exactly zero.
        @variables x1 y1 x2 y2 px1 py1 px2 py2
        @variables m1 m2 q1 q2 cc tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]

        H_corr = zollner_term(
            q_syms,
            p_syms;
            masses = [m1, m2],
            charges = [q1, q2],
            c = cc,
            kappas = [1.0],
            n_particles = 2,
            dims = 2,
        )
        H_simplified = Symbolics.simplify(H_corr)
        @test iszero(H_simplified)
    end

    @testset "kinetic_term yields free-particle EOMs" begin
        # With H = Σ p²/(2m), dq/dt = p/m and dp/dt = 0.
        @variables x1 y1 x2 y2 px1 py1 px2 py2
        @variables m1 m2 tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]
        param_syms = [m1, m2]

        H = kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2)
        sys = HamiltonianSystem(
            H,
            q_syms,
            p_syms;
            param_symbols = param_syms,
            t = tt,
            n_particles = 2,
            dims = 2,
        )

        m1v, m2v = 1.5, 0.25
        q0 = [0.3, -0.4, 1.1, 0.7]
        p0 = [0.2, -0.1, 0.05, 0.4]
        pvals = [m1v, m2v]
        kvals = Float64[]

        out_q = zeros(4)
        sys.dq_dt_compiled(out_q, q0, p0, 0.0, pvals, kvals)
        @test isapprox(
            out_q,
            [p0[1] / m1v, p0[2] / m1v, p0[3] / m2v, p0[4] / m2v];
            atol = 1e-14,
        )

        out_p = zeros(4)
        sys.dp_dt_compiled(out_p, q0, p0, 0.0, pvals, kvals)
        @test all(iszero, out_p)
    end

    @testset "coulomb_term two-body attractive force" begin
        # Opposite charges: equal & opposite forces along the separation axis,
        # magnitude |q₁q₂|/r². Composing kinetic_term + coulomb_term gives a
        # self-contained Hamiltonian with no Weber/Zöllner dependencies.
        @variables x1 y1 x2 y2 px1 py1 px2 py2
        @variables m1 m2 q1 q2 tt
        q_syms = [x1, y1, x2, y2]
        p_syms = [px1, py1, px2, py2]
        param_syms = [m1, m2, q1, q2]

        H =
            kinetic_term(p_syms; masses = [m1, m2], n_particles = 2, dims = 2) +
            coulomb_term(q_syms; charges = [q1, q2], n_particles = 2, dims = 2)
        sys = HamiltonianSystem(
            H,
            q_syms,
            p_syms;
            param_symbols = param_syms,
            t = tt,
            n_particles = 2,
            dims = 2,
        )

        # Place particles at (±1, 0) with r = 2 on x-axis, q₁q₂ = -1.
        q0 = [-1.0, 0.0, 1.0, 0.0]
        p0 = [0.0, 0.3, 0.0, -0.3]
        pvals = [1.0, 1.0, 1.0, -1.0]
        kvals = Float64[]

        out_p = zeros(4)
        sys.dp_dt_compiled(out_p, q0, p0, 0.0, pvals, kvals)

        # Attraction: particle 1 pulled +x, particle 2 pulled −x.
        @test out_p[1] > 0
        @test out_p[3] < 0
        # No y-force (motion along the separation axis).
        @test isapprox(out_p[2], 0.0; atol = 1e-14)
        @test isapprox(out_p[4], 0.0; atol = 1e-14)
        # Equal-and-opposite along x.
        @test isapprox(out_p[1], -out_p[3]; atol = 1e-14)
        # Magnitude |q₁q₂|/r² = 1/4 along the separation axis.
        @test isapprox(abs(out_p[1]), 0.25; atol = 1e-14)
    end

    @testset "kinetic_term + coulomb_term ≡ weber_term(κ=1) at zero momentum" begin
        # Weber's velocity correction is −(q₁q₂/r)·ṙ²/(2c²), which vanishes
        # when p ≡ 0 (so ṙ = 0). Compiled dp/dt from the composed
        # kinetic+coulomb Hamiltonian must match weber_term(κ=1) at p = 0.
        @variables x1 y1 x2 y2 x3 y3 px1 py1 px2 py2 px3 py3
        @variables m1 m2 m3 q1 q2 q3 cc tt
        q_syms = [x1, y1, x2, y2, x3, y3]
        p_syms = [px1, py1, px2, py2, px3, py3]
        m_syms = [m1, m2, m3]
        qc_syms = [q1, q2, q3]

        H_kc =
            kinetic_term(p_syms; masses = m_syms, n_particles = 3, dims = 2) +
            coulomb_term(q_syms; charges = qc_syms, n_particles = 3, dims = 2)
        sys_kc = HamiltonianSystem(
            H_kc,
            q_syms,
            p_syms;
            param_symbols = vcat(m_syms, qc_syms),
            t = tt,
            n_particles = 3,
            dims = 2,
        )

        H_w = weber_term(
            q_syms,
            p_syms;
            masses = m_syms,
            charges = qc_syms,
            c = cc,
            kappas = [1.0, 1.0, 1.0],
            n_particles = 3,
            dims = 2,
        )
        sys_w = HamiltonianSystem(
            H_w,
            q_syms,
            p_syms;
            param_symbols = vcat(m_syms, qc_syms, [cc]),
            t = tt,
            n_particles = 3,
            dims = 2,
        )

        q0 = [0.0, 0.0, 1.0, 0.0, 0.5, 0.9]
        p_zero = zeros(6)
        pvals_kc = [1.0, 1.0, 0.5, 1.0, -1.0, 0.5]
        pvals_w = vcat(pvals_kc, [10.0])  # any finite c works at p = 0
        kvals_empty = Float64[]

        out_kc = zeros(6);
        out_w = zeros(6)
        sys_kc.dp_dt_compiled(out_kc, q0, p_zero, 0.0, pvals_kc, kvals_empty)
        sys_w.dp_dt_compiled(out_w, q0, p_zero, 0.0, pvals_w, kvals_empty)
        @test maximum(abs.(out_kc .- out_w)) < 1e-14

        # dq/dt differs trivially (pᵢ/mᵢ vs weber's identical expression) but
        # should also agree exactly at p = 0 (both identically zero).
        out_kc .= 0;
        out_w .= 0
        sys_kc.dq_dt_compiled(out_kc, q0, p_zero, 0.0, pvals_kc, kvals_empty)
        sys_w.dq_dt_compiled(out_w, q0, p_zero, 0.0, pvals_w, kvals_empty)
        @test maximum(abs.(out_kc .- out_w)) < 1e-14
    end
end
