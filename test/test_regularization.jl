@testset "Regularization" begin
    function make_orbit_problem(dims::Int;
        dt::Float64,
        t_end::Float64,
        v_scale::Float64,
        regularization_enabled::Bool,
        r_on::Float64,
        r_off::Float64,
    )
        m1, m2 = 1.0, 0.1
        q1, q2 = sqrt(0.1), -sqrt(0.1)
        c = 4.0

        r0 = 2.0
        M = m1 + m2
        k = q1 * q2
        v_circ = sqrt(abs(k) * M / (m1 * m2 * r0))

        if dims == 1
            q0 = [-0.2, 0.2]
            p0 = [0.0, 0.0]
            masses = [1.0, 1.0]
            charges = [0.2, -0.2]
        elseif dims == 2
            q0 = [-m2 / M * r0, 0.0, m1 / M * r0, 0.0]
            p0 = [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, m2 * (m1 / M * v_circ * v_scale)]
            masses = [m1, m2]
            charges = [q1, q2]
        else
            q0 = [-m2 / M * r0, 0.0, 0.0, m1 / M * r0, 0.0, 0.0]
            p0 = [0.0, m1 * (-m2 / M * v_circ * v_scale), 0.0, 0.0, m2 * (m1 / M * v_circ * v_scale), 0.0]
            masses = [m1, m2]
            charges = [q1, q2]
        end

        system = WeberSystem(2, dims)
        return WeberProblem(
            system,
            (0.0, t_end),
            q0,
            p0;
            masses = masses,
            charges = charges,
            c = c,
            dt = dt,
            regularization_enabled = regularization_enabled,
            regularization_r_on = r_on,
            regularization_r_off = r_off,
            regularization_max_substeps = 256,
        )
    end

    state_error(sol_a::WeberSolution, sol_b::WeberSolution) =
        norm(sol_a.q[end] - sol_b.q[end]) + norm(sol_a.p[end] - sol_b.p[end])

    @testset "Transform identities" begin
        @testset "Levi-Civita" begin
            rb = WeberElectrodynamics.RegularizationBuffers(2, 2, 0.1, 0.2)
            q_rel = zeros(2)
            p_rel = zeros(2)

            for _ = 1:32
                u1 = randn()
                u2 = randn()
                x = u1^2 - u2^2
                y = 2u1 * u2
                r = sqrt(x^2 + y^2)
                rho2 = u1^2 + u2^2
                @test r ≈ rho2 rtol = 1e-12 atol = 1e-12

                J = [2u1 -2u2; 2u2 2u1]
                @test J * transpose(J) ≈ 4rho2 * I(2) rtol = 1e-12 atol = 1e-12

                rb.rel_q[1] = x
                rb.rel_q[2] = y
                rb.rel_p[1] = randn()
                rb.rel_p[2] = randn()

                WeberElectrodynamics._lc_lift!(rb)
                WeberElectrodynamics._lc_project!(q_rel, p_rel, rb.lc_u, rb.lc_U)

                @test q_rel[1] ≈ x rtol = 1e-10 atol = 1e-10
                @test q_rel[2] ≈ y rtol = 1e-10 atol = 1e-10
            end
        end

        @testset "KS" begin
            rb = WeberElectrodynamics.RegularizationBuffers(2, 3, 0.1, 0.2)
            J = rb.ks_J

            for _ = 1:32
                u = randn(4)
                x1 = u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2
                x2 = 2 * (u[1] * u[2] - u[3] * u[4])
                x3 = 2 * (u[1] * u[3] + u[2] * u[4])
                r = sqrt(x1^2 + x2^2 + x3^2)
                rho2 = sum(abs2, u)
                @test r ≈ rho2 rtol = 1e-12 atol = 1e-12

                rb.ks_u .= u
                WeberElectrodynamics._ks_jacobian!(J, rb.ks_u)
                @test J * transpose(J) ≈ 4rho2 * I(3) rtol = 1e-12 atol = 1e-12

                p = randn(3)
                rb.ks_U[1] = J[1, 1] * p[1] + J[2, 1] * p[2] + J[3, 1] * p[3]
                rb.ks_U[2] = J[1, 2] * p[1] + J[2, 2] * p[2] + J[3, 2] * p[3]
                rb.ks_U[3] = J[1, 3] * p[1] + J[2, 3] * p[2] + J[3, 3] * p[3]
                rb.ks_U[4] = J[1, 4] * p[1] + J[2, 4] * p[2] + J[3, 4] * p[3]

                psi = WeberElectrodynamics._ks_constraint(rb.ks_u, rb.ks_U)
                @test abs(psi) < 1e-10

                rb.ks_U .+= randn(4) * 1e-4
                c_err = WeberElectrodynamics._ks_project_constraint!(rb.ks_U, rb.ks_u, rb.ks_n)
                @test c_err < 1e-10
            end
        end
    end

    @testset "Switching and hysteresis" begin
        rb = WeberElectrodynamics.RegularizationBuffers(3, 2, 0.2, 0.3)

        q_activate = [-0.09, 0.0, 0.09, 0.0, 0.8, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_activate, true)
        @test active
        @test mode == WeberElectrodynamics.REG_MODE_PAIR

        q_between = [-0.13, 0.0, 0.13, 0.0, 0.8, 0.0]
        active, _, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_between, true)
        @test active

        q_off = [-0.17, 0.0, 0.17, 0.0, 0.8, 0.0]
        active, mode, _ = WeberElectrodynamics._detect_regularization_component!(rb, q_off, true)
        @test !active
        @test mode == WeberElectrodynamics.REG_MODE_NONE
    end

    @testset "Pair mode correctness 1D/2D/3D" begin
        @testset "1D" begin
            coarse_reg = make_orbit_problem(1;
                dt = 0.003,
                t_end = 0.8,
                v_scale = 0.0,
                regularization_enabled = true,
                r_on = 0.25,
                r_off = 0.32,
            )
            coarse_plain = make_orbit_problem(1;
                dt = 0.003,
                t_end = 0.8,
                v_scale = 0.0,
                regularization_enabled = false,
                r_on = 0.25,
                r_off = 0.32,
            )
            fine_ref = make_orbit_problem(1;
                dt = 0.001,
                t_end = 0.8,
                v_scale = 0.0,
                regularization_enabled = true,
                r_on = 0.25,
                r_off = 0.32,
            )

            sol_reg = solve(coarse_reg)
            sol_plain = solve(coarse_plain)
            sol_ref = solve(fine_ref)

            @test sol_reg.retcode == :Success
            @test sol_plain.retcode == :Success
            @test all(isfinite, sol_reg.q[end])
            @test all(isfinite, sol_reg.p[end])
            @test sol_reg.regularization.pair_steps > 0

            @test state_error(sol_reg, sol_ref) <= state_error(sol_plain, sol_ref) * 1.05
        end

        for dims in (2, 3)
            coarse_reg = make_orbit_problem(dims;
                dt = 0.003,
                t_end = 3.0,
                v_scale = 0.25,
                regularization_enabled = true,
                r_on = 0.6,
                r_off = 0.9,
            )
            coarse_plain = make_orbit_problem(dims;
                dt = 0.003,
                t_end = 3.0,
                v_scale = 0.25,
                regularization_enabled = false,
                r_on = 0.6,
                r_off = 0.9,
            )
            fine_ref = make_orbit_problem(dims;
                dt = 0.001,
                t_end = 3.0,
                v_scale = 0.25,
                regularization_enabled = true,
                r_on = 0.6,
                r_off = 0.9,
            )

            sol_reg = solve(coarse_reg)
            sol_plain = solve(coarse_plain)
            sol_ref = solve(fine_ref)

            @test sol_reg.retcode == :Success
            @test sol_plain.retcode == :Success
            @test all(isfinite, sol_reg.q[end])
            @test all(isfinite, sol_reg.p[end])
            @test sol_reg.regularization.pair_steps > 0

            @test state_error(sol_reg, sol_ref) <= state_error(sol_plain, sol_ref) * 1.05
        end
    end

    @testset "Chain mode correctness" begin
        sys = WeberSystem(3, 2)
        q0 = [-0.12, 0.0, 0.0, 0.0, 0.12, 0.0]
        p0 = [0.0, 0.02, 0.0, 0.0, 0.0, -0.02]

        prob = WeberProblem(
            sys,
            (0.0, 0.2),
            q0,
            p0;
            masses = [1.0, 1.0, 1.0],
            charges = [0.1, -0.2, 0.1],
            c = 4.0,
            dt = 0.002,
            regularization_enabled = true,
            regularization_r_on = 0.28,
            regularization_r_off = 0.39,
            regularization_chain_enabled = true,
            regularization_max_substeps = 256,
        )

        sol = solve(prob)
        @test sol.retcode == :Success
        @test sol.regularization.chain_steps > 0
    end

    @testset "Backward mode parity" begin
        sys = WeberSystem(2, 2)
        q0 = [-1.0, 0.0, 1.0, 0.0]
        p0 = [0.0, -0.05, 0.0, 0.05]

        kwargs = (
            masses = [1.0, 0.5],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
        )

        prob_off = WeberProblem(
            sys,
            (0.0, 0.2),
            q0,
            p0;
            kwargs...,
            regularization_enabled = false,
        )
        prob_on = WeberProblem(
            sys,
            (0.0, 0.2),
            q0,
            p0;
            kwargs...,
            regularization_enabled = true,
            regularization_r_on = 0.05,
            regularization_r_off = 0.08,
        )

        sol_off = solve(prob_off)
        sol_on = solve(prob_on)

        @test sol_on.regularization.pair_steps == 0
        @test sol_on.regularization.chain_steps == 0
        @test sol_on.regularization.unregularized_steps == length(sol_on) - 1
        @test sol_off.t == sol_on.t
        @test all(sol_off.q .== sol_on.q)
        @test all(sol_off.p .== sol_on.p)
    end

    @testset "Diagnostics validity" begin
        prob = make_orbit_problem(2;
            dt = 0.003,
            t_end = 3.0,
            v_scale = 0.25,
            regularization_enabled = true,
            r_on = 0.6,
            r_off = 0.9,
        )
        sol = solve(prob)
        d = sol.regularization

        @test d.pair_steps + d.chain_steps + d.unregularized_steps == length(sol) - 1
        @test d.total_substeps >= d.pair_steps + d.chain_steps + d.unregularized_steps
        @test d.max_substeps_used >= 1
        @test isfinite(d.min_encounter_distance)
        @test d.max_constraint_violation >= 0
    end

    @testset "Allocation checks" begin
        sys = WeberSystem(2, 2)

        q0 = [-1.0, 0.0, 1.0, 0.0]
        p0 = [0.0, -0.05, 0.0, 0.05]
        prob_unreg = WeberProblem(
            sys,
            (0.0, 0.02),
            q0,
            p0;
            masses = [1.0, 0.5],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
            regularization_enabled = false,
        )
        int_unreg = init(prob_unreg)
        step!(int_unreg)
        alloc_unreg = @allocated step!(int_unreg)
        @test alloc_unreg <= 512

        q1 = [-0.08, 0.01, 0.08, -0.01]
        p1 = [0.0, -0.02, 0.0, 0.02]
        prob_reg = WeberProblem(
            sys,
            (0.0, 0.02),
            q1,
            p1;
            masses = [1.0, 1.0],
            charges = [0.1, -0.1],
            c = 4.0,
            dt = 0.001,
            regularization_enabled = true,
            regularization_r_on = 0.2,
            regularization_r_off = 0.26,
            regularization_max_substeps = 256,
        )
        int_reg = init(prob_reg)
        step!(int_reg)
        @test int_reg.diagnostics.mode_history[1] == WeberElectrodynamics.REG_MODE_PAIR
        alloc_reg = @allocated step!(int_reg)
        @test alloc_reg <= 1024
    end
end
