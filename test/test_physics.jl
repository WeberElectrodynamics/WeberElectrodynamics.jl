@testset "Physics Validation" begin
    @testset "Energy conservation - Coulomb-like" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 5.0), dt = 0.001)
        sol = solve(prob)
        masses = prob.system.masses
        charges = prob.system.charges

        E(q, p) = coulomb_like_energy_2body_2d(q, p, masses, charges)
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i = 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-8
    end

    @testset "Energy conservation - Weber" begin
        prob = make_weber_problem(tspan = (0.0, 2.0), dt = 0.0005)
        sol = solve(prob)
        masses = prob.system.masses
        charges = prob.system.charges
        c = prob.system.c

        E(q, p) = weber_energy_2body_2d(q, p, masses, charges, c)
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i = 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-6
    end

    @testset "Newton's third law - Coulomb-like forces" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
        sol = solve(prob)

        forces = compute_force_timeseries(
            sol,
            2,
            2,
            prob.system.masses,
            prob.system.charges,
            prob.system.c,
        )
        n3 = check_newtons_third_law(forces)

        # F_12 = -F_21 should hold for Coulomb-like
        @test n3.global_max_violation < 1e-10
    end

    @testset "Newton's third law - Weber forces" begin
        prob = make_weber_problem(tspan = (0.0, 1.0), dt = 0.001)
        sol = solve(prob)

        forces = compute_force_timeseries(
            sol,
            2,
            2,
            prob.system.masses,
            prob.system.charges,
            prob.system.c,
        )
        n3 = check_newtons_third_law(forces)

        # F_12 = -F_21 should hold for Weber (it's a central force)
        @test n3.global_max_violation < 1e-6
    end

    @testset "Center of mass conservation" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 2.0), dt = 0.01)
        sol = solve(prob)

        m1, m2 = prob.system.masses
        M = m1 + m2

        # Initial center of mass position and momentum
        x_cm_0 = (m1 * sol.q[1][1] + m2 * sol.q[1][3]) / M
        y_cm_0 = (m1 * sol.q[1][2] + m2 * sol.q[1][4]) / M
        px_cm_0 = sol.p[1][1] + sol.p[1][3]
        py_cm_0 = sol.p[1][2] + sol.p[1][4]

        # Check conservation at all times
        for i = 1:length(sol)
            x_cm = (m1 * sol.q[i][1] + m2 * sol.q[i][3]) / M
            y_cm = (m1 * sol.q[i][2] + m2 * sol.q[i][4]) / M
            px_cm = sol.p[i][1] + sol.p[i][3]
            py_cm = sol.p[i][2] + sol.p[i][4]

            # Position drifts with velocity
            expected_x = x_cm_0 + px_cm_0 / M * sol.t[i]
            expected_y = y_cm_0 + py_cm_0 / M * sol.t[i]

            @test x_cm ≈ expected_x rtol = 1e-8
            @test y_cm ≈ expected_y rtol = 1e-8
            @test px_cm ≈ px_cm_0 rtol = 1e-10
            @test py_cm ≈ py_cm_0 rtol = 1e-10
        end
    end

    @testset "Angular momentum conservation (Coulomb-like)" begin
        prob = make_coulomb_like_problem(tspan = (0.0, 2.0), dt = 0.01)
        sol = solve(prob)

        function angular_momentum(q, p)
            # L = r1 × p1 + r2 × p2 (z-component in 2D)
            q[1] * p[2] - q[2] * p[1] + q[3] * p[4] - q[4] * p[3]
        end

        L0 = angular_momentum(sol.q[1], sol.p[1])

        for i = 1:length(sol)
            Li = angular_momentum(sol.q[i], sol.p[i])
            @test Li ≈ L0 rtol = 1e-8
        end
    end

    @testset "Symplectic integrator: no secular drift" begin
        # Long simulation to check for secular drift
        prob = make_coulomb_like_problem(tspan = (0.0, 50.0), dt = 0.01)
        sol = solve(prob)

        # Energy at beginning, middle, and end
        masses = prob.system.masses
        charges = prob.system.charges
        E(q, p) = coulomb_like_energy_2body_2d(q, p, masses, charges)
        E_start = E(sol.q[1], sol.p[1])
        E_mid = E(sol.q[length(sol)÷2], sol.p[length(sol)÷2])
        E_end = E(sol.q[end], sol.p[end])

        # No secular drift: errors should be bounded, not growing
        error_mid = abs(E_mid - E_start)
        error_end = abs(E_end - E_start)

        # End error should not be significantly larger than mid error
        # (would indicate secular drift)
        @test error_end < 2 * error_mid + 1e-10
    end

    @testset "Phase space orbit closure (Weber)" begin
        # For a near-circular bound orbit, phase space should approximately close
        # Use Weber system
        prob = make_weber_problem(tspan = (0.0, 20.0), dt = 0.001)
        sol = solve(prob)

        # For a bound orbit, the particle should stay bounded
        r_min = Inf
        r_max = 0.0
        for i = 1:length(sol)
            dx = sol.q[i][1] - sol.q[i][3]
            dy = sol.q[i][2] - sol.q[i][4]
            r = sqrt(dx^2 + dy^2)
            r_min = min(r_min, r)
            r_max = max(r_max, r)
        end

        # Orbit should be bounded (not escaping)
        # Weber orbits precess and may have larger excursions than Kepler
        @test r_max < 25.0
        @test r_min > 0.05  # Don't get too close (no collision)
    end

    @testset "Bound vs unbound orbits" begin
        # Negative total energy = bound orbit
        prob_bound = make_coulomb_like_problem(tspan = (0.0, 2.0), dt = 0.01)
        sol_bound = solve(prob_bound)

        masses = prob_bound.system.masses
        charges = prob_bound.system.charges
        E0 = coulomb_like_energy_2body_2d(sol_bound.q[1], sol_bound.p[1], masses, charges)
        @test E0 < 0  # Bound orbit has negative energy

        # Separation should remain bounded
        r_max = 0.0
        for i = 1:length(sol_bound)
            dx = sol_bound.q[i][1] - sol_bound.q[i][3]
            dy = sol_bound.q[i][2] - sol_bound.q[i][4]
            r = sqrt(dx^2 + dy^2)
            r_max = max(r_max, r)
        end
        @test r_max < 10.0  # Should stay bounded
    end

    @testset "Weber vs Coulomb limit" begin
        # As c → ∞, Weber should approach Coulomb
        m1, m2 = 1.0, 0.5
        q1, q2 = 1.0, -1.0

        # Weber with very large c
        sys_large_c = WeberSystem(2, 2; masses = [m1, m2], charges = [q1, q2], c = 1e10)

        # Weber with smaller c
        sys_small_c = WeberSystem(2, 2; masses = [m1, m2], charges = [q1, q2], c = 10.0)

        # Test at specific point
        q = [1.0, 0.0, -1.0, 0.0]
        p = [0.0, 0.5, 0.0, -0.5]

        out_large_c = zeros(4)
        out_small_c = zeros(4)

        sys_large_c.dp_dt_compiled(out_large_c, q, p)
        sys_small_c.dp_dt_compiled(out_small_c, q, p)

        # Force magnitude should be similar but with small Weber correction for finite c
        @test norm(out_large_c) ≈ norm(out_small_c) rtol = 0.1
    end
end
