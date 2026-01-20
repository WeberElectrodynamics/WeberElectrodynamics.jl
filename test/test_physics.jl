@testset "Physics Validation" begin
    @testset "Harmonic oscillator period" begin
        # Analytical period: T = 2π√(m/k)
        m, k = 1.0, 4.0
        T_analytical = 2π * sqrt(m / k)

        H = build_hamiltonian(harmonic_oscillator_H, 1, 1; param_names=[:m, :k])
        prob = WeberProblem(H, (0.0, 2 * T_analytical), [1.0], [0.0]; params=[m, k], dt=0.001)
        sol = solve(prob)

        # Find zero crossings to measure period
        q_vals = getindex.(sol.q, 1)
        crossings = Int[]
        for i in 2:length(q_vals)
            if q_vals[i-1] * q_vals[i] < 0 && q_vals[i-1] > 0
                push!(crossings, i)
            end
        end

        if length(crossings) >= 2
            T_measured = sol.t[crossings[2]] - sol.t[crossings[1]]
            @test T_measured ≈ T_analytical rtol = 0.01
        end
    end

    @testset "Energy conservation - harmonic oscillator" begin
        prob = make_harmonic_problem(tspan=(0.0, 10.0), dt=0.001)
        sol = solve(prob)

        E(q, p) = harmonic_oscillator_H(q, p, [1.0, 1.0])
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i in 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-10
    end

    @testset "Energy conservation - Coulomb" begin
        prob = make_coulomb_problem(tspan=(0.0, 5.0), dt=0.001)
        sol = solve(prob)
        params = [1.0, 0.5, 1.0]

        E(q, p) = coulomb_H(q, p, params)
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i in 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-8
    end

    @testset "Energy conservation - Weber" begin
        prob = make_weber_problem(tspan=(0.0, 2.0), dt=0.0005)
        sol = solve(prob)
        params = [1.0, 0.1, -0.1, 4.0]

        E(q, p) = weber_H(q, p, params)
        E0 = E(sol.q[1], sol.p[1])

        max_error = 0.0
        for i in 1:length(sol)
            Ei = E(sol.q[i], sol.p[i])
            error = abs(Ei - E0) / abs(E0)
            max_error = max(max_error, error)
        end

        @test max_error < 1e-6
    end

    @testset "Newton's third law - Coulomb forces" begin
        prob = make_coulomb_problem(tspan=(0.0, 1.0), dt=0.01)
        sol = solve(prob)

        forces = compute_force_timeseries(sol, 2, 2, [1.0, 0.5], [1.0, -1.0], 1e10)
        n3 = check_newtons_third_law(forces)

        # F_12 = -F_21 should hold exactly for Coulomb
        @test n3.global_max_violation < 1e-10
    end

    @testset "Center of mass conservation" begin
        prob = make_coulomb_problem(tspan=(0.0, 2.0), dt=0.01)
        sol = solve(prob)

        m1, m2 = 1.0, 0.5
        M = m1 + m2

        # Initial center of mass position and momentum
        x_cm_0 = (m1 * sol.q[1][1] + m2 * sol.q[1][3]) / M
        y_cm_0 = (m1 * sol.q[1][2] + m2 * sol.q[1][4]) / M
        px_cm_0 = sol.p[1][1] + sol.p[1][3]
        py_cm_0 = sol.p[1][2] + sol.p[1][4]

        # Check conservation at all times
        for i in 1:length(sol)
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

    @testset "Angular momentum conservation (Coulomb)" begin
        prob = make_coulomb_problem(tspan=(0.0, 2.0), dt=0.01)
        sol = solve(prob)

        function angular_momentum(q, p)
            # L = r1 × p1 + r2 × p2 (z-component in 2D)
            q[1] * p[2] - q[2] * p[1] + q[3] * p[4] - q[4] * p[3]
        end

        L0 = angular_momentum(sol.q[1], sol.p[1])

        for i in 1:length(sol)
            Li = angular_momentum(sol.q[i], sol.p[i])
            @test Li ≈ L0 rtol = 1e-8
        end
    end

    @testset "Symplectic integrator: no secular drift" begin
        # Long simulation to check for secular drift
        prob = make_harmonic_problem(tspan=(0.0, 100.0), dt=0.01)
        sol = solve(prob)

        # Energy at beginning, middle, and end
        E(q, p) = harmonic_oscillator_H(q, p, [1.0, 1.0])
        E_start = E(sol.q[1], sol.p[1])
        E_mid = E(sol.q[length(sol)÷2], sol.p[length(sol)÷2])
        E_end = E(sol.q[end], sol.p[end])

        # No secular drift: errors should be bounded, not growing
        error_mid = abs(E_mid - E_start)
        error_end = abs(E_end - E_start)

        # End error should not be significantly larger than mid error
        # (would indicate secular drift)
        @test error_end < 2 * error_mid + 1e-12
    end

    @testset "Phase space orbit closure" begin
        # For a bound orbit, phase space should approximately close
        m, k = 1.0, 1.0
        T = 2π * sqrt(m / k)  # One period

        prob = make_harmonic_problem(tspan=(0.0, T), dt=0.001, m=m, k=k)
        sol = solve(prob)

        # Final state should be close to initial
        @test sol.q[end][1] ≈ sol.q[1][1] rtol = 0.01
        @test sol.p[end][1] ≈ sol.p[1][1] atol = 0.01
    end

    @testset "Bound vs unbound orbits" begin
        # Negative total energy = bound orbit
        prob_bound = make_coulomb_problem(tspan=(0.0, 2.0), dt=0.01)
        sol_bound = solve(prob_bound)

        params = [1.0, 0.5, 1.0]
        E0 = coulomb_H(sol_bound.q[1], sol_bound.p[1], params)
        @test E0 < 0  # Bound orbit has negative energy

        # Separation should remain bounded
        r_max = 0.0
        for i in 1:length(sol_bound)
            dx = sol_bound.q[i][1] - sol_bound.q[i][3]
            dy = sol_bound.q[i][2] - sol_bound.q[i][4]
            r = sqrt(dx^2 + dy^2)
            r_max = max(r_max, r)
        end
        @test r_max < 10.0  # Should stay bounded
    end
end
