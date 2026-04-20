@testset "Callbacks" begin
    @testset "CollisionBounce basic" begin
        cb = CollisionBounce(0.01)
        @test cb isa HamiltonianCallback
        @test cb.radius == 0.01
        @test_throws AssertionError CollisionBounce(-0.1)
    end

    @testset "default callbacks are empty" begin
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.05), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 100.0,
            dt = 0.01,
        )
        integrator = init(prob, SymmetricProjectionIntegrator())
        @test integrator.callbacks === ()
    end

    @testset "user callback stored on integrator" begin
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.05), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 100.0,
            dt = 0.01,
        )
        cb = CollisionBounce(0.02)
        integrator = init(prob, SymmetricProjectionIntegrator(); callbacks = cb)
        @test length(integrator.callbacks) == 1
        @test integrator.callbacks[1] === cb
    end

    # The CollisionBounce callback synthesised from `collision_bounce_radius`
    # on a `RegularizedIntegrator` must produce trajectories bit-exactly identical
    # to an explicit `CollisionBounce` callback. Sub-critical like-charge head-on
    # oscillation: bounce is active every half-period, so every step exercises
    # the reflection branch.
    @testset "equivalence: collision_bounce_radius synthesises CollisionBounce" begin
        m1 = m2 = 1.0
        q1 = q2 = 1.0
        c = 4.0
        r0 = 0.05
        mu = m1 * m2 / (m1 + m2)
        rho = q1 * q2 / (mu * c^2)
        T = (2π / c) * rho

        sys = HamiltonianSystem(2, 2)
        q_init = [r0 / 2, 0.0, -r0 / 2, 0.0]
        p_init = [0.0, 0.0, 0.0, 0.0]

        prob = HamiltonianProblem(
            sys,
            (0.0, 3 * T),
            q_init,
            p_init;
            masses = [m1, m2],
            charges = [q1, q2],
            c = c,
            dt = T / 500,
        )
        alg_synth = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            r_on = 0.001,
            r_off = 0.002,
            collision_bounce_radius = 0.01,
        )
        alg_plain = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            r_on = 0.001,
            r_off = 0.002,
        )
        sol_synth = solve(prob, alg_synth)
        sol_explicit = solve(prob, alg_plain; callbacks = CollisionBounce(0.01))

        @test sol_synth.retcode === sol_explicit.retcode
        @test length(sol_synth.t) == length(sol_explicit.t)
        max_q = maximum(
            maximum(abs, sol_synth.q[i] .- sol_explicit.q[i]) for
            i in eachindex(sol_synth.q)
        )
        max_p = maximum(
            maximum(abs, sol_synth.p[i] .- sol_explicit.p[i]) for
            i in eachindex(sol_synth.p)
        )
        @test max_q == 0.0
        @test max_p == 0.0
    end

    @testset "RegularizedIntegrator bridge synthesises CollisionBounce" begin
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.05), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 100.0,
            dt = 0.01,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            collision_bounce_radius = 0.05,
        )
        integrator = init(prob, alg)
        @test length(integrator.callbacks) == 1
        @test integrator.callbacks[1] isa CollisionBounce
        @test integrator.callbacks[1].radius == 0.05
    end

    @testset "explicit CollisionBounce overrides bridge synthesis" begin
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.05), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 100.0,
            dt = 0.01,
        )
        alg = RegularizedIntegrator(
            SymmetricProjectionIntegrator();
            backend = :adaptive_cartesian,
            warn_on_fallback = false,
            collision_bounce_radius = 0.05,
        )
        explicit = CollisionBounce(0.07)
        integrator = init(prob, alg; callbacks = explicit)
        @test length(integrator.callbacks) == 1
        @test integrator.callbacks[1] === explicit
    end

    @testset "abstract callback no-op defaults" begin
        struct _NoOpCallback <: HamiltonianCallback end
        cb = _NoOpCallback()
        sys = HamiltonianSystem(2, 2)
        prob = HamiltonianProblem(
            sys, (0.0, 0.03), [1.0, 0.0, -1.0, 0.0], [0.0, 0.1, 0.0, -0.1];
            masses = [1.0, 1.0],
            charges = [1.0, -1.0],
            c = 100.0,
            dt = 0.01,
        )
        sol_a = solve(prob, SymmetricProjectionIntegrator())
        sol_b = solve(prob, SymmetricProjectionIntegrator(); callbacks = cb)
        @test sol_a.q[end] == sol_b.q[end]
        @test sol_a.p[end] == sol_b.p[end]
    end
end
