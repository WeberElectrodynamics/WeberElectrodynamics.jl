# Internal function for testing - access via module
const solve_nonlinear! = WeberElectrodynamics.solve_nonlinear!

@testset "Nonlinear Solvers" begin
    @testset "RelaxedFixedPointSolver convergence - linear" begin
        # Simple linear problem: f(x) = x - 1 (solution x = 1)
        function f_linear!(result, u)
            result[1] = u[1] - 1.0
            return nothing
        end

        u = [0.0]
        fu_buffer = zeros(1)
        u_old_buffer = zeros(1)

        solver = RelaxedFixedPointSolver(relaxation=0.5)
        success, iters, residual = solve_nonlinear!(u, f_linear!, solver, 1e-12, 100;
            fu_buffer=fu_buffer,
            u_old_buffer=u_old_buffer)

        @test success
        @test u[1] ≈ 1.0 atol = 1e-10
        @test residual < 1e-12
        @test iters < 100
    end

    @testset "RelaxedFixedPointSolver convergence - multidimensional" begin
        # 2D problem: f(x) = x - [1, 2]
        function f_2d!(result, u)
            result[1] = u[1] - 1.0
            result[2] = u[2] - 2.0
            return nothing
        end

        u = zeros(2)
        fu_buffer = zeros(2)
        u_old_buffer = zeros(2)

        solver = RelaxedFixedPointSolver(relaxation=0.5)
        success, _, _ = solve_nonlinear!(u, f_2d!, solver, 1e-12, 100;
            fu_buffer=fu_buffer,
            u_old_buffer=u_old_buffer)

        @test success
        @test u ≈ [1.0, 2.0] atol = 1e-10
    end

    @testset "RelaxedFixedPointSolver max iterations" begin
        # Problem that converges slowly - need many iterations
        function f_slow!(result, u)
            # This creates a very slow convergence
            result[1] = 0.99 * (u[1] - 1.0)
            return nothing
        end

        u = [0.0]
        fu_buffer = zeros(1)
        u_old_buffer = zeros(1)

        solver = RelaxedFixedPointSolver(relaxation=0.1)  # Very slow
        success, iters, residual = solve_nonlinear!(u, f_slow!, solver, 1e-12, 5;
            fu_buffer=fu_buffer,
            u_old_buffer=u_old_buffer)

        @test !success  # Should not converge in 5 iterations
        @test iters == 5
        @test residual > 1e-12
    end

    @testset "Relaxation parameter effect" begin
        function f_test!(result, u)
            result[1] = u[1] - 1.0
            return nothing
        end

        # Higher relaxation = faster convergence for well-conditioned problems
        iterations_list = Int[]
        for relaxation in [0.1, 0.25, 0.5, 1.0]
            u = [0.0]
            fu_buffer = zeros(1)
            u_old_buffer = zeros(1)

            solver = RelaxedFixedPointSolver(relaxation=relaxation)
            success, iters, _ = solve_nonlinear!(u, f_test!, solver, 1e-12, 1000;
                fu_buffer=fu_buffer,
                u_old_buffer=u_old_buffer)
            @test success
            push!(iterations_list, iters)
        end

        # Faster relaxation should converge in fewer iterations
        @test iterations_list[end] <= iterations_list[1]
    end

    @testset "Zero initial residual" begin
        # Already at solution
        function f_zero!(result, u)
            result[1] = u[1] - 1.0
            return nothing
        end

        u = [1.0]  # Start at solution
        fu_buffer = zeros(1)
        u_old_buffer = zeros(1)

        solver = RelaxedFixedPointSolver()
        success, iters, residual = solve_nonlinear!(u, f_zero!, solver, 1e-12, 100;
            fu_buffer=fu_buffer,
            u_old_buffer=u_old_buffer)

        @test success
        @test iters == 1  # Should converge immediately
        @test residual < 1e-12
    end

    @testset "High-dimensional problem" begin
        n = 10
        target = collect(1.0:Float64(n))

        function f_nd!(result, u)
            @. result = u - target
            return nothing
        end

        u = zeros(n)
        fu_buffer = zeros(n)
        u_old_buffer = zeros(n)

        solver = RelaxedFixedPointSolver(relaxation=0.5)
        success, _, _ = solve_nonlinear!(u, f_nd!, solver, 1e-12, 200;
            fu_buffer=fu_buffer,
            u_old_buffer=u_old_buffer)

        @test success
        @test u ≈ target atol = 1e-10
    end
end
