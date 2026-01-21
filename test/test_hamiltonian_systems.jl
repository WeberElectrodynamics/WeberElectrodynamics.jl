@testset "Hamiltonian Systems" begin
    @testset "generate_phase_space_symbols" begin
        # 1 particle, 1D
        q, p = generate_phase_space_symbols(1, 1)
        @test q == [:x1]
        @test p == [:px1]

        # 2 particles, 2D
        q, p = generate_phase_space_symbols(2, 2)
        @test q == [:x1, :y1, :x2, :y2]
        @test p == [:px1, :py1, :px2, :py2]

        # 2 particles, 3D
        q, p = generate_phase_space_symbols(2, 3)
        @test q == [:x1, :y1, :z1, :x2, :y2, :z2]
        @test p == [:px1, :py1, :pz1, :px2, :py2, :pz2]

        # 3 particles, 2D
        q, p = generate_phase_space_symbols(3, 2)
        @test length(q) == 6
        @test length(p) == 6
        @test q == [:x1, :y1, :x2, :y2, :x3, :y3]
    end

    @testset "compile_hamiltonian function API" begin
        H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])

        @test H.degrees_of_freedom == 1
        @test H.parameter_names == [:m, :k]

        # Test compiled functions work
        out_q = zeros(1)
        out_p = zeros(1)
        q = [1.0]
        p = [0.5]
        params = [1.0, 2.0]  # m=1, k=2

        H.dq_dt_compiled(out_q, q, p, params)
        H.dp_dt_compiled(out_p, q, p, params)

        # dq/dt = dH/dp = p/m = 0.5/1.0 = 0.5
        @test out_q[1] ≈ 0.5
        # dp/dt = -dH/dq = -k*q = -2.0*1.0 = -2.0
        @test out_p[1] ≈ -2.0
    end

    @testset "@hamiltonian macro API" begin
        H = @hamiltonian 1 1 [:m, :k] (q, p, params) -> begin
            m, k = params
            sum(p .^ 2) / (2m) + k * sum(q .^ 2) / 2
        end

        @test H isa WeberHamiltonian
        @test H.degrees_of_freedom == 1
        @test H.parameter_names == [:m, :k]

        # Test it produces same results as function API
        out_q = zeros(1)
        out_p = zeros(1)
        H.dq_dt_compiled(out_q, [1.0], [0.5], [1.0, 2.0])
        H.dp_dt_compiled(out_p, [1.0], [0.5], [1.0, 2.0])
        @test out_q[1] ≈ 0.5
        @test out_p[1] ≈ -2.0
    end

    @testset "Multi-particle Hamiltonian" begin
        H = compile_hamiltonian(weber_H, 2, 2; parameter_names=[:m1, :m2, :k, :c])

        @test H.degrees_of_freedom == 4  # 2 particles × 2 dims
        @test H.parameter_names == [:m1, :m2, :k, :c]
        @test length(H.dq_dt_symbolic) == 4
        @test length(H.dp_dt_symbolic) == 4
    end

    @testset "Empty parameter_names" begin
        H_free(q, p, params) = sum(p .^ 2) / 2
        H = compile_hamiltonian(H_free, 1, 1; parameter_names=Symbol[])

        @test H.parameter_names == Symbol[]
        @test H.degrees_of_freedom == 1

        # Should still work
        out_q = zeros(1)
        H.dq_dt_compiled(out_q, [0.0], [1.0], Float64[])
        @test out_q[1] ≈ 1.0
    end

    @testset "Display methods" begin
        H = compile_hamiltonian(harmonic_oscillator_H, 1, 1; parameter_names=[:m, :k])

        # show(io, H)
        io = IOBuffer()
        show(io, H)
        str = String(take!(io))
        @test occursin("WeberHamiltonian", str)
        @test occursin("1 DOF", str)

        # show(io, MIME"text/plain", H)
        io = IOBuffer()
        show(io, MIME"text/plain"(), H)
        str = String(take!(io))
        @test occursin("WeberHamiltonian", str)
        @test occursin("DOF: 1", str)
        @test occursin("Parameters:", str)
    end

    @testset "Coulomb Hamiltonian derivatives" begin
        H = compile_hamiltonian(coulomb_H, 2, 2; parameter_names=[:m1, :m2, :k])

        @test H.degrees_of_freedom == 4

        # Test at specific point
        q = [1.0, 0.0, -1.0, 0.0]  # particles at x=1 and x=-1
        p = [0.0, 0.1, 0.0, -0.1]  # small momenta in y
        params = [1.0, 1.0, 1.0]  # equal masses, k=1

        out_q = zeros(4)
        out_p = zeros(4)
        H.dq_dt_compiled(out_q, q, p, params)
        H.dp_dt_compiled(out_p, q, p, params)

        # dq/dt = dH/dp = p/m
        @test out_q[1] ≈ 0.0  # px1/m1
        @test out_q[2] ≈ 0.1  # py1/m1
        @test out_q[3] ≈ 0.0  # px2/m2
        @test out_q[4] ≈ -0.1  # py2/m2

        # dp/dt = -dH/dq should be attractive force along x-axis
        @test out_p[1] < 0  # particle 1 pulled toward particle 2 (negative x)
        @test out_p[3] > 0  # particle 2 pulled toward particle 1 (positive x)
        @test abs(out_p[2]) < 1e-10  # no force in y at this configuration
        @test abs(out_p[4]) < 1e-10
    end
end
