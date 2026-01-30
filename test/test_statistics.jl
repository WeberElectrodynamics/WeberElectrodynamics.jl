@testset "Statistics" begin
    # Setup: run simulations for testing
    prob_coulomb = make_coulomb_like_problem(tspan = (0.0, 1.0), dt = 0.01)
    sol_coulomb = solve(prob_coulomb)

    prob_weber = make_weber_problem(tspan = (0.0, 2.0), dt = 0.001)
    sol_weber = solve(prob_weber)

    @testset "TrajectoryData" begin
        @testset "compute_trajectory_data basic" begin
            traj = compute_trajectory_data(sol_coulomb, 2, 2)

            @test traj isa TrajectoryData
            @test traj.n_particles == 2
            @test traj.dims == 2
            @test length(traj.trajectories) == 2
            @test size(traj.trajectories[1]) == (length(sol_coulomb.t), 2)
            @test length(traj.initial_positions) == 2
            @test length(traj.final_positions) == 2
            @test length(traj.initial_positions[1]) == 2
        end

        @testset "compute_trajectory_data with stride" begin
            traj = compute_trajectory_data(sol_coulomb, 2, 2; stride = 10)

            expected_points = length(1:10:length(sol_coulomb.t))
            @test size(traj.trajectories[1], 1) == expected_points
        end

        @testset "Initial and final positions" begin
            traj = compute_trajectory_data(sol_coulomb, 2, 2)

            # Initial positions should match solution
            for p = 1:2
                for d = 1:2
                    idx = (p - 1) * 2 + d
                    @test traj.initial_positions[p][d] == sol_coulomb.q[1][idx]
                    @test traj.final_positions[p][d] == sol_coulomb.q[end][idx]
                end
            end
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_trajectory_data(
                sol_coulomb,
                2,
                2;
                stride = 0,
            )
            @test_throws ArgumentError compute_trajectory_data(
                sol_coulomb,
                2,
                2;
                stride = -1,
            )
            @test_throws ArgumentError compute_trajectory_data(sol_coulomb, 3, 2)  # Wrong n_particles
            @test_throws ArgumentError compute_trajectory_data(sol_coulomb, 2, 3)  # Wrong dims
        end

        @testset "1D system" begin
            # Create a 1D 2-body system for this test
            sys1d = WeberSystem(2, 1; masses = [1.0, 1.0], charges = [1.0, -1.0], c = 1e10)
            prob1d = WeberProblem(sys1d, (0.0, 0.5), [1.0, -1.0], [0.1, -0.1]; dt = 0.01)
            sol1d = solve(prob1d)
            traj = compute_trajectory_data(sol1d, 2, 1)
            @test traj.n_particles == 2
            @test traj.dims == 1
            @test size(traj.trajectories[1], 2) == 1
        end
    end

    @testset "EnergyData" begin
        # Energy functions for 2-body Weber system
        masses = prob_weber.system.masses
        charges = prob_weber.system.charges
        c = prob_weber.system.c

        total_energy(q, p, params, t) = weber_energy_2body_2d(q, p, masses, charges, c)
        kinetic_energy(q, p, params, t) =
            sum(p[1:2] .^ 2) / (2 * masses[1]) + sum(p[3:4] .^ 2) / (2 * masses[2])
        function potential_energy(q, p, params, t)
            weber_energy_2body_2d(q, p, masses, charges, c) -
            kinetic_energy(q, p, params, t)
        end

        @testset "compute_energy_timeseries basic" begin
            energy = compute_energy_timeseries(
                sol_weber,
                total_energy,
                nothing,
                nothing,
                nothing,
            )

            @test energy isa EnergyData
            @test length(energy.t) == length(sol_weber.t)
            @test length(energy.total_energy) == length(sol_weber.t)
            @test isnothing(energy.kinetic_energy)
            @test isnothing(energy.potential_energy)
        end

        @testset "compute_energy_timeseries with components" begin
            energy = compute_energy_timeseries(
                sol_weber,
                total_energy,
                kinetic_energy,
                potential_energy,
                nothing,
            )

            @test !isnothing(energy.kinetic_energy)
            @test !isnothing(energy.potential_energy)
            @test length(energy.kinetic_energy) == length(sol_weber.t)
            @test length(energy.potential_energy) == length(sol_weber.t)

            # KE + PE = Total (approximately)
            for i = 1:length(energy.t)
                @test energy.kinetic_energy[i] + energy.potential_energy[i] ≈
                      energy.total_energy[i] rtol = 1e-10
            end
        end

        @testset "Energy statistics" begin
            energy = compute_energy_timeseries(
                sol_weber,
                total_energy,
                nothing,
                nothing,
                nothing,
            )

            @test energy.max_local_error >= 0
            # For symplectic integrator, relative range should be small
            if !isnan(energy.relative_energy_range)
                @test energy.relative_energy_range < 1e-6
            end
        end

        @testset "compute_energy_timeseries with stride" begin
            energy = compute_energy_timeseries(
                sol_weber,
                total_energy,
                nothing,
                nothing,
                nothing;
                stride = 10,
            )

            expected_points = length(1:10:length(sol_weber.t))
            @test length(energy.t) == expected_points
            @test length(energy.total_energy) == expected_points
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_energy_timeseries(
                sol_weber,
                total_energy,
                nothing,
                nothing,
                nothing;
                stride = 0,
            )
            @test_throws ArgumentError compute_energy_timeseries(
                sol_weber,
                total_energy,
                nothing,
                nothing,
                nothing;
                stride = -1,
            )
        end
    end

    @testset "ForceData" begin
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        c = 1e6  # Large c for Coulomb limit

        @testset "compute_force_timeseries basic" begin
            forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c)

            @test forces isa ForceData
            @test forces.n_particles == 2
            @test forces.dims == 2
            @test length(forces.t) == length(sol_coulomb.t) - 1  # Forces need acceleration

            # Should have forces for both (1,2) and (2,1)
            @test haskey(forces.forces, (1, 2))
            @test haskey(forces.forces, (2, 1))
        end

        @testset "Force symmetry (Newton's third law)" begin
            forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c)

            F_12 = forces.forces[(1, 2)]
            F_21 = forces.forces[(2, 1)]

            # F_12 = -F_21
            for t = 1:length(forces.t)
                @test F_12[t] ≈ -F_21[t] atol = 1e-10
            end
        end

        @testset "compute_force_timeseries with stride" begin
            forces =
                compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c; stride = 5)

            # Fewer time points due to stride
            @test length(forces.t) < length(sol_coulomb.t) - 1
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_force_timeseries(
                sol_coulomb,
                2,
                2,
                masses,
                charges,
                c;
                stride = 0,
            )
            @test_throws ArgumentError compute_force_timeseries(
                sol_coulomb,
                3,
                2,
                masses,
                charges,
                c,
            )  # Wrong n_particles
            @test_throws ArgumentError compute_force_timeseries(
                sol_coulomb,
                2,
                2,
                [1.0],
                charges,
                c,
            )  # Wrong masses length
            @test_throws ArgumentError compute_force_timeseries(
                sol_coulomb,
                2,
                2,
                masses,
                [1.0],
                c,
            )  # Wrong charges length
            @test_throws ArgumentError compute_force_timeseries(
                sol_coulomb,
                2,
                2,
                masses,
                charges,
                -1.0,
            )  # Negative c
        end
    end

    @testset "NewtonsThirdLawData" begin
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        c = 1e6

        forces = compute_force_timeseries(sol_coulomb, 2, 2, masses, charges, c)

        @testset "check_newtons_third_law" begin
            n3 = check_newtons_third_law(forces)

            @test n3 isa NewtonsThirdLawData
            @test n3.n_pairs == 1  # One pair for 2 particles
            @test length(n3.t) == length(forces.t)
            @test haskey(n3.pair_violations, (1, 2))
            @test haskey(n3.max_violations, (1, 2))
            @test haskey(n3.mean_violations, (1, 2))
            @test haskey(n3.rms_violations, (1, 2))

            # Violations should be small for Coulomb
            @test n3.global_max_violation < 1e-8
        end
    end

    @testset "PhaseSpaceData" begin
        masses = [1.0, 0.5]

        @testset "compute_phase_space_data basic" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses)

            @test ps isa PhaseSpaceData
            @test length(ps.t) == length(sol_coulomb.t)
            @test length(ps.separation_distance) == length(sol_coulomb.t)
            @test length(ps.radial_velocity) == length(sol_coulomb.t)
            @test !isnothing(ps.theta)
            @test !isnothing(ps.angular_momentum)
        end

        @testset "compute_phase_space_data with stride" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses; stride = 10)

            expected_points = length(1:10:length(sol_coulomb.t))
            @test length(ps.t) == expected_points
            @test length(ps.separation_distance) == expected_points
        end

        @testset "Different particle pairs" begin
            ps_12 =
                compute_phase_space_data(sol_coulomb, 2, 2, masses; particle_pair = (1, 2))
            ps_21 =
                compute_phase_space_data(sol_coulomb, 2, 2, masses; particle_pair = (2, 1))

            # separation_distance should be the same (distance is symmetric)
            @test ps_12.separation_distance ≈ ps_21.separation_distance

            # radial_velocity = dot(r_vec, v_vec) / |r|
            # Swapping particles flips both r_vec and v_vec, so radial_velocity is the same
            @test ps_12.radial_velocity ≈ ps_21.radial_velocity
        end

        @testset "Disable optional computations" begin
            ps = compute_phase_space_data(
                sol_coulomb,
                2,
                2,
                masses;
                compute_angle = false,
                compute_angular_momentum = false,
            )

            @test isnothing(ps.theta)
            @test isnothing(ps.angular_momentum)
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_phase_space_data(
                sol_coulomb,
                2,
                2,
                masses;
                stride = 0,
            )
            @test_throws ArgumentError compute_phase_space_data(
                sol_coulomb,
                2,
                2,
                masses;
                particle_pair = (0, 1),
            )
            @test_throws ArgumentError compute_phase_space_data(
                sol_coulomb,
                2,
                2,
                masses;
                particle_pair = (1, 3),
            )
        end

        @testset "Positive separation distance" begin
            ps = compute_phase_space_data(sol_coulomb, 2, 2, masses)
            @test all(ps.separation_distance .> 0)  # Particles should never collide
        end
    end
end
