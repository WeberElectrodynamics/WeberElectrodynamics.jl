using WeberElectrodynamics: compute_total_kinetic_energy, compute_pair_weber_components

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
            sys1d = HamiltonianSystem(2, 1)
            prob1d = HamiltonianProblem(
                sys1d,
                (0.0, 0.5),
                [1.0, -1.0],
                [0.1, -0.1];
                masses = [1.0, 1.0],
                charges = [1.0, -1.0],
                c = 1e10,
                dt = 0.01,
            )
            sol1d = solve(prob1d)
            traj = compute_trajectory_data(sol1d, 2, 1)
            @test traj.n_particles == 2
            @test traj.dims == 1
            @test size(traj.trajectories[1], 2) == 1
        end
    end

    @testset "EnergyData" begin
        @testset "compute_energy_timeseries basic" begin
            energy = compute_energy_timeseries(sol_weber)

            @test energy isa EnergyData
            @test length(energy.t) == length(sol_weber.t)
            @test length(energy.total_energy) == length(sol_weber.t)
            @test length(energy.kinetic_energy) == length(sol_weber.t)
            @test length(energy.total_potential_energy) == length(sol_weber.t)
            @test energy.n_particles == 2
            @test energy.n_pairs == 1
        end

        @testset "Pair energy decomposition" begin
            energy = compute_energy_timeseries(sol_weber)

            @test haskey(energy.pair_energies, (1, 2))
            pair_data = energy.pair_energies[(1, 2)]
            @test pair_data isa PairEnergyData
            @test length(pair_data.coulomb_term) == length(sol_weber.t)
            @test length(pair_data.velocity_term) == length(sol_weber.t)
            @test length(pair_data.radial_velocity) == length(sol_weber.t)

            # Verify pair potential = coulomb + velocity
            for i in eachindex(pair_data.total_pair_potential)
                @test pair_data.total_pair_potential[i] ≈
                      pair_data.coulomb_term[i] + pair_data.velocity_term[i] rtol = 1e-14
            end
        end

        @testset "Energy conservation" begin
            energy = compute_energy_timeseries(sol_weber)

            # KE + PE = Total
            for i in eachindex(energy.total_energy)
                @test energy.kinetic_energy[i] + energy.total_potential_energy[i] ≈
                      energy.total_energy[i] rtol = 1e-14
            end
        end

        @testset "Hamiltonian validation" begin
            energy = compute_energy_timeseries(sol_weber)

            # Our computed energy should match HamiltonianSystem's compiled Hamiltonian
            @test all(energy.hamiltonian_validation_error .< 1e-12)
        end

        @testset "Statistics" begin
            energy = compute_energy_timeseries(sol_weber)
            stats = energy.statistics

            @test stats isa EnergyStatistics
            @test stats.local_error_max >= stats.local_error_min
            @test stats.local_error_max >= 0
            @test stats.local_error_avg >= 0

            # For symplectic integrator, energy should be well-conserved
            @test stats.global_error_percent_max < 1e-2  # < 0.01% error
        end

        @testset "compute_energy_timeseries with stride" begin
            energy = compute_energy_timeseries(sol_weber; stride = 10)

            expected_points = length(1:10:length(sol_weber.t))
            @test length(energy.t) == expected_points
            @test length(energy.total_energy) == expected_points
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_energy_timeseries(sol_weber; stride = 0)
            @test_throws ArgumentError compute_energy_timeseries(sol_weber; stride = -1)
        end

        @testset "3-body system" begin
            # Create a 3-body system to test n_pairs = 3
            sys3 = HamiltonianSystem(3, 2)
            prob3 = HamiltonianProblem(
                sys3,
                (0.0, 0.1),
                [1.0, 0.0, -0.5, 0.866, -0.5, -0.866],  # Triangle
                zeros(6);
                masses = [1.0, 1.0, 1.0],
                charges = [0.1, 0.1, 0.1],
                c = 1e6,
                dt = 0.001,
            )
            sol3 = solve(prob3)
            energy3 = compute_energy_timeseries(sol3)

            @test energy3.n_particles == 3
            @test energy3.n_pairs == 3
            @test haskey(energy3.pair_energies, (1, 2))
            @test haskey(energy3.pair_energies, (1, 3))
            @test haskey(energy3.pair_energies, (2, 3))
        end
    end

    @testset "PairForceData" begin
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        c = 1e6  # Large c for Coulomb limit

        @testset "compute_pair_force_timeseries basic" begin
            forces =
                compute_pair_force_timeseries(sol_coulomb, (1, 2), 2, 2, masses, charges, c)

            @test forces isa PairForceData
            @test forces.pair == (1, 2)
            @test forces.dims == 2
            @test length(forces.t) == length(sol_coulomb.t) - 1  # Forces need acceleration

            # Check all components are present
            @test length(forces.force) == length(forces.t)
            @test length(forces.magnitude) == length(forces.t)
            @test length(forces.coulomb) == length(forces.t)
            @test length(forces.vector_term_vv) == length(forces.t)
            @test length(forces.vector_term_ra) == length(forces.t)
            @test length(forces.vector_term_rv2) == length(forces.t)
            @test length(forces.radial_term_rdot2) == length(forces.t)
            @test length(forces.radial_term_rddot) == length(forces.t)
        end

        @testset "ForceStatistics" begin
            forces =
                compute_pair_force_timeseries(sol_coulomb, (1, 2), 2, 2, masses, charges, c)
            stats = forces.stats

            @test stats isa ForceStatistics
            @test stats.min <= stats.max
            @test stats.mean >= stats.min
            @test stats.mean <= stats.max
            @test stats.range ≈ stats.max - stats.min
        end

        @testset "Weber force decomposition consistency" begin
            forces =
                compute_pair_force_timeseries(sol_weber, (1, 2), 2, 2, masses, charges, 1e6)

            # Vector form: F = Coulomb + v·v term + r·a term + rv² term
            for t = 1:length(forces.t)
                force_from_terms =
                    forces.coulomb[t] .+ forces.vector_term_vv[t] .+
                    forces.vector_term_ra[t] .+ forces.vector_term_rv2[t]
                @test forces.force[t] ≈ force_from_terms rtol = 1e-12
            end

            # Radial form: F = Coulomb + rdot² term + rddot term
            # Should equal vector form
            for t = 1:length(forces.t)
                force_radial =
                    forces.coulomb[t] .+ forces.radial_term_rdot2[t] .+
                    forces.radial_term_rddot[t]
                @test forces.force[t] ≈ force_radial rtol = 1e-10
            end
        end

        @testset "compute_pair_force_timeseries with stride" begin
            forces = compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                2,
                2,
                masses,
                charges,
                c;
                stride = 5,
            )

            # Fewer time points due to stride
            @test length(forces.t) < length(sol_coulomb.t) - 1
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                2,
                2,
                masses,
                charges,
                c;
                stride = 0,
            )
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                3,
                2,
                masses,
                charges,
                c,
            )  # Wrong n_particles
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                2,
                2,
                [1.0],
                charges,
                c,
            )  # Wrong masses length
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                2,
                2,
                masses,
                [1.0],
                c,
            )  # Wrong charges length
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                2,
                2,
                masses,
                charges,
                -1.0,
            )  # Negative c
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 3),
                2,
                2,
                masses,
                charges,
                c,
            )  # Invalid pair index
            @test_throws ArgumentError compute_pair_force_timeseries(
                sol_coulomb,
                (1, 1),
                2,
                2,
                masses,
                charges,
                c,
            )  # Same particle pair
        end
    end

    @testset "PhaseSpaceData (embedded in PairForceData)" begin
        masses = [1.0, 0.5]
        charges = [1.0, -1.0]
        c = 1e6

        @testset "Phase space data embedded in force data" begin
            forces =
                compute_pair_force_timeseries(sol_coulomb, (1, 2), 2, 2, masses, charges, c)
            ps = forces.phase_space

            @test ps isa PhaseSpaceData
            @test length(ps.t) == length(forces.t)
            @test length(ps.separation_distance) == length(forces.t)
            @test length(ps.radial_velocity) == length(forces.t)
            @test !isnothing(ps.theta)
            @test !isnothing(ps.angular_momentum)
        end

        @testset "Phase space with stride" begin
            forces = compute_pair_force_timeseries(
                sol_coulomb,
                (1, 2),
                2,
                2,
                masses,
                charges,
                c;
                stride = 5,
            )
            ps = forces.phase_space

            # Should have fewer points due to stride
            @test length(ps.t) < length(sol_coulomb.t) - 1
            @test length(ps.separation_distance) == length(forces.t)
        end

        @testset "Positive separation distance" begin
            forces =
                compute_pair_force_timeseries(sol_coulomb, (1, 2), 2, 2, masses, charges, c)
            ps = forces.phase_space
            @test all(ps.separation_distance .> 0)  # Particles should never collide
        end

        @testset "Phase space consistency with force data" begin
            forces =
                compute_pair_force_timeseries(sol_weber, (1, 2), 2, 2, masses, charges, 1e6)
            ps = forces.phase_space

            # Same time vector
            @test ps.t == forces.t

            # r and ṙ should be positive for approaching and negative for separating
            # Just verify they're reasonable values
            @test all(isfinite.(ps.separation_distance))
            @test all(isfinite.(ps.radial_velocity))
        end
    end

    @testset "MomentumData" begin
        @testset "compute_momentum_timeseries basic" begin
            momentum = compute_momentum_timeseries(sol_weber)

            @test momentum isa MomentumData
            @test length(momentum.t) == length(sol_weber.t)
            @test length(momentum.linear_momentum) == length(sol_weber.t)
            @test length(momentum.linear_momentum_magnitude) == length(sol_weber.t)
            @test momentum.n_particles == 2
            @test momentum.dims == 2
        end

        @testset "Linear momentum components" begin
            momentum = compute_momentum_timeseries(sol_weber)

            # Verify component matrix dimensions
            @test size(momentum.linear_momentum_components) == (length(sol_weber.t), 2)

            # Verify components match vector data
            for t = 1:length(momentum.t)
                for d = 1:2
                    @test momentum.linear_momentum_components[t, d] ==
                          momentum.linear_momentum[t][d]
                end
            end
        end

        @testset "Angular momentum present for 2D" begin
            momentum = compute_momentum_timeseries(sol_weber)

            @test !isnothing(momentum.angular_momentum)
            @test momentum.angular_momentum isa Vector{Float64}
            @test length(momentum.angular_momentum) == length(sol_weber.t)
        end

        @testset "2D angular momentum magnitude is absolute value" begin
            sys2d = HamiltonianSystem(2, 2)
            prob2d = HamiltonianProblem(
                sys2d,
                (0.0, 0.01),
                [1.0, 0.0, -1.0, 0.0],
                [0.0, -1.0, 0.0, 1.0];
                masses = [1.0, 1.0],
                charges = [0.0, 0.0],
                c = 1e6,
                dt = 0.01,
            )
            sol2d = solve(prob2d)
            momentum2d = compute_momentum_timeseries(sol2d)

            @test momentum2d.angular_momentum[1] < 0
            @test momentum2d.angular_momentum_magnitude[1] ≈
                  abs(momentum2d.angular_momentum[1])
            @test all(momentum2d.angular_momentum_magnitude .>= 0)
        end

        @testset "Momentum conservation (isolated system)" begin
            momentum = compute_momentum_timeseries(sol_weber)

            # For an isolated system, total momentum should be conserved
            # Check that momentum stays approximately constant
            P_initial = momentum.linear_momentum_magnitude[1]
            for t in eachindex(momentum.linear_momentum_magnitude)
                @test momentum.linear_momentum_magnitude[t] ≈ P_initial rtol = 1e-10
            end

            # Angular momentum should also be conserved
            L_initial = momentum.angular_momentum[1]
            for t in eachindex(momentum.angular_momentum)
                @test momentum.angular_momentum[t] ≈ L_initial rtol = 1e-10
            end
        end

        @testset "compute_momentum_timeseries with stride" begin
            momentum = compute_momentum_timeseries(sol_weber; stride = 10)

            expected_points = length(1:10:length(sol_weber.t))
            @test length(momentum.t) == expected_points
            @test length(momentum.linear_momentum) == expected_points
        end

        @testset "Validation errors" begin
            @test_throws ArgumentError compute_momentum_timeseries(sol_weber; stride = 0)
            @test_throws ArgumentError compute_momentum_timeseries(sol_weber; stride = -1)
        end

        @testset "1D system (no angular momentum)" begin
            # Create a 1D 2-body system
            sys1d = HamiltonianSystem(2, 1)
            prob1d = HamiltonianProblem(
                sys1d,
                (0.0, 0.5),
                [1.0, -1.0],
                [0.1, -0.1];
                masses = [1.0, 1.0],
                charges = [1.0, -1.0],
                c = 1e10,
                dt = 0.01,
            )
            sol1d = solve(prob1d)
            momentum = compute_momentum_timeseries(sol1d)

            @test momentum.dims == 1
            @test isnothing(momentum.angular_momentum)
            @test isnothing(momentum.angular_momentum_magnitude)
        end

        @testset "3D system" begin
            # Create a 3D 2-body system
            sys3d = HamiltonianSystem(2, 3)
            prob3d = HamiltonianProblem(
                sys3d,
                (0.0, 0.5),
                [1.0, 0.0, 0.0, -1.0, 0.0, 0.0],
                [0.0, 0.1, 0.0, 0.0, -0.1, 0.0];
                masses = [1.0, 1.0],
                charges = [1.0, -1.0],
                c = 1e10,
                dt = 0.01,
            )
            sol3d = solve(prob3d)
            momentum = compute_momentum_timeseries(sol3d)

            @test momentum.dims == 3
            @test !isnothing(momentum.angular_momentum)
            @test momentum.angular_momentum isa Vector{Vector{Float64}}
            @test length(momentum.angular_momentum[1]) == 3
            @test !isnothing(momentum.angular_momentum_magnitude)
        end
    end

    @testset "compute_total_kinetic_energy" begin
        # Zero momentum → zero KE
        @test compute_total_kinetic_energy([0.0], [1.0], 1) == 0.0
        @test compute_total_kinetic_energy([0.0, 0.0], [1.0, 2.0], 1) == 0.0

        # Single particle 1D: KE = p²/(2m)
        @test compute_total_kinetic_energy([2.0], [1.0], 1) ≈ 2.0   # 4/2
        @test compute_total_kinetic_energy([3.0], [2.0], 1) ≈ 2.25  # 9/4

        # Single particle 2D: KE = (px²+py²)/(2m)
        @test compute_total_kinetic_energy([3.0, 4.0], [1.0], 2) ≈ 12.5  # 25/2

        # Single particle 3D: KE = (px²+py²+pz²)/(2m)
        @test compute_total_kinetic_energy([1.0, 2.0, 2.0], [1.0], 3) ≈ 4.5  # 9/2

        # Two particles 1D: KE = p1²/(2m1) + p2²/(2m2) = 4/2 + 9/6 = 3.5
        @test compute_total_kinetic_energy([2.0, 3.0], [1.0, 3.0], 1) ≈ 3.5

        # Two particles 2D: particle 1 moving (p=[3,4], m=1 → KE=12.5), particle 2 at rest
        @test compute_total_kinetic_energy([3.0, 4.0, 0.0, 0.0], [1.0, 2.0], 2) ≈ 12.5
    end

    @testset "compute_pair_weber_components" begin
        # Two 2D particles along x-axis: particle 1 at origin, particle 2 at x=2, both at rest
        q_rest = [0.0, 0.0, 2.0, 0.0]
        p_rest = [0.0, 0.0, 0.0, 0.0]
        masses = [1.0, 1.0]
        charges = [1.0, -1.0]
        c = 10.0

        coulomb, vel, rdot, zextra =
            compute_pair_weber_components(q_rest, p_rest, 1, 2, masses, charges, c, 2, 1.0)

        # At rest: rdot = 0, velocity term = 0, zollner_extra = 0 (kappa=1)
        @test rdot ≈ 0.0
        @test vel ≈ 0.0
        @test zextra ≈ 0.0
        # Coulomb: kappa*q1*q2/r = 1*(1)(-1)/2 = -0.5
        @test coulomb ≈ -0.5

        # Particle 2 moving right at v=3 → rdot = +3 (separating)
        p_radial = [0.0, 0.0, 3.0, 0.0]
        coulomb2, vel2, rdot2, _ = compute_pair_weber_components(
            q_rest,
            p_radial,
            1,
            2,
            masses,
            charges,
            c,
            2,
            1.0,
        )
        @test rdot2 ≈ 3.0
        # velocity term = -coulomb2 * rdot2²/(2c²)
        @test vel2 ≈ -coulomb2 * rdot2^2 / (2 * c^2)

        # Large-c limit: velocity correction vanishes relative to Coulomb
        coulomb3, vel3, _, _ = compute_pair_weber_components(
            q_rest,
            p_radial,
            1,
            2,
            masses,
            charges,
            1e9,
            2,
            1.0,
        )
        @test abs(vel3) / abs(coulomb3) < 1e-12

        # kappa ≠ 1: Coulomb scales, Zöllner extra = (kappa-1)*q1*q2/r (rdot=0)
        kappa = 1.3
        coul_k, _, _, zextra_k = compute_pair_weber_components(
            q_rest,
            p_rest,
            1,
            2,
            masses,
            charges,
            c,
            2,
            kappa,
        )
        @test coul_k ≈ kappa * charges[1] * charges[2] / 2.0
        @test zextra_k ≈ (kappa - 1.0) * charges[1] * charges[2] / 2.0

        # 1D: particles at x=0 and x=3 converging (particle 1 right, particle 2 left)
        q_1d = [0.0, 3.0]
        p_1d = [1.0, -1.0]
        coulomb_1d, _, rdot_1d, _ =
            compute_pair_weber_components(q_1d, p_1d, 1, 2, [1.0, 1.0], [1.0, 1.0], c, 1)
        @test coulomb_1d ≈ 1.0 / 3.0  # same-sign charges, r=3
        @test rdot_1d < 0.0            # converging
    end
end
