# Example 1: Two-Body Elliptical Orbit with Precession
#
# Weber's velocity-dependent corrections cause apsidal precession in bound orbits.
# Run: julia --project=papers/Computational-Weber-Electrodynamics/examples papers/Computational-Weber-Electrodynamics/examples/example1_precession.jl

using WeberElectrodynamics
using PGFPlotsX
using LaTeXStrings

# Output directory
output_dir = joinpath(@__DIR__, "figures")
mkpath(output_dir)

println("Example 1: Two-Body Elliptical Orbit with Precession")
println("=" ^ 50)

# Parameters (from ExamplesSpecification.md)
m1, m2 = 1.0, 0.1
k = 0.1
r0 = 2.0
c_weber = 4.0

# Setup Hamiltonian using new API
function H_weber(q, p, params)
    m1, m2, k, c = params
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    KE = (px1^2 + py1^2)/(2m1) + (px2^2 + py2^2)/(2m2)

    dx, dy = x1 - x2, y1 - y2
    r = sqrt(dx^2 + dy^2)

    vx1, vy1 = px1/m1, py1/m1
    vx2, vy2 = px2/m2, py2/m2
    rdot = (dx*(vx1-vx2) + dy*(vy1-vy2)) / r

    PE = -k / r * (1 - rdot^2 / (2*c^2))
    return KE + PE
end

# Build Hamiltonian (2 particles, 2D)
H = build_hamiltonian(H_weber, 2, 2; param_names=[:m1, :m2, :k, :c])

# Initial conditions (center-of-mass frame)
M = m1 + m2
v_circ = sqrt(k * M / (m1 * m2 * r0))
v_scale = 0.7  # elliptical orbit

q0 = [-m2/M * r0, 0.0, m1/M * r0, 0.0]
p0 = [0.0, m1 * (-m2/M * v_circ * v_scale), 0.0, m2 * (m1/M * v_circ * v_scale)]

# Integration settings
T_orbit = 2π * sqrt(r0^3 * m1 * m2 / (k * M))
dt = 0.001
tspan = (0.0, 2 * T_orbit)

# Weber simulation
println("\nRunning Weber simulation (c = $c_weber)...")
prob_weber = WeberProblem(H, tspan, q0, p0;
    params=[m1, m2, k, c_weber],
    dt=dt,
    tolerance=1e-12,
    max_iterations=100
)
sol_weber = solve(prob_weber)
println("  $(length(sol_weber.t)) steps")

# Coulomb simulation (c → ∞)
println("Running Coulomb simulation (c → ∞)...")
prob_coulomb = WeberProblem(H, tspan, q0, p0;
    params=[m1, m2, k, 1e6],
    dt=dt,
    tolerance=1e-12,
    max_iterations=100
)
sol_coulomb = solve(prob_coulomb)
println("  $(length(sol_coulomb.t)) steps")

# Extract trajectory data (stride=100 for 4 orbits with dt=0.001)
stride = 100
indices = 1:stride:length(sol_weber.t)

# Weber trajectories
x1_weber = [sol_weber.q[i][1] for i in indices]
y1_weber = [sol_weber.q[i][2] for i in indices]
x2_weber = [sol_weber.q[i][3] for i in indices]
y2_weber = [sol_weber.q[i][4] for i in indices]

# Coulomb trajectories
x1_coulomb = [sol_coulomb.q[i][1] for i in indices]
y1_coulomb = [sol_coulomb.q[i][2] for i in indices]
x2_coulomb = [sol_coulomb.q[i][3] for i in indices]
y2_coulomb = [sol_coulomb.q[i][4] for i in indices]

# Generate figures
println("\nGenerating figures...")

# Figure 1: Weber trajectory (rosette pattern)
fig_weber = @pgf TikzPicture(
    Axis(
        {
            width = "0.8\\textwidth",
            height = "0.8\\textwidth",
            xlabel = L"x",
            ylabel = L"y",
            "axis equal",
            grid = "major",
            "grid style" = "{gray!30, dashed}",
            title = "Weber Electrodynamics (\$c = 4\$)"
        },
        PlotInc({"blue!70!black", "no marks", "line width" = "0.8pt"},
            Coordinates(x1_weber, y1_weber)),
        PlotInc({"red!70!black", "no marks", "line width" = "0.8pt"},
            Coordinates(x2_weber, y2_weber)),
        # Initial positions
        PlotInc({"blue!70!black", "only marks", mark = "o", "mark size" = "2pt"},
            Coordinates([x1_weber[1]], [y1_weber[1]])),
        PlotInc({"red!70!black", "only marks", mark = "o", "mark size" = "2pt"},
            Coordinates([x2_weber[1]], [y2_weber[1]]))
    )
)
pgfsave(joinpath(output_dir, "trajectory_weber.tikz"), fig_weber)
println("  trajectory_weber.tikz")

# Figure 2: Weber vs Coulomb comparison
fig_comparison = @pgf TikzPicture(
    Axis(
        {
            width = "0.8\\textwidth",
            height = "0.8\\textwidth",
            xlabel = L"x",
            ylabel = L"y",
            "axis equal",
            grid = "major",
            "grid style" = "{gray!30, dashed}",
            "legend pos" = "north east"
        },
        PlotInc({"blue!70!black", "no marks", "line width" = "0.8pt"},
            Coordinates(x2_weber, y2_weber)),
        LegendEntry("Weber (\$c=4\$)"),
        PlotInc({"red!70!black", "no marks", "line width" = "0.8pt", dashed},
            Coordinates(x2_coulomb, y2_coulomb)),
        LegendEntry("Coulomb")
    )
)
pgfsave(joinpath(output_dir, "trajectory_comparison.tikz"), fig_comparison)
println("  trajectory_comparison.tikz")

println("\nDone. Figures in: $output_dir")
