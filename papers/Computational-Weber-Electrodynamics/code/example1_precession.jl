# Example 1: Two-Body Elliptical Orbit with Precession
#
# Weber's velocity-dependent potential causes apsidal precession in bound orbits,
# creating a characteristic rosette pattern instead of a closed ellipse.
#
# Run: julia --project=code code/example1_precession.jl

using WeberElectrodynamics
using PGFPlotsX

# Output directory for generated figures
output_dir = joinpath(@__DIR__, "..", "figures")
mkpath(output_dir)

# Physical parameters
m1, m2 = 1.0, 0.1           # masses
q1, q2 = 1.0, -0.1          # charges (opposite signs for attraction)
r0 = 2.0                    # initial separation
c_weber = 4.0               # Weber constant (finite speed of interaction)

# Create Weber system (2 particles in 2D)
system = WeberSystem(2, 2; masses=[m1, m2], charges=[q1, q2], c=c_weber)

# Initial conditions in center-of-mass frame
M = m1 + m2
k = abs(q1 * q2)            # effective coupling strength
v_circ = sqrt(k * M / (m1 * m2 * r0))  # circular orbit velocity
v_scale = 0.7               # scale down for elliptical orbit

q0 = [-m2/M * r0, 0.0, m1/M * r0, 0.0]
p0 = [0.0, m1 * (-m2/M * v_circ * v_scale), 0.0, m2 * (m1/M * v_circ * v_scale)]

# Integration settings
T_orbit = 2π * sqrt(r0^3 * m1 * m2 / (k * M))
dt = 0.001
tspan = (0.0, 2 * T_orbit)

# Run simulation
prob = WeberProblem(system, tspan, q0, p0; dt=dt)
sol = solve(prob)

# Extract trajectory data (subsample for plotting)
stride = 100
indices = 1:stride:length(sol.t)
x1_traj = [sol.q[i][1] for i in indices]
y1_traj = [sol.q[i][2] for i in indices]
x2_traj = [sol.q[i][3] for i in indices]
y2_traj = [sol.q[i][4] for i in indices]

# Generate figure
fig = @pgf TikzPicture(
    Axis(
        {
            width = "0.8\\textwidth",
            height = "0.8\\textwidth",
            xlabel = raw"$x$",
            ylabel = raw"$y$",
            "axis equal",
            grid = "major",
            "grid style" = "{gray!30}"
        },
        PlotInc({"blue!70!black", "no marks", "line width" = "0.6pt"},
            Coordinates(x1_traj, y1_traj)),
        PlotInc({"red!70!black", "no marks", "line width" = "0.6pt"},
            Coordinates(x2_traj, y2_traj)),
        # Mark initial positions
        PlotInc({"blue!70!black", "only marks", mark = "o", "mark size" = "2pt"},
            Coordinates([x1_traj[1]], [y1_traj[1]])),
        PlotInc({"red!70!black", "only marks", mark = "o", "mark size" = "2pt"},
            Coordinates([x2_traj[1]], [y2_traj[1]]))
    )
)

pgfsave(joinpath(output_dir, "example1_trajectory.tikz"), fig)
println("Saved: figures/example1_trajectory.tikz")
