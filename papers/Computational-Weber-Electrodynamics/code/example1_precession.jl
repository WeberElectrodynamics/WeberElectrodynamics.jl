# Example 1: Two-Body Elliptical Orbit with Precession
#
# Weber's velocity-dependent potential causes apsidal precession in bound orbits,
# creating a characteristic rosette pattern instead of a closed ellipse.
#
# Run: julia --project=code code/example1_precession.jl

using Pkg

# Install dependencies if not present
for pkg in ["PGFPlotsX"]
    if !haskey(Pkg.project().dependencies, pkg)
        Pkg.add(pkg)
    end
end

# Add local WeberElectrodynamics package if not present
if !haskey(Pkg.project().dependencies, "WeberElectrodynamics")
    Pkg.develop(path=joinpath(@__DIR__, "..", "..", ".."))
end

using WeberElectrodynamics
using PGFPlotsX

# Output directory for generated figures
output_dir = joinpath(@__DIR__, "..", "figures")
mkpath(output_dir)

# Physical parameters
m1, m2 = 1.0, 0.1    # masses
k = 0.1              # coupling constant
r0 = 2.0             # initial separation
c_weber = 4.0        # Weber constant (finite speed of interaction)

# Weber Hamiltonian: kinetic energy + velocity-dependent potential
function H_weber(q, p, params)
    m1, m2, k, c = params
    x1, y1, x2, y2 = q
    px1, py1, px2, py2 = p

    # Kinetic energy
    KE = (px1^2 + py1^2)/(2m1) + (px2^2 + py2^2)/(2m2)

    # Separation and radial velocity
    dx, dy = x1 - x2, y1 - y2
    r = sqrt(dx^2 + dy^2)
    vx1, vy1 = px1/m1, py1/m1
    vx2, vy2 = px2/m2, py2/m2
    rdot = (dx*(vx1-vx2) + dy*(vy1-vy2)) / r

    # Weber potential with velocity-dependent correction
    PE = -k / r * (1 - rdot^2 / (2*c^2))

    return KE + PE
end

# Compile Hamiltonian (2 particles in 2D)
H = compile_hamiltonian(H_weber, 2, 2; parameter_names=[:m1, :m2, :k, :c])

# Initial conditions in center-of-mass frame
M = m1 + m2
v_circ = sqrt(k * M / (m1 * m2 * r0))  # circular orbit velocity
v_scale = 0.7                           # scale down for elliptical orbit

q0 = [-m2/M * r0, 0.0, m1/M * r0, 0.0]
p0 = [0.0, m1 * (-m2/M * v_circ * v_scale), 0.0, m2 * (m1/M * v_circ * v_scale)]

# Integration settings
T_orbit = 2π * sqrt(r0^3 * m1 * m2 / (k * M))
dt = 0.001
tspan = (0.0, 2 * T_orbit)

# Run simulation
prob = WeberProblem(H, tspan, q0, p0; params=[m1, m2, k, c_weber], dt=dt)
sol = solve(prob)

# Extract trajectory data (subsample for plotting)
stride = 100
indices = 1:stride:length(sol.t)
x1 = [sol.q[i][1] for i in indices]
y1 = [sol.q[i][2] for i in indices]
x2 = [sol.q[i][3] for i in indices]
y2 = [sol.q[i][4] for i in indices]

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
            Coordinates(x1, y1)),
        PlotInc({"red!70!black", "no marks", "line width" = "0.6pt"},
            Coordinates(x2, y2)),
        # Mark initial positions
        PlotInc({"blue!70!black", "only marks", mark = "o", "mark size" = "2pt"},
            Coordinates([x1[1]], [y1[1]])),
        PlotInc({"red!70!black", "only marks", mark = "o", "mark size" = "2pt"},
            Coordinates([x2[1]], [y2[1]]))
    )
)

pgfsave(joinpath(output_dir, "example1_trajectory.tikz"), fig)
println("Saved: figures/example1_trajectory.tikz")
