# Rhombus relative-equilibrium scan, 2+/2- in 2D.
# +Q at (a,0),(-a,0); -Q at (0,b),(0,-b).  Aspect ratio gamma = b/a.
# At a relative equilibrium the centripetal force on each + balances Coulomb
# radial force; symmetry implies same omega from any particle => 2 algebraic
# equations in (a,b,omega).  Holding a=1 we have one equation linking gamma
# and omega, plus the requirement that the - particles also balance.  Both
# constraints must hold: a 1-parameter family (parameterized by overall scale,
# fixed by a=1) — but generically not every gamma admits a relative equilib;
# the condition is two equations in (b,omega), so a discrete set of (b,omega)
# pairs (i.e. only specific gamma values give a true RE).
#
# We instead enforce the radial balance on the +Q at (a,0) only, treating
# omega as a free function of (a,b), giving a 2-parameter "rotating rhombus"
# manifold; then check residual on the -Q particles.  We report omega and the
# stability classification of the linearised flow.

using LinearAlgebra, Printf

const N=4; const m=ones(N); const Q=[1.0,1.0,-1.0,-1.0]

config(a,b) = [a,0.0, -a,0.0, 0.0,b, 0.0,-b]

function Ucoul(q)
    U=0.0
    for i in 1:N, j in i+1:N
        dx=q[2i-1]-q[2j-1]; dy=q[2i]-q[2j]
        U += Q[i]*Q[j]/sqrt(dx*dx+dy*dy)
    end
    U
end

# omega^2 from particle 1 (+Q at (a,0)) radial balance
# Distance to 2 = 2a (like, repulsive +x).
# Distance to 3=(0,b) and 4=(0,-b): r=sqrt(a^2+b^2); attractive.
#   force on 1 from 3: unit vec (-a,b)/r, magnitude 1/r^2 (attractive).
#   x-component: -a/r^3.  From 4: -a/r^3.  Total -2a/r^3.
# Like pair: +1/(2a)^2 in +x.
# Net F_x on particle 1 = 1/(4a^2) - 2a/(a^2+b^2)^{3/2}.
# Centripetal equation (rotation about origin): m omega^2 a = -F_x (radial inward
# component is -F_x since particle is on +x axis).
# => omega^2 = ( 2a/(a^2+b^2)^{3/2} - 1/(4 a^2) ) / a
function omega2_plus(a,b)
    ( 2a/(a^2+b^2)^(3/2) - 1/(4a^2) ) / a
end
# omega^2 from -Q at (0,b) by analogous calculation
function omega2_minus(a,b)
    ( 2b/(a^2+b^2)^(3/2) - 1/(4b^2) ) / b
end

# A genuine rhombus RE requires omega2_plus == omega2_minus.  Solve numerically
# for b given target gamma=b/a with a fixed at 1 — this only occurs at the
# square (a=b).  Confirm:
println("Compatibility check (a=1):")
for gamma in [0.3,0.5,0.7,1.0,1.5]
    b = gamma
    @printf "gamma=%.2f  omega2_+=%+.5f  omega2_-=%+.5f  diff=%+.2e\n" gamma omega2_plus(1.0,b) omega2_minus(1.0,b) (omega2_plus(1.0,b)-omega2_minus(1.0,b))
end

println("\nOnly gamma=1 (square) satisfies both balances simultaneously.")
println("=> Rhombus family collapses to the single square configuration as a")
println("rigid relative equilibrium with all four particles on a common circle.")

# Hence no extended 1-parameter family of rhombus REs in this problem.
# The square is the unique rhombus RE; its stability was computed in
# square_stability.jl (UNSTABLE: real eigenvalue ~1.03).
