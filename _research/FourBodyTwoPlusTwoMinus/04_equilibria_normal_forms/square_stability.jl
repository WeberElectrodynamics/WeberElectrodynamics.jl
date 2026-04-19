# Linear stability of the alternating-square relative equilibrium
# 4-body, 2+/2-, masses=1, charges=±1, R=1, c=10 (Weber correction tiny but
# vanishes anyway at rigid rotation).
#
# Layout (2D): particle 1 = +Q at (R,0), particle 2 = +Q at (-R,0),
#              particle 3 = -Q at (0,R), particle 4 = -Q at (0,-R).
# Adjacent unlike pairs (1-3,1-4,2-3,2-4) at distance s = R*sqrt(2);
# opposite like pairs (1-2,3-4) at distance 2R.
#
# In the rotating frame at angular velocity omega, the effective potential
# in cartesian (q_i) is
#   U_eff(q) = U_Coulomb(q) - (omega^2/2) * sum_i m_i * |q_i|^2
# (centrifugal term).  The Weber velocity-dependent piece is dropped because
# at the rigid rotation point all pair distances are constant => rdot = 0,
# and the linearized correction to it is O(velocity^2). See NOTES.md.
#
# We compute:
#   omega^2 from the radial force balance on particle 1.
#   The 8x8 Hessian of U_eff at the square (positions only).
#   The full 16x16 linearized Hamilton matrix in the rotating frame
#       M = J * Hess(H_rot)
#   where H_rot = sum |p_i + m omega K q_i|^2/(2m) + U_eff and
#   K is the 2x2 rotation-generator (0,-1;1,0) acting on each (x_i,y_i).
#   Equivalently we use  H_rot = T(p) - omega * L_z + U_Coulomb
#   linearised about (q*, p* = m omega K q*).

using LinearAlgebra
using Printf

const N = 4
const dims = 2
const m = ones(N)
const Q = [1.0, 1.0, -1.0, -1.0]
const R = 1.0

# Reference square configuration
const qstar = [ R,0.0,  -R,0.0,  0.0,R,  0.0,-R ]

# Coulomb potential
function Ucoul(q)
    U = 0.0
    for i in 1:N, j in i+1:N
        dx = q[2i-1]-q[2j-1]; dy = q[2i]-q[2j]
        r = sqrt(dx*dx+dy*dy)
        U += Q[i]*Q[j]/r
    end
    U
end

# Force on particle 1 along -x direction at the square, radial component
# Particle 1 at (R,0). Distance to 2 (at -R,0): 2R, like charges => repulsive,
#   force +x (away from origin) magnitude 1/(2R)^2. Radial component (toward
#   origin = -x_hat dot F) = -1/(2R)^2.
# Distance to 3 (0,R): s = R*sqrt(2); unit vector from 1 to 3 = (-1,1)/sqrt(2).
#   Attractive, F_on_1 = +(-1,1)/sqrt(2) * 1/s^2 . Toward origin (-x) component:
#   -F_x = -(-1/sqrt(2))/s^2 = 1/(sqrt(2) s^2).
# Distance to 4 (0,-R): symmetric, -x component: 1/(sqrt(2) s^2).
# Net radial (toward origin) Coulomb force:
#   F_r = 2/(sqrt(2) s^2) - 1/(2R)^2
#       = sqrt(2)/(2 R^2) - 1/(4 R^2)
#       = (2 sqrt(2) - 1) / (4 R^2)
# Centripetal requirement: m omega^2 R = F_r.
#   omega^2 = (2 sqrt(2) - 1) / (4 m R^3)
const omega2 = (2*sqrt(2) - 1) / (4 * 1.0 * R^3)
const omega  = sqrt(omega2)
println("omega^2 = ", omega2, "   omega = ", omega)

# Hessian of U_eff = U_Coulomb - (omega^2/2) sum m_i |q_i|^2  via finite diff
function Ueff(q)
    U = Ucoul(q)
    for i in 1:N
        U -= 0.5*omega2*m[i]*(q[2i-1]^2 + q[2i]^2)
    end
    U
end

function hessian(f, x; h=1e-5)
    n = length(x); H = zeros(n,n)
    for i in 1:n, j in i:n
        xpp=copy(x); xpp[i]+=h; xpp[j]+=h
        xpm=copy(x); xpm[i]+=h; xpm[j]-=h
        xmp=copy(x); xmp[i]-=h; xmp[j]+=h
        xmm=copy(x); xmm[i]-=h; xmm[j]-=h
        H[i,j] = (f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4h*h)
        H[j,i] = H[i,j]
    end
    H
end

Hq = hessian(Ueff, qstar)
println("\nHessian of U_eff (8x8) eigenvalues:")
evq = eigvals(Symmetric(Hq))
for v in evq; @printf "%+.6f  " v; end; println()

# Now the FULL linearised Hamilton matrix in the rotating frame.
# H_rot(q,p) = sum |p_i|^2/(2 m_i) + U_Coulomb(q) - omega * L_z
# with L_z = sum_i (x_i p_yi - y_i p_xi).
# At the relative equilibrium q=qstar, p_i* = m_i omega ( -y_i, x_i ).
# Let z = (q,p) - (qstar,pstar). Then linear dynamics dz/dt = J Hess(H_rot) z.

function Hrot(z)
    q = z[1:8]; p = z[9:16]
    T = 0.0
    for i in 1:N
        T += (p[2i-1]^2 + p[2i]^2)/(2m[i])
    end
    Lz = 0.0
    for i in 1:N
        Lz += q[2i-1]*p[2i] - q[2i]*p[2i-1]
    end
    return T + Ucoul(q) - omega*Lz
end

pstar = zeros(8)
for i in 1:N
    pstar[2i-1] = -m[i]*omega*qstar[2i]
    pstar[2i]   =  m[i]*omega*qstar[2i-1]
end
zstar = vcat(qstar, pstar)

H16 = hessian(Hrot, zstar; h=1e-5)
J = zeros(16,16)
for i in 1:8
    J[i, i+8] = 1.0
    J[i+8, i] = -1.0
end
M = J*H16
ev = eigvals(M)
println("\nLinearised Hamilton-matrix eigenvalues (16):")
for v in ev
    @printf "%+.6f %+.6fi\n" real(v) imag(v)
end

# Classify
re = real.(ev); im_ = imag.(ev)
nzero = count(e -> abs(e) < 1e-6, ev)
nimag = count(i -> abs(re[i]) < 1e-6 && abs(im_[i]) > 1e-6, eachindex(ev))
nreal = count(i -> abs(re[i]) > 1e-6 && abs(im_[i]) < 1e-6, eachindex(ev))
ncplx = 16 - nzero - nimag - nreal
println("\nzero=$nzero  pure-imag=$nimag  pure-real=$nreal  complex=$ncplx")
