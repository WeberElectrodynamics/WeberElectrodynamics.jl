#!/usr/bin/env wolframscript
(* Hello World in Wolfram Language — demonstrates basic capability. *)

Print["=== Hello from Wolfram Language ==="];
Print["Kernel version: ", $Version];
Print["Front end: ", $FrontEnd];
Print[""];

(* --- Symbolic --- *)
Print["Symbolic integration:"];
Print["  Integrate[Sin[x]^2, x] = ", Integrate[Sin[x]^2, x]];

Print["Symbolic sum (closed form):"];
Print["  Sum[k^2, {k,1,n}] = ", Sum[k^2, {k, 1, n}]];

(* --- Numeric --- *)
Print["High-precision numerics:"];
Print["  N[Pi, 50] = ", N[Pi, 50]];

(* --- Linear algebra --- *)
Print["Eigenvalues of a 3x3 matrix:"];
m = {{2, 1, 0}, {1, 2, 1}, {0, 1, 2}};
Print["  M = ", m];
Print["  Eigenvalues[M] = ", Eigenvalues[m]];

(* --- ODE (relevant to this repo's physics work) --- *)
Print["Solve a simple harmonic oscillator ODE:"];
sol = DSolve[{y''[t] + y[t] == 0, y[0] == 1, y'[0] == 0}, y[t], t];
Print["  y''+y=0, y(0)=1, y'(0)=0  =>  y(t) = ", y[t] /. First[sol]];

Print[""];
Print["=== Demo complete ==="];
