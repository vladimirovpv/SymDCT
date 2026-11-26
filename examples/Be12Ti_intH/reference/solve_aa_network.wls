#!/usr/bin/env wolframscript
(* D_aa of Be12Ti-H from the variational (Trinkle / Thomson) form with detailed balance.

   Conductances c_ij = p_i W_ij = p_j W_ji, Z = Sum_i n_i p_i  (n = A1 B16 C8 D8 E8 F2 G8).
   D_aa = a^2 * Cnet/Z  with  Cnet = Emin/288,
   Emin = min over the 8 periodic corrections u of the quadratic form E(u) below.
   The symmetry x -> -x fixes phi = 0 for A0, E0, F0 and pairs the other classes with phi_i = -phi_j.

   Usage:  wolframscript -file examples/Be12Ti_intH/reference/solve_aa_network.wls   (from the public folder)
*)
Eform = 4*cAC*uC0^2 - 24*cAC*uC0 + 36*cAC + 4*cBC*uB0^2 - 8*cBC*uB0*uC0 + 16*cBC*uB0 + 4*cBC*uB1^2 - 8*cBC*uB1*uC0 - 24*cBC*uB1 + 8*cBC*uC0^2 + 8*cBC*uC0 + 52*cBC + 4*cBD*uB0^2 - 8*cBD*uB0*uD0 + 16*cBD*uB0 + 4*cBD*uB1^2 - 8*cBD*uB1*uD1 - 8*cBD*uB1 + 4*cBD*uD0^2 - 16*cBD*uD0 + 4*cBD*uD1^2 + 8*cBD*uD1 + 20*cBD + 4*cBG*uB0^2 - 8*cBG*uB0*uG0 + 16*cBG*uB0 + 4*cBG*uB1^2 + 8*cBG*uB1*uG1 + 16*cBG*uB1 + 4*cBG*uG0^2 - 16*cBG*uG0 + 4*cBG*uG1^2 + 16*cBG*uG1 + 32*cBG + 2*cEF*uE1^2 - 8*cEF*uE1 + 8*cEF + 4*cEG*uE1^2 - 8*cEG*uE1*uG1 + 16*cEG*uE1 + 4*cEG*uG0^2 - 24*cEG*uG0 + 4*cEG*uG1^2 - 16*cEG*uG1 + 52*cEG + 2*cGD*uD0^2 - 4*cGD*uD0*uG1 + 4*cGD*uD0 + 2*cGD*uD1^2 + 4*cGD*uD1*uG0 + 16*cGD*uD1 + 2*cGD*uG0^2 + 16*cGD*uG0 + 2*cGD*uG1^2 - 4*cGD*uG1 + 34*cGD + 8*cBB*uB0^2 - 16*cBB*uB0 + 8*cBB + 2*cDD*uD0^2 - 4*cDD*uD0*uD1 + 8*cDD*uD0 + 2*cDD*uD1^2 - 8*cDD*uD1 + 8*cDD;

uvars = {uB0, uB1, uC0, uD0, uD1, uE1, uG0, uG1};
cs = {cAC, cBC, cBD, cBG, cEF, cEG, cGD, cBB, cDD};

(* Schur complement of the quadratic form: E(u) = E0 + b.u + u.M.u / 2 *)
E0 = Eform /. Thread[uvars -> 0];
bvec = D[Eform, #] & /@ uvars /. Thread[uvars -> 0];
mmat = Table[D[Eform, uvars[[i]], uvars[[j]]], {i, 8}, {j, 8}];
Print["det M = ", Factor[Det[mmat]] // InputForm];
usol = LinearSolve[mmat, -bvec];
Emin = Together[E0 + bvec . usol/2];
Cnet = Together[Emin/288];
Print["Cnet = ", Cnet // InputForm];
Print["Leaf count: ", LeafCount[Cnet]];

(* optimised evaluation sequence *)
opt = Experimental`OptimizeExpression[Cnet, OptimizationLevel -> 2];
Print["Optimised leaf count: ", LeafCount[opt]];
Put[Cnet, "examples/Be12Ti_intH/output/Cnet_aa.m"];
Put[opt, "examples/Be12Ti_intH/output/Cnet_aa_optimized.m"];

(* check against the Mathematica result of solve_dc.wls (p = 1: W_ij = W_ji = c_ij) *)
ref = Get["examples/Be12Ti_intH/output/dc_symbolic_1_0_0_wl.m"];
map = {W0 -> cAC, W16 -> cAC, W1 -> cEG, W6 -> cEG, W2 -> cEF, W18 -> cEF, W5 -> cGD, W15 -> cGD,
  W9 -> cBD, W13 -> cBD, W10 -> cBC, W17 -> cBC, W11 -> cBG, W7 -> cBG, W8 -> cBB, W14 -> cDD};
Do[
  vals = Thread[cs -> RandomInteger[{1, 30}, 9]/RandomInteger[{1, 6}, 9]];
  Print["check: ", (51 (ref /. map)/a^2 /. vals) - (Cnet /. vals)],
  {3}];
