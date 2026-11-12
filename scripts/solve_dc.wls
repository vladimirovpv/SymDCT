#!/usr/bin/env wolframscript
(* ::Package:: *)

(* SymDCT: symbolic diffusion coefficient from a Stage-2 matrix file.

   This script is an alternative to Stage 3 (src.transport_coefficient_solver)
   for reduced matrices that are too large for SymPy.  It evaluates the same
   expression,

       D = -(1/2) d^2 l^2  (d^2/dt^2 Det[A(t)])|_{t=0} / Det[B],

   where A(t) is the reduced transition-rate matrix, B is A(0) with the first
   row replaced by ones, d^2 is the squared interplanar distance and l is the
   minimal jump length.

   Usage:
       uv run python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
       uv run python -m src.matrix_builder examples/hcp_vac/input/hcp_vac.yaml
       wolframscript -file scripts/solve_dc.wls examples/hcp_vac/output/tp_matrix_axis_0_0_1.m

   Input (both written by Stage 2 into the same folder):
       tp_matrix_axis_<h>_<k>_<l>.m   reduced matrix A(t)
       tc_axis_meta_<h>_<k>_<l>.m     association with d^2 and l

   Output:
       dc_symbolic_<h>_<k>_<l>_wl.m   diffusion coefficient in the jump rates W0, W1, ...
*)

fail[message_String] := (WriteString["stderr", "Error: " <> message <> "\n"]; Exit[1]);

args = If[ListQ[$ScriptCommandLine] && Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine], {}];
If[Length[args] < 1,
  fail["missing argument.\nUsage: wolframscript -file scripts/solve_dc.wls <path>/tp_matrix_axis_<h>_<k>_<l>.m"]];

matrixFile = ExpandFileName[First[args]];
If[! FileExistsQ[matrixFile], fail["file not found: " <> matrixFile]];

prefix = "tp_matrix_axis_";
baseName = FileBaseName[matrixFile];
If[! StringStartsQ[baseName, prefix],
  fail["the file name must have the form " <> prefix <> "<h>_<k>_<l>.m"]];
axisLabel = StringDrop[baseName, StringLength[prefix]];
folder = DirectoryName[matrixFile];
metaFile = FileNameJoin[{folder, "tc_axis_meta_" <> axisLabel <> ".m"}];
If[! FileExistsQ[metaFile],
  fail["file not found: " <> metaFile <> "\nRun Stage 2 again: python -m src.matrix_builder <yaml_config>"]];

matrix = Get[matrixFile];
If[! (MatrixQ[matrix] && Equal @@ Dimensions[matrix]), fail["the matrix file does not contain a square matrix."]];
n = Length[matrix];

meta = Get[metaFile];
If[! AssociationQ[meta], fail["the meta file does not contain an association."]];
d2 = meta["interplanarDistanceSquared"];
jumpLength = meta["minJumpLength"];
If[MissingQ[d2] || MissingQ[jumpLength], fail["the meta file lacks interplanarDistanceSquared or minJumpLength."]];

Print["Axis: ", axisLabel, ", matrix size: ", n, " x ", n];

(* Rate symbols merged by the user (tweak_transition_rates) *)
matrix = matrix /. Lookup[meta, "wSubstitutions", {}];

numerator = D[Det[matrix], {t, 2}] /. t -> 0;

b = matrix;
b[[1]] = ConstantArray[1, n];
denominator = Det[b /. t -> 0];
If[PossibleZeroQ[denominator], fail["the denominator is zero; the matrix is singular at t = 0."]];

dc = Simplify[-(1/2) d2 jumpLength^2 numerator/denominator];

outFile = FileNameJoin[{folder, "dc_symbolic_" <> axisLabel <> "_wl.m"}];
Put[dc, outFile];

Print["Jump rates: ", Lookup[meta, "jumpLabels", {}]];
Print["D = ", ToString[dc, InputForm]];
Print["Written: ", outFile];
