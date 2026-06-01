#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Driver: compute the reference Dirichlet-Laplacian spectra on the catalog of
   regions with NDEigensystem and write one CSV (columns n, lambda_n) per
   domain into results/wolfram/. This mirrors the MATLAB experiments drivers:
   the regions and the reportEigenvalues runner live in the core library
   src/wolfram/eigenvalues-various-domains.wls; this file owns the output path
   and the run loop. Run headless with experiments/run_wolfram.sh. *)

ClearAll["Global`*"]

(* Resolve the project root from this driver's location (experiments/). *)
thisDir = If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]];
projectRoot = ParentDirectory[thisDir];

(* Load the core: the region catalog (wolframRegions), the reportEigenvalues
   runner, and the reference eigenvalues (refEigenvaluesLRegion). *)
Get[FileNameJoin[{projectRoot, "src", "wolfram", "eigenvalues-various-domains.wls"}]];

(* Output directory for the generated CSVs (created if missing). *)
outDir = FileNameJoin[{projectRoot, "results", "wolfram"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

nEigenvalues = 40;

(* Compute and export every region in the catalog. *)
KeyValueMap[
	reportEigenvalues[outDir, #1, #2, nEigenvalues] &,
	wolframRegions
];

(* Reference eigenvalues for the L-shaped region (Betcke & Trefethen, 2005). *)
Print["Reference (Betcke & Trefethen) L-shaped eigenvalues:"];
Print[TableForm[
	Table[{n, N@refEigenvaluesLRegion[[n]]}, {n, 1, Length[refEigenvaluesLRegion]}],
	TableHeadings -> {None, {"n", "lambda_n"}}
]];
