#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]
(* Work in the directory holding this file, whether run as a script
   (wolframscript sets $InputFileName) or inside the front end. *)
SetDirectory[If[$InputFileName =!= "", DirectoryName[$InputFileName], NotebookDirectory[]]]


(* ::Title:: *)
(*Eigenvalues and eigenfunctions in special regions*)


(* ::Section:: *)
(*Regions*)


(* ::Subsection:: *)
(*Triangle*)


(* ::Text:: *)
(*Triangular region defined by vertices.*)


v1 = {0, 0};
v2 = {1, 0};
v3 = {0, 1};

tRegion := Triangle[{v1, v2, v3}];
RegionPlot[tRegion]


(* ::Subsection:: *)
(*Rectangle*)


rRegion := Rectangle[{0, 0}, {\[Pi], 2*\[Pi]}]

RegionPlot[rRegion]


(* ::Subsection:: *)
(*L-shaped*)


lRegion := RegionDifference[
  Rectangle[{0, 0}, {2, 2}],
  Rectangle[{0, 0}, {1, 1}]
];

RegionPlot[lRegion] (* This looks strange. *)
Region[lRegion]


(* ::Subsection:: *)
(*Ellipse minus quadrant*)


ellRegion := RegionDifference[
	Ellipsoid[{0, 0}, {2, 1}],
	Rectangle[{0, -1}, {2, 0}]
]


RegionPlot[ellRegion, AspectRatio->Automatic]
Region[ellRegion]


(* ::Subsection:: *)
(*H-shaped region*)


hRegion := RegionUnion[
  (* Left column: 3 stacked unit squares *)
  Rectangle[{-1, -2}, {0, -1}],  (* bottom-left *)
  Rectangle[{-1, -1}, {0,  0}],  (* middle-left *)
  Rectangle[{-1,  0}, {0,  1}],  (* top-left *)
  (* Crossbar: centre square *)
  Rectangle[{0, -1},  {1,  0}],  (* middle-centre *)
  (* Right column: 3 stacked unit squares *)
  Rectangle[{1, -2},  {2, -1}],  (* bottom-right *)
  Rectangle[{1, -1},  {2,  0}],  (* middle-right *)
  Rectangle[{1,  0},  {2,  1}]   (* top-right *)
]

RegionPlot[hRegion]
Region[hRegion]


(* ::Section:: *)
(*Known eigenvalues*)


(* ::Text:: *)
(*Known reference values, L - shaped region, (Betcke & Trefethen, Computed eigenmodes of planar regions, 2005). See also Fox, L., P. Henrici, and C. Moler. "Approximations and bounds for eigenvalues of elliptic operators." SIAM Journal on Numerical Analysis 4, no. 1 (1967): 89-102.*)


refEigenvaluesLRegion = {9.6397238, 15.19725, 2*\[Pi]^2, 29.5214811, 31.9126360, 41.4745099, 44.9484860, 5*\[Pi]^2};


(* ::Section:: *)
(*Eigenvalues computation*)


reportEigenvalues[specRegion_, nEigenvalues_, name_String] := Module[{eigenvalues, eigenfunctions, rows, data, file, table, figures},
	{eigenvalues, eigenfunctions} = NDEigensystem[ {-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], {x, y} \[Element] specRegion, nEigenvalues];
	rows = Table[{n, eigenvalues[[n]]}, {n, 1, nEigenvalues}];
	data = Join[{{"n", "lambda_n"}}, rows];
	file = name <> "-eigenvalues.csv";
	Export[file, data];
	table = TableForm[
		rows,
		TableHeadings -> {None, {"n", "\!\(\*SubscriptBox[\(\[Lambda]\), \(n\)]\)"}}
	];
	figures = Grid[
		Partition[
			Table[Plot3D[eigenfunctions[[n]], {x, y} \[Element] specRegion, ImageSize->Medium, PlotLabel->StringTemplate["\!\(\*SubscriptBox[\(\[Lambda]\), \(`1`\)]\) = `2`"][n, eigenvalues[[n]]]], {n, 1, nEigenvalues}],
			UpTo[4]
		]
	];
	Print[table];
	Print[figures];
	Print["Wrote ", nEigenvalues, " eigenvalues to ", file];
]


(* ::Section:: *)
(*Results*)


(* ::Subsection:: *)
(*Triangle*)


nEigenvalues = 40;
reportEigenvalues[tRegion, nEigenvalues, "triangle"]


(* ::Subsection:: *)
(*Rectangle*)


nEigenvalues = 40;
reportEigenvalues[rRegion, nEigenvalues, "rectangle"]


(* ::Subsection:: *)
(*L-shaped*)


nEigenvalues = 40;
reportEigenvalues[lRegion, nEigenvalues, "L-shaped"]


(* ::Text:: *)
(*These are Trefethen' s eigenvalues.*)


TableForm[Table[{n, N@refEigenvaluesLRegion[[n]]}, {n, 1, Length[refEigenvaluesLRegion]}], TableHeadings -> {None, {"n", "\!\(\*SubscriptBox[\(\[Lambda]\), \(n\)]\)"}}]


(* ::Subsection:: *)
(*Ellipse  minus quadrant*)


nEigenvalues = 40;
reportEigenvalues[ellRegion, nEigenvalues, "ellipse-minus-quadrant"]


(* ::Subsection:: *)
(*H-shaped region*)


nEigenvalues = 40;
reportEigenvalues[hRegion, nEigenvalues, "H-shaped"]
