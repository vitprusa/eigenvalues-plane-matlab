(* ::Package:: *)

(* Core library: domain regions and the reportEigenvalues function. Loaded by
   the driver experiments/compute_wolfram_spectra.wls, which owns the output
   path and the run loop. This file only defines symbols -- it deliberately
   does not ClearAll or SetDirectory, so that Get-ting it does not disturb the
   driver's state. Open it in the front end to use the regions interactively. *)


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


(* ::Subsection:: *)
(*Region catalog*)


(* ::Text:: *)
(*Association mapping output name -> region, iterated over by the driver. The keys become the CSV file stems.*)


wolframRegions = <|
	"triangle"               -> tRegion,
	"rectangle"              -> rRegion,
	"L-shaped"               -> lRegion,
	"ellipse-minus-quadrant" -> ellRegion,
	"H-shaped"               -> hRegion
|>;


(* ::Section:: *)
(*Known eigenvalues*)


(* ::Text:: *)
(*Known reference values, L - shaped region, (Betcke & Trefethen, Computed eigenmodes of planar regions, 2005). See also Fox, L., P. Henrici, and C. Moler. "Approximations and bounds for eigenvalues of elliptic operators." SIAM Journal on Numerical Analysis 4, no. 1 (1967): 89-102.*)


refEigenvaluesLRegion = {9.6397238, 15.19725, 2*\[Pi]^2, 29.5214811, 31.9126360, 41.4745099, 44.9484860, 5*\[Pi]^2};


(* ::Section:: *)
(*Eigenvalues computation*)


reportEigenvalues[outDir_String, name_String, specRegion_, nEigenvalues_] := Module[{eigenvalues, eigenfunctions, rows, data, file, table, figures},
	{eigenvalues, eigenfunctions} = NDEigensystem[ {-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], {x, y} \[Element] specRegion, nEigenvalues];
	rows = Table[{n, eigenvalues[[n]]}, {n, 1, nEigenvalues}];
	data = Join[{{"n", "lambda_n"}}, rows];
	file = FileNameJoin[{outDir, name <> "-eigenvalues.csv"}];
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


(* The run loop (which domains, how many eigenvalues, and where to write the
   CSVs) lives in the driver experiments/compute_wolfram_spectra.wls. *)
