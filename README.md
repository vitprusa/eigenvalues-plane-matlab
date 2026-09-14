# Discrete sine transform (DST) based discretisation of Laplace operator

Discrete sine transform (DST) based method for the discretisation of the Dirichlet Laplacian operator on generic two-dimensional domains.

## Overview

The core idea is to discretise the Laplace operator with homogeneous
Dirichlet boundary conditions using the **discrete sine transform (DST)**,
which diagonalises the 1-D second-derivative operator on an interval. See

> Fusi, Lorenzo, Oliver Křenek, Vít Průša, Casey Rodriguez, Rebecca Tozzi, and Martin Vejvoda. "Discrete versus continuous—Linear lattice models and their exact continuous counterparts." International Journal of Engineering Science 224 (2026): 104530, [10.1016/j.ijengsci.2026.104530](https://doi.org/10.1016/j.ijengsci.2026.104530)

for reference and thorough discussion.


The 2-D Laplacian is obtained by applying the DST-based second derivative along
every grid row and column and summing the contributions.

An arbitrary domain is handled by embedding it in a rectangular **bounding
box**, sampling the box on a uniform grid, and keeping only the interior grid
points selected by an **indicator function** (a mask). Each masked grid line
splits into contiguous interior blocks, and the DST is applied to each block
with zero Dirichlet conditions at its ends.

> **Note.** The domain selected by the indicator function must be fully
> embedded in the rectangular bounding box. This is *not* checked in the
> code — it is the user's responsibility.

The repository also contains finite-difference (FD), finite-element (FEM),
and Chebyshev (Chebfun) implementations under `src/`, used to cross-check the
DST results.

## Requirements

- MATLAB.
  - **Signal Processing Toolbox** — provides `dst`/`idst`.
  - **Parallel Computing Toolbox** — the assembly uses `parfor`; it still runs
    without the toolbox, just serially.
  - **PDE Toolbox** — only for the FEM scripts in `src/fem`.
  - **Chebfun** — only for the Chebyshev script in `src/cheb`; expected in
    `fullfile(userpath, 'chebfun')` (see `startup.m`).
- Wolfram Language (`wolframscript`) — only in case you need yet another software for eigenvalues computation.

## Getting started

Launch MATLAB from the repository root (or run `startup.m`) to put the source
folders on the path:

```matlab
startup
```

### Discrete sine transform --- minimal working example

For a quick, self-contained demonstration on a single domain, the root
scripts `dst_laplace_full_spectrum.m` (full spectrum via `eig`) and
`dst_laplace_partial_spectrum.m` (leading eigenvalues via `eigs`) assemble
and solve the zero Dirichlet Laplace operator eigenvalue problem one domain inline; edit the bounding box, resolution `M`, and
indicator `phi` at the top of each.

### Discrete sine transform -- driver script

For DST-based computation run the experiment driver:

```matlab
compute_spectrum_dst                         % all domains, full + partial
compute_spectrum_dst("L_shaped")             % one domain, both modes
compute_spectrum_dst("L_shaped", "partial")  % one domain, one mode
```

The optional second argument is the **mode**, selecting which part of the
spectrum to compute (default is both, `["full" "partial"]`):

- `"full"` assembles the dense Laplace matrix (`make_dst_laplace_mat_batched`)
  and computes the **entire** spectrum with `eig`, at the coarser resolution
  `M_full` — dense `eig` scales steeply, so the grid is kept small.
- `"partial"` builds the matrix-free operator (`make_dst_laplace_op_batched`)
  and computes only the **leading `k` eigenvalues** with the iterative `eigs`,
  affordable at the much finer resolution `M_partial`.

Each (domain, mode) pair produces one CSV; the two `M` values come from the
catalog row.

The driver `compute_spectrum_dst` computes the eigenvalues for every domain in `domain_catalog_dst` and writes
one CSV per (domain, mode) into `results/eigenvalues/dst/`, each carrying a header that
records the bounding box, resolution, grid spacing, dofs, and indicator
function. To add or change a domain, edit a single row of
`experiments/domain_catalog_dst.m` (name, bounding box, indicator, and the
full/partial resolutions `M`); the number of eigenvalues `k` and the `eigs`
parameters (`subspace_dim`, `tolerance`, `max_iterations`) are set in
`compute_spectrum_dst`. From a shell, `experiments/run_dst.sh` runs the DST
experiments headless. 


### Other methods

1. Script `experiments/run_mps.sh` runs eigenvalue computation for L-shaped
   domain using the **method of particular solutions** (MPS). MPS eigenvalues
   for the L-shaped domain are the "ground truth", MPS computes them with high
   accuracy.
2. Script `experiments/run_wolfram.sh` runs eigenvalue computation for various
   domains using default `NDEigenvalues`/`NDEigensystem` solver in **Wolfram
   Language**. Wolfram Language computed eigenvalues are used for comparison.
3. Script `experiments/run_cheb.sh` runs the **Chebyshev spectral-collocation
   (Chebfun)** computation for the square [0, pi] x [0, pi] and the rectangle
   [0, 2*pi] x [0, pi], sweeping several Chebyshev orders `N` (one CSV per domain
   and `N`) — a cross-check whose leading eigenvalues match the analytic spectra.
   The domains are rows of `experiments/domain_catalog_cheb.m`, mirroring the DST
   catalog; `compute_spectrum_cheb(name, Nvals)` filters by domain and `N` list.
4. Scripts `experiments/run_fem_eig.sh` and `experiments/run_fem_solvepdeeig.sh`
   run **finite-element** (PDE Toolbox) computations on the same nine domains as
   the DST catalog (built as `decsg` geometry, including the ellipse-minus-quadrant
   and the two GWW drums) — two workflows over the same catalog
   `experiments/domain_catalog_fem.m`: `fem_eig` assembles the stiffness
   and mass matrices and solves the dense generalized problem `eig(K, M)`, while
   `fem_solvepdeeig` uses the high-level `solvepdeeig` solver. Each writes one CSV
   per domain into `results/eigenvalues/fem/`; `compute_spectrum_fem_eig(name)` and
   `compute_spectrum_fem_solvepdeeig(name)` filter by domain.
5. Script `experiments/run_fd.sh` runs a **finite-difference** (5-point stencil)
   computation on the same domains, over `experiments/domain_catalog_fd.m` — the
   same domains as the FEM catalog, but expressed as a bounding box plus an
   indicator mask and a grid resolution `M` (the DST representation). It writes
   one CSV per domain into `results/eigenvalues/fd/`; `compute_spectrum_fd(name)`
   filters by domain.

The four discretisations compared in the article — FEM, Cheb, FD and DST — are
described for the paper in
`results_paper/latex_helpers/methods_description_table.tex`, see
[LaTeX helpers](#latex-helpers-results_paperlatex_helpers) below.

## Core API (`src/dst`)

The Laplace operator is built by a small family of factory functions, all
sharing the signature `(x_range, y_range, M, indicator_function)`:

| Function | Returns | Use |
|----------|---------|-----|
| `make_dst_laplace_op` | operator handle `@(v) ...` + `info` | reference matrix-free operator (re-runs the DST sweep per apply) |
| `make_dst_laplace_op_batched` | sparse-backed operator handle + `info` | precompute-once operator, fast per apply — use with `eigs` |
| `make_dst_laplace_mat` | dense matrix `L` + `info` | reference matrix, assembled column by column |
| `make_dst_laplace_mat_batched` | dense matrix `L` + `info` | vectorised matrix assembly (much faster) |
| `make_dst_laplace_sparse` | sparse matrix `L` + `info` | shared assembly core used by the batched builders |

`info` is a struct with the grid/domain data: `h`, `x_vec`, `y_vec`, `X`, `Y`,
`global_mask`, and `dofs`.

Supporting routines: `dst_laplace_grid` (Laplacian of a full grid),
`dst_d2_slice` / `dst_d2_chunk` (1-D DST second derivative), and
`vals_vec_to_vals_grid` / `vals_grid_to_vals_vec` (scatter/gather between the
degrees-of-freedom vector and the masked grid). `dst_laplace_spectrum`
(assemble → solve → CSV via `write_eigs_csv`) is the end-to-end runner for a
single domain, driven over the catalog by `compute_spectrum_dst`.

## Domains (`src/domains`)

`bounding_box(a, b, c, d)` returns the `x_range`/`y_range` of the box. The
`indicator_*` functions return the mask of a particular domain: rectangle,
isosceles triangle, L-shape, H-shape, ellipse-minus-quadrant, and the two
isospectral GWW drums.

The domains are described for the paper in
`results_paper/latex_helpers/domains_description_table.tex`, see
[LaTeX helpers](#latex-helpers-results_paperlatex_helpers) below.

## LaTeX helpers (`results_paper/latex_helpers`)

LaTeX snippets written by hand, kept next to the
generated ones so the article can `\input` them the same way:

| file | content |
| --- | --- |
| `methods_description_table.tex` | a `booktabs` table describing the four discretisations compared in the DOF sweep, one row per method |
| `domains_description_table.tex` | a `booktabs` table describing the computational domains, one row per domain |

The methods table carries the labels the figures and the comparison tables use
(`FEM`, `Cheb`, `FD`, `DST`), and cites MATLAB, the PDE Toolbox and Chebfun. Its
label is `tab:methods_description`. The domains table
(`tab:domains_description`) names the domains by their catalog identifiers of
`experiments/domain_catalog_dst.m`, which are also the names in the grid figure
files. Nothing regenerates either: keep them in step with `src/` and the catalog
by hand.

## Repository layout

```
.                              startup.m + dst_laplace_full/partial_spectrum demo scripts
data/                          published reference eigenvalues (Betcke & Trefethen, SIAM
                               Review 47(3):469-491, 2005) as <domain>.csv, used as the
                               ground truth where no closed form is known
src/dst/                       DST Laplace operator/matrix builders + spectrum runner
src/domains/                   bounding box and domain indicator functions
src/fd/                        finite-difference spectrum runner (fd_laplace_spectrum)
src/fem/                       finite-element (PDE Toolbox) spectrum runner (fem_laplace_spectrum)
src/cheb/                      Chebfun spectral-collocation runner (chebfun_laplace_spectrum)
src/mps/                       method of particular solutions (Betcke & Trefethen)
src/wolfram/                   Wolfram region catalog + reportEigenvalues (NDEigensystem)
experiments/                   spectrum drivers (compute_spectrum_{dst,mps,wolfram,cheb,fem_eig,fem_solvepdeeig,fd}) + catalogs + run_{dst,mps,wolfram,cheb,fem_eig,fem_solvepdeeig,fd}.sh + read_eigs_csv
experiments/eigenvalues_head/  make_eigenvalues_head_comparison_tables.py (per-domain
                               comparison tables)
experiments/eigenvalues_dof_sweep/ DOF-sweep drivers and plots (compute_*_dof_sweep.m,
                               plot_*_dof_sweep.m, load_dof_sweep, run_*_dof_sweep.sh)
experiments/eigenvalues_dof_requirement/ how many DOF the first n eigenvalues cost
                               (compute_eigenvalues_dof_requirement.m,
                               run_eigenvalues_dof_requirement.sh): the largest
                               relative error over the whole block against DOF, for
                               DST/FD/FEM on all seven domains, plus Cheb on the
                               rectangle (the collocation operator is a Kronecker
                               sum, so it needs a domain filling a tensor-product
                               box); each sweep starts just above n dofs, since a
                               coarser run has no n-th eigenvalue at all, and
                               reaches ~15000, which is what a dense eig fits here;
                               the answer per tolerance is read off the sweep where
                               it reaches one and extrapolated along a fitted power
                               law where it does not. The reference is the closed
                               form on the rectangle and the triangle and a DST run
                               at M_ref (~20000 dofs) elsewhere
experiments/grid_visualisation/ DST grid + domain-mask figures (plot_grid_visualisation.m,
                               run_grid_visualisation.sh)
experiments/eigenvalues_weyl/  Weyl-asymptotics visual check (plot_weyl_asymptotics.m,
                               run_weyl_asymptotics.sh)
experiments_paper/             article (paper) variants of the figure/table experiments:
                               same computation, but no plot titles, catalog domain/method
                               names in the labels, and output to results_paper/
experiments_paper/grid_visualisation/ paper grid + domain-mask figures
                               (M and DOF encoded in the file name, no title)
experiments_paper/eigenvalues_dof_sweep/ paper DOF-sweep plots, black-and-white
                               (marker per DOF level) and colour variants, no title
experiments_paper/eigenvalues_indexing/ rectangle eigenvalue match/order tables
                               as a booktabs LaTeX snippet (eigenvalues_rectangle_tables.m)
experiments_paper/eigenvalues_head/ self-computing eigenvalue-comparison tables
                               (make_eigenvalues_head_tables.m): recomputes DST/FD/FEM
                               (and Cheb for the rectangle) at a few DOF resolutions each (four
                               everywhere, except the three FEM meshes of the non-analytic
                               domains, whose fourth DST/FD grid is sized to ~3000 dofs to
                               match the finest FEM mesh), times every run, caches the spectra
                               as CSV under `results_paper/eigenvalues_head/cache/` and writes
                               two booktabs LaTeX tables per domain (upright and transposed).
                               The digits are set per domain in `domain_configs` as `fmt_mode`
                               and `fmt_n`: `'decimals'` prints `fmt_n` digits after the point,
                               `'significant'` prints `fmt_n` significant digits by varying the
                               decimals with the magnitude, and every domain is on five
                               significant digits. Fixed notation is used while it can show
                               exactly that many digits; a cell too large (no room for a
                               decimal) or too small (a long run of leading zeros) for that
                               falls back to scientific notation at two significant digits,
                               per cell rather than per row, so that one wild value cannot
                               widen a whole column. Neither fallback fires on the present
                               spectra
experiments_paper/eigenvalues_convergence/ convergence test for one eigenvalue
                               (make_eigenvalues_convergence_index.m): lambda_k, any k, with
                               DST/FD/FEM over thirteen resolutions from ~100 to ~10000 dofs
                               (fourteen, to ~15000, on every domain but the rectangle and
                               the triangle), for all seven domains, plotted both as the
                               eigenvalue itself against DOF and as the error against the
                               reference on log-log axes with a fitted-slope triangle per
                               method, the rate fitted to the finest four runs; the reference
                               is the closed-form eigenvalue on the rectangle and the
                               triangle, the published MPS ground state elsewhere, and away
                               from the ground state a DST run of its own at M_ref, the
                               finest grid whose dense eig fits in memory (~20000 dofs, about
                               100 s); caches every run's whole spectrum
experiments_paper/eigenvalues_error_distribution/ how one run spends its accuracy
                               across the block it returns
                               (make_eigenvalues_error_distribution.m,
                               run_eigenvalues_error_distribution.sh): the relative error of
                               each of the first n eigenvalues against its index, for
                               DST/FD/FEM at one moderate resolution, the three sized to
                               carry as nearly the same DOF count as the convergence sweep
                               allows so that the figure compares methods and not
                               resolutions; the reference is a DST run at M_ref. The DST
                               series is drawn in two parts, the eigenvalues
                               pi^2*(m^2+n^2) of the unit squares apart from the rest,
                               the sine basis carrying those exactly at any resolution
results/eigenvalues/           per-method spectra CSV, a subdir each (dst, fd, fem,
                               cheb, mps, wolfram)
results/eigenvalues_head/      generated <domain>_comparison.md tables (from
                               experiments/eigenvalues_head/make_eigenvalues_head_comparison_tables.py)
results/eigenvalues_dof_sweep/ L-shaped, rectangle, and isosceles-triangle
                               DOF-sweep spectra CSV + plots;
                               experiments/eigenvalues_dof_sweep/compute_*_dof_sweep.m
                               generates the data (only the missing CSVs),
                               plot_*_dof_sweep.m draws the figure
results/eigenvalues_dof_requirement/ DOF needed for the first n eigenvalues, per method:
                               <domain>_dof_requirement_n<n>_reference-eigenvalues.csv
                               (the reference spectrum: the closed form on the
                               rectangle and the triangle, a DST run at M_ref on the
                               other five, whose own resolution is then reported as a
                               floor),
                               <domain>_dof_requirement_n<n>.csv (one row per run,
                               with the largest relative error over the block),
                               <domain>_dof_requirement_n<n>_answer.csv (the answer
                               table), <domain>_dof_requirement_n<n>.eps (each marker
                               labelled with what that run cost) and the
                               <domain>_dof_requirement_n<n>.tex figure snippet; the
                               CSVs are reused on re-runs, delete them to recompute.
                               Also richardson_extrapolation.tex, a standalone note
                               (pdflatex) on extrapolating the DST sequence to a
                               reference, validated against the published ground state
results/grid_visualisation/    per-domain <domain>_grid.png showing the DST grid and
                               domain mask (from experiments/grid_visualisation/)
results/eigenvalues_weyl/      per-domain <domain>_weyl.png checking lambda_n/n -> 4*pi/A
                               plus weyl_areas.md (from experiments/eigenvalues_weyl/)
results_paper/                 article (paper) figures and LaTeX table/figure snippets
                               produced by experiments_paper/
results_paper/grid_visualisation/ paper grid PNGs (<domain>_grid_M<M>_dofs<DOF>.png) plus
                               domains.tex
results_paper/eigenvalues_dof_sweep/ paper DOF-sweep PNGs (black-and-white and _colour) plus
                               dof_sweep.tex, dof_sweep_colour.tex
results_paper/latex_helpers/   hand-written LaTeX snippets: methods_description_table.tex
                               describing the four discretisations of src/ and
                               domains_description_table.tex describing the domains of
                               the catalog; not regenerated by any experiment
results_paper/eigenvalues_indexing/ rectangle match/order tables snippet (rectangle_tables_M9.tex)
results_paper/eigenvalues_head/ per domain (rectangle, isosceles_triangle, L_shaped,
                               ellipse_minus_quadrant, H_shaped, gww1, gww2), two booktabs
                               tables of the same data, both \scriptsize and five
                               significant digits:
                               <domain>_eigenvalues_head.tex (one row per run, eigenvalues
                               across) and <domain>_eigenvalues_head_transposed.tex (one
                               row per eigenvalue, runs across)
results_paper/eigenvalues_head/cache/ the self-computed comparison spectra
                               (<domain>_<method>_<k>-eigenvalues.csv, with dof + time
                               metadata); a run whose CSV is here is read back rather than
                               recomputed, delete one to compute that run again
results_paper/eigenvalues_convergence/ per domain and eigenvalue index
                               <domain>_convergence_lambda<k>.csv (method, resolution,
                               dofs, lambda_k, timing), the value figure
                               <domain>_convergence_lambda<k>.eps, the error figure
                               <domain>_convergence_lambda<k>_error.eps and the
                               <domain>_convergence_lambda<k>.tex figure snippet, plus the
                               combined <domain>_convergence.tex float holding every index
                               of a domain in one figure; spectra/ holds the whole spectrum
                               of every run of the sweep and of the M_ref reference run,
                               which is the cache the rest is built from, so a further index
                               costs a redraw and nothing else -- delete a domain's spectra
                               to compute it again
results_paper/eigenvalues_error_distribution/ per domain, block size and DOF level
                               <domain>_error_distribution_n<n>_dofs<D>.csv (one row per
                               index: the reference eigenvalue, then each method's
                               eigenvalue and relative error), two figures of it --
                               <domain>_error_distribution_n<n>_dofs<D>.eps with a linear
                               index and ..._log_log.eps with a logarithmic one, the first
                               for the body of the block and the second for its first few
                               indices -- each with its .tex snippet, plus
                               <domain>_error_distribution_n<n>_reference-
                               eigenvalues.csv, which carries no DOF token and is shared by
                               every DOF level of a domain; both CSVs are caches, delete
                               them to recompute
test/mat_batched/              equivalence and timing tests for the builders
test/laplace_action/           Laplace-operator action on a known function
test/manufactured_solution/    BVP solve via the method of manufactured solutions
test/laplace_matrix_symmetry/  symmetry of the assembled DST Laplace matrix on
                               every domain at four resolutions each, and what
                               symmetrising it costs and changes; spectra/ caches
                               both spectra of every matrix
```

The domains (square, rectangle, ellipse minus a quadrant, isosceles triangle,
small rectangle, H-shaped, L-shaped, and the GWW1/GWW2 isospectral drums) are defined
as rows of `experiments/domain_catalog_dst.m`. `compute_spectrum_dst` runs the
`dst_laplace_spectrum` runner over the catalog and writes the results to
`results/eigenvalues/dst/`; the source tree holds no generated per-domain scripts.

The repository also tracks the generated artifacts, so the results are
available without rerunning anything: the per-method spectra `*.csv` under
`results/`, the article (paper) figures and LaTeX table/figure snippets under
`results_paper/`, and the test reports under `test/`. They can be regenerated
with the experiments drivers (or `experiments/run_wolfram.sh` for the reference
spectra), the `experiments_paper/` runners for the paper figures and snippets,
and the `test/*/run_tests.sh` harnesses; MATLAB autosave files (`*.asv`), editor
backups (`*~`), and LaTeX build artifacts (`auto/`, `*.log`, `*.aux`, ...) are
ignored.

## Tests

The tests live under `test/`, one self-contained scenario per subdirectory:

- `test/mat_batched/` — the batched builders and operators are checked for
  equivalence against the reference implementations and benchmarked.
- `test/laplace_action/` — applies the masked DST Laplace operator
  (`make_dst_laplace_op`) on the rectangle to the analytic mode
  `u = sin(mm*x) sin(nn*y)` and compares against its exact Laplacian
  `-(mm^2 + nn^2) u`. A pure sine mode fits the grid, so the match is at
  floating-point round-off.
- `test/manufactured_solution/` — a method-of-manufactured-solutions BVP
  check on the L-shaped domain: a polynomial `u` that vanishes on the
  boundary is chosen, the discrete problem `L u = laplace u` is solved with
  the sparse DST Laplace matrix, and the recovered `u` is compared against
  the manufactured one. Being a polynomial rather than a sine mode, it
  carries a genuine discretisation error, so the test refines the grid and
  verifies the (second-order) convergence.
- `test/laplace_matrix_symmetry/` — the assembled DST Laplace matrix is
  symmetric in exact arithmetic, each one-dimensional block map being
  symmetric and scattered onto the same index set in both directions, but
  the round-off of the FFT leaves it a few units in the last place short of
  it, and `issymmetric` asks for bit-for-bit equality. The test measures how
  far short, on every domain of the catalogue at four resolutions each, and
  what symmetrising with `dst_laplace_symmetrise` is worth: it computes each
  spectrum both ways, times both, and compares them eigenvalue by
  eigenvalue. Both spectra of every matrix are cached under `spectra/`.

Each subdirectory has its own `run_tests.sh` that auto-discovers every
`check_*.m` in that folder, runs each through MATLAB, writes a markdown
report, and exits non-zero if any test errors. For example:

```bash
test/mat_batched/run_tests.sh             # writes test/mat_batched/report.md
test/laplace_action/run_tests.sh          # writes test/laplace_action/report.md
test/manufactured_solution/run_tests.sh   # writes test/manufactured_solution/report.md
test/laplace_matrix_symmetry/run_tests.sh # writes test/laplace_matrix_symmetry/report.md
```

Set `MATLAB_BIN` to override the `matlab` executable, or pass a path to
choose the report location.

## Authors

The scripts were written by Oliver Křenek, Vít Průša
(<vit.prusa@matfyz.cuni.cz>), Rebecca Tozzi and Martin Vejvoda. Vít Průša is
responsible for the conceptualisation of the work.

The core numerical routines — `dst_d2_chunk`, `dst_d2_slice`,
`dst_laplace_grid`, `vals_vec_to_vals_grid`, `vals_grid_to_vals_vec`,
`make_dst_laplace_mat`, and `make_dst_laplace_op` (all in `src/dst`) — were
written by the human authors.

The auxiliary scripts — the domain catalog (`experiments/domain_catalog_dst.m`),
the spectrum runner (`src/dst/dst_laplace_spectrum.m`) and CSV writer
(`src/dst/write_eigs_csv.m`), and the driver
(`experiments/compute_spectrum_dst.m`) — together with the testing scripts, the
batched versions of the core builders, and the documentation strings, were
written by Claude Code (Claude Opus 4.8).

## License

The whole software is distributed under the BSD 3-Clause License.
