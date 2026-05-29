# eigenvalues-plane-matlab

FFT-based method for the discretisation of second-order differential
operators, used to compute the eigenvalues of the Dirichlet Laplacian on
two-dimensional domains.

## Overview

The core idea is to discretise the Laplace operator with homogeneous
Dirichlet boundary conditions using the **discrete sine transform (DST)**,
which diagonalises the 1-D second-derivative operator on an interval. The
2-D Laplacian is obtained by applying the DST-based second derivative along
every grid row and column and summing the contributions.

An arbitrary domain is handled by embedding it in a rectangular **bounding
box**, sampling the box on a uniform grid, and keeping only the interior grid
points selected by an **indicator function** (a mask). Each masked grid line
splits into contiguous interior blocks, and the DST is applied to each block
with zero Dirichlet conditions at its ends.

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

## Getting started

Launch MATLAB from the repository root (or run `startup.m`) to put the source
folders on the path:

```matlab
startup
```

Then run one of the driver scripts:

```matlab
dst_laplace_full_spectrum      % full spectrum via eig
dst_laplace_partial_spectrum   % leading eigenvalues via eigs
```

Edit the bounding box (`a, b, c, d`), the resolution `M`, and the indicator
function `phi` at the top of each script to change the domain and grid. In
`dst_laplace_partial_spectrum` you can also tune the number of eigenvalues
`k` and the `eigs` parameters (`subspace_dim`, `tolerance`, `max_iterations`).

> **Note.** The domain selected by the indicator function must be fully
> embedded in the rectangular bounding box. This is *not* checked in the
> code — it is the user's responsibility.

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
degrees-of-freedom vector and the masked grid).

## Domains (`src/domains`)

`bounding_box(a, b, c, d)` returns the `x_range`/`y_range` of the box. The
`indicator_*` functions return the mask of a particular domain: rectangle,
isosceles triangle, L-shape, H-shape, ellipse-minus-quadrant, and the two
isospectral GWW drums.

## Repository layout

```
.                              driver scripts + startup.m
src/dst/                       DST-based Laplace operator and matrix builders
src/domains/                   bounding box and domain indicator functions
src/fd/                        finite-difference cross-checks
src/fem/                       finite-element cross-checks (PDE Toolbox)
src/cheb/                      Chebyshev (Chebfun) cross-check
test/mat_batched/              equivalence and timing tests for the builders
```

## Tests

The batched builders and operators are checked for equivalence against the
reference implementations and benchmarked under `test/mat_batched`. Run the
whole suite and produce a markdown report with:

```bash
test/mat_batched/run_tests.sh           # writes test/mat_batched/report.md
```

The runner auto-discovers every `check_*.m` in that folder, runs each through
MATLAB, and exits non-zero if any test errors. Set `MATLAB_BIN` to override
the `matlab` executable, or pass a path to choose the report location.

## Authors

The scripts were written by Oliver Křenek, Vít Průša
(<vit.prusa@matfyz.cuni.cz>), Rebecca Tozzi and Martin Vejvoda. Vít Průša is
responsible for the conceptualisation of the work.

The testing scripts, the batched versions of the initial scripts, and the
documentation strings were written by Claude Code (Claude Opus 4.8).

## License

The whole software is distributed under the BSD 3-Clause License.
