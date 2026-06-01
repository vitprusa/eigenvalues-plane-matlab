# L-shaped domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
L-shaped domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not
timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Wolfram | MPS (ground truth) |
|---|----------|-------------|------|---------|-----------------|---------|--------------------|
| **DOFs** | 1776 | 66901 | 1776 | 2051 | 5653 | — ‡ | 51 † |
| **Time (s)** | 1.46 | 170.24 | 1.03 | 2.09 | 2.21 | — | 4.72 |
| 1 |  9.67525 |  9.64298 |  9.66133 |  9.64999 |  9.64532 |  9.65934 |  **9.63972** |
| 2 | 15.19657 | 15.19725 | 15.17674 | 15.19754 | 15.19732 | 15.19878 | **15.19725** |
| 3 | 19.73921 | 19.73921 | 19.71325 | 19.73946 | 19.73925 | 19.74141 | **19.73921** |
| 4 | 29.52119 | 29.52148 | 29.44785 | 29.52240 | 29.52163 | 29.52854 | **29.52148** |
| 5 | 31.99865 | 31.92059 | 31.89650 | 31.93859 | 31.92643 | 31.96703 | **31.91264** |
| 6 | 41.53975 | 41.48049 | 41.34627 | 41.49560 | 41.48512 | 41.52596 | **41.47451** |
| 7 | 44.94498 | 44.94846 | 44.72497 | 44.95245 | 44.94922 | 44.97286 | **44.94849** |
| 8 | 49.34802 | 49.34802 | 49.12767 | 49.35193 | 49.34860 | 49.37750 | **49.34802** |

**†  MPS** is not a mesh method — the size is 51 Fourier–Bessel basis functions
(204 boundary+interior collocation points), not DOFs; it reaches the highest
accuracy with by far the smallest system. **‡  Wolfram** `NDEigensystem` builds
its own internal adaptive mesh, so its DOF count is not exposed. The remaining
DOFs are read from each CSV header: DST/FD = interior masked grid points,
FEM eig = constrained DOFs (`size(K,1)`), FEM solvepdeeig = mesh nodes.

## Sources

| Column | File | Resolution / DOFs |
|--------|------|-------------------|
| DST full | `results/dst/L_shaped_full-eigenvalues.csv` | M = 49 (dense `eig`), 1776 DOFs |
| DST partial | `results/dst/L_shaped_partial-eigenvalues.csv` | M = 299 (`eigs`, leading 10), 66901 DOFs |
| FD | `results/fd/L_shaped_fd-eigenvalues.csv` | M = 49 (5-point stencil, dense `eig`), 1776 DOFs |
| FEM eig | `results/fem/L_shaped_fem_eig-eigenvalues.csv` | Hmax = 0.08, dense `eig(K,M)`, 2051 DOFs |
| FEM solvepdeeig | `results/fem/L_shaped_fem_solvepdeeig-eigenvalues.csv` | Hmax = 0.05, `solvepdeeig`, 5653 nodes |
| Wolfram | `results/wolfram/L-shaped-eigenvalues.csv` | `NDEigensystem`, default adaptive mesh |
| MPS | `results/mps/L_shaped_eigenvalues_MPS.csv` | method of particular solutions, 51 basis fns |

## Reading the results

MPS is the high-accuracy reference (it matches Betcke & Trefethen to 7-8
digits). Relative to it:

- **lambda_1** is the hardest mode (re-entrant-corner singularity; the true
  value is 9.6397238). MPS nails it; DST partial (9.64298) and FEM solvepdeeig
  (9.64532) are the closest of the rest; FD (9.66133) and DST full (9.67525,
  coarse M = 49) overshoot the most. All methods converge from above.
- **lambda_3 = 2*pi^2 ~ 19.73921** and **lambda_8 = 5*pi^2 ~ 49.34802** are
  captured to machine precision by DST (the uniform grid is aligned so these
  analytic modes are exact) and by MPS; FD cannot reproduce them exactly.
- **FD and DST full use the identical 1776 DOFs** (both M = 49 on the same
  masked L-shaped grid) — an apples-to-apples comparison. At equal size DST
  exactly captures lambda_3, lambda_8 that FD misses; for the corner-singular
  lambda_1 the two are comparable.
- DST partial (M = 299) is the most accurate grid method, at ~38x the DOFs of
  DST full. FD is uniformly lowest (2nd-order stencil); Wolfram sits slightly
  high but still within ~0.2%.
- Chebfun is absent: its tensor-product (`diffmat` + Kronecker) construction
  handles rectangles only, so it has no L-shaped entry.

## Cost vs. accuracy

- **MPS is the efficiency winner**: highest accuracy in 4.7 s with just 51
  basis functions — orders of magnitude fewer DOFs than the grid methods.
- **DST partial buys its accuracy dearly**: 170 s (66 901 DOFs), by far the
  most expensive — the cost of the `eigs` solve on the fine M = 299 operator.
- **FD and DST full are the cheapest grid solves** (~1 s at 1776 DOFs); FEM
  sits between (~2 s). So at equal size DST full matches FD's cost while being
  more accurate on the smooth modes.
- Timings are wall-clock for the per-domain solve (assembly + eigensolver),
  measured with `tic`/`toc`; they exclude MATLAB start-up and I/O and vary
  run-to-run (the first DST call also absorbs the parallel-pool start-up).
