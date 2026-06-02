# L-shaped domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
L-shaped domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not
timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Wolfram | MPS (ground truth) |
|---|----------|-------------|------|---------|-----------------|---------|--------------------|
| **DOFs** | 1776 | 66901 | 1776 | 2051 | 5653 § | — ‡ | 51 † |
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

**§  FEM solvepdeeig** reports the number of **mesh nodes**, not the
(constrained) degrees of freedom — unlike FEM eig, which records `size(K,1)`.

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
