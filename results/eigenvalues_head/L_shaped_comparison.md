# L-shaped domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
L-shaped domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Wolfram | MPS (ground truth) |
|---|---|---|---|---|---|---|---|
| **DOFs** | 1776 | 66901 | 1776 | 2051 | 5333 | — ‡ | 51 † |
| **Time (s)** | 1.01 | 170.04 | 0.20 | 0.63 | 0.48 | — | 4.72 |
| 1 | 9.67525 | 9.64298 | 9.66133 | 9.64999 | 9.64532 | 9.65934 | **9.63972** |
| 2 | 15.19657 | 15.19725 | 15.17674 | 15.19754 | 15.19732 | 15.19878 | **15.19725** |
| 3 | 19.73921 | 19.73921 | 19.71325 | 19.73946 | 19.73925 | 19.74141 | **19.73921** |
| 4 | 29.52119 | 29.52148 | 29.44785 | 29.52240 | 29.52163 | 29.52854 | **29.52148** |
| 5 | 31.99865 | 31.92059 | 31.89650 | 31.93859 | 31.92643 | 31.96703 | **31.91264** |
| 6 | 41.53975 | 41.48049 | 41.34627 | 41.49560 | 41.48512 | 41.52596 | **41.47451** |
| 7 | 44.94498 | 44.94846 | 44.72497 | 44.95245 | 44.94922 | 44.97286 | **44.94849** |
| 8 | 49.34802 | 49.34802 | 49.12767 | 49.35193 | 49.34860 | 49.37750 | **49.34802** |

**†  MPS** is not a mesh method — the size is the number of Fourier–Bessel basis functions, not DOFs. **‡  Wolfram** `NDEigensystem` builds its own internal adaptive mesh, so its DOF count is not exposed.

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/L_shaped_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/L_shaped_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/L_shaped_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/L_shaped_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/L_shaped_fem_solvepdeeig-eigenvalues.csv` |
| Wolfram | `results/eigenvalues/wolfram/L-shaped-eigenvalues.csv` |
| MPS (ground truth) | `results/eigenvalues/mps/L_shaped_eigenvalues_MPS.csv` |
