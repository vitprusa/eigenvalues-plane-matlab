# rectangle domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
rectangle domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Cheb (N=40) | Wolfram | Analytic (exact) |
|---|---|---|---|---|---|---|---|---|
| **DOFs** | 1176 | 19701 | 1176 | 3895 | 18351 | 1444 | — ‡ | — |
| **Time (s)** | 0.31 | 14.21 | 0.05 | 5.62 | 3.56 | 0.38 | — | — |
| 1 | 1.25000 | 1.25000 | 1.24860 | 1.25000 | 1.25000 | 1.25000 | 1.25000 | **1.25000** |
| 2 | 2.00000 | 2.00000 | 1.99737 | 2.00000 | 2.00000 | 2.00000 | 2.00001 | **2.00000** |
| 3 | 3.25000 | 3.25000 | 3.24203 | 3.25001 | 3.25000 | 3.25000 | 3.25004 | **3.25000** |
| 4 | 4.25000 | 4.25000 | 4.22891 | 4.25003 | 4.25000 | 4.25000 | 4.25017 | **4.25000** |
| 5 | 5.00000 | 5.00000 | 4.97767 | 5.00005 | 5.00000 | 5.00000 | 5.00017 | **5.00000** |
| 6 | 5.00000 | 5.00000 | 4.97767 | 5.00005 | 5.00000 | 5.00000 | 5.00020 | **5.00000** |
| 7 | 6.25000 | 6.25000 | 6.22234 | 6.25009 | 6.25000 | 6.25000 | 6.25021 | **6.25000** |
| 8 | 7.25000 | 7.25000 | 7.19745 | 7.25015 | 7.25001 | 7.25000 | 7.25074 | **7.25000** |

**‡  Wolfram** `NDEigensystem` builds its own internal adaptive mesh, so its DOF count is not exposed.

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/rectangle_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/rectangle_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/rectangle_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/rectangle_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/rectangle_fem_solvepdeeig-eigenvalues.csv` |
| Cheb (N=40) | `results/eigenvalues/cheb/rectangle_N40-eigenvalues.csv` |
| Wolfram | `results/eigenvalues/wolfram/rectangle-eigenvalues.csv` |
| Analytic (exact) | closed-form spectrum |
