# square domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
square domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Cheb (N=40) | Analytic (exact) |
|---|---|---|---|---|---|---|---|
| **DOFs** | 2500 | 40000 | 2500 | 4301 | 17953 | 1444 | — |
| **Time (s)** | 5.19 | 66.59 | 1.15 | 9.17 | 4.31 | 0.39 | — |
| 1 | 2.00000 | 2.00000 | 1.99937 | 2.00000 | 2.00000 | 2.00000 | **2.00000** |
| 2 | 5.00000 | 5.00000 | 4.99463 | 5.00001 | 5.00000 | 5.00000 | **5.00000** |
| 3 | 5.00000 | 5.00000 | 4.99463 | 5.00001 | 5.00000 | 5.00000 | **5.00000** |
| 4 | 8.00000 | 8.00000 | 7.98989 | 8.00004 | 8.00000 | 8.00000 | **8.00000** |
| 5 | 10.00000 | 10.00000 | 9.97410 | 10.00008 | 10.00000 | 10.00000 | **10.00000** |
| 6 | 10.00000 | 10.00000 | 9.97410 | 10.00008 | 10.00000 | 10.00000 | **10.00000** |
| 7 | 13.00000 | 13.00000 | 12.96936 | 13.00018 | 13.00001 | 13.00000 | **13.00000** |
| 8 | 13.00000 | 13.00000 | 12.96936 | 13.00018 | 13.00001 | 13.00000 | **13.00000** |

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/square_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/square_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/square_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/square_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/square_fem_solvepdeeig-eigenvalues.csv` |
| Cheb (N=40) | `results/eigenvalues/cheb/square_N40-eigenvalues.csv` |
| Analytic (exact) | closed-form spectrum |
