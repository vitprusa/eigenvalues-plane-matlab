# small rectangle domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
small rectangle domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Cheb (N=40) | Analytic (exact) |
|---|---|---|---|---|---|---|---|
| **DOFs** | 253 | 4851 | 253 | 2143 | 6081 | 1444 | — |
| **Time (s)** | 0.02 | 1.91 | 0.01 | 0.66 | 0.35 | 0.46 | — |
| 1 | 20.00000 | 20.00000 | 19.90311 | 20.00004 | 20.00001 | 20.00000 | **20.00000** |
| 2 | 32.00000 | 32.00000 | 31.81765 | 32.00016 | 32.00002 | 32.00000 | **32.00000** |
| 3 | 52.00000 | 52.00000 | 51.44856 | 52.00072 | 52.00009 | 52.00000 | **52.00000** |
| 4 | 68.00000 | 68.00000 | 66.54542 | 68.00154 | 68.00020 | 68.00000 | **68.00000** |
| 5 | 80.00000 | 80.00000 | 78.45996 | 80.00244 | 80.00032 | 80.00000 | **80.00000** |
| 6 | 80.00000 | 80.00000 | 78.45996 | 80.00265 | 80.00033 | 80.00000 | **80.00000** |
| 7 | 100.00000 | 100.00000 | 98.09087 | 100.00480 | 100.00061 | 100.00000 | **100.00000** |
| 8 | 116.00000 | 116.00000 | 112.38967 | 116.00797 | 116.00101 | 116.00000 | **116.00000** |

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/small_rectangle_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/small_rectangle_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/small_rectangle_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/small_rectangle_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/small_rectangle_fem_solvepdeeig-eigenvalues.csv` |
| Cheb (N=40) | `results/eigenvalues/cheb/small_rectangle_N40-eigenvalues.csv` |
| Analytic (exact) | closed-form spectrum |
