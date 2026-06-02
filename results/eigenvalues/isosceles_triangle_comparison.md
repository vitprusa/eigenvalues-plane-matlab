# isosceles triangle domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
isosceles triangle domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Analytic (exact) |
|---|---|---|---|---|---|---|
| **DOFs** | 1225 | 19900 | 1225 | 2115 | 9326 § | — |
| **Time (s)** | 0.38 | 19.01 | 0.06 | 0.53 | 0.80 | — |
| 1 | 4.99996 | 5.00000 | 4.99463 | 5.00001 | 5.00000 | **5.00000** |
| 2 | 9.99974 | 10.00000 | 9.97410 | 10.00008 | 10.00000 | **10.00000** |
| 3 | 12.99990 | 13.00000 | 12.96936 | 13.00018 | 13.00001 | **13.00000** |
| 4 | 16.99908 | 16.99999 | 16.91890 | 17.00040 | 17.00002 | **17.00000** |
| 5 | 19.99941 | 19.99999 | 19.91416 | 20.00066 | 20.00004 | **20.00000** |
| 6 | 24.99980 | 25.00000 | 24.89363 | 25.00128 | 25.00007 | **25.00000** |
| 7 | 25.99763 | 25.99996 | 25.80268 | 26.00142 | 26.00008 | **26.00000** |
| 8 | 28.99818 | 28.99997 | 28.79793 | 29.00198 | 29.00012 | **29.00000** |

**§  FEM solvepdeeig** reports the number of mesh nodes, not the (constrained) degrees of freedom — unlike FEM eig, which records `size(K,1)`.

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/isosceles_triangle_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/isosceles_triangle_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/isosceles_triangle_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/isosceles_triangle_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/isosceles_triangle_fem_solvepdeeig-eigenvalues.csv` |
| Analytic (exact) | closed-form spectrum |
