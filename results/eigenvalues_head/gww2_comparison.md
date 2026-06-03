# GWW2 isospectral drum domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
GWW2 isospectral drum domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig |
|---|---|---|---|---|---|
| **DOFs** | 936 | 15411 | 936 | 2672 | 6214 |
| **Time (s)** | 0.20 | 8.19 | 0.04 | 1.39 | 1.35 |
| 1 | 2.55187 | 2.54022 | 2.54406 | 2.54031 | 2.53927 |
| 2 | 3.66867 | 3.65774 | 3.65668 | 3.65789 | 3.65684 |
| 3 | 5.20264 | 5.18001 | 5.17672 | 5.18021 | 5.17817 |
| 4 | 6.54099 | 6.53840 | 6.50365 | 6.53872 | 6.53814 |
| 5 | 7.27317 | 7.25227 | 7.22841 | 7.25258 | 7.25057 |
| 6 | 9.21919 | 9.21119 | 9.14783 | 9.21171 | 9.21052 |
| 7 | 10.60897 | 10.59947 | 10.49652 | 10.60025 | 10.59861 |
| 8 | 11.54228 | 11.54200 | 11.42287 | 11.54302 | 11.54199 |

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/gww2_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/gww2_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/gww2_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/gww2_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/gww2_fem_solvepdeeig-eigenvalues.csv` |
