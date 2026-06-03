# GWW1 isospectral drum domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
GWW1 isospectral drum domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig |
|---|---|---|---|---|---|
| **DOFs** | 936 | 15411 | 936 | 2677 | 6193 |
| **Time (s)** | 0.37 | 10.07 | 0.03 | 1.48 | 1.56 |
| 1 | 2.55191 | 2.54022 | 2.54406 | 2.54030 | 2.53927 |
| 2 | 3.66866 | 3.65774 | 3.65668 | 3.65789 | 3.65684 |
| 3 | 5.20255 | 5.18001 | 5.17672 | 5.18020 | 5.17818 |
| 4 | 6.54085 | 6.53840 | 6.50365 | 6.53872 | 6.53814 |
| 5 | 7.27323 | 7.25228 | 7.22841 | 7.25257 | 7.25057 |
| 6 | 9.21933 | 9.21120 | 9.14783 | 9.21171 | 9.21052 |
| 7 | 10.61037 | 10.59950 | 10.49652 | 10.60028 | 10.59861 |
| 8 | 11.54139 | 11.54199 | 11.42287 | 11.54298 | 11.54199 |

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/gww1_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/gww1_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/gww1_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/gww1_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/gww1_fem_solvepdeeig-eigenvalues.csv` |
