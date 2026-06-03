# H-shaped domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
H-shaped domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Wolfram |
|---|---|---|---|---|---|---|
| **DOFs** | 1888 | 30888 | 1888 | 4843 | 12549 | — ‡ |
| **Time (s)** | 1.15 | 24.00 | 0.28 | 10.87 | 1.29 | — |
| 1 | 7.81056 | 7.74551 | 7.78327 | 7.74680 | 7.74032 | 7.77338 |
| 2 | 8.62413 | 8.56334 | 8.59776 | 8.56458 | 8.55848 | 8.58911 |
| 3 | 13.96576 | 13.93392 | 13.92665 | 13.93484 | 13.93131 | 13.95073 |
| 4 | 13.96935 | 13.93782 | 13.93050 | 13.93874 | 13.93524 | 13.95465 |
| 5 | 14.30475 | 14.30557 | 14.26203 | 14.30607 | 14.30558 | 14.31257 |
| 6 | 17.71624 | 17.70854 | 17.66211 | 17.70911 | 17.70791 | 17.72061 |
| 7 | 19.73921 | 19.73921 | 19.68310 | 19.73946 | 19.73925 | 19.74799 |
| 8 | 24.84034 | 24.79734 | 24.72015 | 24.79879 | 24.79392 | 24.83486 |

**‡  Wolfram** `NDEigensystem` builds its own internal adaptive mesh, so its DOF count is not exposed.

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/H_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/H_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/H_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/H_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/H_fem_solvepdeeig-eigenvalues.csv` |
| Wolfram | `results/eigenvalues/wolfram/H-shaped-eigenvalues.csv` |
