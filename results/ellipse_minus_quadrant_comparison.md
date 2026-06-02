# ellipse-minus-quadrant domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
ellipse-minus-quadrant domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Wolfram |
|---|---|---|---|---|---|---|
| **DOFs** | 769 | 11866 | 769 | 2143 | 6435 § | — ‡ |
| **Time (s)** | 0.12 | 6.66 | 0.03 | 0.61 | 0.46 | — |
| 1 | 5.40554 | 5.76317 | 5.38900 | 5.87254 | 5.87067 | 5.87631 |
| 2 | 10.19953 | 11.28754 | 10.14887 | 11.53448 | 11.53023 | 11.54380 |
| 3 | 13.91571 | 14.88577 | 13.83892 | 15.14225 | 15.14113 | 15.14708 |
| 4 | 14.65422 | 15.64818 | 14.56923 | 15.92427 | 15.92400 | 15.92704 |
| 5 | 19.24541 | 20.77876 | 19.11228 | 21.12184 | 21.11966 | 21.13200 |
| 6 | 25.01768 | 26.66692 | 24.80676 | 27.08210 | 27.08082 | 27.09676 |
| 7 | 25.70583 | 27.07592 | 25.43114 | 27.50782 | 27.50173 | 27.53449 |
| 8 | 27.64967 | 29.38017 | 27.33256 | 29.92590 | 29.91975 | 29.95455 |

**‡  Wolfram** `NDEigensystem` builds its own internal adaptive mesh, so its DOF count is not exposed. **§  FEM solvepdeeig** reports the number of mesh nodes, not the (constrained) degrees of freedom — unlike FEM eig, which records `size(K,1)`.

## Sources

| Column | File |
|--------|------|
| DST full | `results/dst/ellipse_minus_quadrant_full-eigenvalues.csv` |
| DST partial | `results/dst/ellipse_minus_quadrant_partial-eigenvalues.csv` |
| FD | `results/fd/ellipse_minus_quadrant_fd-eigenvalues.csv` |
| FEM eig | `results/fem/ellipse_minus_quadrant_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/fem/ellipse_minus_quadrant_fem_solvepdeeig-eigenvalues.csv` |
| Wolfram | `results/wolfram/ellipse-minus-quadrant-eigenvalues.csv` |
