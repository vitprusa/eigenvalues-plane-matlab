# ellipse-minus-quadrant domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
ellipse-minus-quadrant domain, computed by each technique. Values to 5 decimals; the **DOFs**
row gives the size of each discrete problem and **Time (s)** the measured
computation time (`tic`/`toc` around each per-domain solve; Wolfram is not timed).

| n | DST full | DST partial | FD | FEM eig | FEM solvepdeeig | Wolfram |
|---|---|---|---|---|---|---|
| **DOFs** | 649 | 11686 | 649 | 2143 | 6095 | — ‡ |
| **Time (s)** | 0.02 | 4.90 | 0.10 | 0.62 | 0.52 | — |
| 1 | 5.66750 | 5.80665 | 5.64838 | 5.87254 | 5.87067 | 5.87631 |
| 2 | 11.20836 | 11.41222 | 11.14025 | 11.53448 | 11.53023 | 11.54380 |
| 3 | 14.63226 | 14.99032 | 14.53816 | 15.14225 | 15.14113 | 15.14708 |
| 4 | 15.13343 | 15.70957 | 15.03352 | 15.92427 | 15.92400 | 15.92704 |
| 5 | 20.46896 | 20.90935 | 20.30016 | 21.12184 | 21.11966 | 21.13200 |
| 6 | 26.19190 | 26.79546 | 25.93766 | 27.08210 | 27.08082 | 27.09676 |
| 7 | 26.80042 | 27.28047 | 26.46824 | 27.50782 | 27.50173 | 27.53449 |
| 8 | 28.46907 | 29.51387 | 28.09660 | 29.92590 | 29.91975 | 29.95455 |

**‡  Wolfram** `NDEigensystem` builds its own internal adaptive mesh, so its DOF count is not exposed.

## Sources

| Column | File |
|--------|------|
| DST full | `results/eigenvalues/dst/ellipse_minus_quadrant_full-eigenvalues.csv` |
| DST partial | `results/eigenvalues/dst/ellipse_minus_quadrant_partial-eigenvalues.csv` |
| FD | `results/eigenvalues/fd/ellipse_minus_quadrant_fd-eigenvalues.csv` |
| FEM eig | `results/eigenvalues/fem/ellipse_minus_quadrant_fem_eig-eigenvalues.csv` |
| FEM solvepdeeig | `results/eigenvalues/fem/ellipse_minus_quadrant_fem_solvepdeeig-eigenvalues.csv` |
| Wolfram | `results/eigenvalues/wolfram/ellipse-minus-quadrant-eigenvalues.csv` |
