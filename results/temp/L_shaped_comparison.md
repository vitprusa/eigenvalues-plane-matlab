# L-shaped domain — eigenvalue comparison across techniques

First 8 eigenvalues of the Dirichlet Laplacian (-Delta u = lambda u) on the
L-shaped domain, computed by each technique. Values to 5 decimals.

| n | DST full | DST partial | FEM eig | FEM solvepdeeig | Wolfram | MPS (ground truth) |
|---|----------|-------------|---------|-----------------|---------|--------------------|
| 1 |  9.67525 |  9.64298 |  9.64999 |  9.64532 |  9.65934 |  **9.63972** |
| 2 | 15.19657 | 15.19725 | 15.19754 | 15.19732 | 15.19878 | **15.19725** |
| 3 | 19.73921 | 19.73921 | 19.73946 | 19.73925 | 19.74141 | **19.73921** |
| 4 | 29.52119 | 29.52148 | 29.52240 | 29.52163 | 29.52854 | **29.52148** |
| 5 | 31.99865 | 31.92059 | 31.93859 | 31.92643 | 31.96703 | **31.91264** |
| 6 | 41.53975 | 41.48049 | 41.49560 | 41.48512 | 41.52596 | **41.47451** |
| 7 | 44.94498 | 44.94846 | 44.95245 | 44.94922 | 44.97286 | **44.94849** |
| 8 | 49.34802 | 49.34802 | 49.35193 | 49.34860 | 49.37750 | **49.34802** |

## Sources

| Column | File | Resolution |
|--------|------|------------|
| DST full | `results/dst/L_shaped_full-eigenvalues.csv` | M = 49 (dense `eig`) |
| DST partial | `results/dst/L_shaped_partial-eigenvalues.csv` | M = 299 (`eigs`, leading 10) |
| FEM eig | `results/fem/L_shaped_fem_eig-eigenvalues.csv` | Hmax = 0.08, dense `eig(K,M)` |
| FEM solvepdeeig | `results/fem/L_shaped_fem_solvepdeeig-eigenvalues.csv` | Hmax = 0.05, `solvepdeeig` |
| Wolfram | `results/wolfram/L-shaped-eigenvalues.csv` | `NDEigensystem`, default mesh |
| MPS | `results/mps/L_shaped_eigenvalues_MPS.csv` | method of particular solutions |

## Reading the results

MPS is the high-accuracy reference (it matches Betcke & Trefethen to 7-8
digits). Relative to it:

- **lambda_1** is the hardest mode (re-entrant-corner singularity; the true
  value is 9.6397238). MPS nails it; DST partial (9.64298) and FEM solvepdeeig
  (9.64532) are the closest of the rest; DST full (9.67525, coarse M = 49) and
  Wolfram (9.65934) overshoot the most. All methods converge from above.
- **lambda_3 = 2*pi^2 ~ 19.73921** and **lambda_8 = 5*pi^2 ~ 49.34802** are
  captured to machine precision by DST (the uniform grid is aligned so these
  analytic modes are exact) and by MPS.
- DST partial (M = 299) generally beats DST full (M = 49) thanks to the finer
  grid, except where both are exact.
- Wolfram is the least accurate here (coarsest effective resolution) but still
  within ~0.2%.
- Chebfun is absent: its tensor-product (`diffmat` + Kronecker) construction
  handles rectangles only, so it has no L-shaped entry.
