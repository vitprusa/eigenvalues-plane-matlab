# Rectangle DST eigenvalue ordering

Full DST-Laplacian spectrum of the rectangle $[0, 2\pi] \times [0, \pi]$ at resolution $M = 9$ (dofs $= 36$). Each cell holds the **position** of the analytic eigenvalue $\lambda_{m,n} = (m\pi/L_x)^2 + (n\pi/L_y)^2 = m^2/4 + n^2$ in the magnitude-ordered sequence of all analytic eigenvalues (1 = smallest).

The grid resolves $M_x = 9$ modes in $x$ and $M_y = 4$ modes in $y$. Cells matched to a computed eigenvalue (to round-off) are shown in **bold**; the unresolved border ($m > M_x$ or $n > M_y$) is in normal weight. Rows are indexed by $m$, columns by $n$. The bold cells are not positions $1 \dots M_x M_y$: some unresolved border modes outrank resolved ones.

| $m \backslash n$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | **1** | **4** | **10** | **19** | 33 | 47 | 66 | 87 | 113 |
| 2 | **2** | **5** | **11** | **21** | 34 | 49 | 67 | 88 | 116 |
| 3 | **3** | **7** | **14** | **24** | 36 | 52 | 70 | 93 | 118 |
| 4 | **6** | **9** | **15** | **25** | 38 | 54 | 74 | 96 | 120 |
| 5 | **8** | **13** | **18** | **29** | 41 | 58 | 76 | 99 | 124 |
| 6 | **12** | **16** | **23** | **31** | 44 | 60 | 80 | 102 | 128 |
| 7 | **17** | **20** | **27** | **37** | 51 | 65 | 86 | 107 | 133 |
| 8 | **22** | **26** | **32** | **43** | 56 | 72 | 89 | 111 | 136 |
| 9 | **28** | **30** | **40** | **48** | 62 | 78 | 98 | 119 | 146 |
| 10 | 35 | 39 | 45 | 57 | 68 | 84 | 105 | 126 | 151 |
| 11 | 42 | 46 | 53 | 63 | 77 | 94 | 110 | 134 | 159 |
| 12 | 50 | 55 | 61 | 73 | 85 | 100 | 121 | 141 | 168 |
| 13 | 59 | 64 | 71 | 82 | 95 | 109 | 130 | 154 | 178 |
| 14 | 69 | 75 | 81 | 90 | 106 | 122 | 139 | 162 | 188 |

Matched 36 of 126 lattice cells to computed eigenvalues. Computed eigenvalues: 36 total, 0 unmatched.
