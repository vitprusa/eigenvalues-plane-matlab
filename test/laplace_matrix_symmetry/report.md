# DST Laplace matrix symmetry — test report

- Tests run: 1 (0 failed)

## Summary

| Test | Status | Duration |
|------|--------|----------|
| `check_laplace_matrix_symmetry` | ✅ PASS | 5.8s |

## Details

### `check_laplace_matrix_symmetry`

- Status: **✅ PASS**
- Duration: 5.8s

```
  square                   M =  30  cached
  rectangle                M =  29  cached
  isosceles_triangle       M =  30  cached
  small_rectangle          M =  31  cached
  L_shaped                 M =  29  cached
  ellipse_minus_quadrant   M =  31  cached
  H_shaped                 M =  29  cached
  gww1                     M =  29  cached
  gww2                     M =  29  cached
  square                   M =  60  cached
  rectangle                M =  59  cached
  isosceles_triangle       M =  60  cached
  small_rectangle          M =  59  cached
  L_shaped                 M =  59  cached
  ellipse_minus_quadrant   M =  59  cached
  H_shaped                 M =  59  cached
  gww1                     M =  59  cached
  gww2                     M =  59  cached
  square                   M =  90  cached
  rectangle                M =  89  cached
  isosceles_triangle       M =  90  cached
  small_rectangle          M =  91  cached
  L_shaped                 M =  89  cached
  ellipse_minus_quadrant   M =  91  cached
  H_shaped                 M =  89  cached
  gww1                     M =  89  cached
  gww2                     M =  89  cached
  square                   M = 120  cached
  rectangle                M = 119  cached
  isosceles_triangle       M = 120  cached
  small_rectangle          M = 119  cached
  L_shaped                 M = 119  cached
  ellipse_minus_quadrant   M = 119  cached
  H_shaped                 M = 119  cached
  gww1                     M = 119  cached
  gww2                     M = 119  cached

domain                      M   dofs    rel asym  max|L-L'| max|Im ev| t general     t sym max|dlambda|
square                     30    900   1.740e-16  6.040e-14  0.000e+00     0.22s     0.03s   2.342e-11
square                     60   3600   3.565e-16  3.578e-13  3.055e-13    11.17s     1.17s   2.583e-10
square                     90   8100   2.359e-16  6.821e-13  6.045e-13   126.76s     7.93s   6.203e-10
square                    120  14400   2.146e-16  1.364e-12  5.025e-12   716.68s    48.29s   4.784e-09
rectangle                  29    406   1.542e-16  1.515e-14  0.000e+00     0.03s     0.00s   4.036e-12
rectangle                  59   1711   1.640e-16  9.104e-14  0.000e+00     0.90s     0.14s   4.025e-11
rectangle                  89   3916   1.908e-16  1.705e-13  0.000e+00    14.50s     1.14s   1.828e-10
rectangle                 119   7021   1.718e-16  3.642e-13  8.897e-13    90.14s     5.50s   3.301e-10
isosceles_triangle         30    435   1.532e-16  1.030e-13  0.000e+00     0.05s     0.01s   2.342e-11
isosceles_triangle         60   1770   2.162e-16  6.679e-13  0.000e+00     1.08s     0.16s   2.910e-10
isosceles_triangle         90   4005   2.493e-16  1.592e-12  0.000e+00    16.46s     1.21s   8.567e-10
isosceles_triangle        120   7140   2.677e-16  2.728e-12  0.000e+00    95.06s     5.76s   2.379e-09
small_rectangle            31    105   6.631e-17  2.842e-14  1.507e-14     0.00s     0.00s   5.457e-12
small_rectangle            59    406   1.542e-16  2.425e-13  0.000e+00     0.02s     0.01s   6.457e-11
small_rectangle            91    990   1.604e-16  6.821e-13  0.000e+00     0.28s     0.04s   2.201e-10
small_rectangle           119   1711   1.640e-16  1.457e-12  0.000e+00     0.90s     0.15s   6.439e-10
L_shaped                   29    616   1.594e-16  1.705e-13  0.000e+00     0.11s     0.01s   3.638e-11
L_shaped                   59   2581   1.648e-16  8.988e-13  1.938e-13     3.78s     0.38s   6.767e-10
L_shaped                   89   5896   1.958e-16  1.819e-12  1.867e-12    47.75s     3.41s   3.572e-09
L_shaped                  119  10561   1.745e-16  3.595e-12  3.522e-12   273.36s    17.88s   9.866e-09
ellipse_minus_quadrant     31    280   1.196e-16  7.105e-14  0.000e+00     0.02s     0.00s   1.319e-11
ellipse_minus_quadrant     59   1024   1.636e-16  2.383e-13  0.000e+00     0.43s     0.04s   7.731e-11
ellipse_minus_quadrant     91   2443   1.882e-16  9.375e-13  0.000e+00     3.12s     0.35s   5.421e-10
ellipse_minus_quadrant    119   4177   2.017e-16  1.245e-12  0.000e+00    18.50s     1.35s   1.177e-09
H_shaped                   29    621   1.446e-16  8.527e-14  0.000e+00     0.11s     0.01s   2.319e-11
H_shaped                   59   2641   1.513e-16  3.997e-13  0.000e+00     4.13s     0.41s   3.274e-10
H_shaped                   89   6061   1.875e-16  7.105e-13  1.364e-12    58.89s     4.06s   1.772e-09
H_shaped                  119  10881   1.685e-16  1.251e-12  0.000e+00   301.35s    20.45s   5.082e-09
gww1                       29    306   1.384e-16  2.842e-14  0.000e+00     0.02s     0.00s   5.059e-12
gww1                       59   1311   1.940e-16  1.772e-13  0.000e+00     0.37s     0.07s   3.729e-11
gww1                       89   3016   2.243e-16  3.988e-13  0.000e+00     6.77s     0.58s   2.569e-10
gww1                      119   5421   2.386e-16  7.090e-13  0.000e+00    36.77s     2.71s   8.276e-10
gww2                       29    306   1.236e-16  2.646e-14  0.000e+00     0.02s     0.00s   5.286e-12
gww2                       59   1311   1.618e-16  1.772e-13  0.000e+00     0.41s     0.07s   6.207e-11
gww2                       89   3016   1.963e-16  3.988e-13  0.000e+00     7.00s     0.60s   2.237e-10
gww2                      119   5421   2.119e-16  7.090e-13  0.000e+00    37.24s     2.68s   7.512e-10

Over the whole table the general EIG costs 1874 s and the symmetric one 127 s,
a factor of 14.8.

Every matrix is symmetric to rounding, its pattern exactly so, and no
eigenvalue moves: the largest relative difference between the two spectra
is 8.14e-11, on square (M = 120).
```

