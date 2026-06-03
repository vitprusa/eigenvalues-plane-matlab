# Method of manufactured solutions — test report

- Generated: 2026-06-03 14:22:48 CEST
- Host: `mroz`
- Project root: `/home/vitek/Documents/Science/Work/Articles_in_progress/Discrete versus continuous/eigenvalues-plane-matlab`
- MATLAB: `matlab`
- Tests run: 1 (0 failed)

## Summary

| Test | Status | Duration |
|------|--------|----------|
| `check_bvp_L_shaped` | ✅ PASS | 65.6s |

## Details

### `check_bvp_L_shaped`

- Status: **✅ PASS**
- Duration: 65.6s

```
L-shaped domain, manufactured solution u = x(x-pi)(x-b) y(y-pi)(y-e)
    M    dofs      max_err      rel_err    order
   23     313    3.666e-02    8.372e-03      NaN
   47    1345    9.848e-03    2.249e-03     1.90
   95    5569    2.556e-03    5.832e-04     1.95
  191   22657    6.510e-04    1.486e-04     1.97
PASS (error decreases under refinement; finest rel error 1.486e-04 <= tol 5.000e-04)
```

