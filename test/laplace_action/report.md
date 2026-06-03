# DST Laplace action — test report

- Generated: 2026-06-03 14:08:56 CEST
- Host: `mroz`
- Project root: `/home/vitek/Documents/Science/Work/Articles_in_progress/Discrete versus continuous/eigenvalues-plane-matlab`
- MATLAB: `matlab`
- Tests run: 1 (0 failed)

## Summary

| Test | Status | Duration |
|------|--------|----------|
| `check_laplace_action_rectangle` | ✅ PASS | 25.6s |

## Details

### `check_laplace_action_rectangle`

- Status: **✅ PASS**
- Duration: 25.6s

```
rectangle [0,pi]^2  M=20  dofs=400  mode (mm,nn)=(5,1)
||laplace_u_exact - laplace_u_num||_2 = 1.025e-12  (rel 3.753e-15)
PASS (rel residual 3.753e-15 <= tol 1.000e-10)
```

