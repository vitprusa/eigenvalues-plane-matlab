# Batched DST Laplace matrix — test report

- Tests run: 4 (0 failed)

## Summary

| Test | Status | Duration |
|------|--------|----------|
| `check_batched_equiv` | ✅ PASS | 27.8s |
| `check_batched_timing` | ✅ PASS | 29.0s |
| `check_op_batched_equiv` | ✅ PASS | 26.4s |
| `check_op_batched_timing` | ✅ PASS | 39.1s |

## Details

### `check_batched_equiv`

- Status: **✅ PASS**
- Duration: 27.8s

```
rectangle    dofs=  400  ||L_loop - L_batch||_F = 0.000e+00  (rel 0.000e+00)
L-shaped     dofs=  675  ||L_loop - L_batch||_F = 2.739e-12  (rel 5.800e-17)
H            dofs=  448  ||L_loop - L_batch||_F = 0.000e+00  (rel 0.000e+00)
```

### `check_batched_timing`

- Status: **✅ PASS**
- Duration: 29.0s

```
dofs=2700  loop=20.15s  batched=0.173s  speedup=116x  err=0.00e+00
```

### `check_op_batched_equiv`

- Status: **✅ PASS**
- Duration: 26.4s

```
rectangle    dofs=  400  ||Lop_ref - Lop_batched||_F = 5.689e-12  (rel 3.033e-16)
L-shaped     dofs=  675  ||Lop_ref - Lop_batched||_F = 4.308e-11  (rel 3.269e-16)
H            dofs=  448  ||Lop_ref - Lop_batched||_F = 9.019e-12  (rel 3.109e-16)
```

### `check_op_batched_timing`

- Status: **✅ PASS**
- Duration: 39.1s

```
M=60   dofs=  2700  ref=6.778e-02s  dense=2.970e-03s (23x)  sparse=2.226e-04s (305x)  err=1.67e-10
M=120  dofs= 10800  ref=7.465e-02s  dense=5.565e-02s (1x)  sparse=5.028e-03s (15x)  err=1.56e-09
```

