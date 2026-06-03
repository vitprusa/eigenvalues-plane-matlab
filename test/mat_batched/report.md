# Batched DST Laplace matrix — test report

- Tests run: 4 (0 failed)

## Summary

| Test | Status | Duration |
|------|--------|----------|
| `check_batched_equiv` | ✅ PASS | 25.8s |
| `check_batched_timing` | ✅ PASS | 28.0s |
| `check_op_batched_equiv` | ✅ PASS | 23.5s |
| `check_op_batched_timing` | ✅ PASS | 31.9s |

## Details

### `check_batched_equiv`

- Status: **✅ PASS**
- Duration: 25.8s

```
rectangle    dofs=  400  ||L_loop - L_batch||_F = 0.000e+00  (rel 0.000e+00)
L-shaped     dofs=  675  ||L_loop - L_batch||_F = 2.739e-12  (rel 5.800e-17)
H            dofs=  448  ||L_loop - L_batch||_F = 0.000e+00  (rel 0.000e+00)
```

### `check_batched_timing`

- Status: **✅ PASS**
- Duration: 28.0s

```
dofs=2700  loop=21.19s  batched=0.169s  speedup=125x  err=0.00e+00
```

### `check_op_batched_equiv`

- Status: **✅ PASS**
- Duration: 23.5s

```
rectangle    dofs=  400  ||Lop_ref - Lop_batched||_F = 5.689e-12  (rel 3.033e-16)
L-shaped     dofs=  675  ||Lop_ref - Lop_batched||_F = 4.308e-11  (rel 3.269e-16)
H            dofs=  448  ||Lop_ref - Lop_batched||_F = 9.019e-12  (rel 3.109e-16)
```

### `check_op_batched_timing`

- Status: **✅ PASS**
- Duration: 31.9s

```
M=60   dofs=  2700  ref=6.235e-02s  dense=2.814e-03s (22x)  sparse=2.078e-04s (300x)  err=1.67e-10
M=120  dofs= 10800  ref=6.293e-02s  dense=4.658e-02s (1x)  sparse=1.741e-03s (36x)  err=1.56e-09
```

