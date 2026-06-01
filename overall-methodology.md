# Overall methodology: launcher → catalog → runner → driver → writer

Each eigenvalue-computation technique (DST, FEM, Chebfun, MPS, Wolfram) is
organised into the same five-part pipeline. Source code lives under `src/` and
`experiments/`; generated CSVs land in `results/<technique>/` and are
gitignored.

## The five roles

1. **Catalog** (`experiments/domain_catalog_<technique>.m`)
   The single source of truth for which domains a technique computes. A
   function returning a struct array, one row per domain, with an optional
   `name` filter that errors on an unknown domain. Each row holds the domain
   `name` plus whatever that technique needs to build the problem (DST: `box`,
   `phi`, `M_full`, `M_partial`; FEM: `decsg` geometry `gd/ns/sf`, `Hmax_eig`,
   `Hmax_solvepdeeig`; Chebfun: `box`). Adding a domain = one new row, no new
   files.

2. **Runner** (`src/<technique>/..._spectrum.m`)
   A pure function `[evals, info] = runner(entry, ...)` taking one catalog
   entry (plus any method/resolution argument) and returning the eigenvalues
   (positive, of -Laplacian, ascending) and an `info` struct of metadata
   (resolution, dofs, ...). No I/O, no path handling. Where a technique has
   variants, the runner takes a `method`/`N` argument and branches internally
   (FEM `"eig"` vs `"solvepdeeig"`; DST `"full"` vs `"partial"`).

3. **Driver** (`experiments/compute_spectrum_<technique>.m`)
   Orchestrates a run: resolves the project root, runs `startup.m` to set the
   path, adds `experiments/` (for the catalog), creates `results/<technique>/`,
   loops over `catalog(name)` x any sweep (modes, N, ...), calls the runner,
   and calls the writer for each result. Exposes a `(name, ...)` signature so a
   single domain or subset can be computed interactively.

4. **Writer** (`src/<technique>/write_*_csv.m`, or inline in the driver)
   Writes one CSV per result with a leading block of `#`-prefixed metadata
   lines (domain, timestamp, resolution, dofs, method) followed by an
   `n,lambda_n` table via `writetable(..., 'WriteMode', 'append')`. Keeps the
   output self-describing and uniform across techniques.

5. **Launcher** (`experiments/run_<technique>.sh`)
   A thin headless entry point: checks the required binary is on `PATH`
   (`matlab` / `wolframscript`), then invokes the driver
   (`matlab -batch "addpath(...); compute_spectrum_<technique>()"`, or
   `wolframscript -file ...`). One launcher per workflow.

## Data flow

```
run_<technique>.sh                   (launcher, shell)
  --> compute_spectrum_<technique>   (driver: sets path + loops)
        |-- domain_catalog_<technique>(name)      (catalog -> entries)
        |-- <technique>_laplace_spectrum(entry,.) (runner -> evals, info)
        \-- write_<technique>_csv(file, evals, .) (writer -> results/<t>/*.csv)
```

## Conventions

- **Naming:** drivers `compute_spectrum_{dst,mps,wolfram,cheb,fem_eig,fem_solvepdeeig}`;
  catalogs `domain_catalog_{dst,cheb,fem}`; launchers `run_<same>.sh`.
- **Output:** `results/<technique>/<domain>_<variant>-eigenvalues.csv`, where
  variant is the mode (`full`/`partial`), Chebyshev order (`N`), FEM workflow,
  etc.; CSVs are gitignored and regenerated via the launchers.
- **Two-resolution pattern:** techniques with a cheap and an expensive solve
  store two resolutions per catalog row (DST `M_full`/`M_partial`, FEM
  `Hmax_eig`/`Hmax_solvepdeeig`).
- **Sign/order:** runners return positive eigenvalues of -Laplacian, ascending.
- **Path setup:** drivers `run` the root `startup.m`; no per-script `addpath`
  of `src/`.

## Per-technique instances

| Technique | Catalog | Runner | Driver(s) | Launcher(s) |
|-----------|---------|--------|-----------|-------------|
| DST | `domain_catalog_dst` | `dst_laplace_spectrum` | `compute_spectrum_dst` | `run_dst.sh` |
| Chebfun | `domain_catalog_cheb` | `chebfun_laplace_spectrum` | `compute_spectrum_cheb` | `run_cheb.sh` |
| FEM | `domain_catalog_fem` | `fem_laplace_spectrum` | `compute_spectrum_fem_eig`, `compute_spectrum_fem_solvepdeeig` | `run_fem_eig.sh`, `run_fem_solvepdeeig.sh` |
| MPS | — (single domain) | `Ldrum_modified` (script) | `compute_spectrum_mps` | `run_mps.sh` |
| Wolfram | `wolframRegions` (in core `.wls`) | `reportEigenvalues` | `compute_spectrum_wolfram.wls` | `run_wolfram.sh` |

Two methods deviate slightly. **MPS** is single-domain: there is no catalog,
and the driver runs the Betcke & Trefethen script `Ldrum_modified` directly.
**Wolfram** packs the catalog (`wolframRegions` association) and the runner
(`reportEigenvalues`) into one core library
`src/wolfram/eigenvalues-various-domains.wls` that the driver `Get`s.
