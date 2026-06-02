#!/usr/bin/env python3
"""Generate per-domain eigenvalue comparison tables.

For every domain in the DST catalog this writes
results/eigenvalues_head/<domain>_comparison.md (reading the per-method CSVs
from results/eigenvalues/),
comparing the first few Dirichlet-Laplacian eigenvalues from each technique that
handles the domain (DST full/partial, FD, FEM eig/solvepdeeig, Chebyshev for the
rectangular domains, Wolfram where it computed the same domain, MPS for the
L-shape), plus the analytic spectrum where a closed form is known. The DOFs and
Time rows are read from each CSV's metadata header.

Run after the experiment drivers have produced the results/eigenvalues/*/ CSVs:
    python3 experiments/make_eigenvalues_head_comparison_tables.py
"""

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results", "eigenvalues")        # per-method CSV inputs
OUT = os.path.join(ROOT, "results", "eigenvalues_head")   # comparison-table outputs

NEIG = 8          # eigenvalues shown per table
CHEB_N = 40       # Chebyshev order used for the Cheb column

DOMAINS = ["square", "rectangle", "isosceles_triangle", "small_rectangle",
           "L_shaped", "ellipse_minus_quadrant", "H", "gww1", "gww2"]

PRETTY = {
    "square": "square", "rectangle": "rectangle",
    "isosceles_triangle": "isosceles triangle", "small_rectangle": "small rectangle",
    "L_shaped": "L-shaped", "ellipse_minus_quadrant": "ellipse-minus-quadrant",
    "H": "H-shaped", "gww1": "GWW1 isospectral drum", "gww2": "GWW2 isospectral drum",
}

# Domains Wolfram computed at the SAME size (its triangle is a different size).
WOLFRAM_FILE = {
    "rectangle": "rectangle", "L_shaped": "L-shaped",
    "ellipse_minus_quadrant": "ellipse-minus-quadrant", "H": "H-shaped",
}

CHEB_DOMAINS = {"square", "rectangle", "small_rectangle"}


def read_csv(relpath):
    """Parse a results CSV -> dict(dofs, time, eigs) or None if missing."""
    path = os.path.join(RES, relpath)
    if not os.path.exists(path):
        return None
    dofs = time = None
    eigs = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                m = re.search(r"dofs\s*=\s*(\d+)", s)
                if m:
                    dofs = int(m.group(1))
                m = re.search(r"Computation time:\s*([0-9.]+)", s)
                if m:
                    time = float(m.group(1))
                m = re.search(r"accuracy N\s*=\s*(\d+)", s)  # MPS basis size
                if m and dofs is None:
                    dofs = int(m.group(1))
            elif s and not s.lower().startswith("n,") and not s.startswith('"n"'):
                parts = s.split(",")
                try:
                    eigs.append(float(parts[1]))
                except (IndexError, ValueError):
                    pass
    return {"dofs": dofs, "time": time, "eigs": eigs}


def analytic(domain):
    """Sorted analytic eigenvalues, or None when no closed form is used."""
    if domain == "square":                 # m^2 + n^2 on [0,pi]^2
        vals = [m * m + n * n for m in range(1, 40) for n in range(1, 40)]
    elif domain == "rectangle":            # m^2/4 + n^2 on [0,2pi]x[0,pi]
        vals = [m * m / 4 + n * n for m in range(1, 60) for n in range(1, 40)]
    elif domain == "small_rectangle":      # 4 m^2 + 16 n^2 on [0,pi/2]x[0,pi/4]
        vals = [4 * m * m + 16 * n * n for m in range(1, 40) for n in range(1, 30)]
    elif domain == "isosceles_triangle":   # m^2 + n^2, m > n >= 1 (legs pi)
        vals = [m * m + n * n for m in range(2, 60) for n in range(1, m)]
    else:
        return None
    return sorted(vals)


def column(header, relpath, kind):
    """Build a column dict from a CSV, or None if the CSV is missing.

    kind: 'grid' (dofs), 'sp' (FEM solvepdeeig -> mesh-node dofs, footnote §),
          'cheb', 'wolfram' (no dofs/time, footnote ‡), 'mps' (basis, footnote †).
    """
    data = read_csv(relpath)
    if data is None:
        return None
    return {"header": header, "path": relpath, "kind": kind, **data}


def fmt_dofs(col):
    if col["kind"] == "wolfram":
        return "— ‡"
    if col["kind"] == "mps":
        return f"{col['dofs']} †"
    if col["kind"] == "sp":
        return f"{col['dofs']} §"
    if col["kind"] == "analytic":
        return "—"
    return str(col["dofs"]) if col["dofs"] is not None else "—"


def fmt_time(col):
    if col["kind"] in ("wolfram", "analytic"):
        return "—"
    return f"{col['time']:.2f}" if col["time"] is not None else "—"


def build_columns(domain):
    cols = []
    for hdr, rel in [("DST full", f"dst/{domain}_full-eigenvalues.csv"),
                     ("DST partial", f"dst/{domain}_partial-eigenvalues.csv"),
                     ("FD", f"fd/{domain}_fd-eigenvalues.csv"),
                     ("FEM eig", f"fem/{domain}_fem_eig-eigenvalues.csv")]:
        c = column(hdr, rel, "grid")
        if c:
            cols.append(c)
    c = column("FEM solvepdeeig", f"fem/{domain}_fem_solvepdeeig-eigenvalues.csv", "sp")
    if c:
        cols.append(c)
    if domain in CHEB_DOMAINS:
        c = column(f"Cheb (N={CHEB_N})", f"cheb/{domain}_N{CHEB_N}-eigenvalues.csv", "cheb")
        if c:
            cols.append(c)
    if domain in WOLFRAM_FILE:
        c = column("Wolfram", f"wolfram/{WOLFRAM_FILE[domain]}-eigenvalues.csv", "wolfram")
        if c:
            cols.append(c)
    if domain == "L_shaped":
        c = column("MPS (ground truth)", "mps/L_shaped_eigenvalues_MPS.csv", "mps")
        if c:
            cols.append(c)
    an = analytic(domain)
    if an is not None:
        cols.append({"header": "Analytic (exact)", "path": None, "kind": "analytic",
                     "dofs": None, "time": None, "eigs": an})
    return cols


def render(domain, cols):
    truth = next((c for c in cols if c["kind"] in ("analytic", "mps")), None)
    lines = []
    lines.append(f"# {PRETTY[domain]} domain — eigenvalue comparison across techniques")
    lines.append("")
    lines.append(f"First {NEIG} eigenvalues of the Dirichlet Laplacian "
                 "(-Delta u = lambda u) on the")
    lines.append(f"{PRETTY[domain]} domain, computed by each technique. Values to 5 "
                 "decimals; the **DOFs**")
    lines.append("row gives the size of each discrete problem and **Time (s)** the "
                 "measured")
    lines.append("computation time (`tic`/`toc` around each per-domain solve; Wolfram "
                 "is not timed).")
    lines.append("")

    headers = ["n"] + [c["header"] for c in cols]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")
    lines.append("| **DOFs** | " + " | ".join(fmt_dofs(c) for c in cols) + " |")
    lines.append("| **Time (s)** | " + " | ".join(fmt_time(c) for c in cols) + " |")
    for i in range(NEIG):
        row = [str(i + 1)]
        for c in cols:
            if i < len(c["eigs"]):
                v = f"{c['eigs'][i]:.5f}"
                if truth is not None and c is truth:
                    v = f"**{v}**"
            else:
                v = ""
            row.append(v)
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")

    # Footnotes (only those that apply).
    notes = []
    if any(c["kind"] == "mps" for c in cols):
        notes.append("**†  MPS** is not a mesh method — the size is the number of "
                     "Fourier–Bessel basis functions, not DOFs.")
    if any(c["kind"] == "wolfram" for c in cols):
        notes.append("**‡  Wolfram** `NDEigensystem` builds its own internal adaptive "
                     "mesh, so its DOF count is not exposed.")
    notes.append("**§  FEM solvepdeeig** reports the number of mesh nodes, not the "
                 "(constrained) degrees of freedom — unlike FEM eig, which records "
                 "`size(K,1)`.")
    lines.append(" ".join(notes))
    lines.append("")

    # Sources.
    lines.append("## Sources")
    lines.append("")
    lines.append("| Column | File |")
    lines.append("|--------|------|")
    for c in cols:
        if c["path"]:
            lines.append(f"| {c['header']} | `results/eigenvalues/{c['path']}` |")
        else:
            lines.append(f"| {c['header']} | closed-form spectrum |")
    lines.append("")
    return "\n".join(lines)


def main():
    os.makedirs(OUT, exist_ok=True)
    for domain in DOMAINS:
        cols = build_columns(domain)
        out = os.path.join(OUT, f"{domain}_comparison.md")
        with open(out, "w") as fh:
            fh.write(render(domain, cols))
        present = ", ".join(c["header"] for c in cols)
        print(f"wrote {os.path.relpath(out, ROOT)}  [{present}]")


if __name__ == "__main__":
    main()
