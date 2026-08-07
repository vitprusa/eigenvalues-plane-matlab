#!/usr/bin/env bash
# Recompute the paper eigenvalue-comparison spectra (DST, FD, FEM, Cheb) at four
# DOF resolutions each for the rectangle, isosceles-triangle and L-shaped
# domains, timing every run, and write the LaTeX tables into
# results_paper/eigenvalues_head/ and the per-run CSVs into its cache/
# subdirectory. Generates its own data; does not touch the existing results/
# CSVs. A run whose cached CSV is present is read back rather than recomputed --
# delete a CSV to compute that run again. Requires the PDE Toolbox and Chebfun.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="${1:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); make_eigenvalues_head_paper_tables('${name}')"
