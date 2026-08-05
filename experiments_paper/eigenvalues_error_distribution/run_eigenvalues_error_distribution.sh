#!/usr/bin/env bash
# Where in the spectrum each method's error sits, at one moderate resolution.
#
# Plots the relative error of every one of the first n eigenvalues against its
# index, for DST, FD and FEM, all three carrying as nearly the same number of
# degrees of freedom as the convergence sweep allows. The run table, the
# reference spectrum and the figure go into
# results_paper/eigenvalues_error_distribution/.
#
#   ./run_eigenvalues_error_distribution.sh                # first 1000, L-shape
#   ./run_eigenvalues_error_distribution.sh 500            # a shorter block
#   ./run_eigenvalues_error_distribution.sh 1000 L_shaped  # naming the domain
#
# Existing CSVs are reused, so a figure-only pass costs nothing; delete them to
# recompute. The reference is a DST run of about 20000 degrees of freedom and
# costs about two minutes the first time.
# Requires the PDE Toolbox.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
count="${1:-1000}"
name="${2:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); make_eigenvalues_error_distribution(${count}, '${name}')"
