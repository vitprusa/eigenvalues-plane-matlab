#!/usr/bin/env bash
# Compute the ground-state eigenvalue lambda_1 with DST, FD and FEM over a sweep
# of resolutions and plot it against the degrees of freedom, writing the run
# table CSV and the figure into results_paper/eigenvalues_convergence/.
# An existing CSV is reused (figure-only regeneration); delete it to recompute.
# Requires the PDE Toolbox.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="${1:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); make_eigenvalues_convergence('${name}')"
