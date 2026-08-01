#!/usr/bin/env bash
# Compute one eigenvalue, given by its index, with DST, FD and FEM over a sweep
# of resolutions and plot it against the degrees of freedom, writing the run
# table CSV and the figure into results_paper/eigenvalues_convergence/.
#
#   ./run_eigenvalues_convergence_index.sh 100             # every domain
#   ./run_eigenvalues_convergence_index.sh 100 L_shaped    # just that one
#
# Domains: rectangle, isosceles_triangle, L_shaped. The first two have a
# closed-form spectrum and are measured against it at every index; the L-shape
# is measured against the published MPS ground state, or against its own finest
# DST run away from it.
#
# An existing CSV is reused (figure-only regeneration); delete it to recompute.
# Requires the PDE Toolbox.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
index="${1:-}"
name="${2:-}"

if [ -z "${index}" ]; then
	echo "Usage: $(basename "$0") <eigenvalue-index> [domain]" >&2
	exit 2
fi

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); make_eigenvalues_convergence_index(${index}, '${name}')"
