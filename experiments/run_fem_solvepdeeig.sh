#!/usr/bin/env bash
# Run the FEM solvepdeeig eigenvalue experiments for every domain in the
# catalog, writing one CSV per domain into results/eigenvalues/fem/. Requires the PDE
# Toolbox.
#
# To run a subset, call the driver directly from MATLAB, e.g.
#   compute_spectrum_fem_solvepdeeig("L_shaped")
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_spectrum_fem_solvepdeeig()"
