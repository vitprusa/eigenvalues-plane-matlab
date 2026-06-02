#!/usr/bin/env bash
# Run the L-shaped DOF sweep (DST, FD, FEM) and write the spectra CSVs and the
# comparison plot into results/eigenvalues_dof_sweep/. Requires the PDE Toolbox.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_L_shaped_dof_sweep()"
