#!/usr/bin/env bash
# Run the MPS L-shaped eigenvalue experiment, writing the eigenvalues to
# results/mps/L_shaped_eigenvalues_MPS.csv.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_L_shaped_spectra_MPS()"
