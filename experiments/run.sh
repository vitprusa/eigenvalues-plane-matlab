#!/usr/bin/env bash
# Run the eigenvalue experiments and write the CSV outputs into results/:
#   - compute_dst_spectra            DST-Laplacian spectra -> results/dst/
#   - compute_L_shaped_spectra_MPS   MPS L-shaped spectrum -> results/mps/
#
# To run a subset, call a driver directly from MATLAB, e.g.
#   compute_dst_spectra("L_shaped")
#   compute_dst_spectra("L_shaped", "partial")
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_dst_spectra(); compute_L_shaped_spectra_MPS()"
