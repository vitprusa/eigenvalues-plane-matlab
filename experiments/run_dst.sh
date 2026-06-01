#!/usr/bin/env bash
# Run the DST-Laplacian eigenvalue experiments for every domain, writing the
# *-eigenvalues.csv files into results/dst/.
#
# To run a subset, call the driver directly from MATLAB, e.g.
#   compute_dst_spectra("L_shaped")
#   compute_dst_spectra("L_shaped", "partial")
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_dst_spectra()"
