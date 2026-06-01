#!/usr/bin/env bash
# Run the Chebfun eigenvalue experiments for every domain in the catalog and a
# sweep of Chebyshev orders N, writing one CSV per (domain, N) into
# results/cheb/. Requires Chebfun (see startup.m).
#
# To run a subset, call the driver directly from MATLAB, e.g.
#   compute_spectrum_cheb("square")
#   compute_spectrum_cheb("", [16 24 32])
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_spectrum_cheb()"
