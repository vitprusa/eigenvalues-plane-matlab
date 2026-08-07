#!/usr/bin/env bash
# Generate the rectangle eigenvalue-indexing LaTeX snippet (values + ordering
# tables) into results_paper/eigenvalues_indexing/. Defaults to M = 9; pass a
# resolution to override, e.g. ./run_eigenvalues_indexing.sh 19
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M="${1:-9}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); eigenvalues_rectangle_tables(${M})"
