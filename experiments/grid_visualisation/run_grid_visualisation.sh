#!/usr/bin/env bash
# Draw the DST grid + domain mask for each DST-catalog domain (at M_full) into
# results/grid_visualisation/. Pass a domain name to restrict the run, e.g.
#   ./run_grid_visualisation.sh L_shaped
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="${1:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); plot_grid_visualisation('${name}')"
