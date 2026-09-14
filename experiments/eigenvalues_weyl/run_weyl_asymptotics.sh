#!/usr/bin/env bash
# Run the Weyl-asymptotics visual check for every DST-catalog domain: compute
# the full spectrum at M_full, the domain area, and write one EPS per domain
# plus the summary weyl_areas.md into results/eigenvalues_weyl/.
#
# Optional first argument restricts the run to a single domain, e.g.
#   run_weyl_asymptotics.sh L_shaped
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
name="${1:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); plot_weyl_asymptotics('${name}')"
